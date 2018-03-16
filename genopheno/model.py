import argparse
import os
import re

import pandas as pd

from genopheno.models.snp_selectors import mutation_difference
from models import elastic_net, decision_tree, random_forest
from util import timed_invoke, expand_path, clean_output

MODELS = {
    'en': elastic_net.build_model,
    'dt': decision_tree.build_model,
    'rf': random_forest.build_model
}


def __read_phenotype_input(input_dir):
    """
    Reads the preprocessed phenotype files from the initialization steps and creates data frames from each file.
    :param input_dir: The directory containing the preprocessed files.
    :return: A map of phenotypes where the key is the phenotype ID and the value is the phenotype data frame.
    """
    # read preprocessed files
    file_prefix = 'preprocessed_'
    file_name_regex = re.compile('^{}.+\.csv.gz$'.format(file_prefix))
    files = os.listdir(input_dir)

    phenotypes = {}
    for f in filter(file_name_regex.match, files):
        # create data frame from preprocessed files
        df = pd.read_csv(os.path.join(input_dir, f), compression='gzip')
        df.set_index('Rsid', inplace=True)

        # add the data frame to the collection of preprocessed phenotypes
        phenotype = f[len(file_prefix):len(f) - len('.csv.gz')]
        phenotypes[phenotype] = df
        print "{} users and {} SNPs for phenotype '{}'".format(len(df.columns)-4, df.shape[0], phenotype)

    if len(phenotypes) == 0:
        raise ValueError('No preprocessed files in directory "{}". '
                         'This directory should contain the output from the preprocess step.'.format(input_dir))

    return phenotypes


def run(preprocessed_dir, invalid_thresh, invalid_user_thresh, abs_diff_thresh, relative_diff_thresh, data_split,
        no_interactions, negative, max_snps, model_id, output_dir):
    """
    Builds a model to predict phenotype
    :param preprocessed_dir: The directory containing the preprocessed data
    :param invalid_thresh: The acceptable percentage of missing data before a SNP is discarded
    :param invalid_user_thresh: The acceptable percentage of missing data before a user is discarded
    :param abs_diff_thresh: The mutation percent difference between phenotypes to be selected as a model feature
    :param relative_diff_thresh: The relative difference in mutation percentage, calculated as a percent of the
                                smaller mutation percent value.
    :param data_split: The percent data used for testing.
    :param no_interactions: If True the model will not contain interactions
    :param negative: The negative phenotype label
    :param model_id: The id for the model to use
    :param output_dir: The directory to write the model in
    """
    # Expand file paths
    preprocessed_dir = expand_path(preprocessed_dir)
    output_dir = expand_path(output_dir)

    # Make sure output directory exists before doing work
    clean_output(output_dir)

    # Get model
    build_model = MODELS.get(model_id)
    if not build_model:
        raise ValueError('Model Id "{}" is not valid'.format(model_id))

    phenotypes = timed_invoke('reading the preprocessed files', lambda: __read_phenotype_input(preprocessed_dir))
    data_set = timed_invoke('creating model data set', lambda: mutation_difference.create_dataset(
                               phenotypes, invalid_thresh, invalid_user_thresh, abs_diff_thresh, relative_diff_thresh)
                            )
    timed_invoke('building model', lambda: build_model(data_set, data_split, no_interactions, negative, max_snps,
                                                       output_dir))
    print 'Output written to "{}"'.format(output_dir)


if __name__ == '__main__':
    # Parse input
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--preprocessed",
        "-p",
        metavar="<directory path>",
        default="resources/data/preprocessed",
        help="The directory containing the output data from the initialization phase."
             "\n\nDefault: resources/data/preprocessed"
    )

    parser.add_argument(
        "--invalid-snp-thresh",
        "-it",
        metavar="percent",
        type=float,
        default=40,
        help="The maximum percentage of missing or invalid user observations a SNP can have before it is not "
             "considered as a feature in the model."
             "\n\nDefault: 40"
    )

    parser.add_argument(
        "--invalid-user-thresh",
        "-iu",
        metavar="percent",
        type=float,
        default=90,
        help="The maximum percentage of missing or invalid SNP observations a user can have before it is not "
             "considered as a valid example in the model."
             "\n\nDefault: 90"
    )

    parser.add_argument(
        "--absolute-diff-thresh",
        "-adt",
        metavar="percent",
        type=float,
        default=20,
        help="The difference threshold required for the SNP to be selected, in percentage points."
             "\n\nDefault: 20"
    )

    parser.add_argument(
        "--relative-diff-thresh",
        "-rdt",
        metavar="percent",
        type=float,
        default=5,
        help="The relative difference threshold. This is the absolute difference in mutation percentage divided "
             "by the minimum mutation value out of the two phenotypes. The purpose of this is to filter out "
             "SNPs where the change meets the difference threshold, but is still a small "
             "magnitude. For example, if the mutation difference is 5 percentage points, but the SNP mutation "
             "levels for each phenotype are 100% and 95% then this is less meaningful than if "
             "the mutation levels were 6% and 1%. The magnitude threshold is meant to filter "
             "out the SNP where the mutations are 100% and 95% and keep the SNP where the "
             "mutations are 6% and 1%."
             "\n\nDefault: 5"
    )

    parser.add_argument(
        "--split",
        "-s",
        metavar="percent",
        type=float,
        default=33,
        help="The percentage of users to use as test data. The remaining users will be used for training the model."
             "\n\nDefault: 33"
    )

    parser.add_argument(
        "--no-interactions",
        "-ni",
        default=False,
        action='store_true',
        help="If set then interactions will not be included in the model. Note, this does not apply for the decision "
             "tree model."
             "\n\nDefault: False"
    )

    parser.add_argument(
        "--negative",
        "-n",
        help="The phenotype value that should be considered as the negative case. "
             "If not set, the alphabetically first phenotype value will be used as the negative value."
    )

    parser.add_argument(
        "--max-snps",
        "-ms",
        default=None,
        type=int,
        help="The maximum number of SNPs to include in the model"
             "\n\nDefault: None"
    )

    parser.add_argument(
        "--model",
        "-m",
        default="en",
        type=str,
        help="The type of model to use."
             "\nen = Elastic net"
             "\ndt = Decision tree"
             "\n\n Default: en"
    )

    parser.add_argument(
        "--output",
        "-o",
        metavar="<directory path>",
        default="resources/data/model",
        help="The directory that the output files should be written to. This will include all files required for the "
             "machine learning input."
             "\n\nDefault: resources/data/model"
    )

    args = parser.parse_args()
    run(args.preprocessed, args.invalid_snp_thresh, args.invalid_user_thresh, args.absolute_diff_thresh,
        args.relative_diff_thresh, args.split, args.no_interactions, args.negative, args.max_snps, args.model, args.output)
