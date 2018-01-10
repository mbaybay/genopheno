import numpy as np
import pandas as pd
import math


"""
The different possible mutation levels.
nm = no mutations
pm = partial mutations
fm = full mutations
"""
MUTATION_LEVELS = ['nm', 'pm', 'fm']


def __remove_missing_data(pheno, snp_data, invalid_thresh):
    """
    Removes missing data from the user data. If a SNP row has a percentage of users with an invalid or missing genotype
    then the SNP row is removed. Missing SNP data is changed from its numeric missing value to NaN.
    :param snp_data: The SNP data for all users
    :param invalid_thresh: The maximum percentage of invalid data for a row or column
    :return: The SNP data for all users with the missing data removed
    """
    non_user_columns = ['Gene_info', 'pct_fm', 'pct_nm', 'pct_pm']
    user_columns = set(snp_data.columns.values).difference(non_user_columns)
    users_count = len(user_columns)

    snp_count = snp_data.shape[0]

    min_required = math.ceil((1 - invalid_thresh / float(100)) * users_count)
    snp_data.dropna(axis=0, thresh=min_required, inplace=True, subset=user_columns)
    print "{} SNPs removed due to too many missing user observations for phenotype '{}'"\
        .format(snp_count - snp_data.shape[0], pheno)


def __filter_snps(row, diff_thresh, magnitude_thresh, selected_snps):
    """
    Determines if there is a significant mutation difference at each SNP between the two phenotype options.
    Each SNP contains the percentage of users with no, partial and full mutations. Each mutation level is compared
    respectively. If the mutation difference is greater than the defined difference and magnitude thresholds at any
    level, then the SNP is selected.
    :param row: The DataFrame row containing the mutation percentages. The index is the RSID.
    :param diff_thresh: The difference threshold required for the SNP to be selected, in percentage points.
    :param magnitude_thresh:  The magnitude of the difference threshold, in percentage points. This is the
                              minimum value for the difference in mutation percentage divided by the minimum
                              mutation value out of the two phenotypes. The purpose of this is to filter out
                              SNPs where the change meets the difference threshold, but is still a small
                              magnitude. For example, if the mutation difference is 5%, but the SNP mutation
                              levels for each phenotype are 100% and 95% then this is less meaningful than if
                              the mutation levels were 6% and 1%. The magnitude threshold is meant to filter
                              out the SNP where the mutations are 100% and 95% and keep the SNP where the
                              mutations are 6% and 1%.
    :param selected_snps: A list containing the selected SNP RSIDs. If the row has a significant mutation difference,
                          then the RSID for this row will be added to the list.
    """
    snp_rsid = row.name

    for mutation_level in MUTATION_LEVELS:
        # Get the mutation percentages at the mutation level for each phenotype
        mutation_pct_a = np.round(row['pct_{}_a'.format(mutation_level)], 3)
        mutation_pct_b = np.round(row['pct_{}_b'.format(mutation_level)], 3)

        # Calculate the mutation difference between phenotypes
        pct_diff = abs(mutation_pct_a - mutation_pct_b)

        # Check if the difference is greater than the defined difference threshold
        if pct_diff >= diff_thresh:
            pct_min = min(mutation_pct_a, mutation_pct_b)
            if pct_min != 0:
                # Check if the difference is greater than the defined magnitude threshold
                if 100 * pct_diff / pct_min >= magnitude_thresh:
                    selected_snps.append(snp_rsid)
                    break
            else:
                selected_snps.append(snp_rsid)
                break


def __identify_mutated_snps(phenotypes, diff_thresh, magnitude_thresh):
    """
    Identifies SNPs of interest based on mutation differences between the two phenotypes.
    :param phenotypes: A dictionary of phenotypes where the key is the phenotype label and the value is the
                       DataFrame with the mutation information for all users with the phenotype.
    :param diff_thresh: The minimum difference between SNP mutations for each phenotype to be selected.
                        See the method __filter_snps for a more complete explanation.
    :param magnitude_thresh: The minimum magnitude in the mutation difference between SNPs for each phenotype.
                             See the method __filter_snps for a more complete explanation.
    :return: A list of selected SNP RSIDs.
    """
    if len(phenotypes) != 2:
        raise ValueError('Using mutation differences to identify SNPs is only valid for two phenotype options.')

    # Merge the data frames top identify common SNPs and be able to do selection in one pass
    pheno_df_a = phenotypes.values()[0]
    pheno_df_b = phenotypes.values()[1]
    mutation_columns = ['pct_nm', 'pct_pm', 'pct_fm']
    merged = pheno_df_a[mutation_columns].merge(pheno_df_b[mutation_columns], left_index=True, right_index=True,
                                                suffixes=('_a', '_b'))

    # Identify selected SNPs
    selected_snps = []
    merged.apply(__filter_snps, args=(diff_thresh, magnitude_thresh, selected_snps), axis=1)

    return selected_snps


def __format_selected_snps(pheno_label, pheno_df, selected_snps):
    """
    Builds the phenotype DataFrame used for the machine learning model based on the selected SNPs
    :param pheno_df: The DataFrame for the phenotype with the user mutation data
    :param selected_snps: A list of selected SNP RSIDs
    :param pheno_label: The phenotype label (i.e. 'Brown' for eye color)
    :return: The DataFrame for the selected SNPs
    """
    # Filter out SNPs that have not been selected
    snp_data = pheno_df.loc[selected_snps]

    # Drop the mutation percentages and gene info columns because they are no longer needed
    snp_data.drop(labels=['Gene_info', 'pct_fm', 'pct_nm', 'pct_pm'], axis=1, inplace=True)

    # Transpose the data and add columns for user Id and phenotype
    transposed_data = snp_data.transpose()
    transposed_data['phenotype'] = pheno_label

    return transposed_data


def create_dataset(phenotypes, invalid_thresh, invalid_user_thresh, diff_thresh, magnitude_thresh):
    """
    Function to return those SNPs that satisfy a criterion to check for differences between blue and brown SNPs
    :param invalid_thresh The percentage of missing user observations a SNP can have before it is removed
    :param invalid_user_thresh The acceptable percentage of missing data before a user is discarded
    :param diff_thresh: Minimum difference in percentage points of mutated alleles between different phenotypes in order
    to include a SNP in set of candidates
    :param magnitude_thresh: # Minimum percent difference in percentage of mutated alleles between different phenotypes
    for a SNP in order to include a SNP in set of candidates
    :return: A DataFrame where each row is a user and each column is a SNP.
    The value is the number of mutations (0,1,2).
    """
    # Filter out SNPs that do not have enough user observations
    for pheno, pheno_df in phenotypes.iteritems():
        __remove_missing_data(pheno, pheno_df, invalid_thresh)

    # Select snps based on mutation differences between phenotypes
    selected_snps = __identify_mutated_snps(phenotypes, diff_thresh, magnitude_thresh)
    print '{} SNPs with mutation differences identified'.format(len(selected_snps))

    # Generate data frame for each phenotype using the selected SNPs
    final_datasets = []
    for pheno_key, pheno_df in phenotypes.items():
        final_datasets.append(__format_selected_snps(pheno_key, pheno_df, selected_snps))

    # Merge and return aggregate data set
    merged = pd.concat(final_datasets)

    # Remove users that do not have enough observations
    user_count = merged.shape[0]
    snp_count = merged.shape[1] - 1
    min_obs = math.ceil((1 - invalid_user_thresh / float(100)) * snp_count)
    merged.dropna(axis=0, thresh=min_obs, inplace=True)
    print '{} users dropped due to too many missing observations. Remaining users: {}'\
        .format(user_count - merged.shape[0], merged.shape[0])

    return merged
