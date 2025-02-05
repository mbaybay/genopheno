import shutil
import gzip
import filecmp
from genopheno import preprocess, model, predict
from os import mkdir
from os.path import join, exists


EXPECTED_DIR = join('resources', 'expected')
EXPECTED_PREPROCESS = join(EXPECTED_DIR, 'preprocess')
EXPECTED_MODEL = join(EXPECTED_DIR, 'model')
EXPECTED_PREDICT = join(EXPECTED_DIR, 'predict')
DATA_DIR = join('..', 'genopheno', 'resources', 'data')
USER_DATA_DIR = join(DATA_DIR, 'users')
SNP_DATA_DIR = join(DATA_DIR, 'snp')
KNOWN_PHENOTYPES_FILE = join(DATA_DIR, 'known_phenotypes.csv')
OUTPUT_DIR_PREPROCESS = 'output_preprocess'
OUTPUT_DIR_MODEL = 'output_model'
OUTPUT_DIR_PREDICT = 'output_predict'


def test_preprocess_eye_color():
    """
    Tests the preprocessing step.
    """
    __create_dir(OUTPUT_DIR_PREPROCESS)
    preprocess.run(USER_DATA_DIR, SNP_DATA_DIR, KNOWN_PHENOTYPES_FILE, OUTPUT_DIR_PREPROCESS)

    for phenotype in ['preprocessed_Blue_Green.csv.gz', 'preprocessed_Blue_Green.csv.gz', 'snp_database.csv.gz']:
        with gzip.open(join(OUTPUT_DIR_PREPROCESS, phenotype)) as exp, \
                gzip.open(join(EXPECTED_PREPROCESS, phenotype)) as act:
            assert exp.readlines() == act.readlines()


def test_model_eye_color():
    """
    Tests the modeling step.
    """
    __create_dir(OUTPUT_DIR_MODEL)
    model.run(EXPECTED_PREPROCESS, 50, 80, 15, 33, False, None, 200, 'en', OUTPUT_DIR_MODEL)
    __assert_files(
        ['confusion_matrix.txt', 'features.csv', 'model_config.pkl', 'roc.png'], EXPECTED_MODEL, OUTPUT_DIR_MODEL
    )


def test_predict_eye_color():
    """
    Tests the prediction step.
    """
    __create_dir(OUTPUT_DIR_PREDICT)
    predict.run(USER_DATA_DIR, EXPECTED_PREPROCESS, EXPECTED_MODEL, OUTPUT_DIR_PREDICT)
    __assert_files(['predictions.csv'], EXPECTED_PREDICT, OUTPUT_DIR_PREDICT)


def __assert_files(files, exp_dir, act_dir):
    """
    Verifies that the files in both directories are the same
    :param files: An array of file names to compare
    :param exp_dir: The directory with the expected files
    :param act_dir: The directory with the actual files generated in the test
    """
    for filename in files:
        assert filecmp.cmp(join(exp_dir, filename),
                           join(act_dir, filename))


def __create_dir(directory):
    """
    Creates a directory. If present, the old directory and all files within the directory will be deleted.
    :param directory: The directory to create.
    """
    if exists(directory):
        shutil.rmtree(directory)
    mkdir(directory)
