import common
from sklearn.ensemble import RandomForestClassifier


def build_model(dataset, data_split, output_dir):
    common.build_model(
        dataset, data_split, RandomForestClassifier(n_estimators=3000, max_features='auto', n_jobs=-1), output_dir
    )

