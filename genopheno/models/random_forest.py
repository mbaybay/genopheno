import common
import pandas as pd

from sklearn.ensemble import RandomForestClassifier


def build_model(dataset, data_split, no_interactions, negative, max_snps, output_dir):
    model_eval = {
        'features': save_features
    }

    param_grid = {
        "criterion": ["gini", "entropy"],
        "n_estimators": [500, 1000, 3000],
        "max_depth": [None, 7, 8],
        "max_features": ["sqrt", .3, .4]        # sqrt is known for use in classification; same as "auto" config
    }

    common.build_model(
        dataset,
        data_split,
        True,
        negative,
        RandomForestClassifier(n_jobs=-1),
        max_snps,
        output_dir,
        param_grid=param_grid,
        model_eval=model_eval
    )


def save_features(model, model_terms, output_dir):
    ftrs = pd.DataFrame()
    ftrs['Feature'] = model_terms
    ftrs['Importance'] = model.feature_importances_
    ftrs.sort_values(by='Importance', ascending=False, inplace=True)
    ftrs.to_csv(output_dir + "rf_features.csv")
