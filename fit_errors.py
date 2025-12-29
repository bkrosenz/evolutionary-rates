import shap
import xgboost as xgb
import sklearn
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import seaborn as sns
from sys import argv

params = {
    "objective": "reg:squarederror",  # Objective function for regression (squared error)
    "n_estimators": 180,  # Number of boosting rounds (trees)
    "learning_rate": 0.2,  # Step size shrinkage to prevent overfitting
    "max_depth": 8,  # Maximum depth of a tree
    "subsample": 1,  # Subsample ratio of the training instance
    "colsample_bytree": 1,  # Subsample ratio of columns when constructing each tree
    "random_state": 42,  # Seed for reproducibility
    "n_jobs": -1,  # Use all available CPU cores
}

# Initialize the XGBoost Regressor model
estimators = {
    "sigma": r"$\hat{\sigma}_{PIC}$",
    "sigma_js": r"$\hat{\sigma}_{JS}$",
    "sigma_cherry": r"$\hat{\sigma}_{cherry}$",
    "sigma_paired": r"$\hat{\sigma}_{paired}$",
    "sigma_ape": r"$\hat{\sigma}^2_{ape}$",
    "sigma_geig": r"$\hat{\sigma}^2_{geiger}$",
}

x = pd.read_csv(
    "/N/project/phyloML/rate_timescaling/data/pic_dendropy_predictions_no_ils.csv.gz",
    index_col=0,
).drop(
    columns=["sigma_js_old", "tree_ix"], errors="ignore"
)  # .reset_index()

query = argv[1]
print(x.columns, "\n", query)
# query = "h>1 & 4<n<1000 & bl_min<1"
# query = 'n>4'

n_samps = 2000
x = x.dropna().query(query)
X_train = x.drop(columns=estimators.keys())

qstr = query.replace(" ", "")

figdir = Path(f"figures/{qstr}")
figdir.mkdir(exist_ok=True, parents=True)
feat_map={
            "bl_min": r"$bl_{min}$",
            "bl_max": r"$bl_{max}$",
            "bl_median": r"$bl_{median}$",
            "bl_mean": r"$bl_{mean}$",
            "bl_std": r"$bl_{std}$",
            "bl_q1": r"$bl_{Q1}$",
            "bl_q3": r"$bl_{Q3}$",
            "tau": r"$\tau$",
            "rho": r"$\rho$",
            "mu": r"$\mu$",
            "lam": r"$\lambda$",
            "h": "H",
            "mu_over_lam": r"$\mu/\lambda$",
            "h_mu_over_lam": r"$H\mu/\lambda$",
            "lam_minus_mu": r"$\lambda-\mu$",
            "h_lam_minus_mu": r"$H(\lambda-\mu)$",
        }

feat_names = list(
    X_train.columns.to_series().replace(
        feat_map
    )
)

with open("feature_importances/" + qstr + ".txt", "w") as logfile:
    logfile.write(f"query: {query}\nN={len(x)}\n")

    for var in estimators:
        model = xgb.XGBRegressor(**params)
        y_train = x[var]  # (x[var] - 1).abs()
        # Train the model on the training data
        model.fit(X_train, y_train)
        logfile.write(var)

        # XGBoost provides a way to get feature importance
        logfile.write(
            "\nFeature Importances (fscore):\n",
        )
        feature_importances = pd.Series(
            model.feature_importances_, 
            index=X_train.columns, 
            name="Feature Importance"
        ).sort_values(ascending=False)
        logfile.write(str(feature_importances))

        feature_importances = feature_importances.rename(index=feat_map)

        plt.figure(figsize=(10, 6))
        sns.barplot(x=feature_importances.values, 
                    y=feature_importances.index) # TODO: fix this so that feat_names actually map to features
        plt.title(estimators[var])
        plt.ylabel("Feature")
        plt.xlabel("Importance (Gain)")
        plt.savefig(figdir / f"xgb_features_gain_{var}.png")
        plt.clf()
        # plt.savefig(f"figures/xgb_features_{var}_mae.png")
        y_pred = model.predict(X_train)

        mse = mean_squared_error(y_train, y_pred)
        rmse = np.sqrt(mse)
        r2 = r2_score(y_train, y_pred)

        logfile.write(f"\nModel Evaluation:")
        logfile.write(f"\nMean Squared Error (MSE): {mse:.2f}")
        logfile.write(f"\nRoot Mean Squared Error (RMSE): {rmse:.2f}")
        logfile.write(f"\nR-squared (R2) Score: {r2:.2f}\n---\n")
        explainer = shap.TreeExplainer(model, data=X_train.sample(n_samps))
        X_sample = X_train.sample(n_samps).rename(columns=feat_map)
        # X_sample
        shap_values = explainer.shap_values(X_sample)
        shap.summary_plot(shap_values, X_sample)# feature_names=feat_names)
        plt.title(estimators[var])
        plt.savefig(figdir / f"xgb_shap_{var}.png")
