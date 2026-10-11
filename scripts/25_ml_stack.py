"""Stacked machine-learning covariate for the step 2 timeliness model (after Bhatt et al. 2017).

Three learners predict the per-contact slip probability of each DHS cluster from the 16 covariates
of 24_ml_covariates.R: random forest, gradient boosting and ridge-penalised logistic regression. Each is
fitted to binomial data (one row for the slipped and one for the timely contacts of each cluster,
weighted by their counts). The stacked covariate is the mean of the three learners' logits. It
enters the INLA-SPDE model as a covariate, alongside the spatial field (26_ml_stack_cv.R).

To avoid leakage the covariate a cluster receives never comes from a learner that saw it:
  full    10-fold out-of-fold predictions for all clusters (for the full stacked INLA fit), and
          learners refitted to all clusters for the northern grid (for step 4)
  lodo    for each northern district d (as 22_district_cv.R): learners fitted to all clusters
          outside d predict d's clusters; the training clusters get 5-fold out-of-fold predictions
          from inside that training set
Also reported: the learners alone (no spatial field) under leave-one-district-out, so that machine
learning can be compared directly with the spatial models.

Outputs: data/derived/ml_oof_full.csv, ml_grid_pred.csv, ml_lodo.csv, ml_learner_cv.csv
"""
import numpy as np
import pandas as pd
from scipy.stats import binom
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SEED = 20261010
COV = ["tt_walk", "tt_motor", "logpopdens", "fac_5km", "chps_5km", "viirs", "dst_road", "dst_roadx",
       "dst_water", "dst_urban", "dst_crop", "dst_tree", "slope", "elevation", "rain_annual", "rain_wetshare"]

cv = pd.read_csv("data/derived/ml_cluster_covariates.csv")
tl = pd.read_csv("data/derived/dhs_clusters_timeliness.csv")[["cluster", "contacts", "slip56"]]
d = cv.merge(tl, on="cluster")
d = d[(d.contacts > 0) & d[COV].notna().all(axis=1)].reset_index(drop=True)
print(f"clusters with contacts and covariates: {len(d)} | northern with a district: {d.district.notna().sum()}")

LEARNERS = {
    "rf": lambda: RandomForestClassifier(n_estimators=500, min_samples_leaf=10, max_features=0.4,
                                         n_jobs=4, random_state=SEED),
    "gbm": lambda: HistGradientBoostingClassifier(learning_rate=0.03, max_iter=200, max_leaf_nodes=8,
                                                  min_samples_leaf=20, l2_regularization=1.0, random_state=SEED),
    "ridge": lambda: make_pipeline(StandardScaler(), LogisticRegression(C=0.5, max_iter=2000)),
}
logit = lambda p: np.log(np.clip(p, 1e-4, 1 - 1e-4) / (1 - np.clip(p, 1e-4, 1 - 1e-4)))


def binom_rows(df):
    """Two weighted rows per cluster: slipped (y=1) and timely (y=0) contacts."""
    X = pd.concat([df[COV], df[COV]], ignore_index=True)
    y = np.r_[np.ones(len(df)), np.zeros(len(df))]
    w = np.r_[df.slip56.values, (df.contacts - df.slip56).values].astype(float)
    keep = w > 0
    return X[keep], y[keep], w[keep]


def fit_predict(train, test):
    """Each learner's logit for the test clusters, plus their mean (the stacked covariate)."""
    X, y, w = binom_rows(train)
    out = {}
    for nm, make in LEARNERS.items():
        m = make()
        if nm == "ridge":
            m.fit(X, y, logisticregression__sample_weight=w)
        else:
            m.fit(X, y, sample_weight=w)
        out[nm] = logit(m.predict_proba(test[COV])[:, 1])
    out["stack"] = np.mean([out[k] for k in LEARNERS], axis=0)
    return out


def oof(df, k, seed):
    res = {nm: np.zeros(len(df)) for nm in list(LEARNERS) + ["stack"]}
    for tr, te in KFold(n_splits=k, shuffle=True, random_state=seed).split(df):
        p = fit_predict(df.iloc[tr], df.iloc[te])
        for nm in res:
            res[nm][te] = p[nm]
    return res


# ---- full: out-of-fold for all clusters; grid from learners on all clusters -------------
o = oof(d, 10, SEED)
full = pd.DataFrame({"cluster": d.cluster, **{f"ml_{k}": v for k, v in o.items()}})
full.round(5).to_csv("data/derived/ml_oof_full.csv", index=False)
g = pd.read_csv("data/derived/ml_grid_covariates.csv")
ok = g[COV].notna().all(axis=1)
gp = fit_predict(d, g[ok])
grid = pd.DataFrame({"cell": g.cell[ok], "ml": gp["stack"]})
grid.round(5).to_csv("data/derived/ml_grid_pred.csv", index=False)
print(f"grid cells predicted: {ok.sum()} of {len(g)}")

# ---- lodo: nested, per northern district ------------------------------------------------
rows, cvrows = [], []
for dist in sorted(d.district.dropna().unique()):
    te = d.district == dist
    train, test = d[~te].reset_index(drop=True), d[te].reset_index(drop=True)
    inner = oof(train, 5, SEED)
    outer = fit_predict(train, test)
    rows.append(pd.DataFrame({"held_out": dist, "cluster": train.cluster, "role": "train", "ml": inner["stack"]}))
    rows.append(pd.DataFrame({"held_out": dist, "cluster": test.cluster, "role": "test", "ml": outer["stack"]}))
    for nm, v in outer.items():                               # learners alone, district level
        p = 1 / (1 + np.exp(-v))
        cvrows.append({"learner": nm, "district": dist, "n": test.contacts.sum(), "y": test.slip56.sum(),
                       "pred": (p * test.contacts).sum() / test.contacts.sum(),
                       "ll": binom.logpmf(test.slip56, test.contacts, p).sum()})
    print("  held out", dist)
pd.concat(rows).round(5).to_csv("data/derived/ml_lodo.csv", index=False)
lc = pd.DataFrame(cvrows)
lc["obs"] = lc.y / lc.n
lc.to_csv("data/derived/ml_learner_cv.csv", index=False)
def wcor(x, y, w):
    w = w / w.sum(); mx, my = (w * x).sum(), (w * y).sum()
    return (w * (x - mx) * (y - my)).sum() / np.sqrt((w * (x - mx) ** 2).sum() * (w * (y - my) ** 2).sum())
s = lc.groupby("learner").apply(lambda r: pd.Series({
    "logscore": -r.ll.sum() / r.n.sum(),
    "mae": (r.n * (r.pred - r.obs).abs()).sum() / r.n.sum(),
    "corr": wcor(r.pred.values, r.obs.values, r.n.values.astype(float))}), include_groups=False)
s.round(4).to_csv("data/derived/ml_learner_cv_summary.csv")
print("\nLearners alone, leave-one-district-out (same scores as 22_district_cv.R):")
print(s.round(4))
