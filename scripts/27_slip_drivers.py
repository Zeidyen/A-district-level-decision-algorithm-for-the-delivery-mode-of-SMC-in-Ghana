"""What predicts a slipped routine contact? Gradient boosting + SHAP on DHS 2022 contacts.

The main finding is that routine contacts slip even close to a CHPS, so distance is not the main
barrier. This asks which child, mother, household and place characteristics predict a slip.

Data: data/derived/slip_contacts_features.csv (27a_slip_features.R); one row per Penta 2/3
contact, all Ghana (the model learns from all regions; SHAP values are summarised for northern
contacts). Outcome: slip (gap > 56 days or dose never given).

Models, compared by 5-fold cluster-grouped cross-validation (a cluster's contacts never split
between training and test):
  travel time only   logistic regression on log2 motorised travel time
  logistic, all      logistic regression on all features (median-imputed, standardised)
  gradient boosting  HistGradientBoostingClassifier on all features (native missing values), its
                     settings chosen by an inner cluster-grouped CV within each training fold
                     (nested, so the outer score is not optimistic)
Reported: AUC and log loss on held-out folds.
SHAP (TreeExplainer) on the gradient-boosting model refitted to all contacts: mean |SHAP| per
feature in northern contacts (log-odds scale) and the direction of each effect (Spearman between
feature value and SHAP value). These are associations, not causes.

Outputs: data/derived/slip_ml_cv.csv, slip_shap_importance.csv, slip_shap_groups.csv,
         slip_shap_by_region.csv, slip_shap_groups_by_region.csv; out/slip_shap.png, slip_shap_by_region.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import shap
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss, roc_auc_score
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SEED = 20261010
d = pd.read_csv("data/derived/slip_contacts_features.csv")
d["log2_tt_motor"] = np.log2(np.maximum(d.tt_motor_any, 5))
d["log2_tt_walk"] = np.log2(np.maximum(d.tt_walk_any, 5))
d["contact3"] = (d.contact == 3).astype(int)
regions = sorted(d.region.unique())
for r in regions:
    d["reg_" + r.replace(" ", "_")] = (d.region == r).astype(int)

LABEL = {"female": "Child is a girl", "birth_order": "Birth order", "birth_month": "Month of birth",
         "twin": "Twin", "mother_age": "Mother's age", "mother_edu": "Mother's education",
         "mother_working": "Mother working", "insured": "Mother has health insurance",
         "media": "Media exposure (number of media)", "problem_permission": "Problem: permission to seek care",
         "problem_money": "Problem: money for care", "problem_distance": "Problem: distance to facility",
         "problem_alone": "Problem: not wanting to go alone", "anc_visits": "Antenatal visits",
         "facility_delivery": "Born in a health facility", "wealth": "Household wealth quintile",
         "hh_size": "Household size", "under5_in_hh": "Children under 5 in household",
         "head_female": "Female household head", "urban": "Urban", "log2_tt_motor": "Motorised travel time (log2)",
         "log2_tt_walk": "Walking travel time (log2)", "logpopdens": "Population density (log)",
         "contact3": "Third dose (vs second)"}
REG = [c for c in d.columns if c.startswith("reg_")]
FEAT = list(LABEL) + REG
for c in REG:
    LABEL[c] = "Region: " + c[4:].replace("_", " ").title()

X, y, g = d[FEAT], d.slip.values, d.cluster.values
cv = GroupKFold(n_splits=5)
folds = list(cv.split(X, y, g))

models = {
    "travel time only": (["log2_tt_motor"], lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))),
    "logistic, all features": (FEAT, lambda: make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                                                           LogisticRegression(max_iter=2000, C=0.5))),
    "gradient boosting": (FEAT, None),
}
GRID = {"learning_rate": [0.02, 0.05], "max_iter": [100, 300], "max_leaf_nodes": [4, 8, 16],
        "min_samples_leaf": [50, 150]}
def tuned_gb(Xtr, ytr, gtr):
    gs = GridSearchCV(HistGradientBoostingClassifier(early_stopping=False, l2_regularization=1.0, random_state=SEED),
                      GRID, scoring="neg_log_loss", cv=GroupKFold(n_splits=4), n_jobs=4)
    gs.fit(Xtr, ytr, groups=gtr)
    return gs.best_estimator_, gs.best_params_
rows = []
for name, (cols, make) in models.items():
    oof = np.zeros(len(d))
    for tr, te in folds:
        if make is None:
            m, _ = tuned_gb(X.iloc[tr][cols], y[tr], g[tr])
        else:
            m = make().fit(X.iloc[tr][cols], y[tr])
        oof[te] = m.predict_proba(X.iloc[te][cols])[:, 1]
    for sub, mask in (("all Ghana", np.ones(len(d), bool)), ("north", d.north.values)):
        rows.append({"model": name, "subset": sub, "auc": roc_auc_score(y[mask], oof[mask]),
                     "logloss": log_loss(y[mask], oof[mask])})
res = pd.DataFrame(rows).round(4)
res.to_csv("data/derived/slip_ml_cv.csv", index=False)
print("Cluster-grouped 5-fold CV (higher AUC, lower log loss is better):")
print(res.to_string(index=False))

# SHAP on the gradient-boosting model refitted to all contacts.
gb, best = tuned_gb(X, y, g)
print("\ngradient boosting settings chosen on all contacts:", best)
ex = shap.TreeExplainer(gb)
sv = ex.shap_values(X)
sv = sv[1] if isinstance(sv, list) else sv
north = d.north.values
imp = []
for j, f in enumerate(FEAT):
    v, s = X[f].values[north], sv[north, j]
    ok = ~np.isnan(v)
    rho = spearmanr(v[ok], s[ok]).correlation if ok.sum() > 10 and np.nanstd(v[ok]) > 0 else np.nan
    imp.append({"feature": f, "label": LABEL[f], "mean_abs_shap_north": np.abs(s).mean(),
                "mean_abs_shap_all": np.abs(sv[:, j]).mean(), "direction_rho": rho})
imp = pd.DataFrame(imp).sort_values("mean_abs_shap_north", ascending=False)
imp["share_north"] = imp.mean_abs_shap_north / imp.mean_abs_shap_north.sum()
imp.round(4).to_csv("data/derived/slip_shap_importance.csv", index=False)
print("\nSHAP importance, northern contacts (mean |SHAP|, log-odds; direction = Spearman of value vs SHAP):")
print(imp[["label", "mean_abs_shap_north", "share_north", "direction_rho"]].round(3).head(15).to_string(index=False))

GROUP = {"Child and schedule": ["contact3", "female", "birth_order", "birth_month", "twin"],
         "Mother": ["mother_age", "mother_edu", "mother_working", "insured", "media", "problem_permission",
                    "problem_money", "problem_distance", "problem_alone"],
         "Care history": ["anc_visits", "facility_delivery"],
         "Household": ["wealth", "hh_size", "under5_in_hh", "head_female"],
         "Place: access and density": ["log2_tt_motor", "log2_tt_walk", "logpopdens", "urban"],
         "Place: region": REG}
grp = pd.DataFrame([{"group": k, "share_north": imp.loc[imp.feature.isin(v), "share_north"].sum()} for k, v in GROUP.items()])
grp.round(3).to_csv("data/derived/slip_shap_groups.csv", index=False)
print("\nSHAP share by group, northern contacts:"); print(grp.round(3).to_string(index=False))

# ---- by region: same model, SHAP summarised within each northern region ------------------
# One model for all Ghana (350-620 contacts per northern region is too few to fit each region
# separately); its SHAP values are summarised region by region. Shares are of the region's total
# mean |SHAP|; "effect" is the mean SHAP difference between the top and bottom of a feature
# (top vs bottom quartile, or 1 vs 0 for binary features), in log-odds.
NREG = ["northern", "savannah", "north east", "upper east", "upper west"]
reg_rows, grp_rows = [], []
for r in NREG:
    m = (d.region == r).values
    tot = np.abs(sv[m]).mean(axis=0).sum()
    for j, f in enumerate(FEAT):
        v, s_ = X[f].values[m], sv[m, j]
        ok = ~np.isnan(v)
        lo, hi = (np.nanquantile(v, 0.25), np.nanquantile(v, 0.75)) if len(np.unique(v[ok])) > 2 else (0, 1)
        eff = s_[ok & (v >= hi)].mean() - s_[ok & (v <= lo)].mean() if (ok & (v >= hi)).any() and (ok & (v <= lo)).any() and hi > lo else np.nan
        reg_rows.append({"region": r, "feature": f, "label": LABEL[f], "share": np.abs(s_).mean() / tot, "effect": eff})
    for k, v in GROUP.items():
        idx_ = [FEAT.index(f) for f in v]
        grp_rows.append({"region": r, "group": k, "share": np.abs(sv[m][:, idx_]).mean(axis=0).sum() / tot})
rr = pd.DataFrame(reg_rows); rg = pd.DataFrame(grp_rows)
rr.round(4).to_csv("data/derived/slip_shap_by_region.csv", index=False)
rg.round(4).to_csv("data/derived/slip_shap_groups_by_region.csv", index=False)
print("\nSHAP share by group and region:")
print(rg.pivot(index="group", columns="region", values="share")[NREG].round(2).to_string())
KEY = ["contact3", "mother_edu", "wealth", "log2_tt_walk", "facility_delivery", "anc_visits", "birth_month", "logpopdens"]
print("\nEffect (log-odds, top vs bottom of feature) by region, key features:")
print(rr[rr.feature.isin(KEY)].pivot(index="label", columns="region", values="effect")[NREG].round(2).to_string())

fig, ax = plt.subplots(figsize=(9, 4.8))
gp = rg.pivot(index="region", columns="group", values="share").loc[NREG]
gp = gp[list(GROUP)]
gp.index = [r.title() for r in NREG]
left = np.zeros(len(gp))
cols = ["#7f8c8d", "#c0392b", "#e67e22", "#8e44ad", "#1b6ca8", "#3a9e5f"]
for c, col in zip(gp.columns, cols):
    ax.barh(gp.index, gp[c], left=left, color=col, label=c)
    for i, v in enumerate(gp[c]):
        if v > 0.06:
            ax.text(left[i] + v / 2, i, f"{100 * v:.0f}%", ha="center", va="center", color="white", fontsize=8)
    left += gp[c].values
ax.invert_yaxis(); ax.set_xlim(0, 1); ax.set_xlabel("Share of SHAP attribution (mean |SHAP|)")
ax.set_title("What predicts a slipped routine contact, by region", fontsize=11)
ax.legend(ncol=3, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.15), frameon=False)
plt.tight_layout(); plt.savefig("out/slip_shap_by_region.png", dpi=170); plt.close(fig)

# Figure: top 12 features, northern contacts, beeswarm.
top = imp.feature.head(12).tolist()
idx = [FEAT.index(f) for f in top]
plt.figure()
shap.summary_plot(sv[north][:, idx], X[top][north].astype(float), feature_names=[LABEL[f] for f in top],
                  show=False, max_display=12, plot_size=(9, 6))
plt.title("What predicts a slipped routine contact? (northern contacts, SHAP)", fontsize=11)
plt.xlabel("SHAP value (change in log-odds of a slip)")
plt.tight_layout()
plt.savefig("out/slip_shap.png", dpi=170)
print("done")
