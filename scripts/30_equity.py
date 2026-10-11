"""Equity: who would lose protection if SMC moved from the campaign to routine delivery?

The door-to-door campaign is assumed to reach children equally (87% per cycle) whatever their
household; routine delivery reaches a child only through the routine system, so its coverage
follows that child's chance of first contact and of a slipped monthly contact. For each northern
DHS child with a dated card (27a_slip_features.R):
  q_i     per-contact slip probability: logistic regression on all features (as 27_slip_drivers.py,
          fitted to all Ghana contacts), averaged over doses 2 and 3
  P1_w    first-contact reach: Penta1 coverage in the child's wealth quintile (north, weighted)
  C_r,i = P1_w x (1 - k q_i), k = 2 (1.4 and 2.8 as sensitivity); C_c = 0.87
  base_i  no-SMC season incidence of the child's district: observed / (1 - PE x C_c), PE = 0.88
  lost_i  = base_i x PE x (C_c - C_r,i): expected cases averted lost per child per season
Summaries use survey weights: by wealth quintile, mother's education and (for the poorest quintile)
region, the mean routine coverage,
cases averted lost per 1,000 children, and each group's share of all lost cases; plus a
concentration index of the loss over the wealth ranking (negative = concentrated among the poor).
These follow the child-level associations, not causal effects, and cover card-holding children
aged 12-35 months, who are likely better off than those without cards.

Outputs: data/derived/equity_by_group.csv, equity_summary.csv, equity_region_ratio_ci.csv, out/equity.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

D = "data/derived/"
CC, PE = 0.87, 0.88
d = pd.read_csv(D + "slip_contacts_features.csv")
d["log2_tt_motor"] = np.log2(np.maximum(d.tt_motor_any, 5))
d["log2_tt_walk"] = np.log2(np.maximum(d.tt_walk_any, 5))
d["contact3"] = (d.contact == 3).astype(int)
for r in sorted(d.region.unique()):
    d["reg_" + r.replace(" ", "_")] = (d.region == r).astype(int)
FEAT = ["female", "birth_order", "birth_month", "twin", "mother_age", "mother_edu", "mother_working", "insured",
        "media", "problem_permission", "problem_money", "problem_distance", "problem_alone", "anc_visits",
        "facility_delivery", "wealth", "hh_size", "under5_in_hh", "head_female", "urban", "log2_tt_motor",
        "log2_tt_walk", "logpopdens", "contact3"] + [c for c in d.columns if c.startswith("reg_")]
model = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), LogisticRegression(max_iter=2000, C=0.5))
model.fit(d[FEAT], d.slip)

# One row per northern child; q averaged over its dose-2 and dose-3 contact.
kids = d[d.north].drop_duplicates("child_id").copy()
q = np.zeros(len(kids))
for c3 in (0, 1):
    X = kids[FEAT].copy(); X["contact3"] = c3
    q += model.predict_proba(X)[:, 1] / 2
kids["q"] = q
p1 = pd.read_csv(D + "penta1_by_wealth_north.csv").set_index("wealth")["penta1"]
kids["p1"] = kids.wealth.map(p1)
cl = pd.read_csv(D + "ml_cluster_covariates.csv")[["cluster", "district"]]
bu = pd.read_csv(D + "district_burden.csv")[["district", "season_inc_u5"]]
kids = kids.merge(cl, on="cluster", how="left").merge(bu, on="district", how="left")
kids["season_inc_u5"] = kids.season_inc_u5.fillna(kids.season_inc_u5.median())
kids["base"] = kids.season_inc_u5 / (1 - PE * CC)


def scenario(k):
    cr = kids.p1 * (1 - np.minimum(1, k * kids.q))
    return cr, kids.base * PE * (CC - cr)


EDU = {0: "No education", 1: "Primary", 2: "Secondary", 3: "Higher"}
rows, summ = [], []
for k in (1.4, 2.0, 2.8):
    cr, lost = scenario(k)
    w = kids.weight
    for var, lab in (("wealth", "Wealth quintile"), ("mother_edu", "Mother's education")):
        for g, sub in kids.groupby(var):
            ww = w[sub.index]
            rows.append({"k": k, "by": lab, "group": int(g) if var == "wealth" else EDU.get(int(g), str(g)),
                         "children": len(sub), "routine_cov": np.average(cr[sub.index], weights=ww),
                         "campaign_cov": CC, "lost_per_1000": 1000 * np.average(lost[sub.index], weights=ww),
                         "share_of_loss": (lost[sub.index] * ww).sum() / (lost * w).sum(),
                         "share_of_children": ww.sum() / w.sum()})
    # Concentration index of the loss over the wealth ranking (fractional rank, ties averaged).
    order = kids.assign(lost=lost, w=w).sort_values("wealth")
    cw = order.w.cumsum() / order.w.sum()
    order["rank"] = order.groupby("wealth")["w"].transform(lambda x: 0) + (cw - order.w / (2 * order.w.sum()))
    order["rank"] = order.groupby("wealth")["rank"].transform("mean")
    mu = np.average(order.lost, weights=order.w)
    ci = 2 * np.average((order.lost - mu) * (order["rank"] - np.average(order["rank"], weights=order.w)),
                        weights=order.w) / mu
    summ.append({"k": k, "routine_cov_mean": np.average(cr, weights=w), "lost_per_1000_mean": 1000 * mu,
                 "concentration_index": ci})
    # By region: does the poorest national wealth quintile bear more than its share of the loss?
    # (Districts hold too few card-holding children for a district-level estimate.)
    for reg, sub in kids.groupby("region"):
        ww = w[sub.index]; poor = sub.wealth == 1
        loss_poor = (lost[sub.index][poor] * ww[poor]).sum() / (lost[sub.index] * ww).sum()
        kids_poor = ww[poor].sum() / ww.sum()
        rows.append({"k": k, "by": "Region: poorest quintile", "group": reg, "children": len(sub),
                     "routine_cov": np.average(cr[sub.index], weights=ww), "campaign_cov": CC,
                     "lost_per_1000": 1000 * np.average(lost[sub.index], weights=ww),
                     "share_of_loss": loss_poor, "share_of_children": kids_poor})
res = pd.DataFrame(rows).round(4)

# Uncertainty in each region's disproportionality ratio (poorest quintile's share of the loss / its
# share of children), k = 2: 500 bootstrap resamples of DHS clusters, refitting the slip model each
# time. The algorithm's equity safeguard uses the lower 95% bound (31_delivery_algorithm.R).
def region_ratios(dd):
    m = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), LogisticRegression(max_iter=2000, C=0.5))
    m.fit(dd[FEAT], dd.slip)
    kk = dd[dd.north].drop_duplicates("child_id").copy()
    qq = sum(m.predict_proba(kk[FEAT].assign(contact3=c3))[:, 1] / 2 for c3 in (0, 1))
    kk["p1"] = kk.wealth.map(p1)
    kk = kk.merge(cl, on="cluster", how="left").merge(bu, on="district", how="left")
    kk["season_inc_u5"] = kk.season_inc_u5.fillna(bu.season_inc_u5.median())
    lost_ = kk.season_inc_u5 / (1 - PE * CC) * PE * (CC - kk.p1 * (1 - np.minimum(1, 2.0 * np.asarray(qq))))
    out = {}
    for reg, sub in kk.assign(lost=lost_.values).groupby("region"):
        poor = sub.wealth == 1
        out[reg] = ((sub.lost[poor] * sub.weight[poor]).sum() / (sub.lost * sub.weight).sum()) / (sub.weight[poor].sum() / sub.weight.sum())
    return out
rng = np.random.default_rng(20261011)
clusters = d.cluster.unique()
byc = {c: g for c, g in d.groupby("cluster")}
boots = []
for b in range(500):
    pick = rng.choice(clusters, len(clusters), replace=True)
    dd = pd.concat([byc[c].assign(child_id=byc[c].child_id + f"_{i}") for i, c in enumerate(pick)], ignore_index=True)
    boots.append(region_ratios(dd))
bt = pd.DataFrame(boots)
est = region_ratios(d)
eqci = pd.DataFrame({"region": list(est), "ratio": [est[r] for r in est],
                     "lo": [bt[r].quantile(0.025) for r in est], "hi": [bt[r].quantile(0.975) for r in est]}).round(3)
eqci.to_csv(D + "equity_region_ratio_ci.csv", index=False)
print("Disproportionality ratio (poorest quintile's loss share / child share), k = 2, 95% bootstrap CI:")
print(eqci.to_string(index=False))
res.to_csv(D + "equity_by_group.csv", index=False)
sm = pd.DataFrame(summ).round(4)
sm.to_csv(D + "equity_summary.csv", index=False)
print(f"northern children: {len(kids)}")
print(sm.to_string(index=False))
print(res[res.k == 2.0].drop(columns="k").to_string(index=False))
print("(Region rows: share_of_loss / share_of_children are the poorest quintile's share of the region's loss / children.)")

# Figure (k = 2): routine coverage by wealth quintile against the campaign, and share of loss.
r2 = res[(res.k == 2.0) & (res.by == "Wealth quintile")]
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
lab = ["Poorest", "Poorer", "Middle", "Richer", "Richest"]
ax[0].bar(lab, r2.routine_cov, color="#c0392b", label="Routine delivery")
ax[0].axhline(CC, color="#1b6ca8", lw=2, label="Campaign (87%)")
ax[0].set_ylim(0, 1); ax[0].set_ylabel("Per-cycle coverage"); ax[0].set_title("Coverage by household wealth", fontsize=10)
ax[0].legend(fontsize=8, loc="lower right", frameon=False)
ax[1].bar(lab, 100 * r2.share_of_loss, color="#7f8c8d")
ax[1].plot(lab, 100 * r2.share_of_children, "o--", color="black", label="share of children")
ax[1].set_ylabel("Share of cases averted lost (%)"); ax[1].set_title("Who loses protection under routine delivery", fontsize=10)
ax[1].legend(fontsize=8, frameon=False)
for a in ax:
    a.tick_params(axis="x", labelsize=8)
fig.suptitle(f"Northern children, k = 2 (concentration index of loss {sm.loc[sm.k == 2.0, 'concentration_index'].iloc[0]:.2f})",
             fontsize=10)
plt.tight_layout(); plt.savefig("out/equity.png", dpi=170)
print("done")
