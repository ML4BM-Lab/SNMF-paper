import sys
import os

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr
from skimage.metrics import structural_similarity as ssim

from scipy.special import gammaln

from scipy.stats import wilcoxon

def pvalue_to_stars(p):
    if p < 1e-4:
        return "****"
    elif p < 1e-3:
        return "***"
    elif p < 1e-2:
        return "**"
    elif p < 0.05:
        return "*"
    return "-"


def add_significance(ax, data, x_col, y_col, order, pairs, alternative="greater"):
    ymin, ymax = ax.get_ylim()
    y_range = ymax - ymin
    if not np.isfinite(y_range) or y_range <= 0:
        y_range = max(abs(data[y_col].max()), 1.0)

    gap = y_range * 0.035
    height = y_range * 0.025
    step = y_range * 0.08
    text_pad = y_range * 0.01

    current_y = ymax + gap
    max_annotation_y = current_y

    for g1, g2 in pairs:
        d1 = data[data[x_col] == g1][y_col].dropna()
        d2 = data[data[x_col] == g2][y_col].dropna()

        if d1.empty or d2.empty:
            continue

        _, p = wilcoxon(d1, d2, alternative=alternative)
        stars = pvalue_to_stars(p)

        x1 = order.index(g1)
        x2 = order.index(g2)
        bracket_top = current_y + height

        ax.plot([x1, x1, x2, x2],
                [current_y, bracket_top, bracket_top, current_y],
                lw=1.5, c="black", zorder=10)

        ax.text((x1 + x2) / 2,
                bracket_top + text_pad,
                stars,
                ha="center",
                va="bottom",
                fontsize=18,
                zorder=11)

        max_annotation_y = max(max_annotation_y, bracket_top + text_pad)
        current_y += step

    ax.set_ylim(ymin, max_annotation_y + gap)

sns.set_theme(style="whitegrid")
sns.set_context(
    "talk",
    font_scale=1.4,
    rc={
        "axes.titlesize": 25,
        "axes.labelsize": 22,
        "xtick.labelsize": 17,
        "ytick.labelsize": 17,
        "legend.fontsize": 17,
        "axes.linewidth": 1.4,
        "grid.linewidth": 0.7,
    },
)

results_path = sys.argv[1]
ground_truth_proportions = pd.read_csv(sys.argv[2], index_col=0)
ground_truth_counts = pd.read_csv(sys.argv[3], index_col=0)


rmses = {}
jsds = {}
pccs = {}
ssims = {}
fit_metrics = {}

name_dict = {
    "spot": r"$\phi = \beta$",
    "gene": r"$\phi = \alpha$",
    "full": r"$\phi = \alpha + \beta$"
}

def nb_loglikelihood(X, Mu, Phi, eps=1e-10):
    X = np.asarray(X, dtype=float)
    Mu = np.maximum(np.asarray(Mu, dtype=float), eps)
    Phi = np.maximum(np.asarray(Phi, dtype=float), eps)

    Theta = 1.0 / Phi

    ll = (
        gammaln(X + Theta)
        - gammaln(Theta)
        - gammaln(X + 1)
        + Theta * np.log(Theta / (Theta + Mu))
        + X * np.log(Mu / (Theta + Mu))
    )

    return np.sum(ll)


def residual_variance(X, Mu, Phi, eps=1e-10):
    X = np.asarray(X, dtype=float)
    Mu = np.maximum(np.asarray(Mu, dtype=float), eps)
    Phi = np.maximum(np.asarray(Phi, dtype=float), eps)

    Var = Mu + Phi * Mu**2
    residuals = (X - Mu) / np.sqrt(np.maximum(Var, eps))

    return np.nanvar(residuals)


def reconstruction_error(X, Mu):
    return np.mean((np.asarray(X) - np.asarray(Mu)) ** 2)

def number_of_parameters(mode, n_spots, n_genes, k):
    # W/H mean parameters
    mean_params = n_genes * k + k * n_spots

    if mode == "spot":
        disp_params = n_spots
    elif mode == "gene":
        disp_params = n_genes
    elif mode == "full":
        # alpha_i + beta_j with mean(beta)=0 constraint
        disp_params = n_genes + n_spots - 1
    else:
        raise ValueError(mode)

    return mean_params + disp_params

def get_S(positions, TAU=0.8):
    from scipy.spatial.distance import cdist
    from scipy.optimize import minimize

    x = positions[:, 0]
    y = positions[:, 1]

    def mean_value(gamma, x, y, tau):
        gamma = gamma[0] if np.ndim(gamma) else gamma

        # Pairwise squared Euclidean distances
        coords = np.column_stack((x, y))
        D2 = cdist(coords, coords, metric='sqeuclidean')

        # Similarity matrix
        S = np.exp(-gamma * D2)
        S[S < 1e-3] = 0

        # Row normalization
        row_sums = S.sum(axis=1, keepdims=True)
        S = S / row_sums

        # Objective function
        return (np.mean(np.diag(S)) - tau) ** 2

    # Optimize gamma
    result = minimize(
        mean_value,
        x0=[1.0],
        args=(x, y, TAU),
        method="BFGS"
    )

    gamma = result.x[0]

    # Compute final normalized similarity matrix
    coords = np.column_stack((x, y))
    D2 = cdist(coords, coords, metric='sqeuclidean')

    S = np.exp(-gamma * D2)
    S[S < 1e-3] = 0

    row_sums = S.sum(axis=1, keepdims=True)
    S = S / row_sums
    return S
    

for dirpath, subdirs, files in os.walk(results_path):
    for f in files:
        if f == "SNMF_proportions.csv":
            mode = dirpath.split("/")[-1]
            value = name_dict.get(mode)
            proportions = pd.read_csv(os.path.join(dirpath, f), index_col=0)

            proportions = proportions.loc[ground_truth_proportions.index]

            # RMSE
            rmse = np.sqrt(((proportions - ground_truth_proportions) ** 2).mean(axis=1))
            rmses[value] = rmse.values

            # JSD
            P = proportions.values
            Q = ground_truth_proportions.values

            jsd = np.array([
                jensenshannon(p, q)**2
                for p, q in zip(P, Q)
            ])

            jsds[value] = jsd

            # PCC
            pcc = np.array([
                pearsonr(p, q)[0]
                if np.std(p) != 0 and np.std(q) != 0
                else np.nan
                for p, q in zip(P, Q)
            ])
            pccs[value] = pcc

            # SSIM
            coords = np.array([
                list(map(int, idx.split("x")))
                for idx in proportions.index
            ])

            rows = coords[:, 0]
            cols = coords[:, 1]

            H = rows.max() + 1
            W = cols.max() + 1

            vals = []

            for celltype in ground_truth_proportions.columns:

                pred_img = np.zeros((H, W), dtype=float)
                true_img = np.zeros((H, W), dtype=float)

                pred_vals = proportions[celltype.replace("-", ".").replace("&", ".")].values
                true_vals = ground_truth_proportions.loc[[f"{r}x{c}" for r,c in zip(rows, cols)], celltype].values

                pred_img[rows, cols] = pred_vals
                true_img[rows, cols] = true_vals

                data_range = max(
                    pred_img.max() - pred_img.min(),
                    true_img.max() - true_img.min(),
                    1e-8
                )

                score = ssim(
                    pred_img,
                    true_img,
                    data_range=data_range
                )

                vals.append(score)

            ssims[value] = np.array(vals)

            # ---------- MODEL-FIT METRICS ----------
            if os.path.exists(os.path.join(dirpath, "signatures.csv")) and os.path.exists(os.path.join(dirpath, "phi.csv")):
                signatures = pd.read_csv(
                    os.path.join(dirpath, "signatures.csv"),
                    index_col=0
                )

                phi = pd.read_csv(
                    os.path.join(dirpath, "phi.csv"),
                    index_col=0
                )

                # Ensure matching order
                X = ground_truth_counts.T.loc[proportions.index, signatures.index].values

                # Mu = reconstructed mean counts
                positions = np.array(ground_truth_counts.columns.str.split("x").tolist()).astype(float)
                S = get_S(positions)
                Mu = (signatures.values @ proportions.values.T @ S).T

                Phi = phi.T.loc[proportions.index, signatures.index].values

                ll = nb_loglikelihood(X, Mu, Phi)
                rec_err = reconstruction_error(X, Mu)
                res_var = residual_variance(X, Mu, Phi)

                n_obs = X.size
                n_params = number_of_parameters(
                    mode=mode,
                    n_spots=X.shape[0],
                    n_genes=X.shape[1],
                    k=proportions.shape[1]
                )

                aic = 2 * n_params - 2 * ll
                bic = np.log(n_obs) * n_params - 2 * ll

                fit_metrics[value] = {
                    "reconstruction_error": rec_err,
                    "nb_loglikelihood": ll,
                    "aic": aic,
                    "bic": bic,
                    "residual_variance": res_var
                }


# Convert to DataFrame (long format)
rmse_df = pd.DataFrame(
    [(val, v) for val, values in rmses.items() for v in values],
    columns=["value", "rmse"]
)

jsd_df = pd.DataFrame(
    [(val, v) for val, values in jsds.items() for v in values],
    columns=["value", "jsd"]
)

pcc_df = pd.DataFrame(
    [(val, v) for val, values in pccs.items() for v in values],
    columns=["value", "pcc"]
)

ssim_df = pd.DataFrame(
    [(val, v) for val, values in ssims.items() for v in values],
    columns=["value", "ssim"]
)

# ---------- FIXED ORDER ----------
order = [
    name_dict.get(i) for i in [
        "spot",
        "gene",
        "full"
    ]
]


palette = dict(zip(
    order,
    sns.color_palette("Set2", n_colors=len(order))
))

# pairwise comparisons for significance, using Negative Binomial as reference
reference_value = name_dict.get("full")
if reference_value not in order:
    raise ValueError(
        f"Reference dispersion mode '{reference_value}' was not found in results: "
        + ", ".join(order)
    )

pairs = [
    (value, reference_value)
    for value in order
    if value != reference_value
]
pairs = sorted(
    pairs,
    key=lambda pair: abs(order.index(pair[0]) - order.index(pair[1]))
)

# ---------- RMSE PLOT ----------
plt.figure(figsize=(10,6))
ax = sns.violinplot(
    data=rmse_df,
    x="value",
    y="rmse",
    hue="value",
    palette=palette,
    order=order,
    cut=0
)

add_significance(
    ax,
    rmse_df,
    "value",
    "rmse",
    order,
    pairs,
    alternative="greater"
)

ax.set_xlabel("")
ax.set_ylabel("RMSE")

plt.tight_layout()

os.makedirs(os.path.join(results_path, "plots"), exist_ok=True)

plt.savefig(os.path.join(results_path, "plots", "rmse_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "rmse_comparison.pdf"), dpi=300)
plt.close()


# ---------- JSD PLOT ----------
plt.figure(figsize=(10,6))
ax = sns.violinplot(
    data=jsd_df,
    x="value",
    y="jsd",
    hue="value",
    palette=palette,
    order=order,
    cut=0
)

add_significance(
    ax,
    jsd_df,
    "value",
    "jsd",
    order,
    pairs,
    alternative="greater"
)

ax.set_xlabel("")
ax.set_ylabel("JSD")

plt.tight_layout()

plt.savefig(os.path.join(results_path, "plots", "jsd_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "jsd_comparison.pdf"), dpi=300)
plt.close()

# ---------- PCC PLOT ----------
plt.figure(figsize=(10,6))
ax = sns.violinplot(
    data=pcc_df,
    x="value",
    y="pcc",
    hue="value",
    palette=palette,
    order=order,
    cut=0
)

add_significance(
    ax,
    pcc_df,
    "value",
    "pcc",
    order,
    pairs,
    alternative="less"
)

ax.set_xlabel("")
ax.set_ylabel("PCC")

plt.tight_layout()

plt.savefig(os.path.join(results_path, "plots", "pcc_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "pcc_comparison.pdf"), dpi=300)
plt.close()

# ---------- SSIM PLOT ----------
plt.figure(figsize=(10,6))
ax = sns.boxplot(
    data=ssim_df,
    x="value",
    y="ssim",
    hue="value",
    palette=palette,
    order=order,
    width=0.65,
    linewidth=1.2,
    fliersize=2
)

add_significance(
    ax,
    ssim_df,
    "value",
    "ssim",
    order,
    pairs,
    alternative="less"
)

ax.set_xlabel("")
ax.set_ylabel("SSIM")

plt.tight_layout()

plt.savefig(os.path.join(results_path, "plots", "ssim_comparison.png"), dpi=300)
plt.savefig(os.path.join(results_path, "plots", "ssim_comparison.pdf"), dpi=300)
plt.close()

# ---------- MODEL-FIT SUMMARY PLOTS ----------
if len(fit_metrics) > 0:

    fit_df = (
        pd.DataFrame(fit_metrics)
        .T
        .reset_index()
        .rename(columns={"index": "value"})
    )

    fit_long = fit_df.melt(
        id_vars="value",
        var_name="metric",
        value_name="score"
    )

    metric_labels = {
        "reconstruction_error": "Reconstruction error",
        "nb_loglikelihood": "NB log-likelihood",
        "aic": "AIC",
        "bic": "BIC",
        "residual_variance": "Residual variance"
    }

    for metric, ylabel in metric_labels.items():

        df = fit_long[fit_long["metric"] == metric].copy()

        plt.figure(figsize=(10, 6))
        ax = sns.barplot(
            data=df,
            x="value",
            y="score",
            hue="value",
            palette=palette,
            order=order
        )

        ax.set_xlabel("")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)

        plt.tight_layout()

        filename = metric + "_comparison"

        plt.savefig(
            os.path.join(results_path, "plots", f"{filename}.png"),
            dpi=300
        )
        plt.savefig(
            os.path.join(results_path, "plots", f"{filename}.pdf"),
            dpi=300
        )
        plt.close()

    fit_df.to_csv(
        os.path.join(results_path, "plots", "model_fit_metrics.csv"),
        index=False
    )
