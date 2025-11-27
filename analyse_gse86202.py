from pathlib import Path
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

FORCE_LOGFC_COL = None
FORCE_PVAL_COL = None
FORCE_GENE_COL = None

try:
    ROOT = Path(_file_).resolve().parents[1]
except NameError:
    ROOT = Path.cwd()

DATA = ROOT / "data" / "GSE86202_Gene_differential_expression.xlsx"
OUT = ROOT / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_excel(DATA, engine="openpyxl")
orig_columns = list(df.columns)
norm_columns = [str(c).strip() for c in orig_columns]
lower_map = {c.lower(): c for c in norm_columns}

def find_column_by_tokens(tokens):
    for lc, orig in lower_map.items():
        if all(tok in lc for tok in tokens):
            return orig
    return None

logfc_col = None
if FORCE_LOGFC_COL and FORCE_LOGFC_COL in orig_columns:
    logfc_col = FORCE_LOGFC_COL
if logfc_col is None:
    for c in [
        find_column_by_tokens(["log", "fc"]),
        find_column_by_tokens(["log2", "fold"]),
        find_column_by_tokens(["fold_change"]),
        find_column_by_tokens(["log2"])
    ]:
        if c:
            logfc_col = c
            break

pval_col = None
if FORCE_PVAL_COL and FORCE_PVAL_COL in orig_columns:
    pval_col = FORCE_PVAL_COL
if pval_col is None:
    for c in [
        find_column_by_tokens(["p_value"]),
        find_column_by_tokens(["p", "value"]),
        find_column_by_tokens(["pval"]),
        find_column_by_tokens(["p", "val"]),
        find_column_by_tokens(["p_"])
    ]:
        if c:
            pval_col = c
            break

gene_col = None
if FORCE_GENE_COL and FORCE_GENE_COL in orig_columns:
    gene_col = FORCE_GENE_COL
if gene_col is None:
    for try_name in ["gene", "gene_id", "ensembl_gene_id", "symbol"]:
        if try_name in lower_map:
            gene_col = lower_map[try_name]
            break
    if gene_col is None:
        for lc, orig in lower_map.items():
            if "gene" in lc:
                gene_col = orig
                break

if logfc_col is None or pval_col is None:
    print("Could not detect required columns.")
    sys.exit(1)

WORK_COL_LOGFC = "logfc"
WORK_COL_PVAL = "pval"
WORK_COL_GENE = "gene_id"

work = df.copy()
work[logfc_col] = pd.to_numeric(work[logfc_col], errors="coerce")
work[pval_col] = pd.to_numeric(work[pval_col], errors="coerce")
work = work.rename(columns={logfc_col: WORK_COL_LOGFC, pval_col: WORK_COL_PVAL})
if gene_col:
    work = work.rename(columns={gene_col: WORK_COL_GENE})
else:
    work = work.rename(columns={orig_columns[0]: WORK_COL_GENE})

work = work.dropna(subset=[WORK_COL_LOGFC, WORK_COL_PVAL]).copy()
clean_path = OUT / "gse86202_clean.csv"
work.to_csv(clean_path, index=False)

work["neglog10_p"] = -np.log10(work[WORK_COL_PVAL].replace(0, 1e-300))
sig_mask = (work[WORK_COL_PVAL] < 0.05) & (work[WORK_COL_LOGFC].abs() >= 1)

plt.figure(figsize=(8,6))
sns.scatterplot(data=work, x=WORK_COL_LOGFC, y="neglog10_p", s=10, alpha=0.6, color="#4C72B0")
if sig_mask.any():
    sns.scatterplot(data=work[sig_mask], x=WORK_COL_LOGFC, y="neglog10_p", s=18, alpha=0.9, color="#DD2C00")
plt.axvline(1, color="grey", linestyle="--")
plt.axvline(-1, color="grey", linestyle="--")
plt.axhline(-np.log10(0.05), color="grey", linestyle="--")
plt.title("Volcano Plot — GSE86202")
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10(p-value)")
plt.tight_layout()
plt.savefig(OUT / "volcano_gse86202.png", dpi=300)
plt.close()

top20 = work.sort_values(by=WORK_COL_PVAL).head(20).copy()
top20.to_csv(OUT / "top20_genes.csv", index=False)
top20 = top20.sort_values(by=WORK_COL_LOGFC)

plt.figure(figsize=(7,10))
plt.hlines(y=top20[WORK_COL_GENE], xmin=0, xmax=top20[WORK_COL_LOGFC], color="lightgray")
colors = ["#DD2C00" if v > 0 else "#1F77B4" for v in top20[WORK_COL_LOGFC]]
plt.scatter(top20[WORK_COL_LOGFC], top20[WORK_COL_GENE], color=colors, s=60)
plt.axvline(0, color="black")
plt.title("Top 20 Significant Genes (Lollipop Plot)")
plt.xlabel("log2 Fold Change")
plt.tight_layout()
plt.savefig(OUT / "lollipop_top20.png", dpi=300)
plt.close()

plt.figure(figsize=(11,5))
plt.subplot(1,2,1)
sns.histplot(work[WORK_COL_LOGFC], bins=40, kde=True, color="#4C72B0")
plt.title("Distribution of log2 Fold Change")
plt.subplot(1,2,2)
sns.histplot(work["neglog10_p"], bins=40, kde=True, color="#DD2C00")
plt.title("Distribution of -log10(p-value)")
plt.tight_layout()
plt.savefig(OUT / "distributions.png", dpi=300)
plt.close()

print("DONE")
print("Saved:", clean_path.name)
print("Saved: volcano_gse86202.png")
print("Saved: lollipop_top20.png")
print("Saved: distributions.png")
print("Saved: top20_genes.csv")