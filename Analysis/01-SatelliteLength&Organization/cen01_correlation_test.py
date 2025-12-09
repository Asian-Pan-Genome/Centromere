import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, kruskal
import statsmodels.api as sm
import statsmodels.formula.api as smf


df = pd.read_csv("apg_chrom_bsat_inversion_info.xls", sep="\t")

df["inversion"] = df["inversion"].astype(int)
df["25_array"] = df["25_array"].astype(int)
df["BSat_array_copy_number"] = df["BSat_array_copy_number"].astype(int)

print("==== 1. inversion vs 25_array====")
cont = pd.crosstab(df["inversion"], df["25_array"])
print("Contingency table:\n", cont)

# 使用卡方检验
chi2, p_chi2, dof, ex = chi2_contingency(cont)
print(f"Chi-square test p-value: {p_chi2}")

# Fisher
if cont.shape == (2, 2):
    _, p_fisher = fisher_exact(cont)
    print(f"Fisher exact test p-value: {p_fisher}")

print("\n==== 2. inversion vs BSat_array_copy_number ====")
grouped = [group["BSat_array_copy_number"].values for name, group in df.groupby("inversion")]
stat, p_kruskal = kruskal(*grouped)
print(f"Kruskal-Wallis test p-value: {p_kruskal}")

print("\n==== 3. 25_array vs BSat_array_copy_number====")
grouped = [group["BSat_array_copy_number"].values for name, group in df.groupby("25_array")]
stat, p_kruskal = kruskal(*grouped)
print(f"Kruskal-Wallis test p-value: {p_kruskal}")
