import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('default')
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr, gmean
library_df = pd.read_csv('data/rubisco_library_uncorrected.csv')
library_df.set_index(library_df.apply(lambda row: row['internal_id'] if row['from_lab'] == 1 else row['identifier'], axis=1), inplace=True)


# Filter variatns measured in both our pipeline and the literature
regression_data = library_df.query("from_lab == 1 & from_literature == 1").copy()

# Apply log10 scaling to rates
regression_data['log_rate_lab'] = np.log10(regression_data['rate_lab'])
regression_data['log_rate_lit'] = np.log10(regression_data['rate_lit_30deg'])

X = regression_data['log_rate_lab'].values
y = regression_data['log_rate_lit'].values

# Fit intercept with slope = 1
b = np.mean(y - X)  # since slope=1, intercept = mean difference
regression_line = X + b

r_squared = 1 - np.sum((y - regression_line)**2) / np.sum((y - np.mean(y))**2)
pearson_coef, _ = pearsonr(X, y)

print(f"log(y) = log(x) + {b:.2f}")
print(f"10^b ≈ {10**b:.2f} multiplicative bias")
print(f"R²: {r_squared:.2f}")
print(f"Pearson correlation coeff: {pearson_coef:.2f}")

# Bootstrap CI for b and 10^b
boot_b = []

for _ in range(10000):
    idx = np.random.choice(len(X), size=len(X), replace=True)
    boot_b.append(np.mean(y[idx] - X[idx]))

boot_b = np.array(boot_b)
ci_b = np.percentile(boot_b, [2.5, 97.5])
ci_factor = 10**ci_b

print(f"95% CI for b: [{ci_b[0]:.2f}, {ci_b[1]:.2f}]")
print(f"95% CI for multiplicative factor: [{ci_factor[0]:.2f}, {ci_factor[1]:.2f}]")

plt.figure(figsize=(2.7, 2.7))
plt.xscale("log")
plt.yscale("log")

plt.errorbar(
    10**X, 10**y,
    xerr=regression_data['std_rate_lab'].values,
    fmt='o', color='black', ecolor='gray', elinewidth=1, capsize=2
)
plt.plot([0, 12], [0, 12 * 10**b], color='black')

plt.xlabel('measured rate [s$^{-1}$]', fontsize=11)
plt.ylabel('published rate [s$^{-1}$]', fontsize=11)
plt.xlim(0.467, 10**X.max() * 1.2)
plt.ylim(1, 10**y.max() * 1.2)

plt.savefig('outputs/rates_lit_our_lab.png', dpi=300, transparent=True, bbox_inches='tight')
plt.show()

library_df.insert(library_df.columns.get_loc('rate_lab') + 3,
                  'rate_lab_corrected',
                  library_df['rate_lab'] * 10**b)

library_df.insert(library_df.columns.get_loc('std_rate_lab') + 3,
                  'std_rate_lab_corrected',
                  library_df['std_rate_lab'] * 10**b)


#create a column that averages corrected ratees from our lab and rates from the lit
library_df['avg_rate_lit_lab'] = library_df[['rate_lab_corrected', 'rate_lit_30deg']].mean(axis=1)

library_df.to_csv('data/rubisco_library.csv', index=False)

Supplementary_Data_1_df = library_df[library_df['from_literature'] == 1]

Supplementary_Data_1_df = Supplementary_Data_1_df.drop(columns=[
    "std_rate_lab", "N_rate_lab", "std_rate_lab_corrected",
    "internal_id", "rbcX", "90pct_clustering", "from_lab",
    "rate_lab", "rate_lab_corrected", "active<0.5s-1", "from_literature",
    "Tdeg_lit", "rate_lit", "avg_rate_lit_lab"
])

Supplementary_Data_1_df.to_csv('outputs/Supplementary_Data_1.csv', index=False)

Supplementary_Data_3_df = library_df[library_df['from_lab'] == 1]

Supplementary_Data_3_df = Supplementary_Data_3_df.rename(columns={"N_rate_lab": "N", "std_rate_lab_corrected": "stdev"})


# Indicate which measured variants are new to this study
Supplementary_Data_3_df["new_to_this_study"] = (
    Supplementary_Data_3_df["form"].isin(["Form 1'", "Form 3"]) |
    Supplementary_Data_3_df["internal_id"].str.startswith("RBC5_") |
    Supplementary_Data_3_df["internal_id"].isin(["RBC_3", "RBC2_44", "RBC2_45", "RBC_63", "RBC2_41", "RBC_67", "RBC_89", "RBC_47", "RBC2_40", "RBCSeed_9", "RBC_51", "RBC_48", "RBC2_39", "RBC2_36", "RPE", "ANC_1", "ANC_3"])
).astype(int)

Supplementary_Data_3_df = Supplementary_Data_3_df.drop(columns=[
    "rbcX", "90pct_clustering", "from_lab", "from_literature", "rate_lab",
    "Tdeg_lit", "rate_lit", "rate_lit_30deg", "reference", "avg_rate_lit_lab", "std_rate_lab"
])

#Reorder
Supplementary_Data_3_df = Supplementary_Data_3_df[["identifier", "internal_id", "species", "prok_euk", "form", "large_subunit", "small_subunit", "new_to_this_study", "rate_lab_corrected", "stdev", "N", "active<0.5s-1"]]

Supplementary_Data_3_df.to_csv('outputs/Supplementary_Data_3.csv', index=False)

