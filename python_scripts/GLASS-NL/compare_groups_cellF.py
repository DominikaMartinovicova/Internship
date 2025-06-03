#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# compare_groups_cellF.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Compare the deconvolved cell fractions among the specified groups based on
# CNA status of IFNE and CDKN2AB
#
#-------------------------------------------------------------------------------
# 0 Import packages and read the data
#-------------------------------------------------------------------------------
# Load packages
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kruskal, shapiro, levene, boxcox, ttest_ind, mannwhitneyu
import scikit_posthocs as sp

# Specify plotting parameters
plt.rcParams.update({
    'font.size': 14,            # Base font size
    'axes.titlesize': 16,       # Title font
    'axes.labelsize': 14,       # X/Y label font
    'xtick.labelsize': 12,      # Tick label font
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
})

#-------------------------------------------------------------------------
# 1 Boxplot of cell fractions in different groups
#-------------------------------------------------------------------------
# Load data
groups = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova/GLASS-NL/GLASS-NL_DNA-seq/GLASS_comparison_groups.csv', index_col=0)
fractions = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova/GLASS-NL/Statescope/Estimated_Fractions_GLASSNL2_fitsDominika.csv', index_col=0)
samples_outside_01_09 = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova/GLASS-NL/GLASS-NL_DNA-seq/GLASS_tumorF_for_Statescope_filtered_outside01and09.csv', index_col=0).index.tolist()     # To possibly highlight the location of these samples in the boxplots

# Add group information by merging on index on index
df = fractions.join(groups)
print(df)

# Exclude the "Malignant" column to get the cell fraction normalized to non-tumor component (Possibly exclude more cell types to get purely stromal/lymphoid/myeloid component)
df_no_malignant = df.drop(columns=['Malignant']) #, 'Endothelial', 'Fibroblast', 'Macrophage', 'Microglia', 'Monocyte', 'Oligodendrocyte', 'Pericyte', 'cDC', 'pDC'])

# Normalize remaining cell type fractions row-wise (i.e., for each sample)
cell_type_cols = df_no_malignant.columns.difference(['group'])  # ensure group column is excluded
df_no_malignant[cell_type_cols] = df_no_malignant[cell_type_cols].div(df_no_malignant[cell_type_cols].sum(axis=1), axis=0)

# Melt for plotting (long-form)
df_melted = df_no_malignant.reset_index().melt(id_vars=['index', 'group'], 
                                                value_vars=cell_type_cols,
                                                var_name='Cell Type',
                                                value_name='Fraction')

print(df_melted)

# Boxplot
plt.figure(figsize=(12, 6))
sns.boxplot(data=df_melted, x='Cell Type', y='Fraction', hue='group', palette='Set2')

# #  Optionally highlight specific samples
# highlight_df = df_melted[df_melted['index'].isin(samples_outside_01_09)]
# # Generate a color palette and map it to samples
# palette = sns.color_palette("husl", len(samples_outside_01_09))
# sample_color_map = dict(zip(samples_outside_01_09, palette))

# # Get consistent order for Cell Types and Groups
# cell_type_order = df_melted['Cell Type'].unique()
# group_order = df_melted['group'].unique()

# # Map Cell Types and Groups to x-axis positions
# n_groups = len(group_order)
# group_offsets = {g: i * 0.2 - 0.3 for i, g in enumerate(group_order)}  # dodge-like shift

# # Overlay each sample point at the correct position
# for sample in samples_outside_01_09:
#     sample_data = highlight_df[highlight_df['index'] == sample]
#     sample_color = sample_color_map[sample]
    
#     for _, row in sample_data.iterrows():
#         x_base = list(cell_type_order).index(row['Cell Type'])
#         offset = group_offsets[row['group']]
#         plt.scatter(
#             x=x_base + offset,
#             y=row['Fraction'],
#             color=sample_color,
#             s=30,
#             marker='x',
#             edgecolor='',
#             label=sample,
#             zorder=10
#         )

# # De-duplicate legend entries
# handles, labels = plt.gca().get_legend_handles_labels()
# unique = dict(zip(labels, handles))

# plt.legend(unique.values(), unique.keys(), bbox_to_anchor=(1.05, 1), loc='upper left', title='Highlighted Samples')

#plt.ylim(0.0, 0.2) 
plt.xticks(rotation=45)
plt.title('Cell Type Fractions by Group')
plt.tight_layout()
plt.show()

# Save the plot
plt.savefig("/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/deconvolution/Statescope/groups_cellfraction_boxplot_GLASSNL2_fitsDominika_Statescope_asterisk.png")

# #-------------------------------------------------------------------------
# # 2 Statistical testing Kruskal-Wallis with Dunn's test for pairwise comparison 
# # and multiple testing correction for all the groups
# #-------------------------------------------------------------------------
# # Store results
# normality_results = {}
# homogeneity_results = {}

## Check if data meets the parametric asssumptions
# print("Checking ANOVA assumptions...\n")
# for cell in cell_type_cols:
#     print(f"--- {cell} ---")
    
#     # 1. Normality test (Shapiro-Wilk) per group
#     group_pvals = []
#     for g in sorted(df['group'].unique()):
#         data = df[df['group'] == g][cell].dropna()
#         if len(data) >= 3:  # Shapiro-Wilk needs at least 3 values
#             stat, p = shapiro(data)
#             group_pvals.append(p)
#             print(f"  Group {g}: Shapiro-Wilk p = {p:.4f}")
#         else:
#             print(f"  Group {g}: Not enough data for normality test.")
    
#     # Flag if any group is not normal
#     normality_results[cell] = all(p > 0.05 for p in group_pvals)
    
#     # 2. Homogeneity of variances (Levene’s test)
#     group_data = [df[df['group'] == g][cell].dropna() for g in sorted(df['group'].unique())]
#     if all(len(gd) > 1 for gd in group_data):  # Levene requires >1 value per group
#         stat, p = levene(*group_data)
#         homogeneity_results[cell] = p > 0.05
#         print(f"  Levene’s test p = {p:.4f}")
#     else:
#         homogeneity_results[cell] = False
#         print("  Not enough data for Levene’s test.")
    
#     print(f"  → Eligible for ANOVA: {normality_results[cell] and homogeneity_results[cell]}\n")

# # Dunn's test post-hoc pairwise comparisons (only for significant Kruskal-Wallis)
# print("\nDunn’s post-hoc test (for significant Kruskal-Wallis):")
# for cell in cell_type_cols:
#     groups_data = [df[df['group'] == g][cell].dropna() for g in sorted(df['group'].unique())]
#     stat, p = kruskal(*groups_data)
#     print(f"{cell}: H = {stat:.3f}, p = {p:.4e}")
#     if p < 0.05:
#         print(f"\n{cell}: Kruskal-Wallis p = {p:.4e} → Running Dunn’s test")
#         dunn_result = sp.posthoc_dunn(df, val_col=cell, group_col='group', p_adjust='bonferroni')
#         print(dunn_result)
#     else:
#         print(f"{cell}: Kruskal-Wallis p = {p:.4e} → Not significant, skipping Dunn’s test")

# #-------------------------------------------------------------------------
# # 3 Statistical testing only between group 3 and 4. Depending on meeting 
# # the parametric assumptions use either Student's t-test or Mann-Whitney U-test.
# #-------------------------------------------------------------------------
print("Statistical testing: Group 3 vs Group 4\n")
# Store results
test_results = {}

# Loop over all the cell types
for cell in cell_type_cols:
    print(f"--- {cell} ---")
    
    # Filter data for groups 3 and 4
    g3 = df[df['group'] == 3][cell].dropna()
    g4 = df[df['group'] == 4][cell].dropna()
    
    # Normality test (Shapiro-Wilk)
    norm_g3 = shapiro(g3)[1] > 0.05 if len(g3) >= 3 else False
    norm_g4 = shapiro(g4)[1] > 0.05 if len(g4) >= 3 else False
    print(f"  Shapiro p-values → Group 3: {shapiro(g3)[1]:.4f}, Group 4: {shapiro(g4)[1]:.4f}")
    
    # Homogeneity of variances (Levene’s test)
    var_equal = False
    if len(g3) > 1 and len(g4) > 1:
        var_equal = levene(g3, g4)[1] > 0.05
        print(f"  Levene’s test p = {levene(g3, g4)[1]:.4f}")
    else:
        print("  Not enough data for Levene’s test.")
    
    # Decide test type
    if norm_g3 and norm_g4 and var_equal:
        # Use independent t-test
        stat, p = ttest_ind(g3, g4, equal_var=True)
        test_used = 't-test'
    else:
        # Use Mann–Whitney U test
        stat, p = mannwhitneyu(g3, g4, alternative='two-sided')
        test_used = 'Mann-Whitney U'
    
    test_results[cell] = {'test': test_used, 'stat': stat, 'p': p}
    print(f"  → {test_used}: stat = {stat:.3f}, p = {p:.4e}\n")


# Filter to only Groups 3 and 4
df_melted_filtered = df_melted[df_melted['group'].isin([3, 4])]

#-------------------------------------------------------------------------
# 4 Boxplot of cell fractions in groups 3 and 4
#-------------------------------------------------------------------------
# # Plot
plt.figure(figsize=(10, 6))
ax = sns.boxplot(
    data=df_melted_filtered,
    x='Cell Type',
    y='Fraction',
    hue='group',
    palette={3:  '#1f77b4', 4: '#ff7f0e'}
)

# Optionally add asterisk for p < 0.1 -------------------------------
for i, cell in enumerate(df_melted_filtered['Cell Type'].unique()):
    res = test_results.get(cell, {})
    pval = res.get('p', 1.0)

    if pval < 0.1:
        # Get max y for this cell type to place line and text
        y_vals = df_melted_filtered[df_melted_filtered['Cell Type'] == cell]['Fraction']
        y_max = y_vals.max()
        y_line = y_max + 0.02
        y_text = y_line + 0.0005

        # Draw horizontal line
        ax.plot([i - 0.2, i + 0.2], [y_line, y_line], color='black', linewidth=1)

        # Draw asterisk
        ax.text(i, y_text, '*', ha='center', va='bottom', color='black', fontsize=16)
# Add asterisk explanation below the plot
plt.figtext(0.11, 0.88, '* p < 0.1', ha='left', fontsize=12)
# -------------------------------------------------------------------

plt.xticks(rotation=45, ha='right')
plt.title('Cell Type Fractions in Group 3 and 4')
plt.tight_layout()
plt.ylim(top=0.3)

# Save
plt.savefig("/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/deconvolution/Statescope/groups34_cellfraction_boxplot_GLASSNL2_fitsDominika_Statescope_asterisk.png")