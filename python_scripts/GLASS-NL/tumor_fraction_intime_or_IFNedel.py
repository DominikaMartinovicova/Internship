import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import shapiro, levene, ttest_ind, mannwhitneyu
import seaborn as sns

plt.rcParams.update({
    'font.size': 14,            # Base font size
    'axes.titlesize': 16,       # Title font
    'axes.labelsize': 14,       # X/Y label font
    'xtick.labelsize': 12,      # Tick label font
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
})


# -------------------------------------------------
# 1. Load your full CSV (adjust filename as needed)
# -------------------------------------------------
df = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_DNA-seq/GLASS_tumorF_for_Statescope_filtered_curated_fitsDominika.csv')

# -------------------------------------------------
# 2. Keep sample + Malignant, parse IDs
# -------------------------------------------------
df = df.iloc[:, [0, df.columns.get_loc('Malignant')]]
df.columns = ['Sample', 'Malignant']

df['Patient']   = df['Sample'].str.extract(r'^(\d{3})')[0]
df['Resection'] = df['Sample'].str.extract(r'_R(\d+)')[0].astype(int)

# -------------------------------------------------
# 3. Sort for correct line order
# -------------------------------------------------
df = df.sort_values(['Patient', 'Resection'])

# -------------------------------------------------
# 4. Plot segment‑wise, colouring each link
# -------------------------------------------------
plt.figure(figsize=(8,10))

for patient, g in df.groupby('Patient'):
    g = g.reset_index(drop=True)
    
    # loop over consecutive pairs
    for i in range(1, len(g)):
        x = g.loc[i-1:i, 'Resection']
        y = g.loc[i-1:i, 'Malignant']
        colour = 'blue' if y.iloc[1] > y.iloc[0] else 'red'
        
        # draw the segment + second point (first point already drawn by previous segment)
        plt.plot(x, y, marker='o', linewidth=2, color=colour)
    
    # draw the first point for this patient if it didn’t get a segment yet
    if len(g) == 1:
        plt.plot(g['Resection'], g['Malignant'], marker='o', color='grey')

# -------------------------------------------------
# 5. Cosmetics
# -------------------------------------------------
plt.xlabel('Resection number')
plt.ylabel('Malignant fraction')
plt.title('Tumor fraction change between consecutive resections\n(blue = rise, red = fall/no change)')
plt.xticks(sorted(df['Resection'].unique()))
plt.tight_layout()
plt.show()

#plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/tumor_fractions/tumor_fractions_intime_updown.png')

group_df = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_DNA-seq/GLASS_comparison_groups.csv')  # Adjust file path

# Assuming group_df has columns 'Sample' and 'Group'
# Merge group information into the main dataframe
df = df.merge(group_df[['Sample', 'group']], on='Sample', how='left')
print(df)

# -------------------------------------------------
# 4. Plot the ∆ values with group information
# -------------------------------------------------
color_map = {1: 'dodgerblue', 2: 'gold', 3: 'darkorange', 4: 'red'}

plt.figure(figsize=(8,10))

# Loop through each patient group
for patient, g in df.groupby('Patient'):
    g = g.reset_index(drop=True)
    
    # Loop over consecutive pairs
    for i in range(1, len(g)):
        x = g.loc[i-1:i, 'Resection']
        y = g.loc[i-1:i, 'Malignant']
        group = g.loc[i, 'group']  # Get the group for the second sample

        colour = color_map.get(group, 'grey')  # Default to grey if group is not in map
        
        # Draw the segment + second point (first point already drawn by previous segment)
        plt.plot(x, y, marker='o', color=colour, linewidth=2)
    
    # Draw the first point for this patient if it didn’t get a segment yet
    if len(g) == 1:
        group = g['group'].iloc[0]  # Get group for the only sample
        colour = color_map.get(group, 'grey')
        plt.plot(g['Resection'], g['Malignant'], marker='o', color=colour)

handles = [plt.Line2D([0], [0], color=color_map[group], lw=4) for group in color_map]
labels = [f'Group {group}' for group in color_map]
plt.legend(handles=handles, labels=labels, title='Group')


plt.xlabel('Resection number')
plt.ylabel('Malignant fraction')
plt.title('Tumor fraction change between consecutive resections\n(Colored by group)')
plt.xticks(sorted(df['Resection'].unique()))
plt.tight_layout()
plt.show()

# Save the plot
#plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/tumor_fractions/tumor_fractions_intime_grouped.png')

# -------------------------------------------------
# 4. Plot the ∆ values with group information (points only)
# -------------------------------------------------
color_map = {1: 'dodgerblue', 2: 'gold', 3: 'darkorange', 4: 'red'}

plt.figure(figsize=(8,10))

# Plot individual points, colored by group
for group, g in df.groupby('group'):
    color = color_map.get(group, 'grey')
    plt.scatter(g['Resection'], g['Malignant'], color=color, label=f'Group {group}', s=60)

# Remove duplicate labels in legend
handles = [plt.Line2D([0], [0], marker='o', color='w', label=f'Group {g}', 
                      markerfacecolor=color_map.get(g, 'grey'), markersize=10) for g in color_map]
plt.legend(handles=handles, title='Group')

plt.xlabel('Resection number')
plt.ylabel('Malignant fraction')
plt.title('Tumor fraction per resection (Colored by group)')
plt.xticks(sorted(df['Resection'].unique()))
plt.tight_layout()
plt.show()

# Save the plot
#plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/tumor_fractions/tumor_fractions_intime_grouped_notconnected.png')

# -------------------------------------------------
# Statistical testing for difference in tumor purity between group 3 and 4
# -------------------------------------------------

# Extract malignant fractions for groups 3 and 4
g3 = df[df['group'] == 3]['Malignant'].dropna()
g4 = df[df['group'] == 4]['Malignant'].dropna()

# Shapiro-Wilk normality test
norm_g3 = shapiro(g3)[1] > 0.05 if len(g3) >= 3 else False
norm_g4 = shapiro(g4)[1] > 0.05 if len(g4) >= 3 else False
print(f"Shapiro-Wilk p-values → Group 3: {shapiro(g3)[1]:.4f}, Group 4: {shapiro(g4)[1]:.4f}")

# Levene's test for equal variances
var_equal = levene(g3, g4)[1] > 0.05 if len(g3) > 1 and len(g4) > 1 else False
print(f"Levene’s test p = {levene(g3, g4)[1]:.4f}" if var_equal is not None else "Not enough data for Levene’s test.")

# Choose statistical test
if norm_g3 and norm_g4 and var_equal:
    stat, p = ttest_ind(g3, g4, equal_var=True)
    test_used = "Independent t-test"
else:
    stat, p = mannwhitneyu(g3, g4, alternative='two-sided')
    test_used = "Mann–Whitney U test"

print(f"\n{test_used} result:")
print(f"  Statistic = {stat:.3f}")
print(f"  p-value = {p:.4e}")


# -------------------------------------------------
# Boxplot for tumor fractions for group 3 and 4 with resection
# -------------------------------------------------
df_filtered = df[df['group'].isin([3, 4])].copy()
print(df_filtered)
df_filtered['group'] = df_filtered['group'].astype(int)

color_map = {3: '#8da0cb', 4: '#e78ac3'}
palette = {str(k): v for k, v in color_map.items()}

# Convert group to string for seaborn's hue
df_filtered['group_str'] = df_filtered['group'].astype(str)

# Plot boxplot
plt.figure(figsize=(8, 10))
sns.boxplot(data=df_filtered, x='Resection', y='Malignant', hue='group_str', palette=palette)

# Stripplot for individual data points
sns.stripplot(data=df_filtered, x='Resection', y='Malignant', hue='group_str', 
              palette=palette, dodge=True, jitter=True, marker='o', alpha=0.6, edgecolor='gray', linewidth=0.5)

# Fix double legend issue from box + strip
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles[:2], ['Group 3', 'Group 4'], title='Group')

# Labels and title
plt.xlabel('Resection number')
plt.ylabel('Malignant fraction')
plt.title('Tumor fraction per resection (Groups 3 & 4 only)')
plt.tight_layout()

# Save and show
#plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/tumor_fractions/groups34_tumor_fractions_boxplot.png')
plt.show()



# -------------------------------------------------
# Boxplot for tumor fractions for group 3 and 4 without resection
# -------------------------------------------------

# Filter to groups 3 and 4 only
df_filtered = df[df['group'].isin([3, 4])].copy()
df_filtered['group'] = df_filtered['group'].astype(int)
df_filtered['group_str'] = df_filtered['group'].astype(str)

# Define color map for just groups 3 and 4
color_map = {3:  '#1f77b4', 4: '#ff7f0e'}
palette = {str(k): v for k, v in color_map.items()}

# Plot
plt.figure(figsize=(6, 6))

# Boxplot by group
sns.boxplot(data=df_filtered, x='group_str', y='Malignant', palette=palette)

# Overlay stripplot for individual points
sns.stripplot(data=df_filtered, x='group_str', y='Malignant', palette=palette, 
              jitter=True, marker='o', alpha=0.6, edgecolor='gray', linewidth=0.5)

# Labels and title
plt.xlabel('Group')
plt.ylabel('Tumor purity')
plt.title('Tumor purity distribution in Groups 3 & 4')
plt.tight_layout()

# Save and show
plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/tumor_fractions/groups34_tumor_fractions_boxplot_woresection.png')
plt.show()

# -------------------------------------------------
# Scatterplot of tumor fractions for group 3 and 4
# -------------------------------------------------
df_filtered = df[df['group'].isin([3, 4])].copy()

palette = {3:  '#1f77b4', 4: '#ff7f0e'}


# Plot
plt.figure(figsize=(max(12, len(df_filtered) * 0.4), 6))
sns.scatterplot(data=df_filtered, x='Sample', y='Malignant', hue='group',palette=palette, s=100)
plt.title('Tumor Purity in Groups 3 and 4')
plt.xlabel('Sample')
plt.ylabel('Tumor Purity')
plt.grid(True)
plt.xticks(rotation=45)
plt.subplots_adjust(bottom=0.3)


plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/tumor_fractions/groups34_tumor_fractions_scatter_woresection.png')





