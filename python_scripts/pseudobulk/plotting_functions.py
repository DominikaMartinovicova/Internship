# Plotting functions to evaluate deconvolution
# 1. Boxplot for mean gene expression and std from signature created with Statescope/OncoBLADE
# 2a. Barplot for PCC betweeen true and estimated cell fractions after deconvolution per cell type
# 2b. Barplot for RMSD between true fractons and deconvolved fractions
# 3. Boxplot for RMSD between true fractons and deconvolved fractions
# 4. Scatterplot for correlation betweeen true and estimated cell fractions after deconvolution per cell type
# 5. Combined scatterplot of correlation for all celltypes
# 6. Scatterplot for malignant cell fraction (of the true bulk)

#-------------------------------------------------------------------------------
# 0 Import packages and prepare variables
#-------------------------------------------------------------------------------
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from scipy.stats import pearsonr
import numpy as np
import colorcet as cc
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

plt.rcParams.update({
    'font.size': 14,            # Base font size
    'axes.titlesize': 16,       # Title font
    'axes.labelsize': 14,       # X/Y label font
    'xtick.labelsize': 12,      # Tick label font
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
})


#------------------------------------------------------------
# 1 Boxplot for mean gene expression and std from signature created with Statescope/OncoBLADE
#------------------------------------------------------------
def check_signature_scExp_scVar(signature, scExp_out, scVar_out):
    print("Checking scExp and scVar from signature...")
    print("Signature:")
    print(signature)

    # Split the dataframe to exp or std
    df_scExp = signature.filter(regex='^scExp')
    df_scVar = signature.filter(regex='^scVar')
    print(df_scExp)
    print(df_scVar)

    # Adjust dfs for seaborn
    df_scExp_melted = df_scExp.reset_index().melt(id_vars='Gene', var_name='CellType', value_name='scExp')
    df_scVar_melted = df_scVar.reset_index().melt(id_vars='Gene', var_name='CellType', value_name='scVar')

    # Plot a boxplot for each cell type
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='CellType', y='scExp', data=df_scExp_melted)
    # Add titles and labels
    plt.title('scExp Distribution for Each Cell Type')
    plt.ylabel('scExp')
    plt.xlabel('Cell Type')
    plt.xticks(rotation=45, ha='right')
    # Save the plot
    plt.savefig(scExp_out)
    print("Plot saved in ", scExp_out)

    # Plot a boxplot for each cell type
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='CellType', y='scVar', data=df_scVar_melted)
    # Add titles and labels
    plt.title('scVar Distribution for Each Cell Type')
    plt.ylabel('scVar')
    plt.xlabel('Cell Type')
    plt.xticks(rotation=45)
    # Save the plot
    plt.savefig(scVar_out)
    print("Plot saved in ", scVar_out)

#------------------------------------------------------------
# 2a Barplot for PCC betweeen true and estimated cell fractions after deconvolution per cell type
#------------------------------------------------------------
def correlation_Tf_Df(true_fractions_file, deconvolved_fractions_file, output_corr, custom_colors):
    true_fractions = pd.read_csv(true_fractions_file, index_col = 0)
    true_fractions.index.name = None
    deconvolved_fractions = pd.read_csv(deconvolved_fractions_file, index_col=0).sort_index()

    print("True fractions:")
    print(true_fractions.head())
    print(true_fractions.shape)
    print("Deconvolved fractions:")
    print(deconvolved_fractions.head())
    print(deconvolved_fractions.shape)

    assert true_fractions.shape == deconvolved_fractions.shape, "Shapes of true and predicted fractions must match!"


    # Sort rows numerically based on the number in the sample name (e.g., Sample 1, Sample 2, ...)
    # deconvolved_fractions.index = deconvolved_fractions.index.map(lambda x: int(x.split()[-1]))
    # deconvolved_fractions = deconvolved_fractions.sort_index()

    # # Sort both rows (index) and columns
    # true_fractions = true_fractions.sort_index(axis=1)
    # deconvolved_fractions = deconvolved_fractions.sort_index(axis=1)

    # # Ensure that the row names in deconvolved_fractions match the row names in true_fractions
    # deconvolved_fractions = deconvolved_fractions.set_index(true_fractions.index)

    # print("True fractions:")
    # print(true_fractions.head())
    # print(true_fractions.shape)
    # print("Deconvolved fractions:")
    # print(deconvolved_fractions.head())
    # print(deconvolved_fractions.shape)
    # Assuming the columns of both dataframes correspond to the same cell types
    correlations = {}

    # Loop through each cell type (each column)
    for column in true_fractions.columns:
        true_values = true_fractions[column]
        deconvolved_values = deconvolved_fractions[column]
        
        # Calculate Pearson correlation
        corr, _ = pearsonr(true_values, deconvolved_values)
        correlations[column] = corr

    # Print the results
    for cell_type, correlation in correlations.items():
        print(f"Pearson correlation for {cell_type}: {correlation:.3f}")
    
    # Calculate and print average Pearson correlation
    avg_corr = sum(correlations.values()) / len(correlations)
    print(f"\nAverage Pearson Correlation across all cell types: {avg_corr:.3f}")

    # colors = sns.color_palette("Paired", len(correlations))

    # Prepare bar colors using custom colormap
    bar_colors = [custom_colors.get(cell_type, '#333333') for cell_type in correlations.keys()]

    # Plot the correlations
    plt.figure(figsize=(10, 6))
    plt.bar(correlations.keys(), correlations.values(), color=bar_colors)
    plt.xlabel('Cell Types')
    plt.ylabel('Pearson Correlation')
    plt.title('Pearson Correlation between True and Deconvolved Cell Fractions')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_corr)
    print("Plot saved in ", output_corr)


#------------------------------------------------------------
# 2b Barplot for RMSD between true fractons and deconvolved fractions
#------------------------------------------------------------
def RMSD_Tf_Df(true_fractions_file, deconvolved_fractions_file, output_path, custom_colors):
    # Assuming 'true_fractions' and 'predicted_fractions' are pandas DataFrames
    true_fractions = pd.read_csv(true_fractions_file, index_col = 0).rename_axis(None)
    deconvolved_fractions = pd.read_csv(deconvolved_fractions_file, index_col = 0)

    # # Sort rows numerically based on the number in the sample name (e.g., Sample 1, Sample 2, ...)
    # deconvolved_fractions.index = deconvolved_fractions.index.map(lambda x: int(x.split()[-1]))
    # deconvolved_fractions = deconvolved_fractions.sort_index()

    # # Sort both rows (index) and columns
    # true_fractions = true_fractions.sort_index(axis=1)
    # deconvolved_fractions = deconvolved_fractions.sort_index(axis=1)

    # # Ensure that the row names in deconvolved_fractions match the row names in true_fractions
    # deconvolved_fractions = deconvolved_fractions.set_index(true_fractions.index)

    # Ensure the true and predicted fractions are aligned
    print("True fractions:")
    print(true_fractions.head())
    print(true_fractions.shape)
    print("Deconvolved fractions:")
    print(deconvolved_fractions.head())
    print(deconvolved_fractions.shape)

    assert true_fractions.shape == deconvolved_fractions.shape, "Shapes of true and predicted fractions must match!"

    # Calculate RMSD for each cell type
    rmsd_values = {}
    for cell_type in true_fractions.columns:
        print(f'Calculating RMSD for {cell_type}...')
        rmsd_values[cell_type] = np.sqrt(((true_fractions[cell_type] - deconvolved_fractions[cell_type]) ** 2).mean())
    
    # Calculate and print average Pearson correlation
    avg_rmsd = sum(rmsd_values.values()) / len(rmsd_values)
    print(f"\nAverage RMSD across all cell types: {avg_rmsd:.3f}")

    # Convert the RMSD values to a pandas Series for easier plotting
    rmsd_series = pd.Series(rmsd_values)
    print(rmsd_series)

    # colors = sns.color_palette("Paired", len(true_fractions.columns))

    bar_colors = [custom_colors.get(cell_type, '#333333') for cell_type in rmsd_values.keys()]
    # Plotting the RMSD for each cell type
    plt.figure(figsize=(10, 6))
    rmsd_series.plot(kind='bar', color=bar_colors)
    plt.title('RMSD for Each Cell Type')
    plt.xlabel('Cell Type')
    plt.ylabel('RMSD')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()
    plt.savefig(output_path)
    print("Plot saved in ", output_path)

#------------------------------------------------------------
# 3 Boxplot for RMSD between true fractons and deconvolved fractions
#------------------------------------------------------------
def box_RMSD_Tf_Df(true_fractions_file, deconvolved_fractions_file, output_path, custom_colors):
    # Assuming 'true_fractions' and 'predicted_fractions' are pandas DataFrames
    true_fractions = pd.read_csv(true_fractions_file, index_col=0).rename_axis(None)
    deconvolved_fractions = pd.read_csv(deconvolved_fractions_file, index_col=0)

    # # Sort rows numerically based on the number in the sample name (e.g., Sample 1, Sample 2, ...)
    # deconvolved_fractions.index = deconvolved_fractions.index.map(lambda x: int(x.split()[-1]))
    # deconvolved_fractions = deconvolved_fractions.sort_index()

    # # Sort both rows (index) and columns
    # true_fractions = true_fractions.sort_index(axis=1)
    # deconvolved_fractions = deconvolved_fractions.sort_index(axis=1)

    # # Ensure that the row names in deconvolved_fractions match the row names in true_fractions
    # deconvolved_fractions = deconvolved_fractions.set_index(true_fractions.index)

    # Ensure the true and predicted fractions are aligned
    print("True fractions:")
    print(true_fractions.head())
    print(true_fractions.shape)
    print("Deconvolved fractions:")
    print(deconvolved_fractions.head())
    print(deconvolved_fractions.shape)

    assert true_fractions.shape == deconvolved_fractions.shape, "Shapes of true and predicted fractions must match!"

    # Calculate RMSD for each sample and store it in a new DataFrame for plotting
    rmsd_values = pd.DataFrame(index=true_fractions.index)
    for cell_type in true_fractions.columns:
        print(f'Calculating RMSD for {cell_type}...')
        rmsd_values[cell_type] = np.sqrt(((true_fractions[cell_type] - deconvolved_fractions[cell_type]) ** 2))

    # Create a DataFrame for plotting that contains RMSD values for all samples and cell types
    rmsd_melted = rmsd_values.melt(var_name='Cell Type', value_name='RMSD')

    # Set up a boxplot for RMSD values for each cell type
    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Cell Type', y='RMSD', data=rmsd_melted, palette=custom_colors, width=0.6)

    # Set plot title and labels
    plt.title('Boxplot of RMSD for Each Cell Type')
    plt.xlabel('Cell Type')
    plt.ylabel('RMSD')

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    # Save the plot as an image
    plt.savefig(output_path)
    print(f"Plot saved in {output_path}")

    # Show the plot
    plt.show()


#------------------------------------------------------------
# 4 Scatterplot for correlation betweeen true and estimated cell fractions after deconvolution per cell type
#------------------------------------------------------------
def scatter_corr_Tf_Df(true_fractions_file, deconvolved_fractions_file, metadata_file, output_path):
    metadata = pd.read_csv(metadata_file)
    true_fractions = pd.read_csv(true_fractions_file, index_col = 0).sort_index()
    deconvolved_fractions = pd.read_csv(deconvolved_fractions_file, index_col=0).sort_index()

    print("True fractions:")
    print(true_fractions.head())
    print(true_fractions.shape)
    print("Deconvolved fractions:")
    print(deconvolved_fractions.head())
    print(deconvolved_fractions.shape)

    # Specify the cell type you're interested in
    cell_type_list = ['pDC'] #, 'CD8_Tcell', 'CD4_Tcell', 'T_reg', 'Malignant']

    # Extract case_id and cohort from the CSV file
    cohort = metadata.set_index('case_id')['Cohort'].sort_index()
    print(cohort)

    # Ensure that true_fractions, predicted_fractions, and cohort are aligned by sample names (index)
    # This assumes the sample names are the index in all dataframes
    assert true_fractions.index.equals(deconvolved_fractions.index), "Sample names must match between true and predicted fractions"
    assert true_fractions.index.equals(cohort.index), "Sample names must match between true fractions and cohort"

    for cell_type in cell_type_list:
        print(f'Processing {cell_type}')
        # Extract the true and predicted fractions for the specified cell type
        true_cell_type = true_fractions[cell_type]
        deconvolved_cell_type = deconvolved_fractions[cell_type]

        # Map cohort categories to colors
        cohort_categories = cohort.unique()
        print(cohort_categories)

        colormap = plt.get_cmap('tab10', len(cohort_categories))

        # Create a color map for each cohort
        color_map = {coh: colormap(i) for i, coh in enumerate(cohort_categories)}

        # Get the colors for each sample based on their cohort
        sample_colors = cohort.map(color_map)
        print(sample_colors)
    # Create the scatter plot with colored dots based on cohort
        plt.figure(figsize=(8, 6))
        scatter = plt.scatter(true_cell_type, deconvolved_cell_type, c=cohort[1], cmap = sample_colors, alpha=0.8)
        plt.title(f'Scatter Plot of True vs Predicted Fractions for {cell_type} by Cohort')
        plt.xlabel('True Fraction')
        plt.ylabel('Predicted Fraction')
        plt.grid(True)
        print(scatter.legend_elements()[0])
        # Create legend handles
        #handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[coh], markersize=10) for coh in cohort_categories]
        plt.legend(handles=scatter.legend_elements()[0], labels=cohort_categories)
        # Add the legend with the correct labels
        #plt.legend(handles=handles, labels=cohort_categories, title='Cohorts')
        
        # Save fig
        plt.savefig(output_path + f'/{cell_type}_scatter_corr.png')
    print("Plots saved in ", output_path)
#------------------------------------------------------------
# 5 Combined scatterplot for correlation betweeen true and estimated cell fractions after deconvolution for all cell types 
#------------------------------------------------------------
def combined_scatter(true_fractions_file, deconvolved_fractions_file, metadata_file, output_path, output_path_zoomedin, custom_colors):
    print('Reading data...')
    metadata = pd.read_csv(metadata_file)
    cohort_df = metadata.set_index('case_id').sort_index()
    cohort_df = cohort_df[['Cohort']]
    print(cohort_df)

    true_fractions = pd.read_csv(true_fractions_file, index_col = 0).sort_index()
    deconvolved_fractions = pd.read_csv(deconvolved_fractions_file, index_col=0).sort_index()

    # Sort rows numerically based on the number in the sample name (e.g., Sample 1, Sample 2, ...)
    # deconvolved_fractions.index = deconvolved_fractions.index.map(lambda x: int(x.split()[-1]))
    # deconvolved_fractions = deconvolved_fractions.sort_index()

    # # Sort both rows (index) and columns
    # true_fractions = true_fractions.sort_index(axis=1)
    # deconvolved_fractions = deconvolved_fractions.sort_index(axis=1)

    # # Ensure that the row names in deconvolved_fractions match the row names in true_fractions
    # deconvolved_fractions = deconvolved_fractions.set_index(true_fractions.index)

    print("True fractions:")
    print(true_fractions.head())
    print(true_fractions.shape)
    print("Deconvolved fractions:")
    print(deconvolved_fractions.head())
    print(deconvolved_fractions.shape)

    # Melt the dataframes to long format
    df_true_long = true_fractions.reset_index().melt(id_vars='case_id', var_name='cell_type', value_name='true_fraction')
    df_pred_long = deconvolved_fractions.reset_index().rename(columns={'index':'case_id'}).melt(id_vars=['case_id'], var_name='cell_type', value_name='predicted_fraction')


    # Merge both dataframes on sample_id and cell_type
    df = pd.merge(df_true_long, df_pred_long, on=['case_id', 'cell_type'])
    print(df)
    # Create the scatter plot
    plt.figure(figsize=(10, 6))
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.subplots_adjust(right=0.7)

    # palette = sns.color_palette(cc.glasbey, n_colors=len(true_fractions.columns))
    sns.scatterplot(data=df, x='true_fraction', y='predicted_fraction', hue='cell_type', palette=custom_colors, alpha=0.8)

    # Labels and title
    plt.plot([0, 1], [0, 1], 'r--', label="Identity Line")
    plt.title('Scatterplot of True vs Predicted Cell Fractions')
    plt.xlabel('True Cell Fractions')
    plt.ylabel('Predicted Cell Fractions')

    plt.legend(title='Cell Type', bbox_to_anchor=(1.4, 0.5), loc='center right')

    plt.savefig(output_path)
    print("Plot saved in ", output_path)

    plt.xlim(0, 0.3)
    plt.ylim(0, 0.3)

    # Title and layout for zoomed-in plot
    plt.title("Zoomed-in Scatterplot (0 to 0.3)")
    plt.legend(title="Cell Type", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    # Save the zoomed-in plot
    plt.savefig(output_path_zoomedin)

#------------------------------------------------------------
# 6 Scatterplot for malignant cell fraction of the actual bulk (GLASS-NL)
#------------------------------------------------------------
def plot_cell_fraction_scatter(true_fraction, deconvolved_fraction, output_path):
    """
    Plots a scatterplot comparing true and estimated cell fractions.

    Parameters:
        true_csv (str): Path to the CSV file with true cell fractions.
        estimated_csv (str): Path to the CSV file with estimated cell fractions.
        index_col (int or str): Column to use as the index when reading the CSVs. Default is 0.
    """

    # Load the CSV files
    true_df = pd.read_csv(true_fraction, index_col=0)
    estimated_df = pd.read_csv(deconvolved_fraction, index_col=0)

    # Ensure the dataframes are aligned
    true_df = true_df.sort_index()
    estimated_df = estimated_df.sort_index()

    # Extract the column of interest
    true_df = true_df['Malignant']
    estimated_df = estimated_df['Malignant']

    # # Flatten the dataframes
    # true_flat = true_df.values.flatten()
    # estimated_flat = estimated_df.values.flatten()

    # Combine into a single DataFrame to handle NaNs
    combined = pd.concat([true_df, estimated_df], axis=1)
    combined.columns = ['True', 'Estimated']

    # Drop rows with NaN in either True or Estimated
    combined = combined.dropna()

    # Extract cleaned values
    true_df = combined['True']
    estimated_df = combined['Estimated']
    

    # Calculate Pearson correlation coefficient
    pcc, _ = pearsonr(true_df, estimated_df)

    # Calculate Root Mean Square Deviation
    rmsd = np.sqrt(mean_squared_error(true_df, estimated_df))

    # Print the results
    print(f"Pearson Correlation Coefficient (PCC): {pcc:.3f}")
    print(f"Root Mean Square Deviation (RMSD): {rmsd:.3f}")

    # Plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=true_df, y=estimated_df, alpha=0.7)
    plt.xlabel('True Cell Fractions')
    plt.ylabel('Estimated Cell Fractions')
    plt.title('Scatterplot of True vs. Predicted Cell Fraction of Malignant cells')
    plt.plot([0, 1], [0, 1], 'r--', label='Identity line')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)

    
def main():
    ## Define a color dictionary
    color_map = {'B_cell':'royalblue',
                'CD4_Tcell':'darkorange',
                'CD8_Tcell':'mediumseagreen',
                'Endothelial':'crimson',
                'Fibroblast':'darkviolet',
                'Malignant':'sienna',
                'Monocyte':'darkkhaki',
                'NK_cell':'palevioletred',
                'Oligodendrocyte':'mediumturquoise',
                'Pericyte':'lightblue',
                'Plasma_B':'dimgrey',
                'TAM':'lightpink',
                'T_reg':'mediumpurple',
                'cDC':'chocolate',
                'pDC':'firebrick',
                'Microglia':'lightcoral',
                'Macrophage':'steelblue'
                }
    # signature = 
    # scExp = 
    # scVar = 
    # check_signature_scExp_scVar(signature, scExp, scVar)

    true_fractions = "/net/beegfs/cfg/tgac/dmartinovicova/Statescope/test/data/BLADE_wo_seed/GT_all.csv"
    deconvolved_fractions = "/net/beegfs/cfg/tgac/dmartinovicova/Statescope/test/data/BLADE_wo_seed/Estimated_fractions_all.csv"
    output_path_corr = "/net/beegfs/cfg/tgac/dmartinovicova/figures/pseudobulk_deconvolution/wo_seed/all_l1_statescope_corr.png"
    output_path_RMSD = "/net/beegfs/cfg/tgac/dmartinovicova/figures/pseudobulk_deconvolution/wo_seed/all_l1_statescope_rmsd.png"
    box_output_path_RMSD = "/net/beegfs/cfg/tgac/dmartinovicova/figures/pseudobulk_deconvolution/wo_seed/all_l1_statescope_box_rmsd.png"
    
    correlation_Tf_Df(true_fractions, deconvolved_fractions, output_path_corr, color_map)
    RMSD_Tf_Df(true_fractions, deconvolved_fractions, output_path_RMSD, color_map)
    box_RMSD_Tf_Df(true_fractions, deconvolved_fractions, box_output_path_RMSD, color_map)

    metadata = "/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/data/metadata/metadata_case_id.csv"
    output_path = "/net/beegfs/cfg/tgac/dmartinovicova/figures/pseudobulk_deconvolution/wo_seed/all_l1_statescope_scatter.png"
    output_path_zoomedin = "/net/beegfs/cfg/tgac/dmartinovicova/figures/pseudobulk_deconvolution/wo_seed/all_l1_scatter_zoomedin.png"
    # scatter_corr_Tf_Df(true_fractions, deconvolved_fractions, metadata, output_path)
    combined_scatter(true_fractions, deconvolved_fractions, metadata, output_path, output_path_zoomedin, color_map)
    #plot_cell_fraction_scatter(true_fractions, deconvolved_fractions, output_path)

    #groups_meta = '/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_DNA-seq/GLASS_comparison_groups.csv'
    # output_boxplot = "/net/beegfs/cfg/tgac/dmartinovicova/figures/pseudobulk_deconvolution/check"
    # boxplot_cell_fractions_in_groups(deconvolved_fractions, groups_meta, output_boxplot)
    
    print("Done")

main()