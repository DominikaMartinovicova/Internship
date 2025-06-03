configfile: "config.yaml"


# Perform normalization and filtering

data_dir = config["all"]["data_dir"]
plot_dir = config["all"]["plot_dir"]

rule all:
  input:
      data_dir + "phenotyping/adata_final.h5ad"
      
 
rule QC_scRNA:
    input:
        adata = data_dir + "adata/{dataset}.h5ad"
    output:
        processed_data = data_dir + "processed/{dataset}_processed.h5ad"   #mcpg = minimal cells per gene
    params:
        max_genes=lambda wildcards: config["filtering"]["specific"].get(wildcards.dataset, {}).get("max_genes", config["filtering"]["default"]["max_genes"]),
        cache_dir = ".cache/{dataset}/",
        plot_dir = plot_dir + "{dataset}/preprocessing/",
        min_genes=config["filtering"]["default"]["min_genes"],
        max_mt_pct=config["filtering"]["default"]["max_mt_pct"],
        min_cells=config["filtering"]["default"]["min_cells"]
    shell: 
        """
        python scripts/QC_scRNA.py \
            -i {input.adata} \
            -o {output.processed_data} \
            -max_genes {params.max_genes} \
            -min_genes {params.min_genes} \
            -max_mt_pct {params.max_mt_pct} \
            -min_cells {params.min_cells} \
            -cache_dir {params.cache_dir} \
            -plot_dir {params.plot_dir}
        """

      
rule concatenate_scRNA:
    input:
        processed_adata_dir = data_dir + "processed/",
        phenotype_map = data_dir + "annotations_markers/mapping.csv"
    output:
        adata_concatenated = data_dir + "combined_all/adata_combined.h5ad"
    params:
        cache_dir = ".cache/"
    shell:
      """
        python3 scripts/concatenate_scRNA.py \
          -i {input.processed_adata_dir} \
          -phen_map {input.phenotype_map} \
          -o {output.adata_concatenated} \
          -cache_dir {params.cache_dir} 
        
      """


rule normalize_scRNA:
    input:
        adata_concatenated = data_dir + "combined_all/adata_combined.h5ad"
    output:
        adata_normalized = data_dir + "combined_all/adata_normalized.h5ad"  #without scaling only normalized and log1p transformed
    params:
        cache_dir = ".cache",
        plot_dir = plot_dir
    shell:
      """
        python3 scripts/normalize_scRNA.py \
          -i {input.adata_concatenated} \
          -o {output.adata_normalized} \
          -plot_dir {params.plot_dir} \
          -cache_dir {params.cache_dir}  
              
      """
      
rule ingest_annotations_scRNA:
    input:
        adata_normalized = data_dir + "combined_all/adata_normalized.h5ad"
    output:
        adata_un_ingested_scaled = data_dir + "combined_all/adata_un_ingested_scaled.h5ad",
        adata_ingested_scaled = data_dir + "combined_all/adata_ingested_scaled.h5ad"
    params:
        cache_dir = ".cache"
    shell:
      """
        python3 scripts/ingest_annotations_scRNA.py \
          -i {input.adata_normalized} \
          -o {output.adata_un_ingested_scaled} \
          -o_comb {output.adata_ingested_scaled}\
          -cache_dir {params.cache_dir}  
              
      """
      

rule clustering_scRNA:
    input:
        adata_ingested_scaled = data_dir + "combined_all/adata_ingested_scaled.h5ad"
    output:
        adata_clustered = data_dir + "phenotyping/adata_clustered.h5ad"
    params:
        cache_dir = ".cache/",
        plot_dir = plot_dir,
        n_pcs = config["clustering"]["n_pcs"],
        n_neighbors = config["clustering"]["n_neighbors"],
        resolution = config["clustering"]["resolution"],
        batch_correction = False
    shell:
      """
        python3 scripts/clustering_scRNA.py \
        -i {input.adata_ingested_scaled} \
        -o {output.adata_clustered} \
        -n_pcs {params.n_pcs} \
        -n_neighbors {params.n_neighbors} \
        -resolution {params.resolution} \
        -cache_dir {params.cache_dir} \
        -plot_dir {params.plot_dir} \
        -bbknn {params.batch_correction}

      """

rule subset_scRNA:
    input:
        adata_ingested_scaled = data_dir + "combined_all/adata_ingested_scaled.h5ad"
    output:
        immune = data_dir + "phenotyping/immune_subset.h5ad",
        myeloid = data_dir + "phenotyping/myeloid_subset.h5ad"
    params:
        cache_dir = ".cache/",
    shell:
     """
       python3 scripts/subset_scRNA.py \
         -i {input.adata_ingested_scaled} \
         -o1 {output.immune} \
         -o2 {output.myeloid} \
         -cache_dir {params.cache_dir} 
     """

rule reclustering_scRNA:
    input:
        adata_ingested_scaled = data_dir + "phenotyping/{phenotype}_subset.h5ad"
    output:
        adata_clustered = data_dir + "phenotyping/{phenotype}_clustered.h5ad"
    params:
        cache_dir = ".cache/",
        plot_dir = plot_dir,
        n_pcs = config["clustering"]["n_pcs"],
        n_neighbors = config["clustering"]["n_neighbors"],
        resolution = config["clustering"]["resolution"],
        batch_correction = False
    shell:
      """
        python3 scripts/clustering_scRNA.py \
        -i {input.adata_ingested_scaled} \
        -o {output.adata_clustered} \
        -n_pcs {params.n_pcs} \
        -n_neighbors {params.n_neighbors} \
        -resolution {params.resolution} \
        -cache_dir {params.cache_dir} \
        -plot_dir {params.plot_dir} \
        -bbknn {params.batch_correction}

      """


rule label_scRNA:
    input:
        adata_clustered = data_dir + "phenotyping/{subset}_clustered.h5ad",
        labels = data_dir + "annotations_markers/labels_{subset}.csv"
    output:
        adata_labeled = data_dir + "phenotyping/{subset}_labeled.h5ad"
    params:
        cache_dir = ".cache/"
    shell:
       """
       python3 scripts/label_scRNA.py \
         -i {input.adata_clustered} \
         -l {input.labels} \
         -o {output.adata_labeled} \
         -cache_dir {params.cache_dir} 
       """

rule combine_labels_scRNA:
    input:
        adata_clustered = data_dir + "phenotyping/adata_labeled.h5ad",
        immune_clustered = data_dir + "phenotyping/immune_labeled.h5ad",
        myeloid_clustered = data_dir + "phenotyping/myeloid_labeled.h5ad"
    output:
        adata_final = data_dir + "phenotyping/adata_final.h5ad"
    params:
        cache_dir = ".cache/"
    shell:
       """
       python3 scripts/combine_labels_scRNA.py \
         -i_all {input.adata_clustered} \
         -i_i {input.immune_clustered} \
         -i_m {input.myeloid_clustered} \
         -o {output.adata_final} \
         -cache_dir {params.cache_dir} 
       """






