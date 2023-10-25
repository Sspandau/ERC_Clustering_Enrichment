# ERC Clustering and Enrichment Analysis

The analysis in this repository groups proteins with correlated evolutionary rates and tests for enriched pathways in those resulting clusters. 
ERCs use a phylogeny-based approach to detect evolutionary signals between pairs of proteins, which may be evidence of shared functions among those pairs.

## Installation

### Environment Requirements

The scripts in this analysis require Python version > 3.9 and R version > 4.1.

### Python Package Dependencies
The external Python dependencies can be installed using pip:
`pip install -r requirements.txt`

### R Package Dependencies
The external R packages can be installed using `install.packages()` function in base R.

*ctc package requires Bioconductor installation with the following commad: `BiocManager::install("ctc")`

### Preparation of Data
Please see ERC-Pipeline repository for instructions on how to calculate the evolutionary rate correlations and the subsequent rank matrix of ERC partners. 
The rank matrix must be in csv format. 

## Preparing the Distance Matrices for Clustering
Execute the following command to generate the distance matrix for clustering.

`python3 matrix.py --input erc_rank_mat.csv --output erc_toppartners_matrix.csv --distance_type partneroverlap --percent 1`

*To specify which distance to use, you must indicate "ranksum" or "partneroverlap" by passing `--distance_type ranksum`. The default distance is ranksum.

*If you are using the partner overlap distance metric, then you need to indicate what percentage of top ranking partners are to be compared between pairs of proteins by pass `--percent 1`. The default percentage is 1.

## Clustering in R
The following analysis is in R.

### Clustering
Execute the following command to generate a newick file, cluster id csv, and cluster counts csv for the ranksum clustering.

`R clustering_ranksum.R path/to/ranksum_matrix.csv path/to/newick.newick path/to/cluster_id.csv path/to/cluster_counts.csv k 150`

Execute the following command to generate a newick file, cluster id csv, and cluster counts csv for the partner overlap clustering.

`R clustering_partneroverlap.R path/to/partneroverlap_matrix.csv path/to/newick.newick path/to/cluster_id.csv path/to/cluster_counts.csv k 150`

*The last two arguments indicate how you want to stop the hierarchical clustering and at what threshold. 

'k' indicates that you want a predetermined number of clusters, and '150' indicates the number of desired clusters in this case. 

'h' is the other option for this argument, which stops the clustering algorithm based on the distance threshold between clusters. '150' would indicate that the algorithm should stop clustering before the distance between a pair of proteins in a cluster exceeds 150.

### Subclustering
Execute the following command to generate subcluster id csv and subcluster counts csv for the ranksum clustering.

`R ranksum_subclustering.R path/to/ranksum_matrix.csv path/to/cluster_id.csv path/to/subcluster_id.csv path/to/subcluster_counts.csv 200`

Execute the following command to generate subcluster id csv and subcluster counts csv for the ranksum clustering.

`R partneroverlap_subclustering.R path/to/partneroverlap_matrix.csv path/to/cluster_id.csv path/to/subcluster_id.csv path/to/subcluster_counts.csv 200`

*The last argument indicates which clusters should be split into subclusters based on that size. '200' would indicate that any cluster with > 200 proteins should be split into subclusters.

## Enrichment with Enrichr interface in Python


