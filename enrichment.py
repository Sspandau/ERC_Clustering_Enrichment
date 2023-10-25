import sys
import click
import gseapy as gp
import pandas as pd
import numpy as np
from functions_cluster import *

@click.command()
@click.option('--Input', help="The path and name of the cluster csv from r scripts")
@click.option('--Output', help="The path and name of the output csv")
@click.option('--clustering_level', default='cluster', help="Indicate whether enrichment is for clusters or subclusters")
@click.option('--Filter', is_flag=True, help="Indicate whether to filter for top result in each cluster")


df1 = pd.read_csv(Input)

df1 = df1.iloc[:, 1:]

gene_list = df1.proteins.to_list()

go_mf = gp.get_library(name='GO_Molecular_Function_2021', organism='Human')
go_bg = gp.get_library(name='GO_Biological_Process_2021', organism='Human')
go_cp = gp.get_library(name='GO_Cellular_Component_2021', organism='Human')
kegg = gp.get_library(name='KEGG_2021_Human', organism='Human')
wp = gp.get_library(name='WikiPathway_2021_Human', organism='Human')

if clustering_level == "cluster":
    enr_bg = enrichment_clusters(df1, gene_list)
    if Filter:
        enr_bg = Filter(enr_bg, "cluster")
elif clustering_level == "subcluster":
    enr_bg = enrichment_subclusters(df1, gene_list)
    if Filter:
        enr_bg = Filter(enr_bg, "subcluster")
        
dfpct2.to_csv(Output)
