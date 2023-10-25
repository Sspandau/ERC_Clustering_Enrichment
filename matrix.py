import sys
import click
import pandas as pd
import numpy as np
from functions_cluster import *

@click.command()
@click.option('--distance_type', default = "ranksum", help="The type of distance metric to use")
@click.option('--Input', help="The rank matrix for clustering")
@click.option('--output', help="Path and name of ouput distance matrix")
@click.option('--percent', default=1, help="What percentage of two proteins' top ranking partners should be compared")

dfpct = pd.read_csv(Input)
dfpct = dfpct.iloc[:, 1:]
proteins = dfpct.columns.to_list()
dfpct.index = proteins

if distance_type == "ranksum":
    dfpct2 = ranksum_mat(dfpct, proteins)
elif distance_type == "partneroverlap":
    dfpct2 = partneroverlap_mat(dfpct, proteins, percent)

dfpct2.to_csv(Output)

