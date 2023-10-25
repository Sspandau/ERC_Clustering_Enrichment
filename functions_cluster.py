import pandas as pd
import numpy as np

def ranksum_mat(df, proteins):
    dct1 = {}
    for col in range(len(dfpct.axes[1])):
        prot = proteins[col]
        if prot not in dct1.keys():
            print(prot)
            dct1[prot] = []
            for row in range(len(dfpct.axes[0])):
                r1 = dfpct.iloc[col, row]
                r2 = dfpct.iloc[row, col]
                rsum = r2 + r1
                dct1[prot].append(rsum)
        else:
            continue
    dfpct2 = pd.DataFrame.from_dict(dct1)
    dfpct2.index = proteins
    return dfpct2

def partneroverlap_mat(df, proteins, pct):
    thresh = round((len(proteins)*pct)/100)
    dct1 = {}
    for prot in proteins:
        #print(prot)
        dfpct = dfpct.sort_values(prot)
        rankl = dfpct['proteins'].to_list()
        dct1[prot] = rankl[0:thresh]
    dfpct2 = pd.DataFrame.from_dict(dct1)
    dct = {}
    for prot in proteins:
        print(prot)
        dct[prot] = []
        list1 = dfpct2[prot]
        for col in range(len(dfpct2.axes[1])):
            list2 = dfpct2.iloc[:, col]
            dct[prot].append(len(list(set(list1) & set(list2))))
    dfpct3 = pd.DataFrame.from_dict(dct)
    return dfpct3

def enrichment_clusters(df1, gene_list):
    for i in range(df1.cluster.unique()): #len(df.subcluster.to_list()) df1.cluster.unique()
        df1sub = df1[df1['cluster'] == i]
        print(i)
        gene_list2 = df1sub.proteins
        if len(df1sub.index) != 1:
            if i == 1:
                enr_bg = gp.enrichr(gene_list=gene_list2, 
                                    gene_sets=['GO_Biological_Process_2021','GO_Molecular_Function_2021','GO_Cellular_Component_2021','KEGG_2021_Human', 'WikiPathway_2021_Human'],
                                    background=gene_list,
                                    outdir=None) # don't write to disk
                enr_bg = pd.DataFrame(enr_bg.results)
                enr_bg['Cluster'] = i
            else:
                enr_bg1 = gp.enrichr(gene_list=gene_list2, 
                                     gene_sets=['GO_Biological_Process_2021','GO_Molecular_Function_2021','GO_Cellular_Component_2021','KEGG_2021_Human', 'WikiPathway_2021_Human'],
                                     background=gene_list,
                                     outdir=None)
                enr_bg1 = pd.DataFrame(enr_bg1.results)
                enr_bg1['Cluster'] = i
                enr_bg = pd.concat([enr_bg, enr_bg1], axis=0)
        else:
            continue
        
    return enr_bg

def enrichment_subclusters(df1, gene_list):
    for i in df.subcluster.unique():
        print(i)
        dfsub = df[df['subcluster'] == i]
        for j in dfsub.cluster.unique():
            
            dfsub2 = dfsub[dfsub['cluster'] == j]
            gene_list2 = dfsub2.proteins
            if i == 1 and j == 0:
                enr_bg = gp.enrichr(gene_list=gene_list2, 
                                    gene_sets=['GO_Biological_Process_2021','GO_Molecular_Function_2021','GO_Cellular_Component_2021','KEGG_2021_Human', 'WikiPathway_2021_Human'],
                                    background=gene_list,
                                    outdir=None) # don't write to disk
                enr_bg = pd.DataFrame(enr_bg.results)
                enr_bg['Cluster'] = i
                enr_bg['Sub Cluster'] = j
            else:
                enr_bg1 = gp.enrichr(gene_list=gene_list2, 
                                     gene_sets=['GO_Biological_Process_2021','GO_Molecular_Function_2021','GO_Cellular_Component_2021','KEGG_2021_Human', 'WikiPathway_2021_Human'],
                                     background=gene_list,
                                     outdir=None)
                enr_bg1 = pd.DataFrame(enr_bg1.results)
                enr_bg1['Cluster'] = i
                enr_bg1['Sub Cluster'] = j
                enr_bg = pd.concat([enr_bg, enr_bg1], axis=0)
        
    return enr_bg

def filter_enrichment(enr_bg, clustertype):
    if clustertype == "subcluster":
        top_results = {}
        for i in enr_bg['Cluster'].unique():
            print(i)
            sub_df1 = enr_bg[enr_bg['Cluster'] == i]
            for j in sub_df1['Sub Cluster'].unique():
                print(j)
                sub_df = enr_bg[(enr_bg['Sub Cluster'] == j) & (enr_bg['Cluster'] == i)]
                top_results[f'{i}'+"."+f'{j}'] = sub_df.iloc[0, :]
    else:
        top_results = {}
        for i in enr_bg['Cluster'].unique():
            print(i)
            sub_df = enr_bg[enr_bg['Cluster'] == i]
            top_results[f'{i}'] = sub_df.iloc[0, :]
            
    top_results = pd.DataFrame.from_dict(top_results, orient='index')
    return top_results    