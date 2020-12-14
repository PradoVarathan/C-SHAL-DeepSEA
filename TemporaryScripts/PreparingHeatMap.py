# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 08:05:04 2020

@author: ppugale
"""

import pandas as pd
import os
import seaborn as sns
from matplotlib import pyplot as plt
import pickle as pkl
import numpy as np
import fastcluster
os.chdir("C:/Users/ppugale/OneDrive - Indiana University/Documents/DeepSEA_Lab")
ss_file = pd.read_csv("SNPs_Thres_0.05",sep=" ")


final_results = pkl.load(open("Final_Output.pkl","rb"))
results = np.array(final_results[i] for i in final_results.keys())
label_names = pd.read_csv("Labels.csv")
Labels = list(label_names['Label_Name'])
results_df = pd.DataFrame.from_dict(final_results, orient="index",columns=Labels)


df = results_df.loc[snps_top_20]
my_palette = dict(zip(df.cyl.unique(), ["orange","yellow","brown"]))
row_colors = df.cyl.map(my_palette)
 
# plot
sns.clustermap(df)


snps_top_p = ['rs35374935','rs10752132','rs4841512','rs12968854','rs373806','rs2823088','rs1178163','rs11074388','rs9860254','rs16979118','rs205144','rs7668273','rs4805143','rs6549972','rs12406816','rs190707322','rs2025440','rs10269719','rs14781','rs4443214']

df = results_df.loc[snps_top_p]
my_palette = dict(zip(df.cyl.unique(), ["orange","yellow","brown"]))
row_colors = df.cyl.map(my_palette)

# plot
sns.clustermap(df)
fastcluster.