import os
import pickle as pl
import pandas as pd
os.chdir("")

final_log_output = pl.load(open("Final_Output.pkl","rb"))
labels = pd.read_csv("Data/Labels.csv")
label_output = {}

for rsid in final_log_output.keys():
    min_temp = max(final_log_output[rsid])
    indx = list(final_log_output[rsid]).index(min_temp)
    print(indx)
    label_output[rsid] = [labels.Label_Name[indx],min_temp]
    
Top_Labels = pd.DataFrame.from_dict(label_output, orient='index')
Top_Labels = Top_Labels.sort_values(by=[1],ascending=False)
Top_Labels.to_csv("Log_Results.csv")

