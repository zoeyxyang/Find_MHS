import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from collections import Counter


'''
Input: a list of NSC IDs

Output: a list of Drug Names
'''

def MatchIdName(IDs, dataDir):
    
    df = pd.read_csv(dataDir +"ID_name.csv", usecols = ["NSC", "DrugName"],index_col=0)
    ID_Name_map = df.to_dict()["DrugName"]
    drugNames = []
    for ID in IDs:
        if ID_Name_map[ID] == "-":
            drugNames.append("(no name)")
        else: 
            drugNames.append(ID_Name_map[ID])
    return drugNames


'''
Input: a list of the hitting sets sizes, the directory to be saved in, the algorithm used

Output: save the bar plot in the directory "plots"
'''
def PlotBar(HS_sizes, plotDir, algo):
    x = list(Counter(HS_sizes).keys()) # equals to list(set(words))
    x = map(str, x)
    x = list(x)
    y = list(Counter(HS_sizes).values())
    fig = plt.figure()
    #  Bar plot
    plt.bar(x, y, color ='blue',
            width = 0.25)
    plt.xlabel("hitting set sizes")
    plt.ylabel("number of occurance")
    plt.title(algo)
    plt.savefig(plotDir+algo+'_barplot.png')

    
    
    
    