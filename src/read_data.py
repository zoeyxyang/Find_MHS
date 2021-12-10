#set the threshold (effective or not effective)
#threshold = 0
#tissues = ["breast", "colon", "prostate"]

import pandas as pd
import numpy as np
import csv
from itertools import chain
import os

'''
#cd the current directory
os.chdir('../')
outDir = os.getcwd()
dataDir = outDir +"/data/"
matrixDir = outDir +"/matrix/"
'''



def get_matrix(threshold, tissues, dataDir, matrixDir ):
    
    with open(dataDir+'simplified.csv') as infile:
        reader = csv.reader(infile)
        headers = next(reader)

    header_indices = [i for i, item in enumerate(headers) if item]
    df = pd.read_csv(dataDir+'simplified.csv', usecols = header_indices)
    df = df.set_index('NSC')
    df = df.replace("na", -100) #if there's no data, we'll assume there is no effect
    df = df. astype(float)
    bin_df, csv_path = specific_zero_one_matrix(df,tissues,threshold, matrixDir)
    return bin_df, csv_path




#Enter tissue type as string parameter of function
def specific_zero_one_matrix(df,tissues,threshold, matrixDir):
    tissuedict = {"breast":range(0,5),
                  "cns":range(5,11),
                  "colon":range(11,18),
                  "leukemia":range(18,24),
                  "melanoma":range(24,34),
                  "non-small cell lung":range(34,43),
                  "ovarian":range(43,50),
                  "prostate":range(50,52),
                  "renal":range(52,60),
                   "all":range(0,60)}
    rng =  []
    name = ""
    for tissue in tissues:
        rng += tissuedict[tissue]
        name += tissue+"_"

    bin_df = df.iloc[:,rng]
    bin_df = bin_df.where(bin_df > threshold, 0)
    bin_df = bin_df.where(bin_df == 0, 1)
    bin_df.to_csv(matrixDir+name + '_zero_one_matrix.csv', index=True) #save the matrix
    csv_path = matrixDir+name + '_zero_one_matrix.csv'
    return bin_df, csv_path

#specific_zero_one_matrix(tissues)











