import argparse
from read_data import get_matrix
from functions import PlotBar
import os



#User Defined Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--tissues", type=str, nargs = "+",
                    help="list of target tissues")
parser.add_argument("-th", "--threshold", type=int,
                    help="the threshold of drug activity")
parser.add_argument("-a", "--algorithms", type=str,
                    help="algorithm to be used")
args = parser.parse_args()

print("set tissues as:", args.tissues)
print("set threshold as", args.threshold)
print("set algorithm as", args.algorithms)


threshold = args.threshold
tissues = args.tissues
algo = args.algorithms


#Set paths for the directory

os.chdir('../')
outDir = os.getcwd() #project directory
dataDir = outDir +"/data/" #project/data
matrixDir = outDir +"/matrix/"#project/matrix
plotDir = outDir + "/plots/"


# get the 0/1 matrix and store it in "matrix"
bin_df, csv_path = get_matrix(threshold, tissues, dataDir, matrixDir)

#execute specified algorithms 
if algo == "SA":
    from SA_algo import SA
    HS_sizes = SA(bin_df, dataDir)
    PlotBar(HS_sizes, plotDir, algo="SA")
    
elif algo == "ILP":
    from ILP_algo import MHS_ILP
    MHS_ILP(csv_path, dataDir)
    
elif algo == "HDF":
    from HDF_algo import HDF
    HS_sizes = HDF(bin_df, dataDir, 1)
    PlotBar(HS_sizes, plotDir, algo="HDF")
    
elif algo == "GD":
    from GD_algo import GD
    GD(bin_df, dataDir)
    
    
else:
    print("algorithm", algo, "doesn't exist")
    
