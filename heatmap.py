
###########################################
# IMPORTS
###########################################

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

import argparse

###########################################
# PARSING
###########################################

def get_file_lines(filepath):
    lineStrings = []
    FILE = open(filepath, "r")
    for line in FILE:
        if (line[-1] == "\n"): lineStrings.append(line[:-1])
        else: lineStrings.append(line)
    FILE.close()
    return lineStrings

def lines_to_dataframe(lines):

    names = []
    dist_dict = dict()
    for line in lines:

        linarr = line.split()
        percents = list(map(lambda x : float(x), linarr[2:]))
        anno = linarr[1].split("|")[1]

        if anno in dist_dict.keys(): print(anno)
        dist_dict[anno] = percents
        names.append(anno)

    DF = pd.DataFrame(columns = names)
    for sample in names: DF.loc[sample] = dist_dict[sample]

    return DF


def pim_to_dataframe(filepath : str) -> pd.DataFrame:
    lines = get_file_lines(filepath)
    return lines_to_dataframe(lines)

###########################################
# VISUALIZATION
###########################################

def heatmap(DF : pd.DataFrame, DATE : str) -> None:

    # create heatmap
    ax = sns.heatmap(data = DF,
                     xticklabels = False, yticklabels = 10)

    # set plot / axes settings
    plt.yticks(fontsize = 5)
    ax.set_xlabel("Sample ID.")

    # show heatmap
    # plt.show()
    plt.savefig("heatmap-" + DATE + ".png")

###########################################
# MAIN METHOD
###########################################

def main():

    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("filepath",
                        help="Path of file containing sampler arguments.")
    args = parser.parse_args()

    # parse file input
    DF = pim_to_dataframe(args.filepath)
    DATE = args.filepath.split("-")[2]
    heatmap(DF, DATE)

main()