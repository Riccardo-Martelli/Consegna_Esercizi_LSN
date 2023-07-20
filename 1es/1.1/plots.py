#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def main():

 plt.rcParams.update({'font.size': 18})

 df = pd.read_csv("test.txt")

 M=100000              # Total number of throws
 N=100                 # Number of blocks
 L=int(M/N)            # Number of throws in each block, please use for M a multiple of N

 x = np.arange(N)
 x*=L

 plt.figure(figsize=(12,8))
 plt.title("Error of the mean")
 plt.errorbar(x,df["mean"]-0.5,yerr=df["error_mean"])
 plt.xlabel('#throws')
 plt.ylabel('<r>-1/2')
 plt.grid(True)
 
 plt.figure(figsize=(12,8))
 plt.title("Error on the std")
 plt.errorbar(x,df["std"]-1/12.,yerr=df["error_std"])
 plt.xlabel('#throws')
 plt.ylabel('<(r-0.5)^2>-1/12')
 plt.grid(True)
 plt.show()

if __name__ == "__main__":
    main()

