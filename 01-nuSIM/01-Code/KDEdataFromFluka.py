'''
    KDEdataFromFluka.py

  Assumes that nuSim code is in python path.

    @file KDEdataFromFluka.py

    @brief   Generate a pion according the KDE derived from the fluka distributions

    @author  Paul Kyberd

    @version     1.0
    @date        08 October 2024

    KDE distribtuion output file as a variable

    @version     1.1
    @date        14 October 2024

    

    @param  

    The main routine runs a set of tests checcking that things work - but no real
    validation
    The code returns the px, py, pz, x, y, xp, yp of the pion
        xp = px/pz
        yp = py/pz
   
    Version 1.0                                 07 October 2024
    @Author:                                     Paul Kyberd
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import traceSpace as tSp
import scipy as sp

class KDEdataFromFluka:

    def __init__(self):
        print("instantiate")
 
    def beamGen_kde(self, data, columns, num_samples):
        kde_data = data[columns]
    # Compute Scott's factor for KDE bandwidth
        scott_factor = kde_data.shape[0]**(-1./(kde_data.shape[1]+4))
        scott_factor /= 10 # Scott's factor is too large for our purposes, so we reduce it

    # Estimate underlying distribution of the data using KDE
        kde = sp.stats.gaussian_kde(kde_data.T,bw_method=scott_factor)

    # Generate samples from the KDE distribution
        samples = kde.resample((num_samples,)).T

    # Create a new dataframe to store the samples
        kde_sample_df = pd.DataFrame(samples, columns=columns)
    
        return kde_sample_df


    def kde_plotting(self, df,new_df):
    # quick plotting to compare the original and the KDE generated data 
        plt.hist(df['x'][:10000], bins=60, range = [-15,15], alpha=0.5, label='Original')
        plt.hist(new_df['x'][:10000], bins=60, range = [-15,15], alpha=0.5, label='KDE')
        plt.legend(loc='upper right')
        plt.show()

        plt.hist(df['xp'][:10000], bins=60, range = [-0.05,0.05], alpha=0.5, label='Original')
        plt.hist(new_df['xp'][:10000], bins=60, range = [-0.05,0.05], alpha=0.5, label='KDE')
        plt.legend(loc='upper right')
        plt.show()

        plt.hist(df['p'][:10000], bins=60, range = [1.7,2.3], alpha=0.5, label='Original')
        plt.hist(new_df['p'][:10000], bins=60, range = [1.7,2.3], alpha=0.5, label='KDE')
        plt.legend(loc='upper right')
        plt.show()

        plt.hist(df['p'][:10000], bins=60, range = [1.7,2.3], alpha=0.5, label='Original')
        plt.hist(new_df['p'][:10000], bins=60, range = [1.7,2.3], alpha=0.5, label='KDE')
        plt.legend(loc='upper right')
        plt.show()

        plt.hist2d(df['x'][:10000],df['xp'][:10000],bins=[50,50], range = [[-15,15],[-0.05,0.05]])
        plt.title('Original data x-px phase space')
        plt.show()
        plt.hist2d(new_df['x'][:10000], new_df['xp'][:10000],bins=[50,50], range = [[-15,15],[-0.05,0.05]])
        plt.title('KDE data x-px phase space')
        plt.show()

    def readFile(self, flukaDataFile):
    #load Fluka data into a pandas dataframe
        df = pd.read_csv(flukaDataFile, sep = "\s+", header=None, names=['x','y','z','px','py','pz','xp','yp','E'])
    #calculate total momentum
        df['p'] = (np.sqrt(df['px']**2+df['py']**2+df['pz']**2)).to_list()

        return df

    def writeFile(self, new_df, fileName):

        new_df.to_csv(fileName, sep=' ', index=False, header=None)

        return
 
    def openGenFile(self, fileName):
        self._genFile = open(fileName, "r")

    def readGenFile(self):

        data = self._genFile.readline()
        data = data.strip()
        data = data.split()
        if not data: return None, None
        #   create a trace space vector, and alter the distance from (assumed) mm to m
        pion = tSp.traceSpace(0.0, float(data[0])/100.0, float(data[1])/100.0, -50.0, float(data[2]), float(data[3]))
        pPion = float(data[5])

        return pion, pPion
 

def main():

    nSamples =1000
    #   Instantiate the KDE class
    flukaDataFile = '/Users/paulkyberd/work/nuSTORMAug24/01-nuSIM/31-Target/v2/selected_pions_5GeV_2.0mm.txt'
    KDEdata = KDEdataFromFluka()

    df = KDEdata.readFile(flukaDataFile)

    print ("df is \n", df)

    #num_samples is used set to the new sample size
    new_df = KDEdata.beamGen_kde(df, columns=['x', 'y', 'xp', 'yp', 'pz', 'p'], num_samples=nSamples)

    print (" new_df is \n", new_df)
    #quick plotting for visual inspection
    KDEdata.kde_plotting(df,new_df)

    #save the new data to a text file
    genDataFile = "testing_kde11.txt"
    KDEdata.writeFile(new_df, genDataFile)

    #   Read the generated data file
    KDEdata.openGenFile(genDataFile)

    #   Read the generated data file
    for i in range(nSamples):
        pion, pPion =  KDEdata.readGenFile()
        if pion == None: break
        print ("main: pion is ", pion)
        print ("main: pPion is ", pPion)


if __name__ == "__main__" :
    main()

22