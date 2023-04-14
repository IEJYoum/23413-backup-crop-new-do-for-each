# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 14:34:40 2023

@author: youm
"""


'''
from cmifA44:
op = ["start over with raw data","log2","scale from -1 to 1", "z-score","elmarScale","trim outliers",
      "make control TMA sample sizes the same","combat",
      "apply TMA combat to other dataset","equalizeBiomLevel","adjust for negative values", 
      "save to csv", "pick subset of data", "manually threshold",
      "cluster by obs catagory","Leiden cluster","GMM cluster","K-means","aggregate",
      "manually celltype random training set","auto-cell-type",
      "convert df to fractions in obs categories","convert to superbiom-only df",
      "remove non-primary biomarkers","calculate biomarker expression in region around each cell",
      "count label fractions in neighborhood","calculate entropy in neighborhood",
      "select ROI","remove cells expressing certain biomarker combinations","pick random subset","clag","clauto"]
fn = [revert,log2,scale1,zscore,elmarScale,outliers,equalizeTMA,combat,TMAcombat,equalizeBiomLevel,remNegatives,save,pick,
      manThresh,obCluster,leiden,gmm,kmeans,aggregate,celltype,autotype,countObs,superBiomDF,
      onlyPrimaries,regionAverage,neighborhoodFractions,neighborhoodEntropy,roi,simulateTherapy,subset,clag,clauto]
'''

import os
import numpy as np
import pandas as pd
import math
import scipy


def main(df,obs,dfxy):
    obs["all"] = 'all'
    dfs = [df,obs,dfxy]
    dfa = []
    ch,uch = obMenu(obs,'repeat analysis on each unique value in:')
    #obcol = obs.columns[ch]
    nn,commands = mainMenu(dfs)
    #print(commands)
    #return()
    for uc in uch:
        key = obs.iloc[:,ch] == uc
        sdfs = []
        for d in dfs:
            sdfs.append(d.loc[key,:])
        #print(sdfs,'sdfs')
        sdfs,nn=mainMenu(sdfs,commands,uc)
        dfa.append(sdfs)
    odfs = []
    for i in range(3):
        bi = []
        for d in dfa:
            bi.append(d[i])
        odfs.append(pd.concat(bi,axis=0))
    return(odfs[0],odfs[1],odfs[2])

'''
main functions
'''


def menu(dfs,options,functions,com=[],cat=''):
    print(com,'com into menu')
    if len(com) == 0:
        coms = []
        while True:
            print("\n")
            for i,op in enumerate(options):
                print(i,op)
            try:
                print("send non-int when done (return to previous menu)")
                ch = int(input("number: "))
            except:
                print(coms,"coms out of menu")
                return([],coms)
            nn,com=functions[ch](dfs,com=[])
            coms.append([ch]+com)

    else:
        for subcom in com:
            if type(subcom) == list:
                ch = subcom[0]
                print('running subcommand:',subcom,options[ch], 'on category',cat)
                dfs,nn = functions[ch](dfs,subcom,cat)
        return(dfs,[])
    
    
    
    
def obMenu(obs,title="choose category:"):
    for i,col in enumerate(obs.columns):
        print(i,col)   
    ch = int(input(title))
    uch = obs[obs.columns[ch]].unique()
    return(ch,uch)


'''
menus  
'''
def mainMenu(dfs,com=[],cat=''):
    print('main menu')
    op = ['general data handling','scaling','clustering','batch-correction',
          'celltyping','neighborhood analysis']
    fn = [selection,scaling,clustering,batchCorrection,celltyping,neighborhoodAnalysis]
    dfs,coms=menu(dfs,op,fn,com,cat)
    #print(coms,'coms out from mainMenu')
    return(dfs,coms)

    


def selection(dfs,com=[],cat=''):
    print('selection (general data handling)')
    op = ['save to csv',] #'save in ram','revert to ram save','pick subset of data'
    fn  = [save,pick]
    dfs,com=menu(dfs,op,fn,com,cat)
    #print(com,'com in selection')
    return(dfs,com)

            
        
def scaling(dfs,com,cat=''):
    print("scaling")
    op = ['zscore across samples (cells)']
    fn=[zscore]
    dfs,com=menu(dfs,op,fn,com,cat)
    return(dfs,com)

def clustering():
    pass

def batchCorrection():
    pass

def celltyping():
    pass

def neighborhoodAnalysis():
    pass


'''
scaling
'''

def zscore(dfs,com=[],cat=''):
    if len(com) == 0:
        return([],[])
    dfs[0] = scipy.stats.zscore(dfs[0])
    return(dfs,[])


'''
selection (general)
'''


def save(dfs,com=[],cat=''):
    print('saving')
    if len(com) == 0:
        filename = input("filename?")
        return([],[filename])
    df,obs,dfxy = dfs[0],dfs[1],dfs[2]
    filename = com[1]+'_'+cat
    print(filename)
    df.to_csv(filename+"_df.csv")
    obs.to_csv(filename+"_obs.csv")
    dfxy.to_csv(filename+"_dfxy.csv")
    return(dfs,[])

def pick(dfs,com,cat=''):
    print('picking subset')








if __name__ == "__main__":
    #'''#"zzz_hta14_tumorneighborhoodcts1"
    folder = ""#r"C:\Users\youm\Desktop\src\zzzzzzzzzzz_current/"
    stem = 'cl56_depth_study_H12'#'89_LC-4_withN'#''96_LC'#'96_LC'#'97_mtma2'#'93_hta14'###'96_hta14_primary'#'97_hta14bx1_primary_celltype'#'99_hta14'#"temp"#"zzz_hta1499"#"zzz14bx1_97"#"hta14bx1 dgram"#folder+"14_both"##"tempHta14_200"#"HTA14f"#"zzzz_hta1498_neighborhoodsOnly"#"hta1415Baf1"#"HTA15f"#"0086 HTA14+15"#"99HTA14"#"z99_ROIs_5bx_HTA1415"#"temp"#"z99_ROIs_5bx_HTA1415"#<-this one has old celltyping no TN #"0084 HTA14+15" #"HTA9-14Bx1-7 only"#"0.93 TNP-TMA-28"#"0.94.2 TNP-TMA-28 primaries"#"1111 96 TNP-28" #'0093 HTA14+15'#"0094.7 manthreshsub primaries HTA14+15"#"0094 HTA14+15" #"096 2021-11-21 px only" #'095.08 primaries only manthreshsub 2021-11-21 px only'#"094 manthreshsub 2021-11-21 px only" #  '095.1 primaries only manthreshsub 2021-11-21 px only'#
    print("axis labels %s on barplots more ticks/lines")
    print(stem)
    df = pd.read_csv(stem+"_df.csv",index_col=0) 
    obs = pd.read_csv(stem+"_obs.csv",index_col=0).astype(str)
    dfxy = pd.read_csv(stem+"_dfxy.csv",index_col=0)
    print(df.shape[0],"cells")
    main(df,obs,dfxy)