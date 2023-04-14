# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 18:08:37 2021

@author: youm
"""


import warnings
warnings.simplefilter(action='ignore')#, category=FutureWarning)

import numpy as np
import pandas as pd
import time
import os
import matplotlib
import matplotlib.pyplot as plt
import copy
import scanpy as sc
import anndata
import math
import seaborn as sns
#import phenograph  #problem with igraph again
#from scipy import sparse
#from sklearn.metrics import adjusted_rand_score
#import sklearn as skl
import re
import scipy as sp
import statistics as stat
import random
#import bokehClusterMap1 as bcm
import allcolors as allc
import sys
import orthogonal7 as ort
import combat1 as combat1
#import recropTma as rec
#import lithresh1 as lithresh
#import RESTORE as RES
#import IFanalysisPackage0 as IF
from skimage import io
import napari7 as NP
import napari as NAPARI
from sklearn.mixture import BayesianGaussianMixture as GMM
import matplotlib.style
import matplotlib as mpl
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
import skimage
import PIL
import tifffile

mpl.style.use('default')

ODF = 9
OOBS = 9
OXY = 9
COMMANDS = []





def main(df,obs,dfxy):
    global ODF
    global OOBS
    global OXY
    global COMMANDS
    if isinstance(ODF,int):
        ODF = copy.deepcopy(df)
        OOBS = copy.deepcopy(obs)
        OXY = copy.deepcopy(dfxy)
    op = ["processing","visualization"]#,"aggregate","split"]
    fn = [processing,visu]
    df,obs,dfxy=menu(op,fn,df,obs,dfxy)
    return(df,obs,dfxy)




def menu(options,functions,df=9,obs=9,dfxy=9,cdf=9):
    while True:
        #try:
            print("\n")
            for i,op in enumerate(options):
                print(i,op)
            try:
                print("send non-int when done (return df)")
                ch = int(input("number: "))
            except:
                return(df,obs,dfxy)
            if isinstance(cdf,int):
                df,obs,dfxy=functions[ch](df,obs,dfxy) 
                print(all(obs.index==df.index),"all index the same")
            else:
                df,obs,dfxy=functions[ch](df,obs,dfxy,cdf)
                print(all(obs.index==df.index),"all index the same")




def processing(df,obs,dfxy):
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
    df,obs,dfxy = menu(op,fn,df,obs,dfxy)
    return(df,obs,dfxy)
    
    
def visu(df,obs,dfxy):
    mpl.style.use('default')
    if input("recluster obs? (y)") == "y":
        df,obs,dfxy=obCluster(df,obs,dfxy)
    op = ["heatmap","cluster-bar plot","bar plot","box plot","correlation matrix","scatterplot","umap","show spatial",
          "scanpy visuals","built in pie chart","silhouette score for clusters","napari visualization","sorted single-cell heatmap","save obs columns as annotation png set"]#,"send cluster df to processing"]
    fn = [heatmap,clusterbar,barplot,boxplot,correlationMatrix,scatterplot,showUmap,spatial,scanpyv,neighPie,silh,nap,sortedMap,
          savepngs]#,cdfToProcessing]
    clusterA,ucl,obs = clusterMeans(df,obs,dfxy)
    cdf = pd.DataFrame(clusterA,index=ucl,columns = df.columns)
    df,obs,dfxy = menu(op,fn,df,obs,dfxy,cdf = cdf)
    return(df,obs,dfxy)



def countObs(df,obs,dfxy,cdf=None):
    if input("this will overwrite your working dataframe and will make some functions crash. Continue? (y)") != "y":
        return(df,obs,dfxy)
    obs = obs.astype(str)
    ch,uch = obMenu(obs,title="sort by:")
    
    cols = []
    for ob in obs.columns:
        uos = obs[ob].unique()
        if len(uos) > 30 or len(uos) < 2:
            continue
        for uo in uos:
            cols.append(str(ob)+"   "+str(uo))
    A = np.zeros((len(uch),len(cols)))
    #O = np.zeros((len(uch),len(cols)),dtype=object)
    for i,cat in enumerate(uch):
        sobs = obs.loc[obs[obs.columns[ch]] == cat,:]
        for j,c in enumerate(cols):
            clis = c.split("   ")
            ob = clis[0]
            uo = clis[1]
            key = sobs[ob] == uo
            A[i,j] = key.sum()/sobs.shape[0]
    print(A)       
    ndf = pd.DataFrame(index=uch,columns=cols,data=A)
    obs = pd.DataFrame(index=uch,columns=[obs.columns[ch]],data=uch)
    xy = pd.DataFrame(index=uch,columns=["x","y"],data=np.zeros((len(uch),2)))
    return(ndf,obs,xy)
            
            

def aggregate(df,obs,dfxy,cdf=None): 
    for i,ob in enumerate(obs.columns):
        print(i,ob)
    print("aggregate by which?")
    ch1 = int(input("enter number:"))
    ch2 = int(input("0:average biomarker values or 1:count fraction positive"))
    color = obs.columns[ch1]
    bins = np.unique(np.array(sorted(obs[obs.columns[ch1]]))) 
    newA = []
    newObs = []
    newXY = []
    for bi in bins:
        #print(bi)
        arow = []
        orow = []
        xrow = []
        key = obs[color] == bi
        sdf = df.loc[key,:]
        sobs = obs.loc[key,:]
        sxy = dfxy.loc[key,:]
        if ch2 == 1:
            for col in sdf.columns:           
                tdf = sdf.loc[sdf[col]>0,:]
                arow.append(tdf.shape[0]/sdf.shape[0])
        else:
            arow = list(sdf.mean(axis=0))
        for col in sobs.columns:
            orow.append(sp.stats.mode(list(sobs[col]))[0][0])
        for col in sxy.columns:
            xrow.append(np.mean(list(sxy[col])))
            #xrow.append(float(sxy[col].mean()))
        newA.append(arow)
        newObs.append(orow)
        newXY.append(xrow)
    newDF = pd.DataFrame(data=newA, columns=df.columns,index=bins)
    newObs = pd.DataFrame(data=newObs,columns = obs.columns,index=bins)
    newXY = pd.DataFrame(data=newXY,columns=dfxy.columns,index=bins)
    print(newObs)
    return(newDF,newObs,newXY)

def clusterMeans(DF,obs,dfxy):
    clusterA = []
    try:
        ucl = np.unique(obs.loc[:,'Cluster'].values)
    except:
        DF,obs,dfxy=obCluster(DF, obs, dfxy)
        ucl = np.unique(obs.loc[:,'Cluster'].values)
    nClusters = len(ucl)
    clusterA = np.zeros((nClusters,DF.shape[1]))
    print("ucl",ucl)
    for i,c in enumerate(ucl):
        cl = DF.loc[obs['Cluster'] == c,:]
        try:
            markerMeans = np.mean(cl.values,axis = 0)
        except:
            markerMeans = np.mean(cl.values.astype(float),axis = 0)
        clusterA[i,:] = markerMeans
    #print(clusterA,"clusterA")
    return(clusterA,ucl,obs)

def silh(df,obs,dfxy,cdf):
    
    obs = obs.fillna("nan")
    cols = []
    scores = []
    for col in sorted(obs.columns):
        ncat = len(obs.loc[:,col].unique())
        if ncat < 2 or ncat > 20:
            continue
        try:
            sco = silhouette_score(df,obs.loc[:,col])
            silhouette_vals = silhouette_samples(df,obs.loc[:,col])
            silhH1(silhouette_vals,obs.loc[:,col])
            print(col,sco)
            cols.append(col)
            scores.append(sco)
        except Exception as e:
            print("could not score",col,"    ",e)
    plt.scatter(cols,scores)
    plt.xticks(rotation = 90)
    plt.show()
    return(df,obs,dfxy)

def silhH1(silhouette_vals,labels):
    fig,ax1 = plt.subplots()
    y_ticks = []
    y_lower, y_upper = 0, 0
    for i, cluster in enumerate(np.unique(labels)):
        cluster_silhouette_vals = silhouette_vals[labels == cluster]
        cluster_silhouette_vals.sort()
        y_upper += len(cluster_silhouette_vals)
        ax1.barh(range(y_lower, y_upper), cluster_silhouette_vals, edgecolor='none', height=1)
        ax1.text(-0.03, (y_lower + y_upper) / 2, str(i + 1))
        y_lower += len(cluster_silhouette_vals)
    plt.show()

'''
PROCESSING
'''

#def cdfToProcessing(df,obs,dfxy,cdf):
    #cobs = 
'''
Endothelial: CD31+ and others.
From others in 1), Tumor: any of Ecad+, CK19+, CK7+, CK8+, CK5+, CK14+, CK17+, MUC1+ and others.
From others in 2), Immune: any of CD20+, CD45+, CD4+, CD8+, CD68+, CD3+ and others.
From others in 3), Active fibroblast: any of aSMA+, VIM+, CAV1 and others.
From others in 4) is Inactive fibroblast.
Tumor in 2), tumor proliferative: Ki67+, PCNA+, pHH3+.
Tumor in 2), luminal: CK19+, CK7+, CK8+, but none of CK5+, CK14+, CK17+.
Tumor in 2), luminal-basal: any of Ecad+, CK19+, CK7+, CK8+, and any of CK5+, CK14+, CK17+.
Tumor in 2), luminal-basal: any of CK5+, CK14+, CK17+, but none of CK19+, CK7+, CK8+.
Tumor in 2), HR+ tumor: any of ER+, PgR+.
Tumor in 2), HER2+ tumor: any of HER2+.
Tumor in 2), none of HER2+, ER+, PgR+, so exclusive 10) and 11).
Immune in 3), immune proliferative: Ki67+, PCNA+, pHH3+.
Immune in 3), T-cell: CD4+, CD8+, CD3+.
Immune in 3), B-cell: CD20+.
Immune in 3), macrophage: CD68+, then others are “other immune”.
'''

#calculate the average distance to 25 nearest neighbors
#for each biomarker, find fraction of X+ cells nearby (X being every other biomarker)
#implement chi-squared testing function


def subset(df,obs,dfxy,n=False):
    if type(n) != int:
        n = int(input("how many cells?"))
    inds = []
    sh = df.shape[0]
    while len(inds) <= n:
        ind = random.randint(0,sh-1)
        if ind not in inds:
            inds.append(ind)
    df= df.iloc[inds,:]
    obs = obs.iloc[inds,:]
    dfxy=dfxy.iloc[inds,:]
    return(df,obs,dfxy)

def simulateTherapy(df,obs,dfxy):    
    maxCombs = 2
    ch,uch = obMenu(obs,title="celltype selector")
    ob = obs.columns[ch]
    print(uch)
    celltype = []
    for i,uc in enumerate(uch):
        print(i,uc)
    ch = int(input("number of celltype to try to eliminate:"))
    celltype.append(uch[ch])
    print(celltype,"type to compare")

    
    zdfq = input("zscore dataframe? (y)")
    thresholdSD = float(input("threshold (in SD if zscored)"))

    #get only non-primary DF (substitute later for non-essential-to-life df): BROKEN?
    pDf,n,nn = onlyPrimaries(df,obs,dfxy)
    nP = []
    for biom in df.columns:
        if biom not in pDf.columns:
            nP.append(biom)
    print(nP,"non primary")
    npDf = df.loc[:,nP]


    #threshA = np.where(npDf.values>thresholdSD,0,1)
    
    
    #find every combination of biomarkers up to maxCombs
    combinations = []
    for i in range(maxCombs):
        if len(combinations) == 0:
            combinations = list(npDf.columns)
        else:        #len(combinations) < maxCombs*len(bioms):
            for comb in combinations:
                if len(comb.split("#")) == maxCombs:
                    continue
                for biom in npDf.columns:
                    if biom not in comb:
                        combinations.append(comb+"#"+biom)

    #make threshold matrix matching dataframe
    if zdfq == "y":
        npDf,n,nn = zscore(npDf,obs,dfxy,ax=0)
    npDf= pd.DataFrame(np.where(npDf.values>thresholdSD,0,1),index=npDf.index,columns=npDf.columns)   
    keptFrames = []
    scores = []
    for comb in combinations:
        frames = [df.copy(),obs.copy(),dfxy.copy()]
        for frame in frames:
            frame.index = pd.Series(frame.index) +"_"+comb
        #print("\n\n",comb)
        key = None
        for biom in comb.split("#"):
            try:
                key = key * npDf.loc[:,biom].values
            except:
                #print("starting key")
                key = npDf.loc[:,biom].values
                
        #print(key)#),~key <-?)
        key = pd.Series(key)==1
        #key.index= npDf.index
        #print(key,key.shape[0])
        #print(df.index,key.index)
        for i,frame in enumerate(frames):
            #print(frame.shape[0],key.shape[0],"frame,key")
            #print(frame.index==key.index)  SUPER WEIRD THE STRS ARE THE SAME BUT  == is [False False False False]
            key.index=frame.index
            frames[i] = frame.loc[key,:]
        totalRem = frames[1].shape[0]/obs.shape[0]*100
        cellt = celltype[0]
        tumKey = obs.loc[:,ob] == cellt
        tumO = obs.loc[tumKey,:]
        otherO = obs.loc[~tumKey,:]
        thTumKey = frames[1].loc[:,ob] == cellt
        thTumRem = frames[1].loc[thTumKey,:]
        thOtherRem = frames[1].loc[~thTumKey,:]
        typeRem = thTumRem.shape[0]/tumO.shape[0]*100
        otherRem = thOtherRem.shape[0]/otherO.shape[0]*100
        #score = (100-typeRem)**2/(100-totalRem)*totalRem**2 #PDPN
        #score = (100-typeRem)/(100-totalRem) #gH2AX
        score = otherRem/typeRem
        scores.append(score)
    Mind = scores.index(max(scores))
    bestComb = combinations[Mind]
    print(bestComb,"Best combination of biomarkers to remove, score:",max(scores))
    for i,biom in enumerate(bestComb.split("#")):
        if i==0:
            key1 = npDf.loc[:,biom].values
        else:
            key1 = key1 * npDf.loc[:,biom].values

            
    key1 = pd.Series(key1)==1
    #print(key,key1,"key key1")
    key1.index=df.index
    try:
        obs["slide"]+="_"+bestComb
        obs["slide_scene"]+="_"+bestComb
    except IndexError:
        print("slide and'or slide_scene missing from obs annotations")
    return(df.loc[key1,:],obs.loc[key1,:],dfxy.loc[key1,:])
    
    
            
        
        
        
def neighborhoodEntropy(df,obs,dfxy):
    
    #calculate sample entropy based on inputs, print sample entropy, ask user if binarize into high-low based on sample total
    radii = []
    while True:
        try:    
            radii.append(float(input("radius (in um) to consider neighbors: (send blank when done)"))/.325)
        except:
            if len(radii) < 1:
                print("invalid radius, using 25 as default")    
                radii = [25/.325]
            break
    ch,uch = obMenu(obs,title="category to add neighborhood fractions to dataframe")
    print(uch)
    goodsts = []
    while True:
        inp = input("string to include in celltypes to consider as centers (blank to include all, blank when done. Leaves NAN values if not all): ")
        if inp == "":
            break
        else:
            goodsts.append(inp)
    #much = len(uch)
    obcol = obs.columns[ch]
    if len(goodsts) == 0:
        goodsts = uch
    for radius in radii:
        ecol = obcol+"_entropy_"+str(radius*.325)
        df[ecol] = 0
        for us in obs["slide_scene"].unique():
            key0 = obs["slide_scene"] == us
            tdfxy = dfxy.loc[key0,:]
            tobs = obs.loc[key0,:]
            for i in range(tdfxy.shape[0]):
                ind = tdfxy.index[i]
                if len(goodsts) > 0:
                    check = False
                    for gs in goodsts:
                        if gs in tobs.loc[ind,obcol]:
                            check = True
                if not check:
                    continue
                    
                if i % 1000 == 1:
                    print(i/tdfxy.shape[0]*100,"% done with",us)
                #neighbors = []
                x,y = tdfxy.iloc[i,0],tdfxy.iloc[i,1]
                nx,ny = tdfxy.iloc[:,0],tdfxy.iloc[:,1]
                distanceV = ((x-nx)**2+(y-ny)**2)**.5
                key = distanceV < radius 
                neighbors = tobs.loc[key,:]
                neighbors = neighbors.drop(pd.Series(tdfxy.index).iloc[i])
                nnei = neighbors.shape[0]
                entropy = 0
                if nnei > 1:   
                    for uc in uch:
                        inNei = neighbors.loc[neighbors.loc[:,obcol] == uc,obcol]
                        Pi = inNei.shape[0]/nnei
                        if Pi != 0:
                            entropy -= Pi * math.log2(Pi)

                        #df.loc[ind,uc+"_"+obcol+"_neighbors_"+str(radius)] = 1/much #add equal fraction of each category
                df.loc[ind,ecol] = entropy
                
    return(df,obs,dfxy)
        



def neighborhoodFractions(df,obs,dfxy):
    tot = False
    if input("count totals instead of fractions? (y)") == "y":
        tot = True
    radii = []
    while True:
        try:    
            radii.append(float(input("radius (in um) to consider neighbors: (send blank when done)"))/.325)
        except:
            if len(radii) < 1:
                print("invalid radius, using 25 as default")    
                radii = [25/.325]
            break
    ch,uch = obMenu(obs,title="category to add neighborhood fractions to dataframe")
    print(uch)
    goodsts = []
    while True:
        inp = input("string to include in celltypes to consider as centers (blank to include all, blank when done. Leaves NAN values if not all): ")
        if inp == "":
            break
        else:
            goodsts.append(inp)
    #much = len(uch)
    obcol = obs.columns[ch]
    if len(goodsts) == 0:
        goodsts = uch
    for radius in radii:
        for uc in uch:
            df[uc+"_"+obcol+"_neighbors_"+str(radius*.325)] = 0
        for us in obs["slide_scene"].unique():
            key0 = obs["slide_scene"] == us
            tdfxy = dfxy.loc[key0,:]
            tobs = obs.loc[key0,:]
            for i in range(tdfxy.shape[0]):
                ind = tdfxy.index[i]
                if len(goodsts) > 0:
                    check = False
                    for gs in goodsts:
                        if gs in tobs.loc[ind,obcol]:
                            check = True
                if not check:
                    continue
                    
                if i % 1000 == 1:
                    print(i/tdfxy.shape[0]*100,"% done with",us)
                #neighbors = []
                x,y = tdfxy.iloc[i,0],tdfxy.iloc[i,1]
                nx,ny = tdfxy.iloc[:,0],tdfxy.iloc[:,1]
                distanceV = ((x-nx)**2+(y-ny)**2)**.5
                key = distanceV < radius 
                neighbors = tobs.loc[key,:]
                neighbors = neighbors.drop(pd.Series(tdfxy.index).iloc[i])
                nnei = neighbors.shape[0]
                if nnei > 1:
                    #nuch = neighbors.loc[:,obcol]
                    for uc in uch:
                        inNei = neighbors.loc[neighbors.loc[:,obcol] == uc,obcol]
                        if tot:
                            df.loc[ind,uc+"_"+obcol+"_neighbors_"+str(radius*.325)] = inNei.shape[0]
                        else:
                            df.loc[ind,uc+"_"+obcol+"_neighbors_"+str(radius*.325)] = inNei.shape[0]/nnei                
                else:
                    for uc in uch:
                        df.loc[ind,uc+"_"+obcol+"_neighbors_"+str(radius*.325)] = 0
                        #df.loc[ind,uc+"_"+obcol+"_neighbors_"+str(radius)] = 1/much #add equal fraction of each category
    return(df,obs,dfxy)


        
        

def regionAverage(df,obs,dfxy): #make only check same slidescene
    try:    
        radius = float(input("radius (um) to consider in average: "))/.325
    except:
        print("invalid radius, using 25u as default")
        radius = 77
    ndf = []
    for us in obs["slide_scene"].unique():
        key0 = obs["slide_scene"] == us
        tdfxy = dfxy.loc[key0,:]
        tobs = obs.loc[key0,:]  
        tdf = df.loc[key0,:]
        for i in range(tdfxy.shape[0]):
            if i % 1000 == 1:
                print(i/tdfxy.shape[0]*100,"% done with",us)
            neighbors = []
            x,y = tdfxy.iloc[i,0],tdfxy.iloc[i,1]
            nx,ny = tdfxy.iloc[:,0],tdfxy.iloc[:,1]
            distanceV = ((x-nx)**2+(y-ny)**2)**.5
            key = distanceV < radius 
            neighbors = tdf.loc[key,:]
            neighbors = neighbors.drop(pd.Series(tdfxy.index).iloc[i])
            if neighbors.shape[0]> 1:
                #print(neighbors)
                #neighbors = pd.DataFrame(neighbors)#,columns = pd.Series(df.columns)+" in range "+str(radius))
                #print(neighbors)
                avg = neighbors.mean(axis=0)
                #print(avg,type(avg))
                #print(avg.values)
                ndf.append(pd.DataFrame(avg.values,index=df.columns,columns = [pd.Series(tdf.index).iloc[i]]).transpose())
            else:
                ndf.append(pd.DataFrame(columns =df.columns ,index=[pd.Series(tdf.index).iloc[i]] ))#
        #print(ndf)
    ndf = pd.concat(ndf,axis=0)
    print("ndf start",ndf,"ndf end")
    avgs = ndf.mean(axis=0)
    for biom in ndf.columns:
        print(biom,"biom")
        ndf.loc[:,biom]=ndf.loc[:,biom].fillna(0)#avgs.loc[biom])
    ndf = ndf.set_axis(pd.Series(ndf.columns)+" in radius "+str(radius*.325),axis=1)
    df= pd.concat([df,ndf],axis=1)
    if input("sort df? (y)") == "y":
        df = df.loc[:,df.columns.sort_values()]
    return(df,obs,dfxy)


def regionAverage1(df,obs,dfxy): #make only check same slidescene
    radius = float(input("radius to consider in average: "))
    ndf = []
    for i in range(dfxy.shape[0]):
        neighbors = []
        x,y = dfxy.iloc[i,0],dfxy.iloc[i,1]
        for j in range(dfxy.shape[0]):
            if i == j:
                continue
            nx,ny = dfxy.iloc[j,0],dfxy.iloc[j,1]
            distance = ((x-nx)**2+(y-ny)**2)**.5
            if distance < radius:
                neighbors.append(df.iloc[j,:])
        if len(neighbors) > 1:
            print(neighbors)
            neighbors = pd.DataFrame(neighbors)#,columns = pd.Series(df.columns)+" in range "+str(radius))
            print(neighbors)
            avg = neighbors.mean(axis=0)
            print(avg,type(avg))
            ndf.append(avg)
        else:
            ndf.append(pd.DataFrame(columns =df.columns ,index=[pd.Series(df.index).iloc[i]] ))#
    print(ndf)
    ndf = pd.DataFrame(ndf)
    avgs = ndf.mean(axis=0)
    for biom in ndf.columns:
        ndf.loc[:,biom].fillna(avgs.loc[biom])
    ndf = ndf.set_axis(pd.Series(ndf.columns)+" in radius "+str(radius))
    df= pd.concat([df,ndf],axis=1)
    return(df,obs,dfxy)
                            


def regionAverage1(df,obs,dfxy): #make only check same slidescene
    radius = float(input("radius to consider in average: "))
    ndf = []
    for i in range(dfxy.shape[0]):
        neighbors = []
        x,y = dfxy.iloc[i,0],dfxy.iloc[i,1]
        for j in range(dfxy.shape[0]):
            if i == j:
                continue
            nx,ny = dfxy.iloc[i,0],dfxy.iloc[i,1]
            distance = ((x-nx)**2+(y-ny)**2)**.5
            if distance < radius:
                neighbors.append(df.iloc[j,:])
        if len(neighbors) > 1:
            print(neighbors)
            neighbors = pd.DataFrame(neighbors)#,columns = pd.Series(df.columns)+" in range "+str(radius))
            print(neighbors)
            avg = neighbors.mean(axis=0)
            print(avg,type(avg))
            ndf.append(avg)
        else:
            ndf.append(pd.DataFrame(columns =df.columns ,index=[pd.Series(df.index).iloc[i]] ))#
    print(ndf)
    ndf = pd.DataFrame(ndf)
    avgs = ndf.mean(axis=0)
    for biom in ndf.columns:
        ndf.loc[:,biom].fillna(avgs.loc[biom])
    ndf = ndf.set_axis(pd.Series(ndf.columns)+" in radius "+str(radius))
    df= pd.concat([df,ndf],axis=1)
    return(df,obs,dfxy)
                                     
            
            
  
def autotype(df,obs,dfxy,chanT=True,name="autoCellType res: ",res=None): #the old version that keeps more information is in cmifAnalysis36
    roundThresh = [1500,1250,1000,750]
    biomRounds = [['CAV1', 'CK17', 'CK5', 'CK7', 'CK8', 'H3K27', 'MUC1', 'PCNA', 'R0c2', 'R6Qc2', 'Vim', 'aSMA', 'pHH3'],
                  ['AR', 'CCND1', 'CD68', 'CD8', 'CK14', 'CoxIV', 'EGFR', 'H3K4', 'HER2', 'PDPN', 'R0c3', 'R6Qc3', 'pS6RP'],
                  ['BCL2', 'CD31', 'CD4', 'CD45', 'ColIV', 'ER', 'Ki67', 'PD1', 'PgR', 'R0c4', 'R6Qc4', 'gH2AX', 'pRB'],
                  ['CD20', 'CD3', 'CD44', 'CK19', 'CSF1R', 'ColI', 'Ecad', 'FoxP3', 'GRNZB', 'LamAC', 'R0c5', 'R6Qc5', 'RAD51']]
  
    odf = df.copy()
    if not res:
        res= float(input("number of standard deviations above mean required to count as +"))
    #chanT = False
    if chanT != False:
        if input("check if threshold above channel threshold (non-z-score) (y)") == 'y':
            chanT = True
            key2 = pd.DataFrame(data=np.zeros_like(odf),columns=odf.columns,index=odf.index)
            means = df.mean(axis=0)
            sds = df.std(axis=0)
            zSer = pd.Series(index = df.columns,data=means+sds*res)
            print(zSer)
            for i,roun in enumerate(biomRounds):
                rawThresh = roundThresh[i]
                for bIm in roun:
                    for bim in df.columns:
                        if bIm in bim:
                            key2.loc[odf.loc[:,bim]>rawThresh,bim] = 1
                            #bimT = zSer.loc[bim]
                            #key2.loc[odf.loc[:,bim]>bimT,bim] = 1
                            #print("ding")
        else:
            chanT = False
                
    #print(list(key2.iloc[:,0]),"ar key")        
    obs[name+str(res)] = " "   
    if input("zscore?") == "y":
    #print(df)
        df,obs,dfxy = zscorev(df,obs,dfxy)  
    #print(df)  
    mapp = {}
    toThresh = []
    for biom in df.columns:
        if "neigh" in biom:
            continue
        cType = fillMap(biom)
        if cType != None:
            mapp[biom]=cType
    toThresh = list(mapp.keys())
    others = ["Ki67", "PCNA", "pHH3","pRB","ER","PgR","AR","HER2","Fox","GRNZB","aSMA","Vim","VIM","ColI","PD1"] #CAV
    for biom in df.columns:
        for o in others:
            if o in biom:
                toThresh.append(biom)
    for biom in toThresh:
        if chanT:
            key3 = key2[biom] == 1
            key = df[biom]>res  
            print(biom,any(key),any(key3),any(key & key3))
            obs.loc[key & key3,name+str(res)] += biom + " "
        else:
            key = df[biom]>res
            obs.loc[key,name+str(res)] += biom + " "
    obs = parseTypes(df,obs,dfxy,column=name+str(res))
    print(obs[name+str(res)].unique(),"uobs")
    obs = parseSecondary(df,obs,dfxy,column=name+str(res))
    #if input("keep z-scoring?") == "y":
        #return(df,obs,dfxy)
    return(odf,obs,dfxy)

def manThresh(df,obs,dfxy):
    ch = input("import manual thresholds from csv? (y)")
    if ch == "y":
        obs["Manual Celltype"] = ""
        chh = input("subtract thresholds? (y)")
        path = input('filepath/name:')
        thresh = pd.read_csv(path)
        for biom in thresh.columns:
            th = float(thresh[biom].iloc[0])
            print(biom,th)
            key = df[biom] > th
            obs.loc[key,"Manual Celltype"] += biom
            if chh == "y":
                df[biom] -= th
                df.loc[df[biom]<0,biom] = 0        
    else:
        if input("manually add annoations for every combination of positivies for custom biomarker set? (y)") == 'y':
            df,obs=ort.tirtiary(df,obs)
            return(df,obs,dfxy)
        chhh = int(input("0:chose obs category or 1:all together"))
        if chhh == 0:
            df,obs = ort.secondary(df,obs)
        else:
            df,obs = ort.primary(df,obs) 
    obs = parseTypes(df,obs,dfxy,column="Manual Celltype")
    obs = parseSecondary(df,obs,dfxy,column="Manual Celltype")
    return(df,obs,dfxy)


def clort(df,obs,dfxy):
    pass

def clauto(df,obs,dfxy):
    obs = obs.astype(str)
    print(obs.shape)
    oobs = obs.copy()
    chs, uchs = [],[]
    while True:
        try:
            ch,uch=obMenu(obs,"obs category to auto-annotate cell types")
            chs.append(ch)
            uchs.append(uch)
        except:
            break
    for i,ch in enumerate(chs):
        uch = uchs[i]
        adf,aobs,axy = clag(df,obs,dfxy,ch,uch)
        x,aobs,xx = autotype(adf,aobs,axy,chanT=False,name=obs.columns[ch]+"cluster autotype",res=float(input('resolution:')))
        print(obs.shape)
        for col in aobs.columns:
            if obs.columns[ch]+"cluster autotype" in col:
                print(col,aobs.loc[:,col].unique(),"!!")
                obs[col] = ""
                for uc in aobs.index:
                    key = obs.iloc[:,ch] == uc
                    obs.loc[key,col] = aobs.loc[uc,col]
        print(obs.shape,df.shape)
        #print(obs,df)
    return(df,obs,dfxy)
                



def clag(df,obs,dfxy,ch=None,uch=None,z=True):
    if not ch:
        ch,uch=obMenu(obs,"obs category to auto-annotate cell types")
    
    if z:
        zdf,zobs,zxy = zscorev(df,obs,dfxy)
    else:
        zdf,zobs,zxy = df,obs,dfxy
    ocol = obs.columns[ch]
    ndf,nobs,nxy = [],[],[]
    for uc in uch:
        key = zobs.loc[:,ocol] == uc
        sdf = zdf.loc[key,:]
        sobs = zobs.loc[key,:]
        sxy = zxy.loc[key,:]
        ndf.append(sdf.mean(axis=0))
        nxy.append(sxy.mean(axis=0))
        #print(sobs.mode(axis=0).iloc[0,:],"/n/n")
        #time.sleep(1)
        nobs.append(sobs.mode(axis=0).iloc[0,:])
    dfs = [ndf,nobs,nxy]
    for i,d in enumerate(dfs):
        dfs[i] =pd.concat(d,axis=1).transpose()
        dfs[i].index = uch.astype(str)
        #print(dfs[i].columns)
        #print(dfs[i].shape)
        #print(dfs[i])         
    return(dfs[0],dfs[1],dfs[2])
    
        
        
'''        
obs.append(getObModes(sobs))     
def getObModes(obs):
    modes = [] #1 per col
    for col in obs.columns:
        mod = obs.loc[:,col].mode(axis=0)
        key = obs.loc[:,col]==mode
        frac = key.sum()/obs.shape[0]
        modes.append(str(mod)+" "+str(frac*100)+"%")
'''        

def parseTypes(df,obs,dfxy,column="none"):
    if column == 'none':
        ch,uch = obMenu("column to apply types to")
        column = obs.columns[ch]
    mapp = {}
    for biom in df.columns:
        cType = fillMap(biom)
        if cType != None:
            mapp[biom]=cType
    #print(mapp)
    mapp = {k: v for k, v in sorted(mapp.items(), key=lambda item: item[1])}
    print(mapp)
    obs["Primary Celltype "+column] = "5 stromal"
    for biom in mapp.keys():
        keyCol = obs[column].str.contains(biom)
        #print(list(keyCol))
        unasKey = obs["Primary Celltype "+column] == "5 stromal"
        obs.loc[keyCol & unasKey,"Primary Celltype "+column] = mapp[biom]    
    return(obs)


def parseSecondary(df,obs,dfxy,column):
    uch = obs[column].unique()
    #print(uch,"uch")
    #obs["Secondary Celltype"+column] = obs["Primary Celltype "+column] #cmif39 has this version
    obs["proliferating "+column] = "no"
    obs["tumor subtype "+column] = np.nan
    obs["receptors "+column] = np.nan
    obs["immune subtype "+column] = np.nan
    obs["immune checkpoints "+column] = np.nan
    obs["cytotoxic "+column] = np.nan
    obs["fibroblast type "+column] = np.nan
    
    
    proL = ["Ki67","PCNA","pHH3","pRB"]#
    lumL = ["CK19","CK7","CK8"]
    basL = ["CK5","CK14","CK17"]
    mesL = ["Vim","VIM","CD44"] #ANY MES MEANS NOT LUM BAS ETC.
    TL4 = ["CD4_"]
    TL8 = ["CD8"]
    TL3 = ["CD3_"] #LOW PRIORITY
    #if all 3 positive, call CD8, otherwise call CD8 CD4 'other T cell' for CD3_ cd4-cd8-
    BL = ["CD20"]
    macL = ["CD68"] #ADD CSF1R?
    Hl = ['ER', 'PgR', 'AR']
    HEl = ["HER2"]
    cpL = ["PD1","Fox"]
    cytL = ["GRNZB"]
    acL = ["aSMA","Vim","VIM","ColI_"]
    for typ in uch:
        key = obs[column] == typ
        typeD = {'pro':0,'Lum':0,'Bas':0,"HR":0,"HER":0,'T-c4':0,'T-c8':0,'T-c3':0,
                 'B-c':0,'Mac':0,"CheckP":0,"CytoT":0,"activeFB":0,"mesen":0}
        if checkL(typ,proL):
            typeD['pro'] = 1        
        if checkL(typ,mesL):
            typeD["mesen"] = 1
        if typeD["mesen"] == 0:
            if checkL(typ,lumL):
                typeD["Lum"] = 1
            if checkL(typ,basL):
                typeD["Bas"] = 1
        if checkL(typ,Hl):
            typeD['HR'] = 1
        if checkL(typ,HEl):
            typeD['HER'] = 1
        if checkL(typ,cpL):
            typeD["CheckP"] = 1
        if checkL(typ,cytL):
            typeD["CytoT"] = 1
        if checkL(typ,TL4):
            typeD['T-c4'] = 1
        elif checkL(typ,TL8):
            typeD['T-c8'] = 1
        elif checkL(typ,TL3):
            typeD['T-c3'] = 1
        elif checkL(typ,BL):
            typeD['B-c'] = 1
        elif checkL(typ,macL):
            typeD['Mac'] = 1
        #if checkL(typ,acL):
            #typeD["activeFB"] = 1
        #print(typ,typeD)
        #print(typ,typeD)
        recSwitch = 0
        tuseSwitch = 0
        for sty in typeD.keys():
            if typeD[sty] == 1: 
                if sty in "mesen":
                    pKey = obs["Primary Celltype "+column] == "3 tumor"
                    obs.loc[key & pKey,"tumor subtype "+column] = sty
                if sty in "Lum Bas":
                    #print('ding')
                    pKey = obs["Primary Celltype "+column] == "3 tumor"
                    if tuseSwitch == 0:
                        obs.loc[key & pKey,"tumor subtype "+column] = sty + " "
                        tuseSwitch = 1
                    else:
                        obs.loc[key & pKey,"tumor subtype "+column] = obs.loc[key & pKey,"tumor subtype "+column]+ sty + " "
                        print(obs.loc[key & pKey,"tumor subtype "+column])
                if sty in "HR HER":
                    pKey = obs["Primary Celltype "+column] == "3 tumor"
                    if recSwitch != 0:
                        obs.loc[key & pKey,"receptors "+column] =obs.loc[key & pKey,"receptors "+column]+" "+ sty
                        #print(any(key&pKey))
                    else:
                        obs.loc[key & pKey,"receptors "+column] = sty
                        recSwitch = 1
                        #print(any(key&pKey),"!")
                if sty in "pro":
                    obs.loc[key,"proliferating "+column] = "yes"
                if sty in "CheckP":
                    pKey = obs["Primary Celltype "+column] == "2 immune"
                    obs.loc[key & pKey,"immune checkpoints "+column] = "yes"
                if sty in "CytoT":
                    pKey = obs["Primary Celltype "+column] == "2 immune"
                    obs.loc[key & pKey,"cytotoxic "+column] = "yes"
                if sty in "T-c4 T-c8 T-c3 B-c Mac":
                    pKey = obs["Primary Celltype "+column] == "2 immune"
                    obs.loc[key & pKey,"immune subtype "+column] = sty
                #if sty in "activeFB":
                    #pKey = obs["Primary Celltype "+column] == "4 stromal"
                    #obs.loc[key & pKey,"fibroblast type "+column] = "active FB"
                #else:
                    #obs.loc[key,"Secondary Celltype"+column] +=" "+ sty
    pKey = obs["Primary Celltype "+column] == "3 tumor"
    key =pd.isna( obs["tumor subtype " + column])
    #print(key.sum(),"number of np.nan")
    obs.loc[key & pKey,"tumor subtype " + column] = "negative"
    
    #print(obs.loc[:,"receptors " + column])
    pKey = obs["Primary Celltype "+column] == "3 tumor"
    key = pd.isna(obs["receptors " + column])
    obs.loc[key & pKey,"receptors " + column] = "TN"
    #print(obs.loc[:,"receptors " + column])
    
    pKey = obs["Primary Celltype "+column] == "2 immune"
    key = pd.isna(obs["cytotoxic " + column])
    obs.loc[key & pKey,"cytotoxic " + column] = "no"

    pKey = obs["Primary Celltype "+column] == "2 immune"
    key = pd.isna(obs["immune subtype " + column])
    obs.loc[key & pKey,"immune subtype " + column] = "unclassified immune"
    
    pKey = obs["Primary Celltype "+column] == "2 immune"
    key = pd.isna(obs["immune checkpoints " + column])
    obs.loc[key & pKey,"immune checkpoints " + column] = "no"

    #pKey = obs["Primary Celltype "+column] == "4 stromal"
    #key = pd.isna(obs["fibroblast type "+column])
    #obs.loc[key & pKey,"fibroblast type " + column] = "support FB"

    return(obs)
        
        

def checkL(biomsS,lis):
    for ent in lis:
        if ent in biomsS:
            #print("donmg")
            return(True)
    return(False)
    

def fillMap(biom):
    #print("NOTE: IMMUNE CURRENTLY RESTRICTED TO CD4* CD8")
    bTypes = [["1 endothelial",["CD31"]],
              ["2 immune",["CD"]],
              ["3 tumor",["CK","Ecad","MUC1",'EGFR',"HER"]],
              ["4 active fibroblast",["aSMA","Vim","VIM","ColI_","CD90"]]] 
    for typeA in bTypes:
        for stem in typeA[-1]:
            if "CD44" in biom or "in radius" in biom or "neighbors" in biom:
                return(None)
            if stem in biom:   
                return(typeA[0])






def superBiomDF(df,obs,dfxy):
    ndf = pd.DataFrame(0,index = df.index, columns = ["1 endothelial","2 immune","3 tumor","4 active fibroblast"])
    for biom in df.columns:
        typ = fillMap(biom)
        if typ != None:
            ndf.loc[:,typ] += df.loc[:,biom]
    return(ndf,obs,dfxy)

def onlyPrimaries1(df):
    print("old function?")
    1/0
    primaries = []
    keys = ["CK","CD","Ecad","Vim","aSMA","Ki67"]#,"Col","PDPN","Fox","PD1"]
    for biom in df.columns:
        for k in keys:
            if k in biom and "CD44" not in biom:
                primaries.append(biom)
    df = df.loc[:,sorted(primaries)]
    return(df)

def onlyPrimaries(df,obs,dfxy):
    primaries = []
    df = doPart(df)
    combd = []
    for biom in df.columns:
        #print(biom)
        if "nuclei_" in biom or "cell_" in biom:
            primaries.append(biom)
        elif "comb" in biom:
            #print(biom,"comb")
            combd.append(biom.split("_")[0])
    #print(combd)
    for biom in df.columns:
        if type(fillMap(biom)) == type(None):
            continue
        cond = 0
        for stem in combd:
            if stem in biom and "comb" not in biom:
                cond = 1
        if cond == 0:
            primaries.append(biom)
    print(primaries)
    df = df.loc[:,primaries]
    return(df,obs,dfxy)
            
            

def doPart(df):   
    dShort=[]
    for col in df.columns:
        shortname = col.split("_")[0]
        if shortname not in dShort:
            dShort.append(shortname)
    for sCol in dShort:
        toComb = []
        for col in df.columns:
            if sCol == col.split("_")[0]:
                toComb.append(col)
        #print("to combine:",sCol,toComb)
        if len(toComb) > 1:
            df=combinePart(df,["_nuc","ellmem","ucadj","ytopla","erinuc","cyto"],toComb,"all",sCol)
            #df=combinePart(df,["uclei","nuc"],toComb,"nuc",sCol)
            #df=combinePart(df,["ellmem"],toComb,"cellmem",sCol)
            #df=combinePart(df,["ucadj","ytopla","erinuc","cyto"],toComb,"cyto",sCol)
    #print(df.columns[df.isna().any()])
    return(df)


def combinePart(df,partitionL,toComb,NName,sCol):
    print(toComb)
    sDF = pd.DataFrame(index=df.index)
    for biomarker in toComb:
        for name in partitionL:
            if name in biomarker:
                sDF[biomarker] = df[biomarker]
                continue
    if sDF.shape[1]>1:
        #print(sDF.columns,"sdf cols\n")
        df.loc[:,sCol+"_"+NName+"_combined"]=sDF.max(axis=1)
    else:
        pass
        #print("sDF only has one entry apparently", sDF.columns)
    return(df)
    
   
def celltype(df,obs,dfxy): #manual 1x1
    odf = df.copy()
    if input("keep only primary biomarkers? (y)") == "y":
        df,obs,dfxy = onlyPrimaries(df,obs,dfxy)
    
    if input("z-score data vertically (y)") == "y":
        df,obs,dfxy = zscorev(df,obs,dfxy)
    if input("z-score data horizontally (y)") == "y":
        df,obs,dfxy = zscoreh(df,obs,dfxy)
    dlen = df.shape[0]
    #dw = df.shape[1]
    try:
        obs.loc[:,"manual cell type"]
        ch = input("manual type exists, start over? (y)")
        if ch == "y" or ch == "Y":
            ch1 = input("keep old cell type? (y)")
            if ch1 == "y" or ch1 == "Y":
                nn = input("name for old cell type?")
                obs[nn] = obs["manual cell type"]
                obs["manual cell type"] = "unidentified"
            else:
                obs["manual cell type"] = "unidentified"                
    except:
        obs["manual cell type"] = "unidentified"
    
    
    while True:
        print("type 'quit' to exit, send blank to skip")
        ind = random.randint(0,dlen)
        if obs["manual cell type"].iloc[ind] == "unidentified":
            plt.bar(df.columns,df.iloc[ind,:].values)
            plt.xticks(rotation = 90)
            plt.grid(axis='x')
            #plt.xlabel(df.columns)
            plt.show()
            ty = input("cell type:")
            if ty == "quit":
                break
            elif ty == "":
                pass
            else:
                obs["manual cell type"].iloc[ind] = ty
    print(obs["manual cell type"].unique())
    return(odf,obs,dfxy)

    

def revert(df,obs,dfxy):
    return(ODF,OOBS,OXY)


def scale1(df,obs,dfxy):
    for col in df.columns:
        mx = max(df[col].max(),-df[col].min())
        df.loc[:,col] = df[col]/mx
    return(df,obs,dfxy)

def elmarScale(df,obs,dfxy):
    if input("jenny scale instead? div by stdev (y)") == "y":
        for col in df.columns:
            lis = list(df[col])
            std = stat.stdev(lis)
            df[col] = df[col]/std
    else:
        for col in df:
            df.loc[df[col]<1,col] = 1
        A = np.log(df.values)
        df= pd.DataFrame(A,index=df.index,columns=df.columns)
        print(df)
        df,obs,dfxy = zscore(df,obs,dfxy,ax=0)
        print(df)
        df += 3
        df = np.exp(df)
        print(df)
    return(df,obs,dfxy)


def log2(df,obs,dfxy):
    npones = np.ones(df.values.shape)/100
    newVals = np.maximum(df.values,npones)
    newVals = np.array(newVals,dtype=float)
    newVals = np.log2(newVals)
    df_calx = pd.DataFrame(data = newVals, index = df.index, columns = df.columns)
    return(df_calx,obs,dfxy)


def zscorev(df,obs,dfxy):
    df,obs,dfxy = zscore(df,obs,dfxy,ax=0)
    return(df,obs,dfxy)

def zscoreh(df,obs,dfxy):
    df,obs,dfxy = zscore(df,obs,dfxy,ax=1)
    return(df,obs,dfxy)

def zscore(df,obs,dfxy,ax=None):
    vals = df.values
    shape = vals.shape
    if ax == None:
        print("0 for vertical (by protein), 1 for horizontal")
        ax = int(input("axis (0/1):"))
    newA = np.zeros(shape)
    if ax == 0:
        for i in range(shape[1]):
            col = vals[:,i].tolist()            
            try:
                zCol = zScoreL(col)
            except:
                zCol = list(np.zeros(len(col)))
            newA[:,i] = zCol
    if ax == 1:
        for i in range(shape[0]):
            col = vals[i,:].tolist()
            zCol = zScoreL(col)
            newA[i,:] = zCol            
    return(pd.DataFrame(data=newA,columns=df.columns,index=df.index),obs,dfxy)

def zScoreL(lis):
    newLis = []
    mean = stat.mean(lis)
    std = stat.stdev(lis)
    if std == 0:
        return(np.zeros(len(lis)))
    for i in lis:
        newLis.append((i-mean)/std)
    return(newLis)

def outliers(df,obs,dfxy):
    df = df.clip(lower=df.quantile(q=0.0013), upper=df.quantile(q=0.9987), axis=1)  # cut by 3sigma = 0.0214 0.9786 or 4sigma = 0.0013 0.9987
    return(df,obs,dfxy)


def equalizeTMA1(df,obs,dfxy):
    try:
        obs.loc[:,"tissue type"]
        tName = "tissue type"
    except:
        tName = "tissue_type"
    utypes = list(obs[tName].unique())
    uSlides = list(obs["slide"].unique())

def equalizeTMA(df,obs,dfxy):
    #0 axis of array will be tissue type, 1 axis will be slide
    try:
        obs.loc[:,"tissue type"]
        tName = "tissue type"
    except:
        tName = "tissue_type"
    key3 = obs["slide_type"] == "JE"
    typesL = list(np.unique(obs.loc[key3,tName]))
    slidesL = list(np.unique(obs.loc[key3,"slide"]))
    print(typesL,slidesL)
    countsA = np.zeros((len(typesL),len(slidesL)))
    for i,ty in enumerate(typesL):
        key1 = obs[tName] == ty
        for j,sl in enumerate(slidesL):            
            key2 = obs["slide"] == sl
            countsA[i,j] = obs.loc[key1&key2,:].shape[0]
    
    shape = countsA.shape
    newdfs = []
    newxys = []
    newobs = []
    print(countsA)
    for i,ty in enumerate(typesL):
        minct = min(countsA[i,:])  
        if minct > 0:
            print(ty)
            maxct = max(countsA[i,:])
            key2 = obs[tName] == ty
            for j,sl in enumerate(slidesL):
                key1 = obs["slide"] == sl           
                sdf = df.loc[key1&key2,:]
                sxy = dfxy.loc[key1&key2,:]
                sobs = obs.loc[key1&key2,:]
                #print(sdf.index)
                if len(sdf.index)>0:
                    inds = random.choices(sdf.index,k=int(maxct))
                    newdfs.append(sdf.loc[inds,:])
                    newxys.append(sxy.loc[inds,:])
                    newobs.append(sobs.loc[inds,:])
    newdfs.append(df.loc[~key3,:])
    newxys.append(dfxy.loc[~key3,:])
    newobs.append(obs.loc[~key3,:])
    df = pd.concat(newdfs)
    df.index=np.arange(df.shape[0]).astype(str)
    dfxy = pd.concat(newxys)
    dfxy.index=np.arange(df.shape[0]).astype(str)
    obs = pd.concat(newobs)
    obs.index=np.arange(df.shape[0]).astype(str)
    return(df,obs,dfxy)


def combat(df,obs,dfxy):
    bayesdata = combat1.main(df.transpose(),obs["batch"])
    df = bayesdata.transpose()
    return(df,obs,dfxy)


def TMAcombat(df,obs,dfxy):
    if input("0: JE-TMA combat  or  1: -KC cell line combat: ") == "0":
        smallDF = df.loc[obs["slide_type"]=="JE",:]
        if smallDF.shape[0] < 10:
            print("ERROR! missing 'JE' 'slide_type' annotation. Must add for JE combat")
        smallOBS = obs.loc[obs["slide_type"]=="JE",:]
    else:
        print("attempting COMBAT with -KC, JEjenum and Tonsil internal cell lines")
        smallDF = df.loc[obs["internal"] == 'True',:]
        if smallDF.shape[0] < 10:
            print("ERROR! missing 'internal' True/False annotation. Must add for KC combat")
        smallOBS = obs.loc[obs["internal"] == 'True',:]
        
    gamma_star, delta_star, stand_mean, var_pooled = combat1.combat_fit(smallDF.transpose(), smallOBS["batch"])
    bayesdata=combat1.combat_transform(df.transpose(), obs["batch"], gamma_star, delta_star, stand_mean, var_pooled)
    df = bayesdata.transpose()
    return(df,obs,dfxy)

def equalizeBiomLevel(df,obs,dfxy):
    print("only makes sense for vertically z-scored data!")
    for i,b in enumerate(pd.Series(df.columns).sort_values()):
        print(i,b)
    housekeep = list(pd.Series(df.columns).sort_values())[int(input('which biom number to use as key:'))]
    print(housekeep,list(pd.Series(df.columns).sort_values()))
    ch,uch = obMenu(obs,title="equalize each (slide/batch/etc.)")
    if input("only consider specific cell type when calculating means? (y)") == 'y':
        ndf,nobs,nxy = pick(df,obs,dfxy)
    else:
        ndf,nobs,nxy = df.copy(),obs.copy(),dfxy.copy()
    means = []
    for bat in uch:
        key = nobs[nobs.columns[ch]] == bat
        ss = ndf.loc[key,housekeep]
        means.append(ss.mean())
    print(means)
    for i,bat in enumerate(uch):
        key = obs[obs.columns[ch]] == bat
        df.loc[key,:] -= means[i]
    return(df,obs,dfxy)
    
    
    
    

def remNegatives(df,obs,dfxy):
    for col in df.columns:
        if min(df[col]) < 0:
            df[col] -= min(df[col])
    return(df,obs,dfxy)


def savepngs(df,obs,dfxy,cdf,savefold=r"\\graylab\BCC_Chin_Lab_RDS\ChinData\Cyclic_Analysis\cmIF_2022-03-xx_mTMA2\figures IY\040423 LC\pngs",nontr = "ask"):
    if type(nontr) != bool:
        if input("make non-transparent for zen/tiff? (y)") == 'y':
            nontr = True
        else:
            nontr = False
    cols = []
    for i,col in enumerate(obs.columns):
        print(i,col)
    while True:
        try:
            imp = input("col save all unique as pngs (send blank to skip) send 'range(0,3)' for [0,1,2])   ")
            #print("imp",imp)
            imp = eval(imp)
            try:
                for im in imp:
                    cols.append(obs.columns[im])
            except:
            #print(imp)
                cols.append(obs.columns[imp])
        except Exception as e:
            i2 = input(str(e)+"\ndone? (y/''):")
            if i2=="y" or i2 == '':
                break
    savefold = checkChange(savefold,"save folder")
    imAL,names = zapari(obs,cols,str(obs.loc[:,"slide_scene"].values[0]))
    alldata = []
    for i,n in enumerate(names):
        print(n)
        #io.imsave(sphelper(savefold,n),imAL[i])
        #with open(sphelper(savefold,n),"w") as file:
            #file.write(imAL[i])
        #print(imAL[i])
        #time.sleep(2)

        #next section from https://www.geeksforgeeks.org/create-transparent-png-image-with-python-pillow/ cause just sending an array of (255,255,255) and (255,255,255,0) doesn't work for some bleeding reason
        if not nontr:
            im0 = PIL.Image.fromarray(imAL[i])
            alldata.append(im0)
            im = im0.convert("RGBA")
            datas = im.getdata() 
            newData = []
            for item in datas:
                if item[0] == 0 and item[1] == 0 and item[2] == 0:  # finding black colour by its RGB value
                    # storing a transparent value when we find a black colour
                    newData.append((255, 255, 255, 0))
                else:
                    newData.append((255, 255, 255))   
                #print(newData)
            im.putdata(newData)
        else:
            im = PIL.Image.fromarray(imAL[i])  
        #alldata.append(im)
        im.save(sphelper(savefold,n),"PNG")
    if False:#nontr:
        with tifffile.TiffWriter(savefold+"/"+input("tiff filename? (without extention)")+".ome.tif") as tif:
            for AD in alldata:
                tif.save(AD)#,photometric="rgb")
    if False:
        tifffile.imwrite(alldata,photometric="rgb")
        
    #io.MultiImage(savefold+"/"+"__".join(cols)+".tif",conserve_memory=False)
    if input("open napari viewer?") == 'y':
        NAPARI.Viewer()
    return(df,obs,dfxy)
        
def sphelper(foln,filn): 
     if not os.path.isdir(foln):
         foln = os.getcwd()
     badchars = [":","/","?",">","<"]
     for bc in badchars:
         if bc in filn:
             filn = filn.replace(bc,".")
     if not os.path.isdir(foln):
         print(foln)
         if input("folder does not exist. Try to create folder?") == 'y':
             os.mkdir(foln) 
         else:
             savefold = input("path to folder:")
             return(sphelper(savefold,filn))
     return(foln+"/"+filn+'.png')   

def zapari(obs,columns,slideScene,cellidn="cellid"):
    OUTLINE = True
    #blank = np.empty((), dtype=object)
    #blank[()] = (255,255,255,0)
    #VAL = np.empty((), dtype=object)
    #VAL[()] = (255,255,255)
    VAL = 1   #at 1, worked for transparent but not nontransp
    segpath = r'\\graylab\BCC_Chin_Lab_RDS\ChinData\Cyclic_Workflow\cmIF_2022-03-25_mTMA2\Segmentation\mTMA2-4_CellposeSegmentation/mTMA2-4_sceneH12_Ecad_nuc30_cell30_matched_CellSegmentationBasins.tif'#r"\\graylab\Chin_Lab\Cyclic_Workflow\cmIF_2021-11-21_HTAN-A\Segmentation\HTA9-14Bx1-7_CellposeSegmentation/HTA9-14Bx1-7_scene001_2200-2000-1500-15001.tif"#HTA9-14Bx1-7_scene001_Ecad_nuc30_cell30_matched_CellSegmentationBasins.tif"
    slideScene=checkChange(slideScene,"slide scene")
    cellidn = checkChange(cellidn,"cell ID category")
    segpath = checkChange(segpath,"segmentation image path (including file and extension- use nuc30_cell30_matched*.tif")
    keyS = obs.loc[:,"slide_scene"] == slideScene
    obs = obs.loc[keyS,:]
    obs = obs.astype(str)
    label = io.imread(segpath)
    iS = label.shape
    #print(label.shape)
    #ucells = np.unique(label)
    if OUTLINE:
        bounds = skimage.segmentation.find_boundaries(label, connectivity=1, background=0)
    j = 0
    imAL = []
    names = []
    for column in columns:
        obs["intcol"] = 0
        uEnts = sorted(list(obs.loc[:,column].unique()))#.astype(str)
        for i,e in enumerate(uEnts):
            print("cluster:",e)
            #chan = np.full_like(label,blank,dtype=object)
            chan = np.zeros_like(label)
            #ochan = np.zeros(iS)
            key = obs.loc[:,column] == e
            sobs = obs.loc[key,:]
            incl = list(sobs.loc[:,cellidn])
            print(len(incl))
            for cell in incl:
                j += 1
                #print(i)
                l = int(cell.split("cell")[-1])
                
                if OUTLINE:
                    chan[np.logical_and(label==l, bounds)] = VAL
                else:
                    chan[label==l] = VAL
            imAL.append(np.copy(chan))
            names.append(column+" "+e)
    return(imAL,names)
            

def checkChange(s,cat):
    if input(cat+":\n"+s+"\nchange? (y):") == 'y':
        return(input(": "))
    else:
        return(s)



def save(df,obs,dfxy):
    tag = input("filename?")
    df.to_csv(tag+"_df.csv")
    dfxy.to_csv(tag+"_dfxy.csv")
    obs.to_csv(tag+"_obs.csv")
    return(df,obs,dfxy)
    
def pick(df,obs,dfxy):
    ch0 = input("0:include,  1:exclude  ")
    for i,ob in enumerate(obs.columns):
        print(i,ob)
    
    print("filter slides by which?")
    ch1 = int(input("enter number:")) 
    categories = np.unique(np.array(sorted(obs[obs.columns[ch1]].astype(str))))
    for i,ob in enumerate(categories):
        print(i,ob)
    print("Enter one number at a time for entries to include/exclude")
    included = []
    if ch0 == "1":
        name = "exclude"
    else:
        name = "include"
    while True:
        try:
            included.append(categories[int(input(name+" #:"))])
        except:
            break
    print(included)
    key = np.empty((obs.shape[0],1),dtype=object)
    for i in included:
        key=np.append(key,np.array([obs[obs.columns[ch1]]==i]).T,axis=1)
    masterkey = []
    for i in range(key.shape[0]):
        row = key[i,:]
        if any(row):
            masterkey.append(True)
        else:
            masterkey.append(False)
    key=pd.Series(masterkey)
    if ch0 == "1":
        key = -key    
    key = list(key)
    df,dfxy,obs = df.loc[key,:],dfxy.loc[key,:],obs.loc[key,:]
    print("df shape",df.shape,obs.shape)
        
    return(df,obs,dfxy)


def roi(df,obs,dfxy):
    usc = sorted(list(obs["slide_scene"].unique()))
    ch1,uch = obMenu(obs,"color by?")
    colors = []
    if len(uch) < 8:
        colors = allc.standcolors
    else:
        colors = allc.colors
    while len(uch) > len(colors):
        colors.append(allc.rgb())
    for scene in usc:
        skey = obs.loc[:,"slide_scene"] == scene
        sxy  = dfxy.loc[skey,:]
        maxInY = max(sxy.iloc[:,1])
        spatialLite(obs,scene,colors,dfxy,uch,ch1)
        try:
            xmin = float(eval(input("x min: ")))
            xmax = float(eval(input("x max: ")))
            ymin = float(eval(input("y min: ")))
            ymax = float(eval(input("y max: ")))
        except:
            print("invalid coordinate, skipping")
            continue
        #Y = -y + max(yseries)
        #-Y + max(yseries) = y
        # = max(sxy.iloc[:,-1]) --- need to add min value
        xkey = dfxy.iloc[:,0] >= xmin
        Xkey = dfxy.iloc[:,0] <= xmax
        ykey = dfxy.iloc[:,1] <= -ymin + maxInY
        Ykey = dfxy.iloc[:,1] >= -ymax + maxInY
        newS = obs.loc[skey & xkey & Xkey & ykey & Ykey,:].copy()
        newD = df.loc[newS.index,:].copy()
        newX = dfxy.loc[newS.index,:].copy()
        print(obs.shape)
        obs = obs.loc[~skey,:]
        print(obs.shape)
        df = df.loc[~skey,:]
        dfxy = dfxy.loc[~skey,:]

        obs = pd.concat([obs,newS],axis=0)
        print(obs.shape)
        df = pd.concat([df,newD],axis=0)
        dfxy = pd.concat([dfxy,newX],axis=0)
        spatialLite(obs,scene,colors,dfxy,uch,ch1,ymin=ymin)
    return(df,obs,dfxy)
        
def spatialLite(nobs,scene,colors,nxy,uch,ch1,ymin=0):
        key=nobs["slide_scene"]==scene
        sobs = nobs.loc[key,:]
        #sdf = ndf.loc[key,:]
        sxy = nxy.loc[key,:]
        #ax.set_aspect('equal')
        #ax.legend(uch,colors,bbox_to_anchor=(1.05, 1), loc='upper left')
        try:
            fig,ax = plt.subplots(figsize=((max(sxy.iloc[:,0])-min(sxy.iloc[:,0]))/500,(max(sxy.iloc[:,1])-min(sxy.iloc[:,1]))/500))
            print(max(sxy.iloc[:,1]),"max Y")
        except Exception as e:
            print(e,"error setting fig and ax",scene)
            print(sxy.isna().any(),"isna")
            fig,ax = plt.subplots()
        for i,ty in enumerate(uch):
            co = colors[i]
            #print(sobs.columns[ch1])
            key1 = sobs[sobs.columns[ch1]]==ty
            #print(key1)
            #tobs = sobs.loc[key,:]
            #tdf = sdf.loc[key,:]
            txy = sxy.loc[key1,:]
            x = []
            y = []
            if txy.shape[0] == 0:
                continue
            for j in range(txy.shape[0]):
                pt = list(txy.iloc[j,:])
                #coords.append((pt[0],pt[1]))
                x.append(pt[0])
                y.append(-pt[1])
            Y = pd.Series(y)
            #sxy = list(sxy.astype(float))
            #print(sxy)
            Y += max(sxy.iloc[:,1])+ymin
            ax.scatter(x,Y,color=co,label=ty,s=1.2)
        lg = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')#, scatterpoints=1, fontsize=10)
        try:
            for k in range(len(uch)):
                lg.legendHandles[k]._sizes = [30]
        except Exception as e:
            print(e,"index k:",k,len(uch),lg.legendHandles)
        plt.title(scene+" "+obs.columns[ch1])    
        plt.show()
    


def obCluster(df,obs,dfxy):
    obs = obs.astype(str)
    for i,ob in enumerate(obs.columns):
        print(i,ob)
    ch = int(input("cluster by?"))
    chob = obs.columns[ch]
    clusters = obs[chob].copy()
    uobs = np.unique(obs[chob].values)
    for i,uo in enumerate(uobs):
        print(i,":",uo)
        clusters.loc[obs[chob]==uo] = uo
    obs["Cluster"] = list(clusters)
    return(df,obs,dfxy)

def kmeans(df,obs,dfxy):
    ncl = int(input("n clusters:"))
    km = KMeans(n_clusters=ncl)
    km.fit(df)
    obs["Kmeans "+str(ncl)] = km.labels_
    return(df,obs,dfxy)


def leiden(df,obs,dfxy):
    print(all(obs.index==df.index),"all index the same")
    res = float(input("recluster with resolution:"))
    adata = anndata.AnnData(df,obs = obs)
    sc.pp.neighbors(adata,use_rep='X')
    sc.tl.leiden(adata, key_added='Cluster', resolution=res)
    cn = "Leiden_"+str(res)
    obs[cn] = adata.obs["Cluster"]
    obs[cn] = obs[cn].astype(str)
    return(df,obs,dfxy)

def gmm(df,obs,dfxy):
    nClusters = int(input("how many different categories? "))
    ctypes = ['full','tied','diag','spherical']
    gmm = GMM(n_components=nClusters).fit(df)
    obs["labels"] = gmm.predict(df)
    #plt.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap='viridis');
    return(df,obs,dfxy)


'''
VISUALIZATION
'''


def nap(df,obs,dfxy,cdf):
    cols = []
    info = []
    for i,col in enumerate(obs.columns):
        print(i,col)
    while True:
        try:
            cols.append(obs.columns[int(input("number of column to include in Napari (send blank when done)"))])
        except:
            break
    print("!!! output will not render until kernel stops running !!!")
    for ss in obs.loc[:,"slide_scene"].unique():
        print(ss)
        fold = input("path to folder with segmentation masks (send blank to skip scene): ")
        if fold == "":
            continue
        i1 = input(".png segmentation file name: ")
        i2 = input(".tiff segmentation file name: ")
        info.append([ss,fold,i1,i2])
    for nfo in info:
        NP.main(slideScene=nfo[0],fold=nfo[1],ip=nfo[2], it=nfo[3],columns=cols)
    return(df,obs,dfxy)
        
        
    

def spatial(df,obs,dfxy,cdf):
    mpl.style.use('default') 
    #ch = input("show all? (y)")
    ch = "y"
    if ch != "y":
        ndf,nobs,nxy = pick(df,obs,dfxy)
    else:
        ndf,nobs,nxy = df.copy(),obs.copy(),dfxy.copy()
    #if  input("use napari? (y)") == "y":
        #napari(ndf,nobs,nxy,cdf)
        #return(df,obs,dfxy)
    SAVE=0
    #if input("save all? (y)") == "y":
        #SAVE = 1
        #location = input("save location? :").replace("\ "[0],"/")+"/"
    usc = sorted(list(nobs["slide_scene"].unique()))
    ch1,uch = obMenu(obs,"color by?")
    if input("pick colors?") == "y":
        colors = []
        for uuu in uch:
            print(uuu)
            colors.append(input("color:"))
    
    #elif len(uch) < 8:
        #colors = allc.standcolors
    else:
        colors = allc.colors
    while len(uch) > len(colors):
        colors.append(allc.rgb())
    for scene in usc:
        key=nobs["slide_scene"]==scene #was obs not nobs, haven't tested change
        sobs = nobs.loc[key,:]
        #sdf = ndf.loc[key,:]
        sxy = nxy.loc[key,:]
        #ax.set_aspect('equal')
        #ax.legend(uch,colors,bbox_to_anchor=(1.05, 1), loc='upper left')
        try:
            fig,ax = plt.subplots(figsize=((max(sxy.iloc[:,0])-min(sxy.iloc[:,0]))/500,(max(sxy.iloc[:,1])-min(sxy.iloc[:,1]))/500))
        except Exception as e:
            print(e,"error setting fig and ax",scene)
            print(sxy.isna().any(),"isna")
            fig,ax = plt.subplots()
        for i,ty in enumerate(uch):
            co = colors[i]
            #print(sobs.columns[ch1])
            key1 = sobs[sobs.columns[ch1]]==ty
            #print(key1)
            #tobs = sobs.loc[key,:]
            #tdf = sdf.loc[key,:]
            txy = sxy.loc[key1,:]
            x = []
            y = []
            if txy.shape[0] == 0:
                continue
            for j in range(txy.shape[0]):
                pt = list(txy.iloc[j,:])
                #coords.append((pt[0],pt[1]))
                x.append(pt[0])
                y.append(-pt[1])
            Y = pd.Series(y)
            #sxy = list(sxy.astype(float))
            #print(sxy)
            Y += max(sxy.iloc[:,-1])
            ax.scatter(x,Y,color=co,label=ty,s=0.8) #was 1.2 for HTA14
        lg = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')#, scatterpoints=1, fontsize=10)
        try:
            for k in range(len(uch)):
                lg.legendHandles[k]._sizes = [30]
        except Exception as e:
            print(e,"index k:",k,len(uch),lg.legendHandles)
        plt.title(scene+" "+obs.columns[ch1])        
        
        if SAVE == 1:
            while True:
                try:
                    plt.savefig(location +scene+" "+obs.columns[ch1]+".png")
                    break
                except Exception as e:
                    print("error with save location",e)
                    location = input("save location? :").replace("\ "[0],"/")+"/"
        plt.show()
    return(df,obs,dfxy)
                


def napari(df,obs,dfxy,cdf):
    ch,uch = obMenu(obs,title="Color by:")
    obCol = obs.iloc[:,ch]
    usc = obs["slide_scene"].unique()
    a = pd.Series(obs.index)
    b = a.apply(lambda n: n.split("cell")[-1])
    c = b.astype(int)
    c = pd.Series(data=c.values, index=obs.index)
    print(a,"a",b,"b",c,"c")
    obs["intID"] = ""
    obs["intID"] = c
    print(obs["intID"],"\n\n\n")
    outlabels = []
    for scene in usc:
        print(scene)
        uf = input("segmentation file folder (or 'skip'):")
        if uf == "skip":
            continue
        
        url = input("png image (CytoProj png): ").replace("\ "[0],"/")
        image = io.imread(uf+"/"+url)     
        url = input("copy in path to segmentation file (Nuclei segmentation basins tif): ").replace("\ "[0],"/")
        segments = io.imread(uf+"/"+url)  
        seg=segments.copy()#.astype(object)
        useg = np.unique(seg)
        
    
        
        print(seg,"seg")
        viewer = NAPARI.Viewer()
        viewer.add_image(image,blending='additive') 
        labels = obs.iloc[:,ch].loc[obs["slide_scene"] == scene]
        print(labels,"labels")
        #print(obCol,"obCol")
        
        
        print(useg,"useg")
        print(c.unique(),"int")
        for se in useg:
            shap = obCol.loc[c==se].shape[0]
            if shap == 0:
                seg[segments==se] = 0
                continue
            uc = obCol.loc[c==se].values[0]
            if uc not in outlabels:
                i = len(labels)
                outlabels.append([uc,i])
            else:
                i = outlabels.index(uc)
            seg[segments==se] = i
        
        viewer.add_labels(seg,blending='additive')
            
    outlabels = pd.DataFrame()
    outlabels.to_csv("Napari labels.csv")   
    return(df,obs,dfxy)
            
        

def obMenu(obs,title="choose category:"):
    for i,col in enumerate(obs.columns):
        print(i,col)   
    ch = int(input(title))
    uch = obs[obs.columns[ch]].unique()
    return(ch,uch)
    

def heatmap(df,obs,dfxy,cdf,title=None):
    print(min(10+cdf.shape[0]/5,2**15/100))
    h = min(10+cdf.shape[0]/5,2**15/100)
    f, ax = plt.subplots(figsize=(20, h))  
    bbox = ax.get_window_extent().transformed(f.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    height *= f.dpi
    print(height,"HEIGHT")
    sns.heatmap(cdf,xticklabels=cdf.columns,yticklabels=cdf.index,center=np.mean(cdf.values))
    if title != None:
        ax.title.set_text(title)    
    plt.show()
    cdf,obs,dfxy = zscore(cdf,obs,dfxy,ax=0)
    f, ax = plt.subplots(figsize=(20, min(10+cdf.shape[0]/5,2**15/100)))
    sns.heatmap(cdf,xticklabels=cdf.columns,yticklabels=cdf.index,center=np.mean(cdf.values))
    if title != None:
        ax.title.set_text(title)
    plt.show()
    return(df,obs,dfxy)

def clusterbar(df,obs,dfxy,cdf):
    if input("single cell heatmap instead?") == "y":
        singHeatmap(df,obs,dfxy)
        return(df,obs,dfxy)
    #fig,ax = plt.subplots(figsize=(15,min(10+cdf.shape[0]/5,2**15/100))) 
    chrc = input("add row colors? (y)")
    if chrc == "y":
        #ch,uch = obMenu(obs,title="color based on?")
        rowColors = {}
        for op in sorted(list(obs["Cluster"].unique())):
            rc = input(op+" color:")
            if rc in matplotlib.colors.BASE_COLORS or rc in matplotlib.colors.CSS4_COLORS:            
                rowColors[op] = rc
            else:
                print("invalid color, using gray")
                rowColors[op] = "gray" 
        rowColors = pd.DataFrame(index=rowColors.values(),data=rowColors.keys()).transpose()
        #rowColors = pd.DataFrame(sorted(list(obs["Cluster"].unique())))[0].map(rowColors)
    if input("z-score first? (y)") == "y":
        cdf,obs,dfxy = zscorev(cdf,obs,dfxy)

    try:
        sns.set(font_scale=1)
        if chrc == "y":
            
            print(rowColors)
            ay = sns.clustermap(cdf,yticklabels=True, xticklabels=True,center=0, cmap='bwr',figsize=(25,min(10+cdf.shape[0]/5,2**15/100)), row_colors =[rowColors])#vmin=-10, vmax=10,
        else:
            ay = sns.clustermap(cdf,yticklabels=True, xticklabels=True,center=0, cmap='bwr',figsize=(25,min(10+cdf.shape[0]/5,2**15/100)), )#,metric="euclidean")#vmin=-10, vmax=10
    except Exception as e:
        print(e,"error")
        return(df,obs,dfxy)
    plt.setp(ay.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    ayax = ay.ax_heatmap.yaxis.get_majorticklabels()
    for i,op in enumerate(obs.columns):
        print(i,op)
    ch = int(input("color by which?"))
    color = list(obs.columns)[ch]
    ucols = obs[color].unique()
    a = np.zeros((len(ayax),len(ucols)))
    bins = []
    for i,Text in enumerate(ayax):
        Text=Text.get_text()
        key1 = obs["Cluster"].astype(str) == str(Text)
        #print(str(Text))
        bins.append(str(Text))
        #print(key1,"1")
        for j,c in enumerate(ucols):
            key2 = obs[color] == c
            #print(key2,"2")
            bi = obs.loc[key1 & key2,:]
            a[i,j] = bi.shape[0]
    fig, ax = plt.subplots(figsize=(10,min(10+cdf.shape[0]/5,2**15/100)))#20))     
    shape = a.shape
    height = np.zeros(shape[0])
    for i in range(shape[1]):     
        try:
            ax.barh(bins,a[:,i],label=ucols[i],left=height,color=allc.colors[i])
        except:
            ax.barh(bins,a[:,i],label=ucols[i],left=height,color=allc.rgb())
        height += a[:,i]
    plt.gca().invert_yaxis()
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',title=color)
    plt.show()    
    mpl.style.use('default')
    return(df,obs,dfxy)


def sortedMap(df,obs,dfxy,vout=10,col_cluster=False):
    mpl.style.use('default')
    sns.set(font_scale=2)
    df,obs,dfxy = zscorev(df,obs,dfxy)
    print(obs.columns)
    #print(df.columns,"df cols")
    #cols = ['Primary Celltype autoCellType res: 1.0','immune checkpoints autoCellType res: 1.0']  
    for col in obs.columns:
        lines = [0]
        print(col)
        nkey = obs[col].isna()
        obs.loc[nkey,col] = "na"
        obs.loc[:,col] = obs.loc[:,col].astype(str)
        utys = sorted(list(obs.loc[:,col].astype(str).unique()))
        print(utys)
        if len(utys) < 2 or len(utys) > 20:
            print('no good')
            continue
        title = col
        data = []
        colors = []
        for i,uo in enumerate(utys):
            data.append(df.loc[obs.loc[:,col]==uo,:])
            colors.append(pd.Series(np.full_like(obs.loc[obs.loc[:,col]==uo,col],allc.colors[i])))
            title += "__"+uo+"-"+allc.colors[i]
            lines.append(lines[-1]+data[-1].shape[0])
        print(colors,'colors')
        data = pd.concat(data,axis=0)
        #print(data,data.index)
        #print(data.shape[0],pd.Series(data.index).unique().shape[0])
        colors = pd.concat(colors,axis=0)
        colors.index = data.index

        cN = col
        if ":" in cN:
            cN = ".".join(cN.split(":"))

        
        ax=sns.clustermap(data, cmap='bwr',row_colors=[colors],
                              yticklabels=False, xticklabels=True,figsize=(25,25), 
                              row_cluster=False, col_cluster=col_cluster,vmin=-10,vmax=10) #, vmin=-vout, vmax=vout,,center=0
        #ax=sns.clustermap(data, vmin=-vout, vmax=vout, cmap='bwr',row_colors=[colors],
                              #yticklabels=False, xticklabels=True,center=0,figsize=(25,25), 
                              #row_cluster=False, col_cluster=col_cluster)
        ay = ax.ax_heatmap
        for l in lines:
            ay.plot([0,df.shape[1]],[l,l])
        plt.title(title)
        #SD = r"Y:\Cyclic_Analysis\20210413_AMTEC_Analysis\data\figures_iy\2021-11-21\single heatmaps\hta14"+BX
        #if SAVEfIGS:
            #plt.savefig(save(0,BX+"_single-heatmap",cN,typ="png"))
        plt.show()
    mpl.style.use('default')
    return(df,obs,dfxy)

def singHeatmap(df,obs,dfxy,pc=False,vout=10):
    if input("z-score first? (y)") == "y":
        df,obs,dfxy = zscorev(df,obs,dfxy)
    if input("pick colors? (y)") == "y":
        pc = True
    for i,col in enumerate(obs.columns):
        print(i,col)
    cols = []
    while True:
        try:
            cta = int(input("number of column to include in colors (blank to escape): "))
            uty = obs.iloc[:,cta].unique()
            colCol = obs.iloc[:,cta].copy()
            colCol[:] = "gray"
            if pc:
                for ut in uty:
                    coch = input(ut+' color?')
                    colCol.loc[obs.iloc[:,cta]==ut] = coch
            else:
                for j,ut in enumerate(uty):
                    colCol.loc[obs.iloc[:,cta]==ut] = allc.colors[j]
            cols.append(colCol)
        except Exception as e:
            print(e)
            break
    stt = time.time()
    try:
        ay = sns.clustermap(df,yticklabels=True, xticklabels=True,center=0, cmap='bwr',figsize=(25,25),vmin=-vout, vmax=vout, row_colors =cols)#
    except Exception as e:
        print("plotting failed!\n",e)
        return(df,obs,dfxy)
    #plt.show()
    dgram = ay.dendrogram_row.dendrogram
    #print(dgram)
    #sortD = df.iloc[dgram["leaves"],:]
    distances = []
    for i in range(df.shape[0]-1):
        c1 = df.iloc[dgram["leaves"][i],:]
        c2 = df.iloc[dgram["leaves"][i+1],:]
        dists = (c1-c2)**2
        #print(dists)
        ds = dists.sum()**.5
        distances.append(ds)
    fz,az = plt.subplots()
    az.plot(distances,np.ones(len(distances))*len(distances)-1-np.arange(len(distances)),marker="o",markersize=.5,linestyle="None")
    md = stat.mean(distances)
    sd = stat.stdev(distances)
    lines = []
    res = int(input("minimum cells in cluster: ")) #must be highest distance for this number of distances each direction 
    sres = float(input("minimum SD of any cluster border: "))
    i = 0
    lastD = None
    dc = 0
    cN = "dgram cluster "+str(res)+"x"+str(sres)
    obs[cN] = 9999
    while i < len(distances):
        #print(i)
        d = distances[i]     
        if d > md+sd*sres:
            nearCells = distances[max(0,i-res):min(i+res,len(distances))]
            m = distances.index(max(nearCells))
            #print(i,m)
            if i == m:
                lines.append(m)
                if lastD != None:
                    lev = dgram['leaves'][lastD+1:i+1]
                    #print(lev)
                    obs.loc[:,cN].iloc[lev] = dc
                    dc += 1
                lastD = i
                i += res+1             
            elif m > i:
                i = m
            else:
                i += 1
        else:
            i += 1
    obs.loc[obs.loc[:,cN]== 9999,cN] = dc
    ax = ay.ax_heatmap
    print(lines)
    for l in lines:
        ax.plot([0,df.shape[1]],[l+1,l+1], 'k-', lw = 1)#ax.get_xlim()
    #ay.hlines(lines, *ay.get_xlim())
    plt.show()
    #heatmap(df,obs,dfxy,df)
    
    print((time.time()-stt)/60," minute runtime")
    mpl.style.use('default')
    return(df,obs,dfxy)

def neighPie(df,obs,dfxy,cdf):
    print("standard pie chart requires consistent radii at 25,50,70,100u")
    if input("custom pie chart? (y)")=="y":
        for Col in df.columns:
           print(Col)
        neighbors = []
        keySL = []
        while True:
            keySt = input("string to be found in all compared entities (any, send blank when done): ")
            if keySt == "":
                break
            else:
                keySL.append(keySt)
        goodSL = []
        while True: 
            inp = input("substring (radius) to include as layer (blank to exit): ")
            if inp == "":
                break
            goodSL.append(inp)
            neighbors.append([])
        
        for col in df.columns:
            for keyS in keySL:
                if keyS in col:
                    for i,st in enumerate(goodSL):
                        if st in col:
                            neighbors[i].append(col)
        #neighPieH1(df,obs,dfxy,neighbors)
        for slide in obs["slide"].unique():
            key = obs["slide"] == slide        
            sdf,sxy,sobs = df.loc[key,:],dfxy.loc[key,:],obs.loc[key,:]
            neighPieH1(sdf,sobs,sxy,neighbors)
            cats = ["Primary Celltype autoCellType res: 1.0",
            'tumor subtype autoCellType res: 1.0',
            'receptors autoCellType res: 1.0',
            'immune subtype autoCellType res: 1.0',
            'immune checkpoints autoCellType res: 1.0',
            'cytotoxic autoCellType res: 1.0',]
            for cat in cats:
                for uty in obs[cat].unique():
                    print(cat,uty)
                    key = obs[cat] == uty
                    tdf,txy,tobs = sdf.loc[key,:],sxy.loc[key,:],sobs.loc[key,:]
                    neighPieH1(tdf,tobs,txy,neighbors,cat=cat)
    else:
        odf,oobs,oxy = df.copy(),obs.copy(),dfxy.copy()
        neighPieH2(df,obs,dfxy)
        for slide in obs["slide"].unique():
            key = obs["slide"] == slide        
            sdf,sxy,sobs = df.loc[key,:],dfxy.loc[key,:],obs.loc[key,:]
            neighPieH2(sdf,sobs,sxy)
            cats = ["Primary Celltype autoCellType res: 1.0",
            'tumor subtype autoCellType res: 1.0',
            'receptors autoCellType res: 1.0',
            'immune subtype autoCellType res: 1.0',
            'immune checkpoints autoCellType res: 1.0',
            'cytotoxic autoCellType res: 1.0',]
            for cat in cats:
                for uty in obs[cat].unique():
                    
                    print(cat,uty)
                    key = obs[cat] == uty
                    tdf,txy,tobs = sdf.loc[key,:],sxy.loc[key,:],sobs.loc[key,:]
                    neighPieH2(tdf,tobs,txy,cat=cat)
    return(df,obs,dfxy)
        
        
        
'''        
        try:
            df,obs,dfxy = pick(odf.copy(),oobs.copy(),oxy.copy())
            neighPieH1(df,obs,dfxy)
        except ValueError:
            return(odf,oobs,oxy)
'''
    
def neighPieH1(df,obs,dfxy,neighbors,cat=None):  
    n = neighbors

    rds = []
    for i in range(len(neighbors)):
        rds.append(df.loc[:,n[i]])
    rds.reverse()
    #print(rds)
    size = 0
    fig, ax = plt.subplots()
    #fig1, ax1 = plt.subplots()
    allinf = []
    sizeinc = .9/len(rds)
    for rd in rds:
        #print(rd.sum(axis=0)) 
        lab = round(rd.sum(axis=0) *100/rd.sum(axis=0).sum(),1)
        #allinf.append(list(lab)+[mixingScore(rd)])
        lab = rd.columns+"_"+lab.astype(str)+"%"
        try:
            ax.pie(rd.sum(axis=0), radius=1-size, colors=plt.colormaps["tab20c"]([0,4,6,8,16]),
                   wedgeprops=dict(width=sizeinc, edgecolor='w'),labels=lab) 
        
        except:
            print("could not print to plot",cat) 
            return(df,obs,dfxy)

        size += sizeinc
    #print(allinf)
    try:
        CL = []
        for LAb in rd.columns:
            CL.append(LAb.split("_")[0])
        tab = ax.table(cellText=allinf,colLabels=CL)#+["tum_N/imm_N"])
        tab.auto_set_font_size(False)
        tab.set_fontsize(10)
        plt.tight_layout()
    except Exception as e:
        print(e,"can't make table")
    if cat == None:         
        ax.set(title=" ".join(list(obs['slide'].unique())+list(obs["Primary Celltype autoCellType res: 1.0"].unique())))
    else:
        ax.set(title=" ".join(list(obs['slide'].unique())+list(obs[cat].unique()))+" "+cat+"  total cells="+str(df.shape[0]))

    
    plt.show()
    return(df,obs,dfxy)

def neighPieH2(df,obs,dfxy,cat=None):  
    #plt.rcParams["figure.autolayout"] = True
    neighbors = []
    while len(neighbors) < 4:
        neighbors.append([])
    
    for i in range(4):
        keyS = "ors_"+str(25*(i+1))
        for col in df.columns:
            if keyS in col:
                neighbors[i].append(col)
    n = neighbors
    #print(n)
    d25,d50,d75,d100 = df.loc[:,n[0]],df.loc[:,n[1]],df.loc[:,n[2]],df.loc[:,n[3]]
    rds = [d25,d50,d75,d100]
    rds.reverse()
    #print(rds)
    size = 0
    fig, ax = plt.subplots(figsize=(8,8))
    #fig1, ax1 = plt.subplots()
    allinf = []
    for rd in rds:
        #print(rd.sum(axis=0)) 
        lab = round(rd.sum(axis=0) *100/rd.sum(axis=0).sum(),1)
        allinf.append(list(lab)+[mixingScore(rd)])
        lab = rd.columns+"_"+lab.astype(str)+"%"
        try:
            ax.pie(rd.sum(axis=0), radius=1-size, colors=plt.colormaps["tab20c"]([0,4,6,8,16]),
                   wedgeprops=dict(width=.2, edgecolor='w'),labels=lab) 
        
        except:
            print("could not print to plot",cat) 
            return(df,obs,dfxy)

        size += .2  
    #print(allinf)
    CL = []
    for LAb in rd.columns:
        CL.append(LAb.split("_")[0])
    tab = ax.table(cellText=allinf,colLabels=CL+["tum_N/imm_N"])
    tab.auto_set_font_size(False)
    tab.set_fontsize(10)
    ax.axis('tight')
    #ax.axis('off')
    if cat == None:         
        ax.set(title=" ".join(list(obs['slide'].unique())+list(obs["Primary Celltype autoCellType res: 1.0"].unique())))
    else:
        ax.set(title=" ".join(list(obs['slide'].unique())+list(obs[cat].unique()))+" "+cat+"  total cells="+str(df.shape[0]))
    
    plt.show()
    return(df,obs,dfxy)

def mixingScore(rd):
    for col in rd.columns:
        if "umor_Primary" in col:
            tcol = col
            continue
        if "mmune_Primary" in col:
            icol = col
            continue
    try:
        score = rd[tcol].sum()/rd[icol].sum()
        print(score,"tumor interactions / immune interactions")
        return(round(score,2))
    except:
        return(None)
            


def barplot(df,obs,dfxy,cdf,piechart = False,biomchart=False):
    oobs = obs.copy()
    obs = obs.astype(str)
    #rch = int(input("show actual number (0) or percentages (1)?"))
    for i,he in enumerate(obs.columns):
        print(i,he)
    ch = int(input("sort x axis by:"))
    binCol = obs.columns[ch]
    obs = obs.sort_values(binCol,axis=0) #THIS MESSES UP THE INDEX
    bins = list(obs.iloc[:,ch].unique())
   
    for i,he in enumerate(obs.columns):
        print(i,he)
    ch2 = int(input("color bars by:"))  
    pc = input("pick colors for barplot? (y)")
    colors = sorted(list(obs.iloc[:,ch2].unique()))
    colCol = obs.columns[ch2] 
    a = np.zeros((len(bins),len(colors)))   
    b = np.zeros((len(bins),len(colors)))   
    c = b.copy()
    d = c.copy()
    #table = np.zeroes(len(bins),len(colors))    
    nClrs = []
    for clr in colors:
        nClrs.append(obs.loc[obs.loc[:,colCol] == clr,:].shape[0])
    
    for i in range(len(bins)):
        Bin = obs.loc[obs[binCol]==bins[i],:] #sobs
        nBin = Bin.shape[0]
        for j in range(len(colors)):
            n = Bin.loc[Bin[colCol]==colors[j],:].shape[0] #tobs
                
            if nBin == 0:
                nBin = 1
            a[i,j] = n
            b[i,j] = n/nBin
            c[i,j] = n/nClrs[j]
    ctots = np.sum(c,axis=1)
    print(ctots,ctots.shape)
    for i in range(c.shape[0]):
        d[i,:] = c[i,:]/ctots[i]*100
    #print(bins,len(bins),len(bins[0]))              
    if biomchart:
        in2 = False
        in3 = False
        if input("use z-scored (by biomarker) data? (y) for errorbarplot") == "y":
            in2 = True
        if input("pick colors for errorbarplot? (y)") == "y":
            in3 = True
        zdf,ttt,tttt = zscorev(df,obs,dfxy)
        #barx = 0
        for binName in bins:
            kEY1 = obs.iloc[:,ch] == binName
            if kEY1.sum() < 2:
                continue
            for coloName in sorted(list(obs.loc[:,colCol].unique())):
                kEY2 = obs.loc[:,colCol] == coloName
                if in2:
                    bdf = zdf.loc[kEY1 & kEY2,:] 
                else:
                    bdf = df.loc[kEY1 & kEY2,:] 
                if bdf.shape[0] > 2:
                    fig1,ax1=plt.subplots(figsize=(bdf.shape[1]/10,6))
                    mean = bdf.mean(axis=0)
                    stdev = bdf.std(axis=0)
                    if in3:
                        keyss = []
                        colorss = []
                        while True:
                            keyS = input("key string for color assignment (sent blank to exit)")
                            if keyS == "":
                                break
                            else:
                                keyss.append(keyS)
                                colorss.append(input("color for bars containing this string"))
                        for i in range(bdf.shape[1]):
                            cn = bdf.columns[i]
                            for j,kS in enumerate(keyss):
                                if kS in cn:
                                    clr = colorss[j]
                                    ax1.errorbar(i,mean[i],yerr=stdev[i],
                                                 marker="o",ls='none',ecolor=clr,c="black") 
                                    break
                            else:
                                ax1.errorbar(i,mean[i],yerr=stdev[i],
                                             marker="o",ls='none',ecolor="black",c="black")                                
                                    
                                    
                                
                            
                    else:
                        ax1.errorbar(np.arange(mean.shape[0]),mean,yerr=stdev,
                                     marker="o",ls='none') 
                    ax1.set(title=binName+" "+coloName)
                    ax1.grid(axis='x')
                    plt.xticks(np.arange(zdf.shape[1]),labels=zdf.columns, rotation='vertical', fontsize=7)
                    plt.show()

                #barx += mean.shape[0]
    if piechart:
        for k in range(a.shape[0]):
            f,ax = plt.subplots()
            ax.pie(a[k,:],autopct='%.0f%%', labels=colors)
            plt.title(bins[k])
            ax.axis('equal')
            plt.show()
    
    colorDict = {}
    if pc == "y":
        for i in range(a.shape[1]):
            print(colors[i])
            colorDict[colors[i]] = input("color:")
    for array in [a,b,c,d]:   
        a = array.copy()
        fig, ax = plt.subplots(figsize=(15+len(bins)/5,10))
        shape = a.shape
        #print("a",a)
        height = np.zeros(shape[0])
        
        for i in range(shape[1]):    
            
            try:
                ax.bar(bins,a[:,i],label=colors[i],bottom=height,color=colorDict[colors[i]]) 
            except Exception as e:
                ax.bar(bins,a[:,i],label=colors[i],bottom=height)
            #cont = ax.BarContainer()
            #labs = np.round(a[:,i],3)
            #print(labs)
            #ax.bar_label(cont, labels=labs, label_type='center')
            height += a[:,i]
        plt.xticks(rotation = 85)
        plt.xlabel(binCol)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',title=colCol)
        plt.show()



    daf = pd.DataFrame(data=b,index=bins,columns=colors).transpose()
    f, ax = plt.subplots(figsize=(10+daf.shape[1]/5, 10+daf.shape[0]/5))
    sns.heatmap(daf,xticklabels=daf.columns,yticklabels=daf.index,center=np.mean(daf.values),annot=True)
    plt.show()    
    return(df,oobs,dfxy)


def boxplot(df,obs,dfxy,cdf):
    fig,ax = plt.subplots(figsize=(len(list(df.columns)),8))
    ax=sns.boxplot(data=df)
    plt.xticks(rotation = 85)
    plt.show()    
    vich = input("violin plot instead? (y): ")
    for i,he in enumerate(obs.columns):
        print(i,he)
    
    ch = int(input("sort x axis by:"))
    binCol = obs.columns[ch]
    bins = sorted(list(obs.iloc[:,ch].unique()))
    for i,he in enumerate(obs.columns):
        print(i,he)
    try:
        print("send non-int to have same color")
        ch = int(input("color boxes by:"))
        colCol = obs[obs.columns[ch]]  
        ucols = colCol.unique()
    except Exception as e:
        print(e)
        colCol = None
        ucols = [1]
    dfo = df.merge(obs,left_index=True,right_index=True).sort_values(binCol)
    #print(dfo)
    for marker in df.columns:
        fig,ax = plt.subplots(figsize=(max(8,int(len(bins)*len(ucols)**.7/2)),8))
        if vich == "y":
            ax=sns.violinplot(hue=colCol, data=dfo, x=binCol, y=marker,showfliers=False)
        else:
            ax=sns.boxplot(hue=colCol, data=dfo, x=binCol, y=marker,showfliers=False)
        try:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',title=colCol)
        except:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xticks(rotation = 85)
        plt.show()
    return(df,obs,dfxy)

def boxplot1(df,obs,dfxy,cdf):
    switch = 0
    if input("violin plot instead? (y)") == "y":
        switch = 1
    for i,he in enumerate(obs.columns):
        print(i,he)
    ch = int(input("sort x axis by:"))
    binCol = obs.columns[ch]
    bins = list(obs.iloc[:,ch].unique())
    for i,he in enumerate(obs.columns):
        print(i,he)
    ch = int(input("color boxes by:"))
    colCol = obs.columns[ch]    
    dfo = df.merge(obs,left_index=True,right_index=True)
    #print(dfo)
    for marker in df.columns:
        fig,ax = plt.subplots(figsize=(24,8))
        if switch == 0:
            ax=sns.boxplot(hue=colCol, data=dfo, x=binCol, y=marker,showfliers=False)
        else:
            sns.violinplot(hue=colCol, data=dfo, x=binCol, y=marker,showfliers=False)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xticks(rotation = 85)
        plt.show()
    return(df,obs,dfxy)



def correlationMatrix(df,obs,dfxy,cdf):
    ordf,oobs,oxy = df.copy(),obs.copy(),dfxy.copy()
    switch = False
    if input("save correlations csv?") == "y":
        switch = True
    if input("z-score data vertically?") == "y":
        df,obs,dfxy=zscorev(df,obs,dfxy)
    col = df.columns
    
    ch,uch = obMenu(obs,"sort correlation matrices by:")
    odfs = []
    df1 = df
    a = df1.values
    shape = a.shape
    out = np.zeros((shape[1],shape[1]))
    for e in range(shape[0]):
        for i in range(shape[1]):
            for j in range(shape[1]):
                out[i,j] += a[e,i] * a[e,j] * 1/shape[0]
    odf = pd.DataFrame(data=out,columns = col,index=col)
    heatmap(df,obs,dfxy,odf,title="all data")
    odfs.append(odf)
    for uc in uch:
        df1 = df.loc[obs.iloc[:,ch]==uc,:]
        a = df1.values
        shape = a.shape
        out = np.zeros((shape[1],shape[1]))
        for e in range(shape[0]):
            for i in range(shape[1]):
                for j in range(shape[1]):
                    out[i,j] += a[e,i] * a[e,j] * 1/shape[0]
        odf = pd.DataFrame(data=out,columns = col,index=col)
        heatmap(df,obs,dfxy,odf,title=uc)
        if switch:
            fn = input("filename?")
            if ".csv" not in fn:
                fn = fn+".csv"
            odf.to_csv(fn)
        odfs.append(odf)
    for i in range(len(odfs)-1):
        for j in range(len(odfs)-1):
            if j > i:
                nodf = odfs[j] - odfs[i]
                heatmap(df,obs,dfxy,nodf,title=uch[j]+" - "+uch[i])
    return(ordf,oobs,oxy)

def scatterplot(df,obs,dfxy,cdf):
    df["celln"]=np.arange(df.shape[0])
    for i,col in enumerate(df.columns):
        print(i,col)
    ch1=int(input("x axis number:"))
    ch2=int(input("y axis number:"))
    x,y =df.iloc[:,ch1],df.iloc[:,ch2]
    #y = x*h
    #x-*y = x- * x * h
    #(xT*x)- * xt*y) = h
    X = pd.concat((x,pd.Series(np.ones_like(x),index=x.index)),axis=1).values
    print(X)
    theta = np.linalg.inv(X.T.dot(X))#.dot(X.T.dot(y))
    print(theta)
    theta = theta.dot(X.T.dot(y))
    print(theta,"theta")
            #np.linalg.inv(X.T.dot(X)).dot(X.T).dot(y)
    #print(theta)
    #y = theta*x
    dumx = np.arange(int(max(x/10)))
    print(dumx,"dx")
    dumy = np.matmul(theta,np.append([np.ones_like(dumx)],[dumx],axis=0))
    print(dumy,"dy")
    xn,yn = df.columns[ch1],df.columns[ch2]
    fig,((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex='col', sharey='row',figsize=(10, 10))
    ax3.scatter(x,y, color='red', marker='x',s=.1)   
    nb = int(200)#df.shape[0]/100)
    X,Y=makeHist(y,nb,orientation="horizontal")
    ax4.plot(X,Y)
    X,Y=makeHist(x,nb)
    ax1.plot(X,Y)
    if input('include line of best fit? (y)') == 'y':
        ax3.plot(dumy)
    ax3.set_xlabel(xn+" level")
    ax3.set_ylabel(yn+" level")
    ax4.set_xlabel(yn+" number of cells (log10)")
    ax1.set_ylabel(xn+" number of cells (log10)")
    plt.show()
    return(df,obs,dfxy)

def makeHist(x,nb,orientation="vertical"):#),log2=False):
    mx,Mx = min(x),max(x)
    rx = Mx-mx
    sx = rx/nb
    binCts=[]
    bins = np.arange(mx,Mx+sx,sx)
    for i in range(1,len(bins)):
        key = x>=bins[i-1]
        key1 = x<bins[i]
        ss = x.loc[key&key1]
        if ss.shape[0] > 1:
            binCts.append(np.log10(ss.shape[0])/np.log10(2))
        else:
            binCts.append(0)
    y = np.arange(len(binCts))*sx
    X = binCts
    #print(X,y)    
    if orientation == "vertical":        
        return(y,X)
    else:
        return(X,y)

def scanpyv(df,obs,dfxy,cdf=None):
    if input("z-score data vertically? (y)") == "y":
        df,obs,dfxy = zscorev(df,obs,dfxy)
    ch,uch = obMenu(obs,"sort x axis by:")
    gb = obs.columns[ch]
    adata = anndata.AnnData(df,obs = obs)
    sc.pl.dotplot(adata, adata.var_names, groupby=gb, show=True)
    sc.pl.stacked_violin(adata, adata.var_names, groupby=gb, show=True)
    sc.pl.matrixplot(adata, adata.var_names, groupby=gb, show=True,dendrogram=True )
    return(df,obs,dfxy)
                       
def showUmap(df,obs,dfxy,cdf=None): #make it so there's an option to re-color without re-genning whole map
    #26obs = obs.astype(str)    
    print(df.index)
    print(obs.index)
    adata = anndata.AnnData(df,obs = obs)
    sc.pp.neighbors(adata,use_rep='X')
    sc.tl.umap(adata)
    plt.rcParams['figure.figsize'] = 8, 8    
    while True:
        try:
            for i,clr in enumerate(obs.columns):
                print(i,clr)
            color = obs.columns[int(input("color by (send non-int to escape)"))]  
            #obs[color] = obs[color].astype(str)
            print("0:show all or 1:select whic or 2:each unique value in own map?")
            ch = input("0, 1, or 2:")
            itm = []
            clr = []
            #shapes = []
            if ch == "1":
                uColors = sorted(list(obs[color].unique()))
                for c in uColors:
                    print("show? ",c)
                    try:
                        cho = input("0 for no, 1 for yes:")
                    except:
                        cho = "0"
                    if cho == "0":
                        #obs[color].loc[obs[color]==c] = float('nan')
                        itm.append(c)
                        clr.append('lightgray')
                        #shapes.append("s")
                    else:
                        itm.append(c)
                        clr.append(input("which color?"))
                        #shapes.append(input("which marker? . v p s * ><^ "))
        
            print("umap color:",color)
            cd = dict(zip(itm,clr))  
            if ch=="1":
                sc.pl.umap(adata,color = color,palette=cd)#,na_color='lightgray'
            elif ch == "2":
                uColors = sorted(list(obs[color].unique()))
                for i,clr in enumerate(uColors):
                    itm = []
                    clr = []
                    for j in range(len(uColors)):
                        itm.append(uColors[j])
                        if i == j:
                            clr.append("blue")
                        else:
                            clr.append("lightgray")
                    cd = dict(zip(itm,clr)) 
                    sc.pl.umap(adata,color = color,palette=cd)
                            
            else:
               sc.pl.umap(adata,color = color,na_color='lightgray') 
            plt.show()
        except Exception as e:
            print(type(e),e)
            break
    ch = input("show umap of biomarkers? (y/n)")
    if ch == "":
        pass
    elif ch[0] == "y" or ch[0] == "Y":
        for biom in df.columns:
            sc.pl.umap(adata,color=biom,vmin=np.mean(df[biom])-np.std(df[biom]),vmax=np.mean(df[biom])+2*np.std(df[biom]),color_map='viridis')     
            plt.show()
    return(df,obs,dfxy)


if __name__ == "__main__":
    #'''#"zzz_hta14_tumorneighborhoodcts1"
    folder = ""#r"C:\Users\youm\Desktop\src\zzzzzzzzzzz_current/"
    stem = '93_hta14'#87_LC-4'##'89_LC-4_withN'#''96_LC'#cl56_depth_study_H12'#'96_LC'#'97_mtma2'#'93_hta14'###'96_hta14_primary'#'97_hta14bx1_primary_celltype'#'99_hta14'#"temp"#"zzz_hta1499"#"zzz14bx1_97"#"hta14bx1 dgram"#folder+"14_both"##"tempHta14_200"#"HTA14f"#"zzzz_hta1498_neighborhoodsOnly"#"hta1415Baf1"#"HTA15f"#"0086 HTA14+15"#"99HTA14"#"z99_ROIs_5bx_HTA1415"#"temp"#"z99_ROIs_5bx_HTA1415"#<-this one has old celltyping no TN #"0084 HTA14+15" #"HTA9-14Bx1-7 only"#"0.93 TNP-TMA-28"#"0.94.2 TNP-TMA-28 primaries"#"1111 96 TNP-28" #'0093 HTA14+15'#"0094.7 manthreshsub primaries HTA14+15"#"0094 HTA14+15" #"096 2021-11-21 px only" #'095.08 primaries only manthreshsub 2021-11-21 px only'#"094 manthreshsub 2021-11-21 px only" #  '095.1 primaries only manthreshsub 2021-11-21 px only'#
    print("axis labels %s on barplots more ticks/lines")
    print(stem)
    df = pd.read_csv(stem+"_df.csv",index_col=0) 
    obs = pd.read_csv(stem+"_obs.csv",index_col=0).astype(str)
    dfxy = pd.read_csv(stem+"_dfxy.csv",index_col=0)
    print(df.shape[0],"cells")
    #print(manual_type_tma-35_dfxy"\n\n THIS DATA IS ALREADY Z-SCORED \n\n")
    '''
    df = pd.read_csv("zz_C_L.1_AMTEC_df.csv",index_col=0)
    obs = pd.read_csv("zz_C_L.1_AMTEC_obs.csv",index_col=0)
    dfxy = pd.read_csv("zz_C_L.1_AMTEC_dfxy.csv",index_col=0)#zL.1_final_amtec_dfxy
    '''
    
    
    
    
    main(df,obs,dfxy)
    
'''
tissue-to-tissue biomarker accessability varies
standard procedure for RNA-seq analysis (GAPDH gene)
normalization method
    housekeeping proteins      - every cell should have baseline level for basic cell homeostasis
        lamAC c750 (expression varies between tumor/nontumor)
        CD4 c647 (only immune cells, should be similar in all immune cells)
        CD31 647 (endothelial specific)
    
    continuous antigen   - expression level changes between cell activities/age
        CK7 c555 (PDAC)
        CK7 488  (Br)
        CK8 488
        CK5 488
        CK19 750
        aSMA 488
        Vim 488
        
        
cd31+ aSMA- CollIV- imature vessle/angeogenesis 
all+ mature blood vessle


...
ck7 ck5 ki67 cd45
'''        


'''
def a-utotype(df,obs,dfxy):
    odf = df.copy()
    if input("zscore?") == "y":
        df,obs,dfxy = zscorev(df,obs,dfxy)  
    df = onlyPrimaries(df)
    df["tumor"] = 0    
    df["stromal"] = 0
    df["immune"] = 0
    df["Ki67+"] = 0
    df["CD31+"] = 0
    df[" fibroblast"] = 0
    tums = 0
    stroms = 0
    ims = 0
    for biom in df.columns:
        ser = df[biom]
        #ser.loc[df[biom]<0] = 0
        if "CK" in biom or "Ecad" in biom or "Muc1" in biom:
            df["tumor"] += ser
            tums += .5
        if "CD" in biom and "31" not in biom and "44" not in biom:
            df["immune"] += ser
            ims += .5
        if "aSMA" in biom or "Vim" in biom or "Coll" in biom or "CAV1" in biom:
            df["stromal"] += ser
            stroms += .5
        if "CD31" in biom:
            df["CD31+"] += ser
        if "Ki67" in biom:
            df["Ki67+"] += ser
    df["tumor"] = df["tumor"]/tums
    df["immune"] = df["immune"]/ims
    df["stromal"] = df["stromal"]/stroms
            
    res = float(input("number of standard deviations (or raw intensity if not zscored) above mean required to count as +"))
    obs["autoCellType res: "+str(res)] = " "
    superB = ["CD31+","immune"," fibroblast","tumor","Ki67+"]#,"CD44+"]
    for sb in superB:
        key = df[sb]>res
        obs.loc[key,"autoCellType res: "+str(res)] += sb + " "
    if input("keep superbioms and z-scoring?") == "y":
        return(df,obs,dfxy)
    return(odf,obs,dfxy)
'''
        
        
        