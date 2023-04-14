# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 11:16:50 2021

@author: youm
"""

import os
import numpy as np
import pandas as pd
import math
import IFprocessing1 as ifp
import IFvisualization0 as ifv
import cmifAnalysis45 as cm
import SVM5 as sv
import time
import scipy as sp
import matplotlib as mpl
mpl.style.use('default')

    
DATAFOLDER = r'W:\ChinData\Cyclic_Analysis\cmIF_2022-03-xx_mTMA2\GM'#r'Y:\Cyclic_Workflow'#r"C:\Users\youm\Desktop\src"#r"Y:\Cyclic_Analysis\cmIF_2019-10-04_MB-TMA\Chin_2019-10" #r'Y:\Cyclic_Workflow\cmIF_2019-12-02_JP-TMA1\Segmentation'#r"Y:\Cyclic_Workflow\cmIF_2021-11-21_HTAN-A\Segmentation\HTA9-14Bx1-7_CellposeSegmentation"
SAVEFOLDER = r"C:\Users\youm\Desktop\src"
DATAFOLDER = DATAFOLDER.replace("\ "[0],"/")
SAVEFOLDER = SAVEFOLDER.replace("\ "[0],"/")


def main(dataFolder=DATAFOLDER,saveFolder=SAVEFOLDER):
    options = ["Import and clean data","load prepared data"]
    functions = [buildDataFrame,load]
    df,obs,dfxy = menu(options,functions)
    print(df.shape,obs.shape,dfxy.shape)     
    obs = obs.astype(str)
    op = ["data editing","processing","visualization","Support Vector Machine","old analysis tool"]
    fn = [loadingMenu,ifp.main,ifv.main,sv.main,cm.main]
    while True:
        df,obs,dfxy=menu(op,fn,df,obs,dfxy)
        inp = input("quit?")
        if inp == "y" or inp == "":
            return(df,obs,dfxy)

    
def menu(options,functions,df=9,obs=9,dfxy=9): #MANUAL MENU
    print("menu")
    while True:
        print("\n")
        for i,op in enumerate(options):
            print(i,op)
        try:
            print("send non-int when done (return df)")
            ch = int(input("number: "))
        except:
            return(df,obs,dfxy)
        df,obs,dfxy=functions[ch](df,obs,dfxy)  
        print(all(obs.index==df.index),"all index the same")
        try:
            print(obs.columns)
        except:
            print("no obs")
        try:
            print(df.shape,"df shape")
        except:
            pass
    
    
def loadingMenu(df=9,obs=9,dfxy=9): 
    op = ["save","drop cells with less than N% of data",
          "edit observations","drop columns based on key string","import biomarkers/observations","combine obs",
          "save unique list of obs for making import table","rename column",
          "combine another prepared dataset (and handle mixed partitions)","handle mixed partitons in existing data",
          "scale data (z-score, etc)","autoclean NA values","edit observation labels"]
    fn = [save,dropCells,editObs,dropCols,importBiom,
          combineObs,saveObs,renCol,combineData,doPart,scale,autoClean,editLabels] 
    df,obs,dfxy = menu(op,fn,df,obs,dfxy)
    print(df.shape,obs.shape,dfxy.shape)
    return(df,obs,dfxy)

def scale(df,obs,dfxy):
    op = ["z-score","log2","scaled raw (div by stdev)"]
    fn = [zscore, log2, scaledRaw]
    df,obs,dfxy = menu(op,fn,df,obs,dfxy)
    return(df,obs,dfxy)

def zscore(df,obs,dfxy):
    df = sp.stats.zscore(df,axis=0,nan_policy="omit")
    return(df,obs,dfxy)

def log2(df,obs,dfxy):
    df = np.log2(df)
    return(df,obs,dfxy)

def scaledRaw(df,obs,dfxy):
    stds = np.std(df,axis=0)
    df = df/stds
    return(df,obs,dfxy)
        


def doPart(df,obs,dfxy):   
    dShort=[]
    for col in df.columns:
        if df[col].isna().any():
            shortname = col.split("_")[0]
            if shortname not in dShort:
                dShort.append(shortname)

    for sCol in dShort:
        toComb = []
        for col in df.columns:
            if sCol == col.split("_")[0]:
                toComb.append(col)
        print("to combine:",sCol,toComb)
        if len(toComb) > 1:
            df=combinePart(df,["uclei","nuc"],toComb,"nuc",sCol)
            df=combinePart(df,["ellmem"],toComb,"cellmem",sCol)
            df=combinePart(df,["ucadj","ytopla","erinuc","cyto"],toComb,"cyto",sCol)
    print(df.columns[df.isna().any()])
    return(df,obs,dfxy)


def combinePart(df,partitionL,toComb,NName,sCol):
    print(toComb)
    sDF = pd.DataFrame(index=df.index)
    for biomarker in toComb:
        for name in partitionL:
            if name in biomarker:
                sDF[biomarker] = df[biomarker]
                continue
    if sDF.shape[1]>1:
        print(sDF.columns,"sdf cols\n")
        df.loc[:,sCol+"_"+NName+"_combined"]=sDF.max(axis=1)
    else:
        print("sDF only has one entry apparently", sDF.columns)
    return(df)


def doPartHierarchical(df,obs,dfxy):
    #KC:
    #Take this order to pick one partition per marker, except pERK.
    #cellmem2p25 > Cytoplasm > exp5 > perinuc5 for Ecad, HER2, EGFR, pARK, which values are in this order it they express.
    order = "cellmem > Cytoplasm > exp5 > perinuc".split(" > ")
    order.append("nuclei")
    print("order:",order)
    
    dShort=[]

    for col in df.columns:
        if df[col].isna().any():
            shortname = col.split("_")[0]
            if shortname not in dShort:
                dShort.append(shortname)

    for sCol in dShort:
        toComb = []
        for col in df.columns:
            if sCol == col.split("_")[0]:
                toComb.append(col)
        print("to combine:",sCol,toComb)
        if len(toComb) > 1:
            df[sCol+"_combined"]=np.nan
            for partition in order:
                print(partition)
                for biomarker in toComb:
                    if partition in biomarker:                  
                        c = toComb.pop(toComb.index(biomarker))
                        df.loc[df[sCol+"_combined"].isna(),sCol+"_combined"]=df[c]
                        print(toComb)
        if len(toComb) >0:
            print("ERROR: ",toComb," partition not included in order, adding as lowest priority")     
            for biomarker in toComb:
                df.loc[df[sCol+"_combined"].isna(),sCol]=df[biomarker]

    return(df,obs,dfxy)
        
    

def combineData(df,obs,dfxy):
    df1,obs2,xy3 = load(9,9,9)
    D = [df,obs,dfxy]
    D1 = [df1,obs2,xy3]
    if input("handle mixed partitions? (y)") == "y":
        dShort = []
        d1Short = []
        for col in df.columns:
            dShort.append(col.split("_")[0])
        for col in df1.columns:
            d1Short.append(col.split("_")[0])  
        
        for col in df.columns:
            if col not in df1.columns:
                #print(col,"tnp28 only")
                sCol = col.split("_")[0]
                if sCol in d1Short:
                    print(sCol)
                    df[sCol] = df[col]

        for col in df1.columns:
            if col not in df.columns:
                #print(col,"tnp28 only")
                sCol = col.split("_")[0]
                if sCol in dShort:
                    print(sCol)
                    df1[sCol] = df1[col]
                    
                    
    for i in range(3):
        d = D[i]
        print(d.shape)
        d1 = D1[i]
        d  = pd.concat([d,d1],axis=0)
        D[i] = d
        print(D[i].shape)
    print(D[0].columns[D[0].isna().any()])
    for co in D[0].columns:
        print(co,D[0][co].isna().sum())
    return(D[0],D[1],D[2])


def combineObs(df,obs,dfxy):
    tdf,nobs,txy = load(9,9,9)
    tdf,txy = 9,9
    if input("combine all? (does not overwrite same names) (y)") == "y":
        for col in nobs:
            if col not in obs:
                obs[col] = ""
                obs.loc[nobs.index,col] = nobs[col]
    else:
        for col in nobs.columns:
            if input("include (/overwrite)"+col+"? (y)") == "y":
                obs[col] = ""
                obs.loc[nobs.index,col] = nobs[col]                        
    return(df,obs,dfxy)

def importBiom(df,obs,dfxy):
    if input("import annotations from csv with matching column? (y)") == "y":
        return(importObs(df,obs,dfxy))
    return(impB(df,obs,dfxy))

    

def impB(df,obs,dfxy):
    ndf,tobs,txy = load(9,9,9)
    tobs,txy = 9,9
    for col in ndf:
        if input("import "+col+" ? (y)") == "y":
            df[col] = ""
            df.loc[ndf.index,col] = ndf[col]
    return(df,obs,dfxy)            

def importObs(df,obs,dfxy):
    while True:
        file = input("path to file including extension")
        try:
            nobs = pd.read_csv(file)
            break
        except Exception as e:
            print(e,"couldn't read file")
            if input("return to main menu? (y)") == 'y':
                return(df,obs,dfxy)
    ch,uch = obMenu(obs,title="obs column with matching values to new file")
    ch1,uch1 = obMenu(nobs,title="new column with matching values")
    nobs.index = nobs.iloc[:,ch1]
    toch = []
    for col in nobs.columns:
        print(col,nobs.loc[:,col].unique())
        if input("include column? (y)") == 'y': 
            if col not in obs.columns:
                obs[col] = ""
                toch.append(col)
            else:
                obs[col+"999"]  = ""
                toch.append(col+"999")
    for uc in uch:
        if uc not in uch1:
            continue
        key = obs.iloc[:,ch] == uc
        for ncol in toch:
            
            val = nobs.loc[uc,ncol]
            print("               ",uc,ncol,val)
            obs.loc[key,ncol] = val
    return(df,obs,dfxy)
            
        
        

def renCol(df,obs,dfxy):
    ds = [df,obs,dfxy]
    for d in ds:
        while True:
            print(d.columns)
            ip = input("column to rename?")
            nn = input("new name")
            try:
                d[nn] = d.pop(ip)
            except Exception as e:
                print("invalid renames",e)
            ch = input("edit another in dataframe? (y)")
            if ch == "":
                break
            elif ch[0] != "y" and ch[0] != "Y":
                break
    ch = input("sort dfs? (y)")
    if ch == "":
        return(ds[0],ds[1],ds[2])
    if ch[0] == "Y" or ch[0] == "y":
        for i,d in enumerate(ds):
            ds[i] = d.loc[:,d.columns.sort_values()] #this apparently doesn't work
        return(ds[0],ds[1],ds[2])
    return(ds[0],ds[1],ds[2])


def editLabels(df,obs,dfxy):
    ch,uch = obMenu(obs,"categorty to edit labels")
    obCol = obs.columns[ch]
    print("send blank to skip")
    for uc in uch:
        print(uc)
        nn = input("new label:")
        if nn != "":
            obs.loc[obs[obCol]==uc,obCol] = nn
    return(df,obs,dfxy)
        
        


def importMissing(df,obs,dfxy):
    d1,o1,x1 = buildDataFrame(9,9,9)
    for i in range(df.shape[0]):
        if df.iloc[i,:].isnull().sum()>0:
            ind = obs["index"].iloc[i]
            key = o1["index"] == ind
            newData = d1.loc[key,:]
            
            if newData.shape[0] > 0:
                print(newData,newData.shape)
                for col in df.columns:
                    if math.isnan(float(df[col].iloc[i]))>0:
                        if col in newData.columns:
                            #try:
                                df[col].iloc[i] = newData[col].values[0]
                            #except:
                                #print(ind,col,"no values")
    return(df,obs,dfxy)


def fillNA(df,obs,dfxy):  
    if input("Use average based on mean? (y):") == "y":
        counts = df.isnull().sum(axis=1)
        #mdf = df.loc[counts>0,:]
        mobs = obs.loc[counts>0,:]
        ch,uch = obMenu(mobs,title="import average means based on")
        for c in uch:
            key = obs.iloc[:,ch]==c
            sdf = df.loc[key,:]
            means = np.nanmean(sdf.values,axis=0)
            for i,col in enumerate(sdf.columns):
                sdf[col] = sdf[col].fillna(means[i])
            df.loc[key,:] = sdf            
            return(df,obs,dfxy)
    if input("fill with numerical value? (y):") == "y":
        val = float(input("value to fill for all missing values:"))
        df = df.fillna(val)
        return(df,obs,dfxy)
    return(df,obs,dfxy)

def flexMenu(title="String to include in list"):
    lis = []
    while True:
        ch=input(title+" (send blank when done): ")
        if ch == "":
            return(lis)
        lis.append(ch)
        
 
    
def obMenu(obs,title="choose category:"):
    for i,col in enumerate(obs.columns):
        print(i,col)   
    ch = int(input(title))
    uch = obs[obs.columns[ch]].unique()
    return(ch,uch)


def dropCols(df,obs,dfxy):
    if input("include based on keystring instead? (y)") == 'y':
        df, obs, dfxy = inclCols(df,obs,dfxy)
        return(df,obs,dfxy)
    lis = [df,obs,dfxy]
    for k,d in enumerate(lis):
        print(list(d.columns))
        toRem = flexMenu(title="remove all columns containing these strings (end with ! for exact string)")
        print("\nbefore:\n",lis[k].columns)
        dr = []
        for col in d.columns:
            for t in toRem:
                if t in col:
                    dr.append(col)
                elif '!' in t:
                    if t[:-1] == col:
                        dr.append(col)
        lis[k] = tryDrop(d,dr)
        print("after:\n",lis[k].columns)
    for d in lis:
        print("\n",d.columns)      
    return(lis[0],lis[1],lis[2])

def inclCols(df,obs,dfxy):
    lis = [df,obs,dfxy]
    for k,d in enumerate(lis):
        incl = []
        print(list(d.columns))
        toRem = flexMenu(title="include all columns containing these strings")
        if len(toRem) == 0:
            print("skipping!")
            continue
        print("\nbefore:\n",lis[k].columns)
        for col in d.columns:   
            for t in toRem:
                if t in col and col not in incl:
                    incl.append(col)
        print("\n\n",incl,"\n\n")
        lis[k] = d.loc[:,incl]
    return(lis[0],lis[1],lis[2])                    
    
    


def countNA(df,obs,dfxy):
    uSlide = obs['slide'].unique()
    uBiom = df.columns.unique()
    outD = pd.DataFrame(index=uSlide,columns = uBiom)
    for s in uSlide:
        key = obs["slide"] == s
        sdf = df.loc[key,:]
        denom = sdf.shape[0]
        for biom in uBiom:
            try:
                num = sdf[biom].isna().sum()
                outD.loc[s,biom] = num/denom
            except:
                print(s,biom)
                #pass
    outD.to_csv("NAN counts.csv")
    return(df,obs,dfxy)
    
def autoClean(df,obs,dfxy):
    #df,obs,dfxy,cdf = DFs[0],DFs[1],DFs[2],DFs[3]
    cho = 1 #drop 0:cells   1:columns(bioms)
    ch = 90 #"max missing % threshold integer (0 to drop all cells with missing values, 100 to keep all
    while ch > 0:
        if cho == 0:
            counts = df.isnull().sum(axis=1)
            #print(counts,counts.shape,df.shape)
            Mx = df.shape[1]
            pts = (np.ones(counts.shape[0]) - counts/Mx)*100
            pts = pd.Series(pts)
            key = pts >= 100-ch
            df = df.loc[key,:]
            obs = obs.loc[key,:]
            dfxy = dfxy.loc[key,:]
            cho = 1
            ch -= 10
        else:
            counts = df.isnull().sum(axis=0)
            #print(counts,counts.shape,df.shape)
            Mx = df.shape[0]
            pts = (np.ones(counts.shape[0]) - counts/Mx)*100
            pts = pd.Series(pts)
            key = pts >= 100-ch
            df = df.loc[:,key]
            cho = 0
            ch -= 10
    #for co in df.columns:
    #    print(co,df[co].isna().sum()/df.shape[0]*100)
    #print(counts,counts.shape,df.shape)
    return(df,obs,dfxy)          
        

def dropCells(df,obs,dfxy):
    if input("countNA to .csv?") == "y":
        df,obs,dfxy = countNA(df,obs,dfxy)
    for co in df.columns:
        print(co,df[co].isna().sum()/df.shape[0]*100)
    print("missing percents)")
    while True:
        try:
            cho = int(input("drop 0:cells   1:columns(bioms) : "))
            ch = int(input("max missing % threshold integer (0 to drop all cells with missing values, 100 to keep all):"))
            break
        except Exception as e:
            print(e)
    if cho == 0:
        counts = df.isnull().sum(axis=1)
        #print(counts,counts.shape,df.shape)
        Mx = df.shape[1]
        pts = (np.ones(counts.shape[0]) - counts/Mx)*100
        pts = pd.Series(pts)
        key = pts >= 100-ch
        df = df.loc[key,:]
        obs = obs.loc[key,:]
        dfxy = dfxy.loc[key,:]
    else:
        counts = df.isnull().sum(axis=0)
        #print(counts,counts.shape,df.shape)
        Mx = df.shape[0]
        pts = (np.ones(counts.shape[0]) - counts/Mx)*100
        pts = pd.Series(pts)
        key = pts >= 100-ch
        df = df.loc[:,key]
    for co in df.columns:
        print(co,df[co].isna().sum()/df.shape[0]*100)
    print(counts,counts.shape,df.shape)
    if input("label edge distance?") == "y":
        if "edge_distance" not in obs.columns:
            print("add edge_distance column to obs")
        else:
            minD = float(input("minimum distance from edge to keep: "))
            key = obs["edge_distance"].astype(float) > minD
            print(key)
            print(obs.shape)
            obs["edge distance >"+str(minD)] = "false"
            obs.loc[key,"edge distance >"+str(minD)] = "true"
            #obs=obs.loc[key,:]
            #df = df.loc[key,:]
            #dfxy = dfxy.loc[key,:]
            print(obs.shape)
    return(df,obs,dfxy)
    

def buildDataFrame(bl1,bl2,bl3,unp=True):
    print("build")
    DFs,names,goodStrs = getDFs()
    df = sortDFs(DFs,names,goodStrs)
    df,obs,dfxy = makeObs(df)
    if unp:
        obs = unpackObs(obs)
    return(df,obs,dfxy)

def saveObs(df,obs,dfxy):
    ch,uch = obMenu(obs,title="save unique entries in category:")
    pd.Series(uch).to_csv(input("filename? ")+".csv")
    return(df,obs,dfxy)

def importObs1(df,obs,dfxy):
    file = navigate("C:/Users/youm/.spyder-py3/src")
    file = pd.read_csv(file,dtype=object,index_col=0)
    file=file.values
    print(file)
    keys = file[:,0]
    print("keys:",keys)
    ch,uob = obMenu(obs,"obs to apply keys to")
    names = []
    #newEnts1 = []
    for i in range(file.shape[1]-1):
        newEnts = file[:,i+1]
        print(newEnts)
        #newEnts1+=list(newEnts)
        names.append(input("name for above set of obs: "))
    #di = pd.concat([pd.Series(keys),pd.Series(newEnts1)],axis=1)
    #print(di)
    for j,name in enumerate(names):
        obs[name] = "other"
        #for ob in uob:
        for i,okey in enumerate(keys):
            try:
                key = obs[obs.columns[ch]] == okey
                if not any(key):
                    print("no result for",okey)
                obs.loc[key,name] = file[i,j+1]
            except:
                print(name,okey,"does not have entry")
    return(df,obs,dfxy)
        
        
        


def unpackObs(obs):
    obs["index"] = obs.index
    while True:
        for i,col in enumerate(obs.columns):
            print(i,col)
        try:
            ch = int(input("split column number:"))
            print("example:",obs.iloc[0,ch])
            c = input("char to split with:")
        except:
            break
        newCols = obs.iloc[0,ch].split(c)
        for j,col in enumerate(newCols):
            print(col)
            name=input("type name or hit enter to discard: ")
            if name == "":
                continue
            else:                
                try:
                    obs[name] = obs.iloc[:,ch].apply(lambda n: n.split(c)[j])
                except:
                    print(name,"could not convert")
    ch = input("manually add obs? (1/y for yes)")
    if ch == "1" or ch == "y" or ch == "Y":
        while True:
            print("\n\nAdding new category...")
            name = input("name of category (send blank to exit): ")            
            if name == "":
                break
            for i,col in enumerate(obs.columns):
                print(i,col)
            try:
                ch = int(input("apply based on category number:"))
            except:
                break
            uobs = obs.iloc[:,ch].unique()
            try:
                obs.loc[:,name]
            except:
                obs[name] = ""
            for uo in uobs:
                print("\n",uo)
                try:
                    print("current:",obs.loc[obs.iloc[:,ch]==uo,name].iloc[0])
                except:
                    print("could not display current name")
                annot = input("enter annotation:")
                if annot == "":
                    continue
                obs.loc[obs.iloc[:,ch]==uo,name] = annot
    while True:
        try:
            ch1,uch = obMenu(obs,"observation 1 to combine")    
            ch2,uch = obMenu(obs,"observation 2 to combine")  
            obs[obs.columns[ch1]+"_"+obs.columns[ch2]] = obs.iloc[:,ch1]+"_"+obs.iloc[:,ch2]
        except:
            break
    
    obs = obs.astype(str)
    return(obs)
        
        
    
    
def load(bl1,bl2,bl3,path = "none"):
    if path == "none":
        path = SAVEFOLDER
    print(path)
    #"C:/Users/youm/.spyder-py3/src"
    while True:
        path = navigate(path,text="select dataframe to load")
        if path is list:
            print("please select specific file")
            continue
        if not os.path.isdir(path):
            break

    name = "_".join(path.split("_")[:-1])
    #print(name)
    name=name.split("/")[-1]+"_"
    #print(name,path)
    path = "/".join(path.split("/")[:-1])
    for file in os.listdir(path):
        fn = "_".join(file.split("_")[:-1])
        fn=fn.split("/")[-1]+"_"
        if fn == name:
            print(file)
            if "dfxy" in file:
                dfxy = pd.read_csv(path+"/"+file,index_col=0)
            elif "df" in file:
                df = pd.read_csv(path+"/"+file,index_col=0)
            elif "obs" in file:
                obs = pd.read_csv(path+"/"+file,index_col=0)    
    return(df,obs,dfxy)


def editObs(df,obs,dfxy):
    ch = input("re-organize all obs? (1/y) :")
    if ch == "1" or ch == "y":
        df,obs,dfxy=makeObs(pd.concat([df,obs,dfxy]))
    obs=unpackObs(obs)
    return(df,obs,dfxy)

def makeObs(df):
    df,dfxy = splitDF(df,"X and Y coordinate columns")
    df,obs = splitDF(df,"Observation columns")
    df = makeDtype(df,dtype=float)
    obs = makeDtype(obs,dtype=str)
    dfxy = makeDtype(dfxy,dtype=float)
    return(df,obs,dfxy)


def makeDtype(df,dtype=str):
    print("\n",dtype)
    for i in range(df.shape[1]):
        try:
            df.iloc[:,i] = df.iloc[:,i].astype(dtype)
        except:
            u = df.iloc[:,i].unique()
            if len(u) < 20:
                new = np.arange(len(u)).astype(dtype)
                print(u,new)
                for j,un in enumerate(u):
                    key = df.iloc[:,i] == un
                    df.iloc[:,i].loc[key] = new[j]
                df.iloc[:,i] = df.iloc[:,i].astype(dtype)
            else:
                df.iloc[:,i] = np.zeros(df.shape[0]).astype(dtype)                
    return(df)


def splitDF(df,titleStr="new"):
    print(titleStr)
    for i,col in enumerate(df.columns):
        print(i,col)
    newDF = []
    dropList = []
    while True:
        try:
            ch = input("column to split off into "+titleStr+": ")
            try:
                ch = int(ch)
            except:
                try:
                    ch = eval(ch)
                except:
                    chl = ch.split(":")
                    ch = range(int(chl[0]),int(chl[1]))
        except ValueError:
            break
        print(ch)
        newDF.append(df.loc[:,pd.Series(df.columns).iloc[ch]])
        nl =  df.columns[ch]
        if type(nl) == str:
            nl = [df.columns[ch]]
        else: 
            nl = list(df.columns[ch])
   
        dropList+=nl
    try:
        newDF = pd.concat(newDF,axis=1)
        df = tryDrop(df,dropList)
    except ValueError as e:
        print(e)
        newDF = pd.DataFrame(index=df.index,data=df.index)
    return(df,newDF)

def tryDrop(df,dropList):
    for colName in dropList:
        try:
            df = df.drop([colName],axis = 1)            
        except Exception as e:
            print(e,colName)
            #print(colName,'not in dataframe')
    return(df)            
 
def sortDFs(DFs,names,goodStrs):
    DFs,names = pd.Series(DFs),pd.Series(names)
    i = 0
    print(names)
    for ch in goodStrs:
        print("\n",ch)
        key = names.str.contains(ch)
        print(key)
        sDFs = DFs.loc[key]
        #print(np.array(sDFs).shape)
        if sDFs.shape[0] == 0:
            continue
        if i == 0:
            DF = pd.concat(list(sDFs),axis=0)
            print("first",DF.shape)
            i = 1          
        else:
            DF = DF.merge(pd.concat(list(sDFs),axis=0),how="outer",left_index=True,right_index=True)
            print(DF.shape)
    print("DF Final:",DF.shape)
    return(DF)
        
    
        
        
        
def getDFs():  
    global STARTFOLDER      
    DFs = []
    names = []
    goodStrs = []
    paths = []
    while True:
        path = getFile(STARTFOLDER)
        #print(type(path))
        if path == 'done':#isinstance(path,type(None)):
            break
        elif type(path) == list:
            path = path[0]
            goodStrs1 = flexMenu(title="add string that file name must contain along axis (e.g. centroid), send blank when done :")
            goodStrs += goodStrs1
            for f in os.listdir(path):
                cond=0
                for gs in goodStrs1:
                    if gs in f:
                        cond=1  
                if ".csv" in f and cond!=0:
                    try:
                        ndf = pd.read_csv(path+"/"+f,index_col=0)
                        #print(ndf.shape[1],"shape1")
                        if ndf.shape[1] < 2:
                           ndf = pd.read_csv(path+"/"+f,dtype=object,sep=" ") 
                           print(ndf.shape,"shape!")
                        DFs.append(ndf)
                        names.append(f)  
                    except:
                        print("error processing",f)
            return(DFs,names,goodStrs)   
        else:
            #goodStrs.append(path.split("/")[-1])
            paths.append(path)      
    for path in paths:
        print(path)
        names.append(path.split("/")[-1]) 
        ndf = pd.read_csv(path,index_col=0)
        if ndf.shape[1] < 2:
           ndf = pd.read_csv(path,dtype=object,sep=" ") 
           print(ndf.shape,"shape!")
        DFs.append(ndf)
    goodStrs1 = flexMenu(title="add string that file name must contain along axis (e.g. centroid), send blank when done :")
    goodStrs += goodStrs1
    print(names,"loaded files")
    return(DFs,names,goodStrs)    


def getFile(folder):   
    path = folder
    if not os.path.isdir(path):
        print("error: ",path," is invalid path")
        path = input("manually select path")
        return(getFile(path))
    try:
        while os.path.isdir(path):
            try:
                path=navigate(path)
            except Exception as e:
                print(e,"error 2")
                path = folder
            if type(path)==list:
                return(path)
    except TypeError as e:
        print(e,"error 1, probably no connection to folder")
    return(path)
            
def navigate1(path,text="send blank to go back to parent directory, send 'all' to return entire folder, 'done' to return what's loaded"): #auto
    global COMPOS    
    folder = sorted(os.listdir(path))
    for i,thing in enumerate(folder):
        print(i,thing)
    print("\n"+text)
    #ch = input("access which number?")
    ch = COMMANDLIST[COMPOS]
    COMPOS += 1
    if ch == "":
        print("going to parent directory")
        plis = path.split("/")
        plis = plis[:-1]
        path = "/".join(plis)
    elif ch == "all":
        return([path])
    elif ch == "quit" or ch[0] == "q" or ch == "done":
        return('done')
    else:
        ch = int(ch)
        path = path+"/"+folder[ch]
    print(path)
    return(path)
        
def navigate(path,text="send blank to go back to parent directory, send 'all' to return entire folder, 'done' to return what's loaded"): #manual
    folder = sorted(os.listdir(path))
    for i,thing in enumerate(folder):
        print(i,thing)
    print("\n"+text)
    ch = input("access which number?")
    if ch == "":
        print("going to parent directory")
        plis = path.split("/")
        plis = plis[:-1]
        path = "/".join(plis)
    elif ch == "all":
        return([path])
    elif ch == "quit" or ch[0] == "q" or ch == "done":
        return('done')
    else:
        try:
            ch = int(ch)
            path = path+"/"+folder[ch]
        except:
            for i,opt in enumerate(folder):
                if ch in opt and ".csv" in opt:
                    print(i,opt)
            ch = input("access which number?")
            path = path+"/"+folder[int(ch)]
    print(path)
    return(path)

def multisave(df,obs,dfxy):
    filename = input("prefix: ")
    ch,uch = obMenu(obs,title="category to divide along to save as .csvs")
    if len(uch) > 50:
        print("error, trying to save more than 50 .csvs")
        return(df,obs,dfxy)
    for uc in uch:
        key = obs.iloc[:,ch] == uc
        sdf = df.loc[key,:]
        sobs = obs.loc[key,:]
        sxy = dfxy.loc[key,:]
        sdf.to_csv(filename+"_"+uc+"_df.csv")
        sobs.to_csv(filename+"_"+uc+"_obs.csv")
        sxy.to_csv(filename+"_"+uc+"_dfxy.csv")
    return(df,obs,dfxy)

def save(df,obs,dfxy):
    filename = input("filename: ")
    df.to_csv(filename+"_df.csv")
    obs.to_csv(filename+"_obs.csv")
    dfxy.to_csv(filename+"_dfxy.csv")
    return(df,obs,dfxy)


if __name__ == "__main__":    
    main()#r"Y:\ChinData\Cyclic_Analysis\20210413_AMTEC_Analysis")
    
    
    '''
                
            #df[sCol+"_combined"]=np.nan
            sDF = pd.DataFrame(index=df.index)
            for biomarker in toComb:
                if "uclei" in biomarker or "perinuc" in biomarker or "nucadj" in biomarker:  
                    sDF[biomarker] = df[biomarker]
            print(sDF.columns)
            df.loc[:,sCol+"_nuc_combined"]=sDF.max(axis=1)
            sDF = pd.DataFrame(index=df.index)
            for biomarker in toComb:
                if "uclei" in biomarker or "perinuc" in biomarker or "nucadj" in biomarker:  
                    sDF[biomarker] = df[biomarker]
            print(sDF.columns)
            df.loc[:,sCol+"_nuc_combined"]=sDF.max(axis=1)
            
            
            
            sDF1 = pd.DataFrame(index=df.index)
            for biomarker in toComb:
                if "uclei" in biomarker or "perinuc" in biomarker or "nucadj" in biomarker or "exp5" in biomarker or "ellmem" in biomarker: 
                    continue
                else:
                    sDF1[biomarker] = df[biomarker]
            print(sDF1.columns)
            df.loc[:,sCol+"_cyto_combined"]=sDF1.max(axis=1)
            
            
            
    '''
    