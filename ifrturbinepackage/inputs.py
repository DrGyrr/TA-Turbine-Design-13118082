import os
import inspect
from tkinter import N
import pandas as pd
import numpy as np
from ifrturbinepackage.definitions import *



def whichcycle(k):
    global T_1,T_5,P_1,P_5,fluid,mflow
    dfcycle=pd.read_csv(os.path.join(ROOT_DIR,"Inputs\cyclelist.csv"),skiprows=1,header=0,index_col=0)
    # T_1 = dfcycle.iloc[k-1]['T_1']
    # P_1 = dfcycle.iloc[k-1]['P_1']
    # T_5 = dfcycle.iloc[k-1]['T_5']
    # P_5 = dfcycle.iloc[k-1]['P_5']
    # fluid = dfcycle.iloc[k-1]['fluid']
    if inspect.stack()[1].function == 'ComputeR1':
        mflow = dfcycle.iloc[k-1]['mflow']
    cycledict= {'T_1'   :np.float64(dfcycle.iloc[k-1]['T_1']),
                'P_1'   :np.float64(dfcycle.iloc[k-1]['P_1']),
                'T_5'   :np.float64(dfcycle.iloc[k-1]['T_5']),
                'P_5'   :np.float64(dfcycle.iloc[k-1]['P_5']),
                'fluid' :dfcycle.iloc[k-1]['fluid'],
                'mflow' :np.float64(dfcycle.iloc[k-1]['mflow'])}
    return cycledict

def whichgparamset(l):
    global Rr5r4,Rb5b4,Rb4r4,RZrr4,NR
    gparamdict=dict()
    dfgparams=pd.read_csv(os.path.join(ROOT_DIR,"Inputs\gparamslist.csv"),skiprows=1,header=0,index_col=0)
    for gparams in list(dfgparams):
        gparamdict[gparams]=np.float64(dfgparams.iloc[l-1][gparams])
        # for x in range(np.shape(gparamlist)[0]):
    #     # assign gparams ke valuenya masing-masing
    #     gparamlist[x][0]=dfgparams.iloc[l-1][gparamlist[x][1]]
    #     # print(gparamlist[x][1],gparamlist[x][0])
    
    return gparamdict

def whatrpm(m):
    global rpm
    dfrpms  = pd.read_csv(os.path.join(ROOT_DIR,"Inputs","rpmlist.csv"),skiprows=0,header=0,index_col=0)
    rpm     = np.float64(dfrpms.iloc[m-1]['rpm'])
    return rpm

def whichfitfun(z):
    global p
    n   = int(z)
    p=np.zeros(shape=(6, 6))
    dffitfun = pd.read_excel(os.path.join(ROOT_DIR,"Inputs","fittedpolycoeffs",f"{n}.xlsx"),sheet_name='COEFFS',engine='openpyxl',skiprows=1,usecols=range(1,7),header=0,index_col=False)
    for i in range(0,6):
        for j in range(0,6):
            p[i][j]=np.float64(dffitfun.iloc[i][j])
    return p


