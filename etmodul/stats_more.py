# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 08:26:45 2019

@author: matevz
"""
from scripts.statistics import RMSE, nanp, R2 

import numpy as np
import pandas as pd
import os
import glob
def stats1(obs, compare):
    Rmse = []
    Nnp = []
    r2 = []
    for y in compare:
        rmse = round(RMSE(obs, compare[y]),2)
        Rmse.append("Data: " + y + "- RMSE= " + str(rmse))
        nnp = round(nanp(obs, compare[y]),2)
        Nnp.append("Data: " + y + "- nanp= " + str(nnp))
        r = round(R2(obs, compare[y]),2)
        r2.append("Data: " + y + "- R2= " + str(r))
    RmSe= "    ".join(Rmse)
    NNP = "     ".join(Nnp)
    r_2 = "         ".join(r2)
    
    return RmSe, NNP, r_2  

def stats(observed, compare):
    dbRMSE = []
    dbNnp = []
    for y in observed:
        Rmse = []
        Nnp = []
        for z in compare:
            rmse = RMSE(compare[z], observed[y]["ET"])
            Rmse.append("Data: " + y + "- RMSE(" + z + ")= " + str(rmse))
            nnp = nanp(compare[z], observed[y]["ET"])
            Nnp.append("Data: " + y + "- nanp(" + z + ")= " + str(nnp))
        dbRMSE.append(Rmse)
        dbNnp.append(Nnp)
    return dbRMSE, dbNnp

#def stats(observed, compare):
#    dbRMSE = pd.DataFrame()
#    dbNnp = pd.DataFrame()
#    for y in observed:
#        Rmse = []
#        Nnp = []
#        index = []
#        for z in compare:
#            rmse = RMSE(compare[z], observed[y]["ET"])
#            Rmse.append("Data: " + y + "- RMSE(" + z + ")= " + str(rmse))
#            nnp = nanp(compare[z], observed[y]["ET"])
#            Nnp.append("Data: " + y + "- nanp(" + z + ")= " + str(nnp))
#            index.append(y)
#        dbRMSE[y] = Rmse
#
#        dbNnp[y] = Nnp
#
#    return dbRMSE, dbNnp