# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 08:28:24 2019

@author: matevz
"""
import numpy as np
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
from scripts.stats_more import stats, stats1 

def plot_sums(obs):
    import matplotlib.pyplot as plt
    dates = (pd.date_range("2015-04-01","2015-10-31"),
             pd.date_range("2016-04-01","2016-10-31"),
             pd.date_range("2017-04-01","2017-10-31"),
             pd.date_range("2018-04-01","2018-10-31"))
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(15,9))
    fig.suptitle("Lysimeter sums", fontsize=16)
    for x,y in zip(dates,np.arange(4)):
        for z in obs:
            axes[y].plot(x, obs[z][x], label = z)
            axes[y].set_ylabel('ET [mm]')
            axes[y].set_xlabel('Time (d)')
            axes[y].grid(True)
            axes[y].legend(loc = "upper left", fontsize=10)
        
    



def plot1(obs, compare):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(15,9))
    fig.suptitle("PM vs FAO-PM vs PT", fontsize=16)
    for x,y in zip(compare,np.arange(3)):
        axes[y].plot(compare.index, compare[x], label = x)
        axes[y].plot(obs.index, obs, label = "Lysimeter data GS-2")
        axes[y].set_ylabel('PET [mm]')
        axes[y].set_xlabel('Time (d)')
        axes[y].grid(True)
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        axes[y].text(0.01, 0.98, rmse[y] + "\t" + nnp[y]+ "\t" + r2[y], transform=axes[y].transAxes, fontsize=14,
                  verticalalignment='top', bbox=props)
        axes[y].legend(loc = "upper right", fontsize=14)
#    ax_twin = axes[y].twinx()
#    ax_twin.plot(daily_swc[start:end].index, daily_swc['SWC_020 (111) [Vol.%]'][start:end],
#                 color="lightgray", label = "SWC",dashes=[6, 2])
#    ax_twin.fill_between(d, daily_swc['SWC_020 (111) [Vol.%]'][start:end], 0.2,facecolor='gray', alpha=0.1)
#    ax_twin.set_ylim(bottom=0.2, top=0.4)
    
def plot2(obs, compare):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(15,9))
    fig.suptitle("PM vs FAO-PM vs PT", fontsize=16)
    for x,y in zip(compare, np.arange(3)):
        axes.plot(compare.index, compare[x], label = x)
    axes.plot(obs.index, obs, label = "Lysimeter data GS-2")
    axes.set_ylabel('PET [mm]')
    axes.set_xlabel('Time (d)')
    axes.grid(True)
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    axes.text(0.01, 0.98, rmse + "\n" + nnp + "\n" + r2, transform=axes.transAxes, fontsize=14,
                  verticalalignment='top', bbox=props)
    axes.legend(loc = "lower right", fontsize=14)
  
def plot5(observed, compare):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(15,9))
    fig.suptitle(str(compare.name) +" vs observed", fontsize=16)
    start = ("2015-04-01","2016-04-01","2017-04-01","2018-04-01")
    end = ("2015-10-31","2016-10-31","2017-10-31","2018-10-31")
    for x,y,z in zip(start,end,np.arange(4)):
        dates= pd.date_range(x, y, freq='d')
        d = dates.to_pydatetime()
        axes[z].plot(compare[x:y].index, compare[x:y], label = str(compare.name))
        axes[z].plot(observed[x:y].index, observed["GS-2"][x:y], label = "C0T0")
        ax_twin = axes[z].twinx()
        ax_twin.plot(daily_swc[x:y].index, daily_swc['SWC_020 (111) [Vol.%]'][x:y],
                         color="lightgray", label = "SWC",dashes=[6, 2])
        ax_twin.fill_between(d, daily_swc['SWC_020 (111) [Vol.%]'][x:y], 0.2,facecolor='gray', alpha=0.1)
        axes[z].set_ylabel('ET [mm]')
        axes[z].set_xlabel('Time (d)')
        axes[z].grid(True)
        axes[z].legend(loc = "upper right", fontsize=14)  
        
def plot4(observed, compare):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(15,9))
    fig.suptitle(str(compare.name) +" vs observed", fontsize=16)
    start = ("2015-04-01","2016-04-01","2017-04-01","2018-04-01")
    end = ("2015-10-31","2016-10-31","2017-10-31","2018-10-31")
    for x,y,z in zip(start,end,np.arange(4)):
        dates= pd.date_range(x, y, freq='d')
        d = dates.to_pydatetime()
        compare1 = compare[compare.index.isin(dates)]
        observed1 = observed[observed.index.isin(dates)]
        
        compare2 = compare1.cumsum()
        observed2 = observed1["GS-2"].cumsum()
        
        axes[z].plot(compare2.index, compare2, label = str(compare.name))
        axes[z].plot(observed2.index, observed2, label = "C0T0")
        ax_twin = axes[z].twinx()
        ax_twin.plot(daily_swc[x:y].index, daily_swc['SWC_020 (111) [Vol.%]'][x:y],
                         color="lightgray", label = "SWC",dashes=[6, 2])
        ax_twin.fill_between(d, daily_swc['SWC_020 (111) [Vol.%]'][x:y], 0.2,facecolor='gray', alpha=0.1)
        axes[z].set_ylabel('ET [mm]')
        axes[z].set_xlabel('Time (d)')
        axes[z].grid(True)
        axes[z].legend(loc = "upper right", fontsize=14)            



#dates_2015 = pd.date_range("2015-04-01","2015-10-31")
#dates_2016 = pd.date_range("2016-04-01","2016-10-31")
#dates_2017 = pd.date_range("2017-04-01","2017-10-31")
#dates_2018 = pd.date_range("2018-04-01","2018-10-31")
#
#
#compare_2015 = lysimeter_daily["GS-2"][lysimeter_daily.index.isin(dates_2015)]
#compare_2016 = lysimeter_daily["GS-2"][lysimeter_daily.index.isin(dates_2016)]
#compare_2017 = lysimeter_daily["GS-2"][lysimeter_daily.index.isin(dates_2017)]
#compare_2018 = lysimeter_daily["GS-2"][lysimeter_daily.index.isin(dates_2018)]
#
#compare_2015["Fao_PM_sum"] = compare["Fao_PM"][daily2015].cumsum()
#compare_2016["Fao_PM_sum"] = lysimeter_daily["GS-2"][lysimeter_daily.index.isin(dates_2016)]
#compare_2017["Fao_PM_sum"] = lysimeter_daily["GS-2"][lysimeter_daily.index.isin(dates_2017)]
#compare_2018["Fao_PM_sum"] = lysimeter_daily["GS-2"][lysimeter_daily.index.isin(dates_2018)]

def plot_obs_pm(obs, compare):
    import matplotlib.pyplot as plt
    dates = (pd.date_range("2015-04-01","2015-10-31"),
             pd.date_range("2016-04-01","2016-10-31"),
             pd.date_range("2017-04-01","2017-10-31"),
             pd.date_range("2018-04-01","2018-10-31"))
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(15,9))
    fig.suptitle("PM vs FAO-PM vs PT", fontsize=16)
    for x,y in zip(dates,np.arange(4)):
        rmse, nnp, r2 = stats1(obs["GS-2"].loc[x], compare.loc[x])
        axes[y].plot(x, obs["GS-2"][x], label = "Lysimeter data GS-2")
        for z in compare:
            axes[y].plot(x, compare[z][x], label = z)  
            axes[y].set_ylabel('ET [mm]')
            axes[y].set_xlabel('Time (d)')
            axes[y].grid(True)       
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            axes[y].text(0.01, 0.98, rmse + "\n" + nnp+ "\n" + r2, transform=axes[y].transAxes, fontsize=7,
                  verticalalignment='top', bbox=props)            
            axes[y].legend(loc = "upper right", fontsize=8)   
def plot_kc(obs,label):
    from sklearn.linear_model import LinearRegression
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.dates import date2num
    dates = (pd.date_range("2015-04-01","2015-10-31"),
             pd.date_range("2016-04-01","2016-10-31"),
             pd.date_range("2017-04-01","2017-10-31"),
             pd.date_range("2018-04-01","2018-10-31"))
    cuts = (["2015-05-27","2015-07-29","2015-10-11"],
            ["2016-05-30","2016-07-27","2016-10-04"],
            ["2017-05-29","2017-07-25","2017-10-03"],
            ["2018-05-28","2018-07-23","2018-10-01"])
    
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(15,9))
    fig.suptitle(label, fontsize=16)
    for x,y,c in zip(dates,np.arange(4),cuts):
        for z in obs:
            axes[y].plot(x, obs[z][x], label = z)
            axes[y].set_ylabel('kc')
            axes[y].set_xlabel('Time (d)')
            axes[y].grid(True)
            axes[y].legend(loc = "upper left", fontsize=10)
#            this_data = pd.DataFrame()
#            this_data["date"]= x.map(lambda a: date2num(a))
#            print(this_data)
#            sns.regplot('Date_Ord', 'High', data=this_data.query('Date > "2017-04-01" and Date < "2017-09-01"'),
#               truncate=True, color='C{}'.format(ii))
            # Plot cutting days

            i0 = dates[y][0]
            for i in range(3):
                # Plot cutting days
                axes[y].axvline(c[i])
                #plot regression line
                this_data = pd.DataFrame()
                this_data["Date"]= pd.date_range(i0,c[i])
                this_data["Date_Ord"]=this_data["Date"].map(lambda a: date2num(a))
                this_data["kc"]= obs["kc"][pd.date_range(i0,c[i])].values
                this_data = this_data.dropna()
                from sklearn.linear_model import LinearRegression
                model = LinearRegression(fit_intercept=True)
                X=this_data["Date"]
                x1=this_data["Date_Ord"]
                y1=this_data["kc"]
                model.fit(x1[:, np.newaxis], y1)
                xfit = this_data["Date_Ord"]
                yfit = model.predict(xfit[:, np.newaxis])
                i0 = c[i]
                axes[y].plot(X, yfit);
                
def plot_obs_pm_sum(obs, compare):
    import matplotlib.pyplot as plt
    dates = (pd.date_range("2015-04-01","2015-10-31"),
             pd.date_range("2016-04-01","2016-10-31"),
             pd.date_range("2017-04-01","2017-10-31"),
             pd.date_range("2018-04-01","2018-10-31"))
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(15,9))
    fig.suptitle("PM vs FAO-PM vs PT", fontsize=16)
    for x,y in zip(dates,np.arange(4)):
        rmse, nnp, r2 = stats1(obs["GS-2"].loc[x], compare.loc[x])
        axes[y].plot(x, obs["GS-2"][x], label = "Lysimeter data GS-2")
        for z in compare:
            axes[y].plot(x, compare[z][x], label = z)  
            axes[y].set_ylabel('ET [mm]')
            axes[y].set_xlabel('Time (d)')
            axes[y].grid(True)       
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            axes[y].text(0.01, 0.98, rmse + "\n" + nnp+ "\n" + r2, transform=axes[y].transAxes, fontsize=7,
                  verticalalignment='top', bbox=props)            
            axes[y].legend(loc = "lower right", fontsize=8) 

        
def plotl(observed, compare):
    import matplotlib.pyplot as plt
    for x,z in zip(observed,np.arange(6)):
        fig, axes = plt.subplots(nrows=6, ncols=3, figsize=(15,9))
        fig.suptitle(x + "- optimal", fontsize=16)
        axes[0].plot(compare.index, compare["Fao_pm"], label = "Fao_pm")
        axes[1].plot(compare.index, compare["pm"], label = "pm")
        axes[2].plot(compare.index, compare["pt"], label = "pt")
        axes[3].plot(compare.index, compare["pm_co2"], label = "pm_co2")
        axes[4].plot(compare.index, compare["fao_pm_py"], label = "fao_pm_py")
        axes[5].plot(compare.index, compare["har_py"], label = "har_py")
        for y in range (6):
            axes[y].plot(observed[x].index, observed[x]["ET"], label = x)
            ax_twin = axes[y].twinx()
            ax_twin.plot(daily_swc[start:end].index, daily_swc['SWC_020 (111) [Vol.%]'][start:end],
                         color="lightgray", label = "SWC",dashes=[6, 2])
            ax_twin.fill_between(d, daily_swc['SWC_020 (111) [Vol.%]'][start:end], 0.2,facecolor='gray', alpha=0.1)
            ax_twin.set_ylim(bottom=0.2, top=0.4)
            axes[y].set_ylabel('ET [mm]')
            axes[y].set_xlabel('Time (d)')
            axes[y].grid(True)
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            axes[y].text(0.01, 0.98, rmse[z][y] + "\t" + nnp[z][y], transform=axes[y].transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
            axes[y].legend(loc = "upper right", fontsize=14,alpha=0)
            
def plot_co2(observed, compare):
    import matplotlib.pyplot as plt
    for x,z in zip(observed,np.arange(6)):
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(15,9))
        fig.suptitle(x + "- optimal", fontsize=16)
        axes[0].plot(compare.index, compare["Fao_pm"], label = "fao_pm")
        axes[0].plot(compare.index, compare["fao_pm_co2"], label = "fao_pm_co2")

        axes[1].plot(compare.index, compare["pm"], label = "pm")
        axes[1].plot(compare.index, compare["pm_co2"], label = "pm_co2")
        for y in range (2):
            axes[y].plot(observed[x].index, observed[x]["ET"], label = x)
            ax_twin = axes[y].twinx()
            ax_twin.plot(daily_swc[start:end].index, daily_swc['SWC_020 (111) [Vol.%]'][start:end],
                         color="lightgray", label = "SWC",dashes=[6, 2])
            ax_twin.fill_between(d, daily_swc['SWC_020 (111) [Vol.%]'][start:end], 0.2,facecolor='gray', alpha=0.1)
            ax_twin.set_ylim(bottom=0.2, top=0.4)
            axes[y].set_ylabel('ET [mm]')
            axes[y].set_xlabel('Time (d)')
            axes[y].grid(True)
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            axes[y].text(0.01, 0.98, rmse[z][y*2] + " -- " + rmse[z][y*2+1] + "\n" + \
                nnp[z][y*2] + " -- " + nnp[z][y*2+1], transform=axes[y].transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
            axes[y].legend(loc = "center right", fontsize=14)
            
def plot_t3(observed, compare):
    import matplotlib.pyplot as plt
    for x,z in zip(observed,np.arange(6)):
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(15,9))
        fig.suptitle(x + "- optimal", fontsize=16)
        axes[0].plot(compare.index, compare["Fao_pm"], label = "fao_pm")
        axes[0].plot(compare.index, compare["fao_pm_t3"], label = "fao_pm_t3")
        axes[1].plot(compare.index, compare["pm"], label = "pm")
        axes[1].plot(compare.index, compare["pm_t3"], label = "pm_t3")
        for y in range (2):
            axes[y].plot(observed[x].index, observed[x]["ET"], label = x)
            ax_twin = axes[y].twinx()
            ax_twin.plot(daily_swc[start:end].index, daily_swc['SWC_020 (111) [Vol.%]'][start:end],
                         color="lightgray", label = "SWC",dashes=[6, 2])
            ax_twin.fill_between(d, daily_swc['SWC_020 (111) [Vol.%]'][start:end], 0.2,facecolor='gray', alpha=0.1)
            ax_twin.set_ylim(bottom=0.2, top=0.4)
            axes[y].set_ylabel('ET [mm]')
            axes[y].set_xlabel('Time (d)')
            axes[y].grid(True)
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            axes[y].text(0.01, 0.98, rmse[z][y*2] + " -- " + rmse[z][y*2+1] + "\n" + \
                nnp[z][y*2] + " -- " + nnp[z][y*2+1], transform=axes[y].transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
            axes[y].legend(loc = "center right", fontsize=14)

def plot_rs(observed, compare):
    import matplotlib.pyplot as plt
    for x,z in zip(observed,np.arange(6)):
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(15,9))
        fig.suptitle(x + "- optimal", fontsize=16)
        axes.plot(compare.index, compare["pm_co2"], label = "pm_co2")
        axes.plot(compare.index, compare["pm_70"], label = "pm_70")
        axes.plot(observed[x].index, observed[x]["ET"], label = x)
        ax_twin = axes.twinx()
        ax_twin.plot(daily_swc[start:end].index, daily_swc['SWC_020 (111) [Vol.%]'][start:end],
                         color="lightgray", label = "SWC",dashes=[6, 2])
        ax_twin.fill_between(d, daily_swc['SWC_020 (111) [Vol.%]'][start:end], 0.2,facecolor='gray', alpha=0.1)
        ax_twin.set_ylim(bottom=0.2, top=0.4)
        axes.set_ylabel('ET [mm]')
        axes.set_xlabel('Time (d)')
        axes.grid(True)
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
        axes.text(0.01, 0.98, rmse[z][1] + " -- " + nnp[z][1] + "\n" + 
                  rmse[z][2] + " -- " + nnp[z][2] + "\n"+ 
                  rmse[z][3] + " -- " + nnp[z][3], 
                  transform=axes.transAxes, fontsize=14,
                verticalalignment='top', bbox=props)
        axes.legend(loc = "center right", fontsize=14)
