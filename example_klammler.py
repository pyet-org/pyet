#compare
import pandas as pd
from etmodule.et0 import et0, et0_daily 
#-----------------------------------------------------------------------------
#Type of solar radiation input (global/net)
radiation = "globalr"
elevation = 266
latitude = 46.7
longitude = 15.56
#10min
min10 = pd.read_csv("data\\meteo_10min.csv", index_col = 0, parse_dates = True)
#radiation in MJ/m2s
#-----------------------------------------------------------------------------
#hourly

#data_use["global"][data_use["global"]<0] = 0
data_h_mean = min10.resample("h").mean()
data_h_sum = min10.resample("h").sum()
data_h_max = min10.resample("h").max()
data_h_min = min10.resample("h").min()
data_h = pd.DataFrame(data_h_max["T"].values,data_h_mean.index, columns=['T_max'])
data_h["T_min"] = data_h_min["T"]
#data_h["rh_min"] = data_h_min["rh"]
#data_h["rh_max"] = data_h_max["rh"]
data_h["wind"] = data_h_mean["wind"]
data_h["solar"] = data_h_sum["global"] * 0.0006
data_h["T"] = data_h_mean["T"]
data_h["rh"] = data_h_mean["rh"]
data_h = data_h.drop(columns=['T_min', 'T_max'])


fao_hour = et0(data_h, radiation, "hour", elevation, latitude,longitude)
#-----------------------------------------------------------------------------
#daily

#test daily


data_day_mean = min10.resample("d").mean()
data_day_sum = min10.resample("d").sum()
data_day_max = min10.resample("d").max()
data_day_min = min10.resample("d").min()
data_day = pd.DataFrame(data_day_max["T"].values,data_day_mean.index, columns=['T_max'])
data_day["T_min"] = data_day_min["T"]
data_day["rh_min"] = data_day_min["rh"]
data_day["rh_max"] = data_day_max["rh"]
data_day["wind"] = data_day_mean["wind"]
data_day["solar"] = data_day_sum["global"] * 0.0006   
#data_day["T"] = data_day_mean["T"]
data_day["T"] = (data_day["T_min"] +data_day["T_max"])/2

fao_daily = et0_daily(data_day, radiation, "day", elevation, latitude,longitude)



#-----------------------------------------------------------------------------
#gernot data
start = "2006-06-01 00:00:00"
end = "2008-03-31 23:50:00"
gdata = pd.read_csv("data\\evap_fao.csv", index_col = 0, parse_dates = True)
g_use = gdata[start:end]

g_use.plot()
fao_daily["ET"].plot()


fao_daily["ET_sum"] = fao_daily["ET"].cumsum()
g_use["ET_sum"] = g_use["Data"].cumsum()
g_use["ET_sum"].plot()
fao_daily["ET_sum"].plot()
