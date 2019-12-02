import pandas as pd
import numpy as np
def et_loop(glob, rad, svp, Cn, G_d, psi, Cd_d, Cd_n, G_n, es, ea,  wind, T): 
    numerator2 = psi * Cn /(T + 273) * wind * (es - ea)
    if glob > 0:
        numerator1 = 0.408 * svp * (rad - G_d) 
        denominator = svp + psi * (1 + Cd_d * wind)            
    else:
        numerator1 = 0.408 * svp * (rad - G_n) 
        denominator = svp + psi * (1 + Cd_n * wind)
    fao = (numerator1 + numerator2)/denominator
    
    return fao

def get_omega_s(date, J, lm , lz):
    hour = pd.to_numeric(date.strftime('%H'))+0.5
    b = 2*np.pi *(J-81)/364
    Sc = 0.1645*np.sin(2*b)-0.1255*np.cos(b) - 0.025*np.sin(b)
    omega_s = np.pi/12*((hour+0.06667*(lm-lz)+Sc)-12)
    return omega_s

def get_fcd(rs, solar_rat, beta):
    night_solar = 0.8
    fcd = np.where(rs<0, 1.35 * solar_rat - 0.35,night_solar)
    fcd = np.maximum(fcd, 0.05)
    fcd = np.minimum(fcd, 1.0)
    return fcd
def get_sol_rat(globalr, r_so):            
    solar_rat = np.where(r_so>0.05,globalr /r_so,0.8)
    solar_rat = globalr/r_so
    solar_rat = np.maximum(solar_rat, 0.3)
    solar_rat = np.minimum(solar_rat, 1.0)
    return solar_rat

def get_Ra(omegas1, omegas2, dr, sol_dec, lat):
    Gsc = 0.082
    Ra = 12*60/ np.pi * Gsc * dr * ((omegas2-omegas1) * np.sin(sol_dec) * np.sin(lat) + np.cos(sol_dec) *
                                              np.cos(lat) * (np.sin(omegas2)-np.sin(omegas1)))
    return Ra
def get_extraterrestrial_radiation(lista_copy, time, lat, longitude):
        """
        Calculates the solar radiation we would receive if there were no
        atmosphere. This is a function of date, time and location.

        If adatetime is a datetime object, it merely returns the
        extraterrestrial radiation R_a; if it is a date object, it returns a
        tuple, (R_a, N), where N is the daylight hours.
        """
        import pandas as pd
        import numpy as np
        import datetime as dt
        j = pd.to_numeric(lista_copy.index.strftime('%j'))

        # Inverse relative distance Earth-Sun, eq. 23, p. 46.
        dr = 1 + 0.033 * np.cos(2 * np.pi * j / 365)

        # Solar declination, eq. 24, p. 46.
        decl = 0.409 * np.sin(2 * np.pi * j / 365 - 1.39)

        if time == "day":  # Daily?
            phi = lat / 180.0 * np.pi
            omega_s = np.arccos(-np.tan(phi) * np.tan(decl))  # Eq. 25 p. 46

            r_a = (
                24
                * 60
                / np.pi
                * 0.0820
                * dr
                * (
                    omega_s * np.sin(phi) * np.sin(decl)
                    + np.cos(phi) * np.cos(decl) * np.sin(omega_s)
                )
            )  # Eq. 21 p. 46
#            n = 24 / np.pi * omega_s  # Eq. 34 p. 48
            return r_a
        else:

            #longitude
            lm = 360 - longitude
            #center of time zone
            lz = 360 - 15.5
            # Seasonal correction for solar time, eq. 32, p. 48.
            b = 2 * np.pi * (j - 81) / 364
            sc = 0.1645 * np.sin(2 * b) - 0.1255 * np.cos(b) - 0.025 * np.sin(b)
            
            # Longitude at the centre of the local time zone
#            utc_offset = adatetime.utcoffset()
#            utc_offset_hours = utc_offset.days * 24 + utc_offset.seconds / 3600.0
#            lz = -utc_offset_hours * 15
    
            # Solar time angle at midpoint of the time period, eq. 31, p. 48.
            t = pd.to_numeric(lista_copy.index.strftime('%H'))+0.5
            omega = np.pi / 12 * ((t + 0.06667 * (lz - lm) + sc) - 12)
    
            # Solar time angles at beginning and end of the period, eqs. 29 and 30,
            # p. 48.
            t1 = 3600 / 3600.0
            omega1 = omega - np.pi * t1 / 24
            omega2 = omega + np.pi * t1 / 24
    
            # Result: eq. 28, p. 47.
            phi = lat / 180.0 * np.pi
            omega_s = np.arccos(-np.tan(phi) * np.tan(decl))
            r_a = (12 * 60 / np.pi * 0.0820 * dr * 
                   ((omega2 - omega1) * np.sin(phi) * np.sin(decl)
                    + np.cos(phi) * np.cos(decl) * (np.sin(omega2) - np.sin(omega1))))
            return r_a
        
def et0(lista, solar, time, elev, latitude,longitude):
    import pandas as pd
    import numpy as np
    lista_copy = lista.copy()
    #air temperature [°C]
    data = pd.DataFrame(lista_copy["T"], lista_copy.index, columns = ["T"])
    data["wind"] = lista_copy["wind"]
    #day of year
    J = pd.to_numeric(lista_copy.index.strftime('%j'))
    #atmospheric pressure [kPa]. Based on equation 7, page 31 in Allen et al (1998).
    P = 101.3 * ((293-0.0065 * elev)/293)**(5.26)    
    #psychrometric constant [kPa °C-1]. Based on equation 8, page 32 in Allen et al (1998).
    psi = 0.665 *10**(-3)*P
    #saturation vapour pressure at the air temperature T [kPa].
    #Based on equation 11, page 36 in Allen et al (1998).
    e0 = 0.6108 * np.exp((17.27 * data["T"]) / (data["T"] + 237.3))
    #slope of saturation vapour pressure curve at air temperature T [kPa °C-1]. 
    #Based on equation 13, page 37 in Allen et al (1998).
    data["svp"] = 4098 * e0 / (data["T"] + 237.3)**2
    if time == "day":
        print("day")
        T_min = lista_copy["T_min"]
        T_max = lista_copy["T_max"] 
        #saturation vapour pressure at daily minimum temperature [kPa]
        emin = 0.6108 * np.exp((17.27 * T_min) / (T_min + 237.3))
        #saturation vapour pressure at daily maximum temperature [kPa]
        emax = 0.6108 * np.exp((17.27 * T_max) / (T_max + 237.3))
        #mean saturation vapour pressure for a day, week...
        #Based on equation 12, page 36 in Allen et al (1998).
        data["es"] = (emax + emin) / 2
        #minimum relative humidity [%]
        rh_min = lista_copy["rh_min"]
        #maximum relative humidity [%]
        rh_max = lista_copy["rh_max"]
        #actual vapour pressure [kPa]
        #Based on equation 17, page 38 in Allen et al (1998).
        data["ea"] = ((emin * rh_max/100) + (emax * rh_min/100)) / 2
    else:
        #mean saturation vapour pressure hour...
        data["es"] = e0
        rh = lista_copy["rh"]
        #actual vapour pressure [kPa]
        #Based on equation 19, page 39 in Allen et al (1998).
        data["ea"] = (rh/100) * data["es"]
        
    if solar == "globalr":
        print("globalr")
        #inverse relative distance Earth-Sun
        #Based on equation 23, page 46 in Allen et al (1998).
        dr = 1 + 0.033 * np.cos(2 * np.pi / 365 * J)
        #solar declination [rad]
        #Based on equation 24, page 46 in Allen et al (1998).
        sol_dec = 0.409 * np.sin(2 * np.pi / 365 * J - 1.39)
        #Calculate sunset hour angle (*Ws*) from latitude and solar
        # declination.
        lat = np.pi / 180 *latitude
        #solar constant = 0.0820 MJ m-2 min-1
        Gsc = 0.082
        #extraterrestrial radiation [MJ m-2 tunit-1]
        #Based on equation 21, page 46 in Allen et al (1998).
        
        if time == "day":
            omega = np.arccos(-np.tan(lat) * np.tan(sol_dec))
            Ra = 24 * 60 / np.pi * Gsc * dr * (omega * np.sin(sol_dec) * 
                np.sin(lat) + np.cos(sol_dec) * np.cos(lat) * np.sin(omega))
        else:
            lm = 15.5
            lz = 22.5 
            omega = get_omega_s(lista_copy.index, J, lm , lz)
            omega_s = np.arccos(-np.tan(lat) * np.tan(sol_dec))

            os1 = omega -(np.pi*1/24)
            os2 = omega +(np.pi*1/24)
            os11 = np.where(os1<-omega_s, -omega_s, os1)
            os21 = np.where(os2<-omega_s, -omega_s, os2)
            os12 = np.where(os11>omega_s, omega_s, os11)
            omega_s2 = np.where(os21>omega_s, omega_s, os21)
            omega_s1 = np.where(os12>omega_s2,omega_s2, os1)
            
            Ra = get_Ra(omega_s1, omega_s2, dr, sol_dec, lat)
        #Ra = get_extraterrestrial_radiation(lista_copy,time, latitude, longitude)
        beta1 = np.sin(lat) * np.sin(sol_dec)
        beta2 = np.cos(lat)*np.cos(sol_dec)*np.cos(omega)
        beta = np.arcsin(beta1+beta2)                         
        #clear-sky solar radiation [MJ m-2 tunit-1]
        #Based on equation 37, page 51 in Allen et al (1998).
#        fk = pd.DataFrame((lista_copy["solar"] /r_so), lista_copy.index,columns = ["k1"])
        r_so = (0.75 + (2 * 10**-5) * elev) * Ra
        #net solar or shortwave radiation [MJ m-2 tunit-1]
        #Based on equation 38, page 51 in Allen et al (1998).
        Rns = (1-0.23) * lista_copy["solar"]
        
        if time == "day":
            #Based on equation 39, page 52 in Allen et al (1998).
            stef = 4.903e-09
            tmp1 = stef * ((T_max+273.2)**4 + (T_min+273.2)**4)/2
            solar_rat = lista_copy["solar"] / r_so
#        elif time== "10min":
#            solar_rat = np.where(r_so>0.05,lista_copy["solar"] / Ra,0.8)
#            solar_rat = np.maximum(solar_rat, 0.3)
#            solar_rat = np.minimum(solar_rat, 1.0)
#            stef = 2.043e-10
#            tmp1 = stef * ((data["T"]+273.32)**4)
            tmp3 = 1.35 * solar_rat - 0.35
        else:
            sol_rat = get_sol_rat(lista_copy["solar"], r_so)
            #Based on equation 53, page 74 in Allen et al (1998).
            stef = 2.043e-10
            tmp1 = stef * ((data["T"]+273.32)**4)
            tmp3 = get_fcd(lista_copy["solar"],sol_rat, beta) 
        tmp2 = 0.34 - 0.14 * (data["ea"]**0.5)

        #net outgoing longwave radiation [MJ m-2 tunit-1]
        Rnl = tmp1*tmp2*tmp3
        #net radiation
        Rn = Rns - Rnl
    else:
        Rn = lista_copy["solar"]                
    data["Rn"] = Rn
    data["globalr"] = lista_copy["solar"] 
    #determined Cn coefficient concerning calculation time step, 
    #Cd coefficient concerning surface resistance   
    #G soil heatflux density [MJ m−2Δt−1] determined based on table 4 in Klammer and Fnk (2014)
    if time == "day":
        Cn = 900
        Cd_d = 0.34
        Cd_n = 0.34
        data["G_d"] = 0
        data["G_n"] = 0
        print("tunit = day")
    else:
        data["es"] = e0  
        data["G_d"] = 0.1 * Rn
        data["G_n"] = 0.5 * Rn     
        if time == "houra":
            Cn = 37
            Cd_d = 0.34
            Cd_n = 0.34
            print("tunit = houra")
        elif time == "hourb":
            Cn = 37
            Cd_d = 0.24
            Cd_n = 0.96         
        elif time == "10min":
            Cn = 37
            Cd_d = 0.24
            Cd_n = 0.96  
            print("tunit = 10min")

    data["ET"] = data.apply(lambda row: et_loop(row["globalr"], row["Rn"], row["svp"], Cn, row["G_d"],
        psi, Cd_d, Cd_n, row["G_n"], row["es"], row["ea"], row["wind"], row["T"]), axis = 1)
    return data

def et0_daily(lista, solar, time, elev, latitude,longitude):
    import pandas as pd
    import numpy as np
    lista_copy = lista.copy()
    #air temperature [°C]
    data = pd.DataFrame(lista_copy["T"], lista_copy.index, columns = ["T"])
    data["wind"] = lista_copy["wind"]
    #day of year
    J = pd.to_numeric(lista_copy.index.strftime('%j'))
    #atmospheric pressure [kPa]. Based on equation 7, page 31 in Allen et al (1998).
    P = 101.3 * ((293-0.0065 * elev)/293)**(5.26)    
    #psychrometric constant [kPa °C-1]. Based on equation 8, page 32 in Allen et al (1998).
    psi = 0.665 *10**(-3)*P
    #saturation vapour pressure at the air temperature T [kPa].
    #Based on equation 11, page 36 in Allen et al (1998).
    e0 = 0.6108 * np.exp((17.27 * data["T"]) / (data["T"] + 237.3))
    #slope of saturation vapour pressure curve at air temperature T [kPa °C-1]. 
    #Based on equation 13, page 37 in Allen et al (1998).
    data["svp"] = 4098 * e0 / (data["T"] + 237.3)**2

    T_min = lista_copy["T_min"]
    T_max = lista_copy["T_max"] 
    #saturation vapour pressure at daily minimum temperature [kPa]
    emin = 0.6108 * np.exp((17.27 * T_min) / (T_min + 237.3))
    #saturation vapour pressure at daily maximum temperature [kPa]
    emax = 0.6108 * np.exp((17.27 * T_max) / (T_max + 237.3))
    #mean saturation vapour pressure for a day, week...
    #Based on equation 12, page 36 in Allen et al (1998).
    data["es"] = (emax + emin) / 2
    #minimum relative humidity [%]
    rh_min = lista_copy["rh_min"]
    #maximum relative humidity [%]
    rh_max = lista_copy["rh_max"]
    #actual vapour pressure [kPa]
    #Based on equation 17, page 38 in Allen et al (1998).
    data["ea"] = ((emin * rh_max/100) + (emax * rh_min/100)) / 2
        
    if solar == "globalr":
        print("globalr")
        #inverse relative distance Earth-Sun
        #Based on equation 23, page 46 in Allen et al (1998).
        dr = 1 + 0.033 * np.cos(2 * np.pi / 365 * J)
        #solar declination [rad]
        #Based on equation 24, page 46 in Allen et al (1998).
        sol_dec = 0.409 * np.sin(2 * np.pi / 365 * J - 1.39)
        #Calculate sunset hour angle (*Ws*) from latitude and solar
        # declination.
        lat = np.pi / 180 *latitude
        #solar constant = 0.0820 MJ m-2 min-1
        Gsc = 0.082
        #extraterrestrial radiation [MJ m-2 tunit-1]
        #Based on equation 21, page 46 in Allen et al (1998).
        

        omega = np.arccos(-np.tan(lat) * np.tan(sol_dec))
        Ra = 24 * 60 / np.pi * Gsc * dr * (omega * np.sin(sol_dec) * 
                np.sin(lat) + np.cos(sol_dec) * np.cos(lat) * np.sin(omega))

        #Ra = get_extraterrestrial_radiation(lista_copy,time, latitude, longitude)
        beta1 = np.sin(lat) * np.sin(sol_dec)
        beta2 = np.cos(lat)*np.cos(sol_dec)*np.cos(omega)
        beta = np.arcsin(beta1+beta2)                         
        #clear-sky solar radiation [MJ m-2 tunit-1]
        #Based on equation 37, page 51 in Allen et al (1998).
#        fk = pd.DataFrame((lista_copy["solar"] /r_so), lista_copy.index,columns = ["k1"])
        r_so = (0.75 + (2 * 10**-5) * elev) * Ra

        #net solar or shortwave radiation [MJ m-2 tunit-1]
        #Based on equation 38, page 51 in Allen et al (1998).
        Rns = (1-0.23) * lista_copy["solar"]
        
        if time == "day":
            #Based on equation 39, page 52 in Allen et al (1998).
            stef = 4.903e-09
            tmp1 = stef * ((T_max+273.2)**4 + (T_min+273.2)**4)/2
            solar_rat = lista_copy["solar"] / r_so
#        elif time== "10min":
#            solar_rat = np.where(r_so>0.05,lista_copy["solar"] / Ra,0.8)
#            solar_rat = np.maximum(solar_rat, 0.3)
#            solar_rat = np.minimum(solar_rat, 1.0)
#            stef = 2.043e-10
#            tmp1 = stef * ((data["T"]+273.32)**4)
            tmp3 = 1.35 * solar_rat - 0.35
        else:
            sol_rat = get_sol_rat(lista_copy["solar"], r_so)
            #Based on equation 53, page 74 in Allen et al (1998).
            stef = 2.043e-10
            tmp1 = stef * ((data["T"]+273.32)**4)
            tmp3 = get_fcd(lista_copy["solar"],sol_rat, beta) 
        tmp2 = 0.34 - 0.14 * (data["ea"]**0.5)

        #net outgoing longwave radiation [MJ m-2 tunit-1]
        Rnl = tmp1*tmp2*tmp3
        #net radiation
        Rn = Rns - Rnl
    else:
        Rn = lista_copy["solar"]                
    data["Rn"] = Rn
    data["globalr"] = lista_copy["solar"] 
    #determined Cn coefficient concerning calculation time step, 
    #Cd coefficient concerning surface resistance   
    #G soil heatflux density [MJ m−2Δt−1] determined based on table 4 in Klammer and Fnk (2014)
    
    Cn = 900
    Cd_d = 0.34
    Cd_n = 0.34
    G_d = 0
    numerator1 = 0.408 * data["svp"] * (data["Rn"] - G_d)
    numerator2 = psi * Cn /(data["T"] + 273) * data["wind"] * (data["es"] - data["ea"])
    denominator = data["svp"] + psi * (1 + Cd_d * data["wind"])            

    fao = (numerator1 + numerator2)/denominator
    data["ET"] = fao
    return data