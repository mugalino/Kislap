# author: mugalino 2022 -- light_curves.py
# This is a code that calculates the nebular phase light curve of type Ia supernovae
# using the model developed by Valenti et al (2008a)

# Computes for the luminosity contribution of key isotopes produced in a type Ia explosion
# Rates were obtained from Seitenzahl et al 2014 and Jacobson Galan et al 2021
# Model was obtained from Jacobson Galan et al 2021, where they consistently consider contributions due to
# gamma ray and positrons.

# Imports packages for calculations and plotting
import numpy as np
import matplotlib.pyplot as plt

# 11fe photometric data from Tucker et al (2021)
times11fe, lum11fe = np.loadtxt('luminosity_data11fe.dat', unpack=True)[:,0], np.loadtxt('luminosity_data11fe.dat', unpack=True)[:,1]
# 19ehk photometric data from Jacobson Galan et al (2021)
times19ehk, lum19ehk = np.loadtxt('luminosity_data19ehk.dat', unpack=True)[:,0], np.loadtxt('luminosity_data19ehk.dat', unpack=True)[:,1]

times_obs = {"11fe" : times11fe, "19ehk" : times19ehk}
lum_obs = {"11fe" : lum11fe, "19ehk" : lum19ehk}
label_obs = {"11fe" : "$L_\mathrm{bol}$ (Tucker et al 2021)", "19ehk" : "$L_\mathrm{bol}$ (Jacobson Galan et al 2021)"}

# This function converts the decay energy rates from keV to specific rates erg s-1 g-1

def conversion_to_erg (energy_decay_rate, decay_time, mass_number): 
    return (decay_time) * (1 / 86400) * (energy_decay_rate) * (1.60e-9) * (1 / mass_number) * 6.022e23

# Computes the luminosity contribution of a particular isotope for an event at some time
# Either considers positron kinetic energy or not

def luminosity(time, isotope, event, pos = True): # total luminosity (erg s-1, input time in days)
    """
    This gives the total luminosity due to the following sources:
    1. Ni-56 decay to Co-56 
    2. Co-56 decay to Fe-56, formation of gamma rays
    3. Production of gamma-rays due to positron annihilation
    4. Kinetic energy of positrons
    5. Co-57 decay to Fe-57
    _. Other decay mechanisms (defined using Bateman equation) Ti-44, Fe-55
    """
    
    # specific energy decay rates (erg s-1 g-1) JG2021
    qNi56 = 3.9e10
    qCo56 = 6.8e9
    qNi57 = 8.9e6
    qTi44 = conversion_to_erg(3.222e-5, 11.351, 44)
    qFe55 = conversion_to_erg(6.916e-4, 3.973, 55)
    
    # decay timescales (days) JG2021
    tauNi56 = 8.77
    tauCo56 = 111.3
    tauCo57 = 392.11
    tauFe55 = 1445.92
    tauTi44 = 31036.62
    
    
    if isotope == "ni56":
        massNi56 = abundances[event][isotope] * 1.989e+33
        L = massNi56 * qNi56 * np.exp(- time / tauNi56)

    elif isotope == "co56":
        massNi56 = abundances[event][isotope] * 1.989e+33
        E = massNi56 * qCo56 * (np.exp(- time / tauCo56) - np.exp(- time / tauNi56)) 
        L = 0.81 * E * (1.0 - np.exp(-(tgamma[event][isotope] / time)**2))

        if pos:
            L += 0.164 * E * (1 - np.exp(-(tgamma[event][isotope] / time)**2)) * (1 - np.exp(-(tpos[isotope] / time)**2)) # Contribution from gamma ray annihilation
            #L += 0.036 * E * (1 - np.exp(-(tpos[isotope] / time)**2)) # Kinetic energy contribution from positrons
        else:
            L += 0.19 * E * (1 - np.exp(-(tgamma[event][isotope] / time)**2)) * (1 - np.exp(-(tpos[isotope] / time)**2))

    elif isotope == "co57":
        massCo57 = abundances[event][isotope] * 1.989e+33
        L = massCo57 * qNi57 * np.exp(- time / tauCo57)

    elif isotope == "fe55":
        massFe55 = abundances[event][isotope] * 1.989e+33
        L = massFe55 * qFe55 * np.exp(- time / tauFe55)

    elif isotope == "ti44":
        number = len(abundances[event][isotope])
        L = []

        for i in range(number):
            massTi44 = np.array(abundances[event][isotope][i]) * 1.989e+33
            L.append(massTi44 * qTi44 * np.exp(- time / tauTi44))

    return np.array(L)

# Just a simple function to return labels for the isotopes included
# in the light curve calculation, returns color and style as well
# Removes the need for a separate array for style inside the function
# Makes everything uniform too

def labels_colors(event, isotope, range_number):
    if isotope == "ni56":
        return ["$^{56}\mathrm{Ni}$ = %.2e $M_\odot$" % abundances[event][isotope], "plum", "-"]
    elif isotope == "co56":
        return ["$^{56}\mathrm{Co}$ = %.2e $M_\odot$" % abundances[event][isotope], "brown", "-."]
    elif isotope == "co57":
        return ["$^{57}\mathrm{Co}$ = %.2e $M_\odot$" % abundances[event][isotope], "purple", "-."]
    elif isotope == "fe55":
        return ["$^{55}\mathrm{Fe}$ = %.2e $M_\odot$" % abundances[event][isotope], "orange", "-."]
    elif isotope == "ti44":
        return [None, ["red", "blue", "green", "pink"][range_number], "dotted"]
    elif isotope == None:
        return [None, ["red", "blue", "green", "pink"][range_number], "-"]

def curve_plotter(time_interval, isotope_list, event, file_name):
    lum_list = []
    lum_sum = []

    for isotope in isotope_list:
        if isotope != "ti44":
            lum_list.append([isotope, luminosity(time_interval, isotope, event, pos = True)])
        elif isotope == "ti44":
            total = np.sum(np.array(lum_list),0)
            for lum44 in luminosity(time_interval, isotope, event, pos = True):
                lum_sum.append([None, lum44 + total[1]])
                lum_list.append([isotope,lum44])
    
    lum_list += lum_sum
   
    range_number = 0
    count = 0

    for i, array in enumerate(lum_list):
        label_name = labels_colors(event, array[0], range_number)[0]
        line_color = labels_colors(event, array[0], range_number)[1]
        line_style = labels_colors(event, array[0], range_number)[2]

        plt.plot(time_interval, array[1], label = label_name, color = line_color, linestyle = line_style)
        
        if (array[0] == "ti44" or array[0] == None) and (count%2 != 0):
            range_number += 1
            count += 1
        else:
            count +=1
            
        if (array[0] == None) and (count%2 != 0):
            j = i - ([i[0] for i in lum_list].count("ti44") - 1)
            y1 = lum_list[i][1]
            y2 = lum_list[j][1]
            plt.fill_between(time_interval, y1, y2, where = y2 > y1, facecolor = line_color, alpha = 0.6)

        if range_number == len([i[0] for i in lum_list if (i[0] != "ti44") and (i[0] != None)]) -1:
            range_number = 0
            
    plt.scatter(times_obs[event], lum_obs[event], c='red', marker= "*", label = label_obs[event])
    plt.yscale("log")
    plt.ylim(1.e30, 1.e40)
    plt.legend()

# Rates and timescales related to radioactive decay, isotopic abundances in different events
# tgamma - gamma ray annihilation
# tpos - positron trapping 

tgamma_11fe = { \
          "ni56" : 0.0,  \
          "co56" : 35.0, \
          "co57" : 160.0,\
          "fe55" : 0.0,\
          "ti44" : 0.0}

tgamma_19ehk = { \
          "ni56" : 0.0,  \
          "co56" : 59.0, \
          "co57" : 0.0,\
          "fe55" : 0.0,\
          "ti44" : 0.0}

tpos = { \
          "ni56" : 0.0,  \
          "co56" : 1200.0, \
          "co57" : 0.0,\
          "fe55" : 0.0,\
          "ti44" : 0.0}

abundances_11fe = {\
             "ni56" : 0.5,\
             "co56" : 0.5,\
             "co57" : 0.019 * 0.5,\
             "fe55" : 0.019 * 0.5 * 0.245,\
             "ti44" : [1.e-4, 1.e-3, 1.15e-6, 4.25e-5, 3.27e-5, 1.4e-3]}
             # ti44 : Roy 2022, Pakmor 2022, Leung/Nomoto 2018 (2), Leung/Nomoto 2020 (2)

abundances_19ehk = {\
             "ni56" : 3.1e-2,\
             "co56" : 2.8e-2,\
             "co57" : 8.3e-4,\
             "fe55" : None,\
             "ti44" : [7.34e-3, 1.46e-2]}
             # ni56 : photospheric phase fitting Jacobson Galan
             # co56 : nebular phase fitting Jacobson Galan
             # ti44 : Zenati et al 2022 (2)

abundances = {"11fe" : abundances_11fe, "19ehk" : abundances_19ehk}
tgamma = {"11fe" : tgamma_11fe, "19ehk" : tgamma_19ehk}

