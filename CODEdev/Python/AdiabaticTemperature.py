#!/usr/bin/env python
#
# This program calculates adiabatic temperature with variable equivalence ratio
# by using Cantera.
#
import numpy as np
from Cantera import *
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# phi variables
phimin		= 0.05
phimax		= 3.0	# if phiRes is 1, phimax will be picked as a referenced phi
phiRes		= 10	# Phi variable resolution (# of points)

# include stoichimetry mixture calculation
Stoicswitch	= 1

# plot switch
PltSwt_TvsPhi	= 0
PltSwt_TvsZ	= 1

# file write switch
fwr_TvsZ	= 1
if (fwr_TvsZ == 1):
    fwrTZ = open('Tad_vs_Z.dat','w')
    #fwrTZ.write('Z vs Tad\n')

#Chemical Mechanism
rxnmech		= 'KERO_BFER.cti'
gas		 = IdealGasMix(rxnmech)

#
# LDI parameter values: T and P will not change fuel mass flow rate
# But it is needed to set a state for Cantera running!
#
Tchamber        = 750           # [K]
Pchamber        = 1.2269e+6     # [Pa]
#
# Inlet air and fuel composition
#
YfInlet		= 1.0		# Pure fuel
YOxInlet	= 0.233		# Pure AIR

#
# Stoichimetric condition
# MSR: Mass Stoichiometric Ratio = (Yo/Yf)st
#
compStoic	= 'KERO:1.0, O2:15.0, N2:56.4' # Stoichimetry composition
gas.setState_TPX(Tchamber,Pchamber,compStoic)
MSR		= gas.massFraction(1)/gas.massFraction(0)
print "######################################"
print "Mass Stoichiometric Ratio:",MSR
print "######################################"

#
# Set equivalence ratio variable
#
if (phiRes != 1):
    dphi	= abs(phimax - phimin) / (phiRes - 1)
    phi		= phimin	# phi initialization
else:
    dphi	= 0
    phi		= phimax

#
# Creates array for storing data
#
phiArray	= np.array([])
TadArray	= np.array([])
ZArray		= np.array([])

#
# Initialize phiArray and ZArray
#
for i in range(phiRes):
    phi		= phi + dphi * i  
    if (Stoicswitch == 1 and phi != 1):
        phiArray	= np.append(phiArray,phi)
    elif (Stoicswitch != 1):
        phiArray	= np.append(phiArray,phi)

if (Stoicswitch == 1):
    phiArray	= np.append(phiArray,1)

sortindex	= np.argsort(phiArray)
phiArray	= phiArray[sortindex]

print phiArray       

#
# Loop for variable phi (Equivalence ratio)
#
for i in range(len(phiArray)):
    phi		= phiArray[i]
    compMix	= 'KERO:%s, O2:15.0, N2:56.4' % (phi)
    print "========================="
    print "working on phi = %4.4s" % (phi)
    gas.setState_TPX(Tchamber,Pchamber,compMix)
    #
    # mixture fraction (Z): fuel/oxidizer ratio 
    # see Poinsot's book p.85
    #
    Yf		= gas.massFraction(0)	# fuel mass fraction in mixture
    YO2		= gas.massFraction(1)	# O2 mass fraction in mixture
    Z		= ((MSR * Yf) - YO2 + YOxInlet) / ((MSR * YfInlet) + YOxInlet)
    ZArray	= np.append(ZArray,Z)
    print "Mixture fraction(Z):",Z

    #for n in range(gas.nSpecies()):
    #    print "Y_",gas.speciesName(n),":",gas.massFraction(n)
    gas.equilibrate('HP', solver=1)
    Tad         = gas.temperature()
    TadArray	= np.append(TadArray,Tad)
    print "Tad = ",Tad,"[K]"

    if (fwr_TvsZ == 1):
        fwrTZ.write('%6.6s\t%7.7s\n' % (Z,Tad))
print "========================="
print "Collecting data DONE"
print "========================="

if (fwr_TvsZ == 1):
    fwrTZ.close()

#
# Plotting data
#
if( PltSwt_TvsPhi == 1 ):
    MinX	= min(phiArray)
    MaxX	= max(phiArray)
    MinY	= 0.98*min(TadArray)
    MaxY	= 1.1*max(TadArray)
    
    plt.axis([MinX,MaxX,MinY,MaxY])
    plt.plot(phiArray, TadArray, 'k--o', linewidth=2)
    plt.xlabel('Equivalence ratio', fontsize=15)
    plt.ylabel('Temperature [K]', fontsize=15)
    plt.grid(True)

    ax		= plt.gca ()
    ylabels	= plt.getp(ax, 'yticklabels')
    plt.setp(ylabels, fontsize=10)
    #plt.legend(
    #           loc='upper right',
    #           borderpad=0.25,
    #           handletextpad=0.25,
    #           borderaxespad=0.25,
    #           labelspacing=0.0,
    #           handlelength=2.0,
    #           numpoints=1)
    #legendText = plt.gca().get_legend().get_texts()
    #plt.setp(legendText, fontsize=15)
    #legend = plt.gca().get_legend()
    #legend.draw_frame(False)
    pltFile = 'Tad_vs_phi.png'
    plt.savefig(pltFile, format='png')
    plt.close()

if( PltSwt_TvsZ == 1 ):
    MinX        = min(ZArray)
    MaxX        = max(ZArray)
    MinY        = 0.98*min(TadArray)
    MaxY        = 1.1*max(TadArray)

    plt.axis([MinX,MaxX,MinY,MaxY])
    plt.plot(ZArray, TadArray, 'k--o', linewidth=2)
    plt.annotate('Adiabatic temperature',  xy=(0.01,2600), fontsize=15, color='r')
    plt.xlabel('Mixture fraction (Z)', fontsize=15)
    plt.ylabel('Temperature [K]', fontsize=15)
    plt.grid(True)

    ax          = plt.gca ()
    ylabels     = plt.getp(ax, 'yticklabels')
    plt.setp(ylabels, fontsize=10)
    #plt.legend(
    #           loc='upper right',
    #           borderpad=0.25,
    #           handletextpad=0.25,
    #           borderaxespad=0.25,
    #           labelspacing=0.0,
    #           handlelength=2.0,
    #           numpoints=1)
    #legendText = plt.gca().get_legend().get_texts()
    #plt.setp(legendText, fontsize=15)
    #legend = plt.gca().get_legend()
    #legend.draw_frame(False)
    pltFile = 'Tad_vs_Z.png'
    plt.savefig(pltFile, format='png')
    plt.close()

