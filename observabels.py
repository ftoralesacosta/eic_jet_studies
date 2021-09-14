import ROOT
import numpy as np
from ROOT import TGraphErrors
from ROOT import TVectorT
import matplotlib.pyplot as plt
from matplotlib import ticker
from root_np_functions import *
from plotting_functions import *

from resolutions import *

filename_14T = "DeltaR_Histograms_Jet_Callibration_1.400000T.root"
filename = "DeltaR_Histograms_Jet_Callibration_3.000000T.root"

file = ROOT.TFile(filename)
file_14T = ROOT.TFile(filename_14T)

#Constants (EIC Luminosity, Pythia Cross Section, N Events)
EIC_Luminosity = 1E10 #10 inverse femtobarn in mililbarn units
E_Jet_Pythia_CrossSection = 9.27E-5 #milibarns
N_Events_1p4Tesla = 2209703 #number of events in 1.4T NTuple
N_Events_3Tesla = 2295430 #number of events in 3.0T NTuple

#Account for cuts in MyJetAnalysis before events are written to NTuple
pre_selection_scale = 1000/993

N_Events_1p4Tesla = N_Events_1p4Tesla*pre_selection_scale
N_Events_3Tesla = N_Events_3Tesla*pre_selection_scale

#Obtain GENERATED Luminosities
Lum_1p4T = N_Events_1p4Tesla/E_Jet_Pythia_CrossSection
Lum_3T = N_Events_3Tesla/E_Jet_Pythia_CrossSection

#Ratio of Luminosities is the SCALE
Scale_1p4T = EIC_Luminosity/Lum_1p4T
Scale_3T = EIC_Luminosity/Lum_3T

print(Scale_3T)
print(Scale_1p4T)

dPhi_bins,dPhi_centers,dPhi_widths = get_th1_binning_np(file,"dPhi_e_TrueJet")
RJ_dPhi_M, RJ_dPhi_errors_M = TH1_to_numpy_wErrors(file,"sigmamin_dPhi_e_RecoJet",True,False)
RJ_dPhi_P, RJ_dPhi_errors_P = TH1_to_numpy_wErrors(file,"sigmaplus_dPhi_e_RecoJet",True,False)
#print(RJ_dPhi_M)

Tesla14_RecodPhi_M,Tesla14_RecodPhi_errors_M = TH1_to_numpy_wErrors(file_14T,"sigmamin_dPhi_e_RecoJet",True,False)
Tesla14_RecodPhi_P,Tesla14_RecodPhi_errors_P = TH1_to_numpy_wErrors(file_14T,"sigmaplus_dPhi_e_RecoJet",True,False)
#print(Tesla14_RecodPhi_P)

RJ_dPhi_M = RJ_dPhi_M/dPhi_widths
RJ_dPhi_errors_M = RJ_dPhi_errors_M/dPhi_widths
RJ_dPhi_P = RJ_dPhi_P/dPhi_widths
RJ_dPhi_errors_P = RJ_dPhi_errors_P/dPhi_widths

Tesla14_RecodPhi_M = Tesla14_RecodPhi_M/dPhi_widths
Tesla14_RecodPhi_errors_M = Tesla14_RecodPhi_errors_M/dPhi_widths
Tesla14_RecodPhi_P = Tesla14_RecodPhi_P/dPhi_widths
Tesla14_RecodPhi_errors_P = Tesla14_RecodPhi_errors_P/dPhi_widths

cool = get_colors(plt.cm.winter,4,False)
fig = plt.figure(figsize=(14,12))

#Sigma +/- plots
plt.errorbar(dPhi_centers,RJ_dPhi_M,yerr=RJ_dPhi_errors_M,
        fmt='-',color="skyblue",linewidth=3,alpha=0.8)
#label=r"$|\Delta\varphi_\mathrm{Reco}|~B = 3.0$ T -" )
plt.errorbar(dPhi_centers,RJ_dPhi_P,yerr=RJ_dPhi_errors_P,
        fmt='-',color="skyblue",linewidth=3,alpha=0.8)
#label=r"$|\Delta\varphi_\mathrm{Reco}|~B = 3.0$ T +" )

plt.errorbar(dPhi_centers,Tesla14_RecodPhi_M,yerr=Tesla14_RecodPhi_errors_M,
        fmt='-',color="palegreen",linewidth=3,alpha=0.8)
#label=r"$|\Delta\varphi_\mathrm{Reco}|~B = 1.4$ T -")

plt.errorbar(dPhi_centers,Tesla14_RecodPhi_P,yerr=Tesla14_RecodPhi_errors_P,
        fmt='-',color="palegreen",linewidth=3,alpha=0.8)
#label=r"$|\Delta\varphi_\mathrm{Reco}|~B = 1.4$ T +")


#Pretty Bands    
plt.fill_between(dPhi_centers,RJ_dPhi_M,RJ_dPhi_P,color="skyblue",alpha=0.8,label="3.0 T Resolution")
plt.fill_between(dPhi_centers,Tesla14_RecodPhi_M,Tesla14_RecodPhi_P,color="palegreen",alpha=1.0,label="1.4 T Resolution")


#Original Plots
#plt.errorbar(dPhi_centers,TJ_dPhi,yerr=TJ_dPhi_errors,
#                         fmt='-',color="darkgrey",fillstyle='none',markersize=7,linewidth=3,
#                        label=r"$|\Delta\varphi_\mathrm{Truth}|~B = 3.0$ T" )
plt.errorbar(dPhi_centers,RJ_dPhi,yerr=RJ_dPhi_errors,
        fmt='--',color=cool[1],linewidth=3,alpha=0.8,
        label=r"$|\Delta\varphi_\mathrm{Reco}|~B = 3.0$ T" )

        #plt.errorbar(dPhi_centers,Tesla14_TruthdPhi,yerr=Tesla14_TruthdPhi_errors,
        #                         fmt='-',color=cool[2],fillstyle='none',markersize=7,linewidth=3,
        #                         label=r"$|\Delta\varphi_\mathrm{Truth}|~B = 1.4$ T" )
        #label=r"$|\varphi_{truth}^{jet} - \varphi^{e}-\pi|$")
        plt.errorbar(dPhi_centers,Tesla14_RecodPhi,yerr=Tesla14_RecodPhi_errors,
                fmt='--',color="limegreen",linewidth=3,alpha=0.8,
                label=r"$|\Delta\varphi_\mathrm{Reco}|~B = 1.4$ T")
                #label=r"$|\varphi_{reco}^{jet} - \varphi^{e}-\pi|$")

plt.text(.01,0.1,r'$L_\mathrm{EIC} = 10\ \mathrm{fb}^{-1}$',ha="left",va="top",size=25,alpha=1.)

plt.tight_layout()
plt.xlim(0,.5)
plt.ylabel("Normalized Counts $\mathrm{d}N/\mathrm{d}\Delta\varphi$",fontsize=25,y=0.5)
plt.xlabel(r"$|\Delta\varphi|=|\varphi^{jet} - \varphi^{e}-\pi|$",fontsize=25,x=0.5)
plt.tick_params(which='both',direction='in',right=True,top=True,bottom=True,length=10,labelsize=20)
plt.legend(fontsize=25,loc='upper right')
plt.savefig("DeltaR_Check_azimuthal_correlations.pdf")
