import ROOT
import numpy as np
from ROOT import TGraphErrors
from ROOT import TVectorT
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.optimize import curve_fit

# import warnings
# warnings.filterwarnings("ignore")

from root_np_functions import *
from plotting_functions import *

############################################ For Gaussian Fitting ############################################

def gauss(x, a, sigma):
    return a * np.exp(-(x) ** 2 / (2 * sigma ** 2))

def double_gauss(x,a,sigma0,b,sigma1):
    return (a * np.exp(-(x) ** 2 / (2 * sigma0 ** 2)) +
            b * np.exp(-(x) ** 2 / (2 * sigma1 ** 2)) )

    #Important: scipy returns the covariance matrix of the fit, from
    #which the standard error of each parameter can be obtained. If
    #a parameter is extremeley small (i.e. close to zero), it should
    #be omitted, due to float point overflow/underflow. If not 
    #handled properly, it could result in the matrix with INF vals.



def fit_range(fit_range,arr):
    #resizes the shape of the array to the range used for gauss fit

    fit_min_value = fit_range[0]
    fit_max_value = fit_range[1]
    min_array = np.absolute(arr-fit_min_value)
    max_array = np.absolute(arr-fit_max_value)

    min_index = np.argmin(min_array)
    max_index = np.argmin(max_array)+1

    return min_index,max_index

############################ Initialize Filename, binning, B Fields, and strings ###############################

#class resolutions():
#    def __init__(
#            self,tweak_string,fit_type,args=settings
#            ):
                
#        self.tweak_string = tweak_string
#        self.fit_type = fit_type
#        self.args = args

#    import resolutions
#    resolutions.cook_meat(minutes=10) #fuck the order if you call the name of the arg correctly (name defined in function definition)

#    def get_th1f_dictionary(self):

#     for B in B_Fields:
#         file = ROOT.TFile.Open(file_dict["%1.1f"%(B)])
    
#         for s in res_IDs:
             
#             for eta in range(N_eta):
    
#                  for p in range(N_mom):
    
#                     h1=file.Get(settings[s+"dir"]+"/%s_et_%i_p_%i_bin"%(s,eta,p))
#                     identifier = ("%s_B_%1.1f_eta_%i_p_%i"%(s,B,eta,p))
    
#                     h1_dictionary[identifier+"_vals"],h1_dictionary[identifier+"_errors"] = TH1_to_numpy_wErrors(h1,False,True)
#                     h1_dictionary[identifier+"_avg"] = h1.GetMean()
#                     h1_dictionary[identifier+"_sigma"] = h1.GetStdDev()
#                     h1_dictionary[identifier+"_sigma_error"] = h1.GetStdDev()/np.sqrt(2*h1.GetEntries() - 2)
#                     # error on stdev = sigma/sqrt{2*(n-1)}
     
#                      if (p==0 and eta==0):
#                         h1_dictionary[s+"edges"],h1_dictionary[s+"centers"],h1_dictionary[s+"widths"] = get_th1_binning_np(h1)
#                     #print(h1.GetMean())
    
#     for s,title in zip(res_IDs,titles):
#         h1_dictionary[s+"title"] = title
#     np.save("np_arrays/%s_indv_histos.npy"%(tweak_string),h1_dictionary)

#Main Parameters
B_Fields = [1.4,3.0]
res_IDs = ["h1_dpp","h1_dph","h1_dth","h1_eDelta_dph"] #Types of Resolutions
fit_type = "double"

#Files
tweak_string = "DeltaR_"
# tweak_string = ""
file_string = tweak_string+"histograms_reco_NoCuts_output_mom_res_sigma_eta_5_p_6_B"
file_dict = {}

for B in B_Fields:
    file_dict["%1.1f"%(B)] = "%s_%1.1f.root"%(file_string,B)

#Settings for different types of resolutions
settings = {}
TDirectories = ["dpp_histos","dph_histos","dth_histos","eD_dph_histos"]

xlims=[[-0.1,0.3],[-0.02,0.02],[-0.02,0.02], [-0.04,0.04]]
# ylims=[ [0.,0.25], [0.0,0.35], [0.0,0.15] ]
if (tweak_string == "DeltaR_"):
    ylims=[ [0.0,15], [0.,2.25], [0.0,1.], [0.0,2.25] ]
else:
    ylims=[ [0.0,15], [0.,2.25], [0.0,1.], [0.0,2.25] ]
    # ylims=[ [5.0,17.5], [0.,2.0], [0.0,1.] ]

fit_lims = [[-0.02,0.02], [-0.005,0.005], [-0.002,0.002], [-0.04,0.04]]

titles=[r"$(P_\mathrm{Truth} - P_\mathrm{Reco})/P_\mathrm{Truth}$",
        r"$\varphi_\mathrm{Reco}-\varphi_\mathrm{Truth}$",
        r"$\theta_\mathrm{Reco}-\theta_\mathrm{Truth}$",
        r"$\Delta\varphi_\mathrm{Reco}-\Delta\varphi_\mathrm{Truth}$"]

for i,s in enumerate(res_IDs):
    settings[s+"dir"] = TDirectories[i] 
    settings[s+"plotlim"]   =  xlims[i] 
    settings[s+"fitlim"]  = fit_lims[i] 
    settings[s+"title"]   =   titles[i] 
    settings[s+"ylim"]    =    ylims[i]


#Binning
bin_file =  ROOT.TFile.Open(next(iter(file_dict.items()))[1]) #grabs first dict val using iter object

eta_binning = {}
mom_binning = {}

eta_binning["edges"],N_eta= TVT_to_numpy(bin_file,"TVT_eta_bin")
mom_binning["edges"],N_mom = TVT_to_numpy(bin_file,"TVT_mom_bin")
eta_binning["centers"], eta_binning["widths"] = get_binning_from_edges(eta_binning["edges"])
mom_binning["centers"], mom_binning["widths"] = get_binning_from_edges(mom_binning["edges"])
N_eta = len(eta_binning["centers"])
N_mom = len(mom_binning["centers"])

eta_colors = get_colors(plt.cm.winter, N_eta)#would be nice to move to plotting py
mom_colors = get_colors(plt.cm.autumn, N_mom)

rad_to_mrad = 1000
to_percent = 100
max_uncertainty = 0.5 #Fits with large uncertainty omitted

h1_dictionary = {}

def get_th1f_dictionary():

    for B in B_Fields:
        file = ROOT.TFile.Open(file_dict["%1.1f"%(B)])

        for s in res_IDs:
            
            for eta in range(N_eta):

                for p in range(N_mom):

                    h1=file.Get(settings[s+"dir"]+"/%s_et_%i_p_%i_bin"%(s,eta,p))
                    identifier = ("%s_B_%1.1f_eta_%i_p_%i"%(s,B,eta,p))

                    h1_dictionary[identifier+"_vals"],h1_dictionary[identifier+"_errors"] = TH1_to_numpy_wErrors(h1,False,True)
                    h1_dictionary[identifier+"_N"] = h1.GetEntries()
                    h1_dictionary[identifier+"_avg"] = h1.GetMean()
                    h1_dictionary[identifier+"_sigma"] = h1.GetStdDev()
                    h1_dictionary[identifier+"_sigma_error"] = h1.GetStdDev()/np.sqrt(2*h1.GetEntries() - 2)
                    # error on stdev = sigma/sqrt{2*(n-1)}

                    if (p==0 and eta==0):
                        h1_dictionary[s+"edges"],h1_dictionary[s+"centers"],h1_dictionary[s+"widths"] = get_th1_binning_np(h1)
                    #print(h1.GetMean())

    for s,title in zip(res_IDs,titles):
        h1_dictionary[s+"title"] = title
    np.save("np_arrays/%s_indv_histos.npy"%(tweak_string),h1_dictionary)




res_dict = {}

def double_gauss_resolutions():
    #Obtains the sigma's from the narrow gaussian
    for s in res_IDs:

        for B in B_Fields:

            for eta in range(N_eta):

                for p in range(N_mom):

                    identifier = ("%s_B_%1.1f_eta_%i_p_%i"%(s,B,eta,p))

                    ydata = h1_dictionary[identifier+"_vals"]
                    xdata=h1_dictionary[s+"centers"]

                    xdata=xdata[~np.isnan(ydata)]#removes NaN
                    ydata=ydata[~np.isnan(ydata)]

                    min_i,max_i = fit_range(settings[s+"fitlim"],xdata)
                    xdata=xdata[min_i:max_i]
                    ydata=ydata[min_i:max_i]

                    #Single Gauss Implementation
                    if (s=="h1_dpp"):
                        popt, pcov = curve_fit(gauss, xdata, ydata, p0=[ydata.max(), xdata.std()/4])
                        perr = np.sqrt(np.diag(pcov))

                        res_dict[identifier] = h1_dictionary[identifier+"_sigma"] * to_percent 
                        res_dict[identifier+"Error"] = h1_dictionary[identifier+"_sigma_error"] * to_percent
                        res_dict[identifier+"Params"] = popt
                        # res_dict[identifier] = popt[2]*to_percent
                        # res_dict[identifier+"Error"] = perr[2]*to_percent

                    #Double Gauss Implementation
                    else:
                        popt, pcov = curve_fit(double_gauss, xdata, ydata, p0=[ydata.max(), xdata.std()/4,ydata.max()/5,xdata.std()])
                        perr = np.sqrt(np.diag(pcov))

                        res_dict[identifier] = popt[1]*rad_to_mrad 
                        res_dict[identifier+"Error"] = perr[1]*rad_to_mrad
                        res_dict[identifier+"Params"] = popt

        np.save("np_arrays/%s_resolutions.npy"%(tweak_string),res_dict)

def single_gauss_resolutions():
    #Obtains the sigma's from the narrow gaussian
    for s in res_IDs:

        for B in B_Fields:

            for eta in range(N_eta):

                for p in range(N_mom):

                    identifier = ("%s_B_%1.1f_eta_%i_p_%i"%(s,B,eta,p))

                    ydata = h1_dictionary[identifier+"_vals"]
                    xdata=h1_dictionary[s+"centers"]

                    xdata=xdata[~np.isnan(ydata)]#removes NaN
                    ydata=ydata[~np.isnan(ydata)]

                    min_i,max_i = fit_range(settings[s+"fitlim"],xdata)
                    xdata=xdata[min_i:max_i]
                    ydata=ydata[min_i:max_i]

                    #Single Gauss Implementation
                    if (s=="h1_dpp"):
                        popt, pcov = curve_fit(gauss, xdata, ydata, p0=[ydata.max(), xdata.std()/4])
                        perr = np.sqrt(np.diag(pcov))

                        res_dict[identifier] = h1_dictionary[identifier+"_sigma"]*to_percent 
                        res_dict[identifier+"Error"] = h1_dictionary[identifier+"_sigma_error"]*to_percent
                        res_dict[identifier+"Params"] = popt

                    else:
                        popt, pcov = curve_fit(gauss, xdata, ydata, p0=[ydata.max(), xdata.std()/4])
                        perr = np.sqrt(np.diag(pcov))

                        res_dict[identifier] = popt[1]*to_percent
                        res_dict[identifier+"Error"] = perr[1]*to_percent
                        res_dict[identifier+"Params"] = popt

        np.save("np_arrays/%s_resolutions.npy"%(tweak_string),res_dict)


def plot_indv_resolutions(fit_type="double"):

    #Load from File
    res_dict = np.load('np_arrays/%s_resolutions.npy'%(tweak_string),allow_pickle=True)[()]
    h1_dict = np.load('np_arrays/%s_indv_histos.npy'%(tweak_string),allow_pickle=True)[()]

    for s in res_IDs:

        for B in B_Fields:

            fig = plt.figure(figsize=(18,12))
            plt.xticks([])
            plt.yticks([])
            plt.title(h1_dictionary[s+"title"]+" ($B = %1.1f$)"%(B),fontsize=40,y=1.04,verticalalignment="bottom")

            for eta in range(N_eta):

                for p in range(N_mom):
                    ax = fig.add_subplot(N_eta,N_mom,p+eta*N_mom+1)

                    #Labels for eta and p ranges
                    if eta==0:
                        p_range_string = r"$%1.1f < p < %1.1f$"%(mom_binning["edges"][p],mom_binning["edges"][p+1])
                        plt.text(0.5,1.18,p_range_string,ha="center",va="top",size=20,transform=ax.transAxes)
                    if p==0:
                        eta_range_string = r"$%1.1f < \eta < %1.1f$"%(eta_binning["edges"][eta],eta_binning["edges"][eta+1])
                        plt.ylabel(eta_range_string,fontsize=16)

                    identifier = ("%s_B_%1.1f_eta_%i_p_%i"%(s,B,eta,p))
                    plt.errorbar(h1_dictionary[s+"centers"],h1_dictionary[identifier+"_vals"],yerr=h1_dictionary[identifier+"_errors"],fmt=".",color="midnightblue",alpha=0.6)

                    #Get X-Y data
                    xdata=h1_dictionary[s+"centers"]
                    ydata = h1_dictionary[identifier+"_vals"]

                    #remove NaN
                    xdata=xdata[~np.isnan(ydata)]
                    ydata=ydata[~np.isnan(ydata)]

                    #Apply Fit Range to X-Y data
                    min_i,max_i = fit_range(settings[s+"fitlim"],xdata)
                    xdata=xdata[min_i:max_i]
                    ydata=ydata[min_i:max_i]

                    #Single Gauss Implementation
                    if (s=="h1_dpp"):
                        popt = res_dict[identifier+"Params"]
                        plt.plot(xdata, gauss(xdata, *popt), 'r-', label='fit')
                        plt.text(0.05,0.9,r"$\sigma = %1.3f$"%(popt[1]),fontsize=8,transform=ax.transAxes)

                    #Double Gauss Implementation
                    else:
                        popt = res_dict[identifier+"Params"]
                        if (fit_type == "double"):
                            plt.plot(xdata, double_gauss(xdata, *popt), 'r--', label='fit')
                        elif (fit_type == "single"):
                            plt.plot(xdata, gauss(xdata, *popt), 'r--', label='fit')
                        plt.text(0.05,0.9,r"$\sigma = %1.3f$"%(popt[1]),fontsize=8,transform=ax.transAxes)

                    #aesthetics
                    plt.xlim(settings[s+"plotlim"])
                    plt.xticks(fontsize=8)
                    plt.yticks(fontsize=8)

                    plt.text(0.98,0.95,"B = %1.1f T"%(B),ha="right",va="top",size=15,alpha=0.7,transform=ax.transAxes)
                    plt.tight_layout()

        plt.savefig(tweak_string+s+"_"+"%1.1f"%(B)+"_doubleGauss_Fits.pdf")

def get_indv_resolutions(fit_type="double"):
    if (fit_type == "double"):
        double_gauss_resolutions() 
    elif (fit_type == "single"):
        single_gauss_resolutions()
    else: print("Fit type must be 'single' or 'double'")
        
res_array = np.zeros((len(res_IDs),len(B_Fields),N_eta,N_mom))
res_error_array = np.zeros((len(res_IDs),len(B_Fields),N_eta,N_mom))
def get_resolution_array():

    for i,s in enumerate(res_IDs):

        for b,B in enumerate(B_Fields):

            for eta in range(N_eta):

                for p in range(N_mom):

                    identifier = ("%s_B_%1.1f_eta_%i_p_%i"%(s,B,eta,p))
                    res_array[i][b][eta][p] = res_dict[identifier]
                    res_error_array[i][b][eta][p] = res_dict[identifier+"Error"]

def plot_resolutions_p_eta():

    for b,B in enumerate(B_Fields):

        fig = plt.figure(figsize=(24,12))
        # N_divs = (int(len(res_array[0,:,:])/N_mom)) #1/2 factor from Errors in dict, 2 factor from eta panels
        N_divs = 6
        ylabels = ["dp/p[%]",r"d$\varphi$ [mrad]",r"d$\theta$ [mrad]"]
        fmt = "o-"
        tick_leg_size = 20
        label_size = 25

        for i in range(N_divs):
        # for si,s in enumerate(res_IDs):
        #     i = si + len(B_Fields) 
            s = i%(len(res_IDs)-1) #for dpp,dth,dph switching
            #FIXME: put s inside p and eta loops

            if (s=="h1_eDelta_dph"): continue #not plotting DeltaPhi Res
            #aesthetics
            plt.tight_layout()
            ax = fig.add_subplot(2,int(N_divs/2),i+1)
                                    
            if (i < N_divs/2):

                for p in range(N_mom):

                    p_array,p_error,binning = truncate_uncertain(res_array[s,b,:,p],res_error_array[s,b,:,p],eta_binning["centers"])
                    # plt.errorbar(eta_binning["centers"],res_array[s,b,:,p],yerr=res_error_array[s,b,:,p],
                    plt.errorbar(binning,p_array,yerr=p_error,
                            fmt=fmt,color = mom_colors[p],label="$%1.0f<p<%1.0f~\mathrm{GeV}/c$"
                            %(mom_binning["edges"][p],mom_binning["edges"][p+1]))
                plt.xlabel(r'$\eta$',fontsize=25,x=0.5)
                # ylims = get_ylims(p_array, p_error)
                # plt.ylim(ylims[0],ylims[1])
                plt.ylim(settings[res_IDs[s]+"ylim"][0],settings[res_IDs[s]+"ylim"][1])
                plt.xlim(np.min(eta_binning["edges"]),np.max(eta_binning["edges"]))
                                                                                                    
            else:
                for eta in range(N_eta):
                    eta_array,eta_error,binning = truncate_uncertain(res_array[s,b,eta,:],res_error_array[s,b,eta,:],mom_binning["centers"])
                    plt.errorbar(binning,eta_array,eta_error,
                    # plt.errorbar(mom_binning["centers"],res_array[s,b,eta,:],yerr=res_error_array[s,b,eta,:],
                            fmt=fmt,color = eta_colors[eta],label="$%1.1f<\eta<%1.1f$"
                            %(eta_binning["edges"][eta],eta_binning["edges"][eta+1]))#confusing...
                plt.xlabel(r'$p~\mathrm{[GeV}/c]$',fontsize=25,x=0.5)
                # ylims = get_ylims(eta_array, eta_error)
                # plt.ylim(ylims[0],ylims[1])
                plt.ylim(settings[res_IDs[s]+"ylim"][0],settings[res_IDs[s]+"ylim"][1])
                plt.xlim(np.min(mom_binning["edges"]),np.max(mom_binning["edges"]))

            if (i==0): plt.text(0.05,.96,r'$B = %1.1f$ T'%(B),ha="left",va="top",size=25,alpha=0.7,transform=ax.transAxes)
            # if (i==0): plt.text(0.05,0.88,r'$%s$'%(cut_label),ha="left",va="top",size=25,alpha=0.7,transform=ax.transAxes)
            if (i==1): plt.legend(fontsize=tick_leg_size-3)
            if (i==5): plt.legend(fontsize=tick_leg_size-1)
            plt.tick_params(which='both',direction='in',right=True,top=True,bottom=True,length=10,labelsize=20)
            yticks = ticker.MaxNLocator(6)
            plt.ylabel(ylabels[s],fontsize=label_size,y=0.5)
            ax.yaxis.set_major_locator(yticks)
        plt.savefig(tweak_string+"B_%1.1f_resolutions_eta_mom.pdf"%(B))
        plt.show()


def truncate_uncertain(array, error, xbins):
    mask = np.logical_and(error/array < 0.33, array > 0.) 
    array = array[mask]
    error = error[mask]
    xbins = xbins[mask]
    return array,error,xbins

def get_ylims(array, error):

    a = 0.8
    y_min = np.nanmin(array-error) * a

    A = 1.2
    y_max =  np.nanmax(array+error) * A

    return y_min, y_max 

def weighted_avg_resolutions():
    for s in res_IDs:
        for B in B_Fields:
            weighted_avg = 0
            sum_N = 0
            for eta in range(N_eta):            
                for p in range(N_mom):
                    identifier = ("%s_B_%1.1f_eta_%i_p_%i"%(s,B,eta,p))
                    #print(B,res_dict[identifier+"_dDeltaPhi_res"], h1dict[identifier+"_N"])
                    weighted_avg += res_dict[identifier]*h1_dictionary[identifier+"_N"]
                    sum_N +=h1_dictionary[identifier+"_N"]
            print(s,B,weighted_avg/sum_N)
