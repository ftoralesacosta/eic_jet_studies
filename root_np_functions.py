import ROOT
import numpy as np
from ROOT import TGraphErrors
from ROOT import TVector
import math

def TVT_to_numpy(file,TVectorT_name):
    TVectorT = file.Get(TVectorT_name)
    N = TVectorT.GetNoElements()
    np_bins = np.zeros(N)
    for i in range(N):
        np_bins[i] = TVectorT[i]
    return np_bins,N-1 #return number of bin centers

def TH1_to_numpy_wErrors(th1, normalize=False, zero_suppress=False, scale_factor=1.):    

    Nx = th1.GetNbinsX()
    th1_array = np.zeros(Nx)
    th1_errors = np.zeros(Nx)
    
    for i in range(Nx):
        val = th1.GetBinContent(i+1)
        error = th1.GetBinError(i+1) * math.sqrt(scale_factor)
        if zero_suppress and val == 0:
            val = np.nan
        if not(zero_suppress) and val!=0 and error/val >= 0.5:
            val=np.nan #skip converged fits with poor statistics
        th1_array[i] = val
        th1_errors[i] = error
    if normalize:
        N = np.nansum(th1_array)*scale_factor
        th1_array = th1_array/N # Equivalent to (1/sqrt(scale_factor))
        th1_errors = th1_errors/N
    return th1_array, th1_errors

def get_th1_binning_np(th1):
    N = th1.GetNbinsX()
    min = th1.GetXaxis().GetXmin()
    max = th1.GetXaxis().GetXmax()
    bins = np.linspace(min,max,N+1)
    centers = (bins[1:] + bins[:-1]) /2
    widths = [(j-i)/2 for i, j in zip(bins[:-1], bins[1:])]
    if (len(centers)!=N):
        print("get_th1_binning_np: something went wrong")
    return bins,centers,widths


def get_binning_from_edges(bin_edges):
    centers = (bin_edges[1:] + bin_edges[:-1]) /2
    widths = [(j-i)/2 for i, j in zip(bin_edges[:-1], bin_edges[1:])]
    return centers,widths
