import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np

def get_colors(color_map, N,test=False):
    evenly_spaced_interval = np.linspace(0, 1, N)
    colors = [color_map(x) for x in evenly_spaced_interval]
    if test:
        fig = plt.figure()
        x = np.linspace(0, 10, 100)
        for i in range(N):
            if i%2 == 0:
                plt.plot(np.sin(x+4*i),color=colors[i])
            else:
                plt.plot(np.cos(x+4*i),color=colors[i])
    return colors

def get_ranges_from_dict(dict,keys,key_base,n_subkeys=0,Errors=True):
    max_scale = 1.1
    min_scale = 0.001
    ranges = {}
    for s in keys:
        min = 999.
        max = 0.
        if (n_subkeys !=0):
            for i in range(n_subkeys):
                if Errors:
                    temp_min = np.nanmin(dict[key_base+"_%s_%i"%(s,i)] - dict[key_base+"_%s_%i_Errors"%(s,i)])
                    temp_max = np.nanmax(dict[key_base+"_%s_%i"%(s,i)] + dict[key_base+"_%s_%i_Errors"%(s,i)])
                else:
                    temp_min = np.nanmin(dict[key_base+"_%s_%i"%(s,i)])
                    temp_max = np.nanmax(dict[key_base+"_%s_%i"%(s,i)])

                if (temp_min < min):
                    min = temp_min
                if (temp_max > max):
                    max = temp_max
                
        else:
            if Errors:
                min = np.nanmin(dict[key_base+"_%s"%(s)] - dict[key_base+"_%s_Errors"%(s)])
                max = np.nanmax(dict[key_base+"_%s"%(s)] + dict[key_base+"_%s_Errors"%(s)])
            else:
                min = np.nanmin(dict[key_base+"_%s"%(s)])
                max = np.nanmax(dict[key_base+"_%s"%(s)])

        ranges["%s_min"%(s)] = min*min_scale
        ranges["%s_max"%(s)] = max*max_scale
 
    return ranges

def get_element_from_val(array,val):
    array = np.asarray(array) #for list inputs
    diff_array = np.absolute(array-val)
    return np.argmin(diff_array) #returns element (np.amin returns value)