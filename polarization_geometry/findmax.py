import numpy as np
import scipy.interpolate as interp
import scipy.signal as sgnl

def findmax(w, sigma, high):
    
    # Interpolation
    w_int = np.logspace(np.log10(np.amin(w)),
                        np.log10(np.amax(w)),1000)
    sig_int = interp.interp1d(w, sigma, kind='cubic', fill_value='extrapolate')
    
    sig_int_arr = sig_int(w_int)
    
    peaks, _ = sgnl.find_peaks(sig_int_arr)
    if len(peaks)==0:
        sig_max = np.nan
        tau = np.nan
    else:
        if high == 0:
            index = peaks[0]
        elif high == 1:
            index = peaks[len(peaks)-1]
        w_max = w_int[index]
        tau = 1/w_max
    
        sig_max = sig_int_arr[index]
    
    return [sig_max, tau]