# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 17:04:51 2021
@author: Adi Pall
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sin
from numpy.fft import fft, ifft, fftfreq, fftshift
    
def dualplot(freq, F, name):
    """ simple plot of real and imaginary parts """
    plt.figure(figsize=(6.5, 2.5))
    plt.suptitle(name)
    plt.subplot(1, 2, 1)
    plt.plot(freq, np.real(F), '.k')
    plt.ylabel('Re(F)')
    plt.subplot(1, 2, 2)
    plt.plot(freq, np.imag(F), '.k')
    plt.ylabel('Im(F)')
    plt.subplots_adjust(wspace=0.5)
    
def low_pass(rate, data, length, cutoff):
    """ basic low pass filter for data
        Arguments:
            rate - the sample rate (Hz) (wav default: 44100 Hz)
            data - the samples
            length - the length (s) of the sample window
            cutoff - remove frequencies above this value
        Returns:
            filtered - the filtered data
    """
    n = data.shape[0]
    freq = fftfreq(n, 1/rate)  # freqs in Hz, (0 ... n/2-1, -n/2, .. -1)/L
    df = fft(data)/n

    # find index for cutoff freq.
    # then cut off frequences k+1 to n/2 and -n/2 to -(k+1)
    k = np.searchsorted(freq[0:n//2], cutoff)
    df[k+1:n//2] = 0
    df[n//2:-k] = 0

    filtered = n*ifft(df)
    return filtered

def plot_spectrum(data, rate, plotname):
    """ plots the log(|F|) for the DFT F vs. frequency (0 to n/2)
        on a semilog x plot. Note that log(magnitude) is proportional
        to decibels, the familiar measure for sound volume here.
    """
    n = data.shape[0]

    transform = fft(data)/n
    log_magnitude = np.log(np.abs(transform))
    freq = fftfreq(n, 1/rate)
    skip = round(freq.shape[0]/2048)  # thin out the plot
    plt.semilogx(freq[:n//2:skip], log_magnitude[:n//2:skip])
    plt.xlabel('freq (Hz)')
    plt.ylabel('log(|F|)')
    plt.title(plotname)
    
##############################Don't worry about above yet##############################################
    
def read_piezo(fname, start = 25, samples = 8192):
    """ using readlines to get all lines first"""
    fname = 'piezos\\' + fname + '.csv'
    fp = open(fname, 'r')
    lines = fp.readlines()
    n = samples
    x = [0]*n
    y = [0]*n
    for k in range(start, start+8192):
        words = lines[k].split(',')
        x[k-start] = float(words[0])
        y[k-start] = float(words[1])
    fp.close()
    return x, y
    
if __name__ == "__main__":
    n = 8192 # num samples
    rate = 500 # 500 Hz sample rate
    sets = 4 # testing 4 data collections (piezo pickups)
    i = 0 # counter
    t = np.zeros([n,sets]) # pre-allocate matrix
    V = np.zeros([n,sets]) # Voltage amplitudes
    total_t = np.zeros(sets) # will need total t-span for fourier domain analysis
    for k in range(16,20): # just the way files named
        t[:,i],V[:,i]=read_piezo(f'acq00{k}')
        shift = t[0,i] # changes relative time to a more comprehensible start from 0.0s
        t[:,i] = t[:,i] + np.abs(shift)
        cur_ts = t[:,i]
        cur_Vs= V[:,i]
        plt.figure(num=i)
        plt.plot(cur_ts,cur_Vs)
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude (V)')
        tit = "Trial {}".format(k)
        plt.title(tit)
        total_t[i] = np.max(cur_ts)
        plt.figure(num=i+sets)
        tit = tit + ' Freq Spectrum'
        plot_spectrum(cur_Vs,rate,tit)
        
        i += 1

# still to do: more detailed fourier domain analysis and convert Voltage to Strain or Accel using constant from
# datasheet
    
    
        