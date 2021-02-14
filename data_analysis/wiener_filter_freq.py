# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 13:16:10 2020
Wiener Filter in Frequency Domain using Variance for Noise Estimate
@author: Adi Pall
"""
import numpy as np
from matplotlib import pyplot as plt
from numpy import fft
from scipy.io import wavfile

def def_sig(N_samp,sig_type,params):
    """ Define impulse or sine signal given input parameters 
        Inputs-
            N_samp: number of samples (int)
            sig_type: 'sin' or 'imp'
            params: dictionary of coefficients
                for impulse: 'A'*exp(-'scale'*(t-'shift')**2)
                for sine: 'A'*sin(2*pi*'f0'*t)
                'dt' is the timestep
        
        Returns-
            t: numpy array of time vals
            signal: numpy array of corresponding signal vals
    """
    t = params['dt'] * np.arange(N_samp)
    if(sig_type=='imp'):
        signal = params['A']*np.exp(-params['scale'] * (t - params['shift']) ** 2)
    elif(sig_type=='sin'):
        signal= params['A']*np.sin(2*np.pi*params['f0']*t)
    return t, signal

def add_white(signal,scale=0.1):
    """ Add Gaussian White Noise to signal
        Inputs-
            signal: numpy array of signal vals
            scale: magnitude of noise (default = 0.1)
        
        Returns-
            signal + noise: numpy array of noisy signal
    """
    np.random.seed(10)
    noise = np.random.normal(0, scale, size=signal.shape) # white noise   
    return signal + noise

def add_tv_white(signal,tv_scale=0.1):
    """ Add Time-Variant Gaussian White Noise to signal
        Inputs-
            signal: numpy array of signal vals
            tv_scale: magnitude of noise (default = 0.1)
        
        Returns-
            signal + noise: numpy array of noisy signal
    """
    np.random.seed(10)
    slope = np.arange(signal.shape[0])/signal.shape[0]  # goes up to 1 (linear increase)
    noise = slope * np.random.normal(0, tv_scale, size=signal.shape)
    return signal + noise

def get_stats(noisy_signal,w_size):
    """ Get mean and variance, for all windows (w/ overlap) of noisy signal in time domain.
        Inputs-
            noisy_signal: numpy array of noisy signal vals
            w_size: number of samples to consider per FFT
        
        Returns-
            mu: numpy array of each window average
            sigma: numpy array of each window variance
    """
    # pad signal so that while loop doesn't miss any data (will introduce some edge error)
    signal=np.concatenate((np.zeros(w_size//2),noisy_signal,np.zeros(w_size//2)),axis=0)
    L=signal.shape[0]
    start=w_size//2
    mu=np.zeros(noisy_signal.shape)
    sigma=np.zeros(noisy_signal.shape)
    while((start+w_size/2)<L):
        mu[start-w_size//2]=np.mean(signal[start-(w_size//2):start+(w_size//2)])
        sigma[start-w_size//2]=np.std(signal[start-(w_size//2):start+(w_size//2)])**2
        start=start+1
    return mu,sigma

def wien_f(noisy_signal,w_size,slide,w_type='rect'):
    """ Wiener Filter in frequency domain using windows and overlapping.
        Inputs-
            noisy_signal: numpy array of noisy signal vals
            w_size: number of samples to consider per FFT
            slide: number of samples to slide window over per FFT
            w_type: type of window to use ['hanning or 'rect'] (default = 'rect')
        
        Returns-
            filt_sig: numpy array of filtered signal
    """
    L=noisy_signal.shape[0]
    filt_sig=np.zeros(L)
    start=0 # index to hold start position
    mu,sigma=get_stats(noisy_signal,w_size)  # assume noise in PSD constant = sigma**2
    # this is the reason for decreased performance in time-varying noise
    noise_est=np.mean(sigma) # sigma contains all windows sigma, so average them
    overlap_fact=slide/w_size # ratio of movement to window size (indicative of overlap)
    while((start+w_size)<L):
        if w_type == "rect":
            cur=noisy_signal[start:start+w_size]
        elif w_type == "hanning":
            # hanning has amplitude 1/2
            cur=noisy_signal[start:start+w_size]*2*np.hanning(w_size)
        X = fft.fft(cur)
        X_psd= abs(X)**2/w_size
        X_est=X_psd-noise_est
        X_est[X_est<0]=0  # magnitude must be non negative
        phi=(X_est)/(X_psd) # 'Optimal Wiener Filter'
        filt_X = phi*X
        # real avoids complex to real warning and ratio compensates for overlap.
        filt_sig[start:start+w_size]=filt_sig[start:start+w_size]+np.real(fft.ifft(filt_X))*overlap_fact 
        start=start+slide
    # if not a whole multiple of w_size, process rest also
    cur=noisy_signal[start:]
    X = fft.fft(cur)
    X_psd = abs(X)**2/w_size
    X_est = X_psd-noise_est
    X_est[X_est<0] = 0  # magnitude must be non-negative
    phi=(X_est)/(X_psd)
    filt_X = phi*X
    filt_sig[start:]=filt_sig[start:]+np.real(fft.ifft(filt_X))*overlap_fact
    return filt_sig

if __name__ == "__main__":
    
    fs=40. # sampling rate
    params=dict() # use dictionary for better user experience
    params['dt']=1/fs
    
    # sine params
    params['A']=1 # amplitude
    params['f0']=0.08 # natural frequency
    N=4000 # number of samples
    
    # impulse params
    params['shift']=30
    params['scale']=0.2
    
    # define noise
    noise_type = 'white' # or tv_white
    scale=0.1
    tv_scale=0.1
    
    # Wiener filter window parameters
    w_size=500
    slide=60
    
    # Testing with generated signals
    for sig_type in ['imp','sin']:
        for noise_type in ['white','tv-white']:
            # define signal
            t,signal=def_sig(N,sig_type,params);
            if(noise_type=='white'):
                noisy_signal = add_white(signal,scale)
            else:
                noisy_signal = add_tv_white(signal,tv_scale)
            
            # filter signal
            for w_type in ['rect','hanning']:
                filtered_sig=wien_f(noisy_signal,w_size,slide,w_type)
                
                plt.figure(figsize = (6.5,4.5))
                plt.plot(t,noisy_signal)
                plt.plot(t,filtered_sig)
                plt.title(sig_type+'_'+noise_type+'_'+w_type)
                plt.legend(['unfiltered','filtered'])
                plt.xlabel('Time (s)')
                plt.ylabel('Amplitude')
                plt.savefig('./demo/'+sig_type+'_'+noise_type+'_'+w_type+'.png')
                plt.show()
                
    # Testing with audio files
    w_size = 1000 # need larger window for the audio because of higher sampling rate
    slide  = 120
    noise_type = 'white'
    scale=0.01
    tv_scale=0.01
    # .wav example soundbites taken from:
        # https://www2.cs.uic.edu/~i101/SoundFiles/
    
    for fname in ['CantinaBand3','StarWars3']:
        fs,data=wavfile.read('./demo/'+fname+'.wav')
        dt = 1/fs
        data=data/32000.0 # data stored in 16 bit signed, scale to between -1 and 1
        tot_t = data.shape[0]*dt
        t = np.arange(0,tot_t,dt)
        for noise_type in ['white','tv-white']:
            if(noise_type=='white'):
                noisy_audio = add_white(data,scale)
            else:
                noisy_audio = add_tv_white(data,tv_scale)
            wavfile.write('./demo/'+fname+'_'+noise_type+'_noisy.wav',fs,noisy_audio)
            for w_type in ['rect','hanning']:
                filt_audio = wien_f(noisy_audio,w_size,slide,w_type)
                wavfile.write('./demo/'+fname+'_'+noise_type+'_filt.wav',fs,filt_audio)
                plt.figure(figsize = (6.5,4.5))
                plt.plot(t,noisy_audio)
                plt.plot(t,filt_audio)
                plt.title(fname+'_'+noise_type+'_'+w_type)
                plt.legend(['unfiltered','filtered'])
                plt.xlabel('Time (s)')
                plt.ylabel('Amplitude')
                plt.savefig('./demo/'+fname+'_'+noise_type+'_'+w_type+'.png')
                plt.show()
