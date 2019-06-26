import cmath
import math

import numpy as np

#Find the delay between too signals
def find_delay(sign, reference, fs=44100):
    #Safe Range
    sr = 16

    ##Cross-Correlation
    R = toFreq(sign) * np.conj(toFreq(reference))
    R = R / np.abs(R)

    #Time-Domain
    cc = toTime(R, n=(sr * n))

    max_samples = int(sr * n / 2)
    cc = np.concatenate((cc[-max_samples:], cc[:max_samples+1]))

    # Find Delay
    samples = np.argmax(np.abs(cc)) - max_samples
    time_diff = samples / (sr * fs)
    
    return time_diff

def toFreq(t,n=None):
    return np.fft.rfft(t,n=n)
    
def toTime(T,n=None):
    return np.fft.irfft(T,n=n)

#Transform a signal in its array vector representation given the array manifold vector
def DeslocSig(sig, v_manifold):
    return sig*np.conj(v_manifold)
