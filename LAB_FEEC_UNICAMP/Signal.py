# -*- coding: utf-8 -*-
from .Utils import toTime,toFreq
import numpy as np

#Encapsulates the signal with some functions
#Store its original signal and the weights applied
class Signal:
    def __init__(self,Spectrum):
        self.Spectrum = Spectrum
        self.Weights = np.ones(len(Spectrum))
        
    def freq(self):
        return self.Spectrum
    
    def time(self):
        return toTime(self.Spectrum)
    
    def add_noise(self,noise):
        self.Spectrum = self.Spectrum + np.squeeze(noise)
    
    def weight(self,weights):
        self.Weights = weights
    
    def weighted(self):
        return self.Spectrum * self.Weights
    
    def weighted_time(self):
        return toTime(self.weighted())