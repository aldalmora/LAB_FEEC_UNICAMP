import cmath
import math
import os
import sys
from collections import deque

import IPython
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
import scipy.io
import scipy.signal
from scipy.io import wavfile
from scipy.ndimage.filters import gaussian_filter1d

from .Beamforming import (applyDAS, applyMVDR, calcManifoldVector,
                          calcWeightDelay)
from .Microphone import Microphone
from .Signal import Signal
from .Utils import find_delay, toFreq, toTime

#Encapsulates the concept of an linear array between a source and a sample, based on its distances.
class Array:
        
    def __init__(self, c):
        self.Microphones = [] #Microphones positions from the closer to the further from the source
        self.c = c #Speed of sound
        self.dis_mic_amos = 0 #Distance between the closest microphone to the sample and the sample
        self.dis_mic_autf = 0 #Distance between the closest microphone to the source and the sample

    #Read ITA Files
    #Calculate the time difference between signals by GCC-PHAT
    @classmethod
    def fromIta(self, c, fs, arq, dis_mic_amos, dis_mic_autf):
        ret = Array(c)
        ret.fs = fs
        ret.dis_mic_amos = dis_mic_amos
        ret.dis_mic_autf = dis_mic_autf

        file_matlab = sp.io.loadmat(arq)
        TF = file_matlab['ITA_TOOLBOX_AUDIO_OBJECT'][0][0][14].T
        #T[3] - Closer mic to the source
        #T[0] - Further mic
        
        #Took the reference measurement to remove the source eigenfunction
        file_ref = sp.io.loadmat('mic_ref.ita')
        TF_ref = file_ref['ITA_TOOLBOX_AUDIO_OBJECT'][0][0][14].T
        
        totalS = len(TF_ref[0])
        sigTref = Signal(TF_ref[0])
        w = calcWeightDelay(ret.fs,totalS,-2/ret.c)
        sigTref.weight(w)
        ret.TF_ref = sp.array(sigTref.weighted())
        signals = [Signal(T/ret.TF_ref) for T in TF]

        #Calculate the distance between each pair of microphones
        dis_1_2 = find_delay(toTime(signals[0].Spectrum),toTime(signals[1].Spectrum),fs)[0]*c
        dis_2_3 = find_delay(toTime(signals[1].Spectrum),toTime(signals[2].Spectrum),fs)[0]*c
        dis_3_4 = find_delay(toTime(signals[2].Spectrum),toTime(signals[3].Spectrum),fs)[0]*c

        ret.Microphones = []
        ret.Microphones.append(Microphone('Mic3',signals[3],dis_mic_autf))
        ret.Microphones.append(Microphone('Mic2',signals[2],dis_mic_autf+dis_3_4))
        ret.Microphones.append(Microphone('Mic1',signals[1],dis_mic_autf+dis_3_4+dis_2_3))
        ret.Microphones.append(Microphone('Mic0',signals[0],dis_mic_autf+dis_3_4+dis_2_3+dis_1_2))

        return ret

    #Read WAV file and simulate it reaching the array and reflecting on the sample
    def simulatefromWav(self, arq, pos_mics, dis_mic_amos, coefInc=1 ,coefRef=0.2):
        self.fs, signal = wavfile.read(arq)
        self.dis_mic_amos = dis_mic_amos
        self.wav_time_signal = signal

        TF = toFreq(signal)
        lenSpec = len(TF)
        pos_mics = sorted(pos_mics)
        pos_amos = pos_mics[-1] + dis_mic_amos
        self.Microphones = []

        for i in range(0,len(pos_mics)):
            sinalInc = Signal(TF)
            w = calcWeightDelay(self.fs,lenSpec,pos_mics[i]/self.c)
            sinalInc.weight(w)

            sinalRef = Signal(sinalInc.weighted())
            w = calcWeightDelay(self.fs,lenSpec,2*(pos_amos-pos_mics[i])/self.c)
            sinalRef.weight(w)

            sinalTotal = Signal(coefInc*sinalInc.weighted() + coefRef*sinalRef.weighted())

            mic = Microphone('',sinalTotal,pos_mics[i])
            self.Microphones.append(mic)
