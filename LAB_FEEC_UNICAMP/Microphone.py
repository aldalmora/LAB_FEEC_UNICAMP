# -*- coding: utf-8 -*-
from .Signal import Signal

class Microphone():
    def __init__(self,name,signal,dist_source):
        self.name = name
        self.signal = signal
        self.dist_source = dist_source