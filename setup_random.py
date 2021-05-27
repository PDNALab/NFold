#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from meld import system
from meld import parse

def setup_system():
    # load the sequence
    sequence = parse.get_sequence_from_AA1(filename='sequence.dat')
    n_res = len(sequence.split())
    
    # build the system
    p = system.subsystem.SubSystemFromSequence(sequence)           
    b = system.builder.SystemBuilder(forcefield="ff14sbside")
    s = b.build_system([p])                  

setup_system()
