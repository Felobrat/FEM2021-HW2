import numpy as np
from sympy import *
import matplotlib.pylab as plt

def sigma(x):
    fid = open(f"{x}","r")

    while True:
        line = fid.readline()
        if line.find("Sigma_x") >= 0:
            break
        
    while True:
        line = fid.readline() 
        s1 = line.split()
        if np.float(s1[0]) > 2:
            break
    while True:
        line = fid.readline() 
        s1 = line.split()
        if np.float(s1[0]) > 0:
            break
    Nelems = np.int32(fid.readline())
    sigmamax=0
    sigmamin=0
    nodo=0
    for i in range(Nelems):
        line = fid.readline()
        sl = line.split()
        if sigmamin > float(sl[1]):
            sigmamin= float(sl[1])
            #print(sigmamin)
        if sigmamax < float(sl[1]):
            sigmamax= float(sl[1])
            #print(sigmamax)
    return sigmamax



def desp(x):
    fid = open(f"{x}","r")

    while True:
        line = fid.readline()
        if line.find("Despl") >= 0:
            break
        
    while True:
        line = fid.readline() 
        s1 = line.split()
        if np.float(s1[0]) > 2:
            break
    while True:
        line = fid.readline() 
        s1 = line.split()
        if np.float(s1[0]) > 0:
            break
    Nnodes = np.int32(fid.readline())
    defo=0
    nodo=0
    for i in range(Nnodes):
        line = fid.readline()
        sl = line.split()
        if defo > float(sl[1]):
            defo= float(sl[1])
    return defo    