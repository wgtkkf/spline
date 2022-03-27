# 3rd order spline interporation & integral calculation
# change input files name & input files number
# read gap & flux data files
# Coded by Takuro TOKUNAGA
# Last modified: May 15 2018

import math
import numpy as np
import cmath
import time
import pandas as pd
import scipy.linalg as linalg
import matplotlib.pyplot as plt
from scipy.integrate import trapz, simps, quad, quadrature, romberg
start = time.time()

# function: begin
def begin():
    print ("begin")

# function: end
def end():
    print ("end")

# parameters
#number = 370 # number of data points x & fx, change here, for Gold
#number = 361 # number of data points x & fx, change here, for Si
number = 61 # number of data points x & fx, change here, for SiC, number=100
integral = 0 # initialization
sum = 0      # initialization
filenum = 46 # number of flux txt files, change here

#
sx=np.zeros((number), dtype='float64')
fx=np.zeros((number), dtype='float64')

#  matrix and vectors
ma=np.zeros((number,number), dtype='float64')
vx=np.zeros((number), dtype='float64')
vb=np.zeros((number), dtype='float64')

# vector for gap
vg=np.zeros((filenum), dtype='float64') # vector for gap

# main
begin()

# output file
f = open('power.txt', 'w')

# file read x & fx
for j in range(0, filenum):
    #data = pd.read_csv("flux"+str(j)+".txt", sep="\t", header=None) # change file name here tab:\t for local PC
    data = pd.read_csv("../fe/results/sic/flux"+str(j)+".txt", sep=" ", header=None) # change file name here tab:\t for CHPC
    data.columns = ["omega", "prop", "evan", "prop+evan"]

    gdata = pd.read_csv("gap.txt", sep=" ", header=None) # gap data
    gdata.columns = ["gap"]
    vg[j] = gdata.iat[j,0] # 0th line

    # input data into tables
    for i in range(0, number):
        sx[i] = data.iat[i,0] # x line
        fx[i] = data.iat[i,3] # fx line, 1:prop, 2:evan, 3:total

    # initialization of matrix & vectors
    # matrix A
    for i in range(0, number): # from 0 to 13, total 14
        # sub-diagonal components
        if i>0 and i<number-1: # 1 to 12, number=14
            ma[i][i-1] = sx[i] - sx[i-1] # left
            ma[i][i+1] = sx[i+1] - sx[i] # right

        #diagonal components
        if i==0:
            ma[i][i] = 1
        elif i>0 and i<number-1: # 1 to 12
            ma[i][i] = 2*(ma[i][i-1]+ma[i][i+1])
        elif i==number-1: # i=13
            ma[i][i] = 1

    # vector b
    for i in range(0, number): # from 0 to 13, total 14
        if i>0 and i<number-1: # from 1 to 12
            h1 = sx[i+1]-sx[i]
            h0 = sx[i]-sx[i-1]
            a2 = fx[i+1]
            a1 = fx[i]
            a0 = fx[i-1]

            vb[i] = (3/h1)*(a2-a1)-(3/h0)*(a1-a0)

        # reset parameters
        h1 = 0
        h0 = 0
        a2 = 0
        a1 = 0
        a0 = 0

        #print(str(vb[i]))

        # calculation of vector x, LU method
        #vx = np.linalg.solve(ma, vb) originally written like this

    # calculation of vector x, LU method
    vx = np.linalg.solve(ma, vb) # maybe this is correct

    # integral calculation
    for i in range(0, number-1):

        # reset coefficients
        ipfx = 0
        sa = 0
        sh = 0
        sb = 0
        sd = 0

        # interporation function
        def ipfunction(x):
            tempx1 = sx[i]
            tempx2 = sx[i+1]

            # coefficients
            sa = fx[i] # small a
            sh = sx[i+1]-sx[i] # small h
            sb = (1/sh)*(fx[i+1]-fx[i]) - (sh/3)*(2*vx[i]+vx[i+1]) # small b
            sc = vx[i] # small c
            sd = (vx[i+1]-vx[i])/(3*sh) # small d

            # define of interporation function
            ipfx = sa + sb*(x-sx[i]) + sc*np.power((x-sx[i]),2)\
            + sd*np.power((x-sx[i]),3)

            return ipfx

        integral = quad(ipfunction, sx[i], sx[i+1])[0]
        sum = sum + integral

    # output for file
    f.write(str(vg[j]))
    f.write(str(' '))
    f.write(str(sum)) # [W/m2]
    f.write('\n')

    # reset parameters for next files
    sum = 0

# file close
f.close()

end()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
