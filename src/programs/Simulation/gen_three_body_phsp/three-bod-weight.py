#!/usr/bin/python3
#
from math import *
import numpy as np
import vector
#import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import sklearn.linear_model
from sklearn.metrics import mean_squared_error, r2_score
#
#
#
#   Define particle masses
#
Masspi0 = 0.1349770
Masspiq = 0.13957039
Masskq  = 0.493677
Masseta = 0.547862
#
#
MassP   = 0.938272081
#
MinEgam = 8.0
MaxEgam = 9.0
#
# Set up random number generator
#
rng = np.random.default_rng(12345)
#print(rng)
#print(rng.random())
#
NumX = 100
MXStep = 1.0/NumX
#
NumPts=10000
MassStep = 1.0/NumPts
#
# ###################################
# Set up the three daughter particles
#
m1 = Masspiq
m2 = Masspiq
m3 = Masseta
#
# ###################################
#
# Generate a random mass between MinMass and MaxMass
#
def Uniform(Min,Max):
    a = rng.random()
    return a*(Max-Min)+Min
#
#----------------
#
# Define the Kallen Function
#
def Lambda(a,b,c):
#
    lam = a*a - 2.*a*(b+c) + (b-c)*(b-c) 
    if lam < 0.:
        lam = 0.
#
    return lam
#
#
# ######################################
#  Start of code
# ######################################
#
#
# Compute the mass squared.
#
m1sq = m1*m1
m2sq = m2*m2
m3sq = m3*m3
#
# Set up arrays to contain the final answer. 
#  MassX[i] is the mass of the eta pi pi system
#  MassXa[i] is the 2 pi mass that maximizes that weight
#  WeightX[i] is the maximum weight associated the this eta pi pi mass
#
MassX   = np.zeros(NumX)
MassXa  = np.zeros(NumX)
WeightX = np.zeros(NumX)
dhat = np.zeros(NumX)
#
#  Set up the mass range for MassX
#
MinMassX = m1 + m2 + m3
MaxMassX = 3.000
#
#
# Loop over values of MassX
#
i=0
while i<NumX:
#
# Calculate the mass of X and its square. 
#   Make sure it is not exacly the sum of the three final state mesons.
#
    M = MinMassX + MXStep * i * ( MaxMassX - MinMassX )
    if i <= 0 :
        M=M+0.0001
#
    Msq  = M*M
#
# Compute the normalization factor
#
    m12sq = (m1+m2)*(m1+m2)
    lam1 = Lambda(Msq,m12sq,m3sq)
#
    Mm3sq = (M-m3)*(M-m3)
    lam2 = Lambda(Mm3sq,m1sq,m2sq)
#
    if lam1<=0.:
        print("Error: lam1 =",lam1,"M=",M)
    if lam2<=0.:
        print("Error: lam2 =",lam2,"M=",M)
#
    wn = (M-m3)/np.sqrt(lam1*lam2)
#
# Loop over the possible masses (ma) for the m1/m2 system
#
    MaxVal = 0.
    MaxMa  = 0.
#
    im=0
    while im<NumPts:
#
#   Get the mass of the resonance.
#
        ma = (m1+m2) + MassStep * im * (M-m1-m2-m3)
        masq = ma*ma
#
        lam1 = Lambda(masq,m1sq,m2sq)
        lam2 = Lambda(Msq,masq,m3sq)
#
        w = wn * np.sqrt(lam1 * lam2) / ma
#
        if w>MaxVal:
            MaxVal = w
            MaxMa  = ma
#
#    print("Mass= {:4f}, weight = {:4f}".format(ma,w))
#
        im+=1
#
#    print("Maximum: MassX= {:4f} Massa= {:4f}, weight= {:4f}".format(M,MaxMa,MaxVal))
#
    MassX[i]   = M
    MassXa[i]  = MaxMa
    WeightX[i] = MaxVal
#
    i+=1 
#
# Data have now been obtained.
# Fit a linear regression
#
# Create a LinearRegression object to fit (Mass a) = A*(Mass X) + B
# 
reg = sklearn.linear_model.LinearRegression(fit_intercept=True)
#
# Fit using the  MassX and MassXa date.
reg.fit(MassX.reshape(-1,1),MassXa) 
#
# Use the linear regression to predict the MassXa values from the MassX
y_hat = reg.predict(MassX.reshape(-1,1))
#
# Compute the difference between the true values of MassXa and the predicted ones.
#
i=0
while i<NumX:
    dhat[i] = MassXa[i]-y_hat[i]
    i+=1
#
# Print the best-fit linear regression expression
#
print('Fit: Slope:     {%.8f}'%reg.coef_)
print('Fit: Intercept: {%.8f}'%reg.intercept_)
# The mean squared error
#
print('Mean squared error: %.4f'%mean_squared_error(MassXa,y_hat))
#
# The coefficient of determination: 1 is perfect prediction
print("Coefficient of determination: %.4f" %r2_score(MassXa, y_hat))
#
# Plot the original data
plt.plot(MassX, MassXa, '.')
plt.xlabel('$M(\eta\pi\pi)$')
plt.ylabel('$M(\pi\pi)_{max}$')
plt.title('$M(\pi\pi)_{max}$ versus $M(\eta\pi\pi)$')
#
#
# Plot the regression-predicted datapoints as a line
plt.plot(MassX, y_hat, '-')
#
plt.savefig("MvsM_fit.png")
