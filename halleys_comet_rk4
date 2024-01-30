import numpy as np 
import matplotlib.pyplot as plt

#starting conditions
M = 1.989e30
G = 6.67e-11/((1.498e11)**3)
ri = 5.2e12/(1.498e11)
rz = 5.2e12/(1.498e11)
thetai = 0
vi = 0

#starting angular velocity
ai = (880/((1.498e11)))/ri
ti = 0
tmax = 40000000000
k = 200000

#arrays for time and the co-ordinates of the planet
theta = []
r = []
t = []
y = []
x = [] 

#radial acceleration
def dv(r, a):
    dv = (r*a**2) - G*M/r**2
    return dv
    
#tangential acceleration
def da(r, a, v):
    da = -(2*a*v)/r
    return da

#individual rkstep
def rkstep(ri, vi, ai, h):
    
    kd1 = h*dv(ri, ai)
    ka1 = h*da(ri, ai, vi)
    kr1 = h*vi
    kt1 = h*ai 
    
    return kd1, ka1, kr1, kt1

#rk program with 4 rksteps
def rk4(ri, vi, ai, h, thetai):
    
    kd1, ka1, kr1, kt1 = rkstep(ri, vi, ai, h)

    kd2, ka2, kr2, kt2 = rkstep(ri + 0.5*kr1, vi + 0.5*kd1, ai + 0.5*ka1, h)

    kd3, ka3, kr3, kt3 = rkstep(ri + 0.5*kr2, vi + 0.5*kd2, ai + 0.5*ka2, h)

    kd4, ka4, kr4, kt4 = rkstep(ri + kr3, vi + kd3, ai + ka3, h)

    rh = ri + (1/6)*(kr1 + 2*kr2 + 2*kr3 + kr4)
    theta = thetai + (1/6)*(kt1 + 2*kt2 + 2*kt3 + kt4)
    dv = vi + (1/6)*(kd1 + 2*kd2 + 2*kd3 + kd4)
    da = ai + (1/6)*(ka1 + 2*ka2 + 2*ka3 + ka4)
    
    return rh, theta, dv, da

#function to work out orbit of the body
def orbit(ri, thetai, vi, ai, ti, tmax, k):
    
    #working out step size
    h = (tmax - ti)/k
    tpoints = np.arange(ti, tmax, h)
    
    for t in tpoints:
        
        #appending data arrays
        y.append(ri * np.sin(thetai))
        x.append(ri * np.cos(thetai))
        r.append(ri)
        theta.append(thetai)
        
        ri, thetai, vi, ai, = rk4(ri, vi, ai, h, thetai)
        
    #graphical output
    plt.plot(x, y)
    plt.xlabel("x/ au")
    plt.ylabel("y/ au")
    return plt.show

orbit(ri, thetai, vi, ai, ti, tmax, k)
