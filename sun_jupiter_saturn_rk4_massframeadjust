import numpy as np 
import matplotlib.pyplot as plt

#intial conditions
g = 6.67e-11/((1.498e11)**3)

m1 = 1.989e30
m2 = 1.898e27
m3 = 5.683e26

x1 = 0
x2 = 5.2
x3 = 9.555                         

y1 = 0
y2 = 0
y3 = 0

#function to work our actual distance between bodies from x and y co-ords
def r(x, y):
    r = (x**2 + y**2)**0.5
    return r

r1 = r(x1, y1)
r2 = r(x2, y2)
r3 = r(x3, y3)

dx1 = 0
dx2 = 0
dx3 = 0
    
#starting velocities in y
dy1 = 0
dy2 = ((g*m1)/(r2))**0.5
dy3 = ((g*m1)/(r3))**0.5

t1 = 0
tmax = 10000000000
h = 100000

#arrays for time and the co-ordinates of the planet
t = []
xs = []
ys = []
xj = []
yj = []
xsa = []
ysa = []

#function acceleration in different directions
def dv(m1, m2, z1, z2, z3, r1, r2):
    dv = ((g*m1)/(r1**3))*(z2-z1) + ((g*m2)/(r2**3))*(z3-z1)
    return dv

#function to work out the centre of mass
def com(a1, a2, a3):
    com = (m1/(m1+m2+m3))*(a1) + (m2/(m1+m2+m3))*(a2) + (m3/(m1+m2+m3))*(a3)
    return com

#individual rkstep
def rkstep(x1, x2, x3, y1, y2, y3, dx1, dx2, dx3, dy1, dy2, dy3, h):
    #working out distances between bodies
    r12 = r((x2-x1), (y2-y1))
    r13 = r((x3-x1), (y3-y1))
    r23 = r((x3-x2), (y3-y2))
    
    kdx1 = h*dv(m2, m3, x1, x2, x3, r12, r13)
    kdy1 = h*dv(m2, m3, y1, y2, y3, r12, r13)
    kx1 = h*dx1
    ky1 = h*dy1
    
    kdx2 = h*dv(m1, m3, x2, x1, x3, r12, r23)
    kdy2 = h*dv(m1, m3, y2, y1, y3, r12, r23)
    kx2 = h*dx2
    ky2 = h*dy2
    
    kdx3 = h*dv(m1, m2, x3, x1, x2, r13, r23)
    kdy3 = h*dv(m1, m2, y3, y1, y2, r13, r13)
    kx3 = h*dx3
    ky3 = h*dy3
    
    return kdx1, kdy1, kx1, ky1, kdx2, kdy2, kx2, ky2, kdx3, kdy3, kx3, ky3

#program to work out the orbit of the bodies
def orbit(t1, x1, x2, x3, y1, y2, y3, dx1, dx2, dx3, dy1, dy2, dy3, h):
    #while loop so rk4 iterates until tmax
    while t1<tmax:
        
        #adjusting the position and velocities of the bodies to be in the centre of mass frame
        x1 = x1 - com(x1, x2, x3)
        x2 = x2 - com(x1, x2, x3)
        x3 = x3 - com(x1, x2, x3)
        
        y1 = y1 - com(y1, y2, y3)
        y2 = y2 - com(y1, y2, y3)
        y3 = y3 - com(y1, y2, y3)
        
        dx1 = dx1 - com(dx1, dx2, dx3)
        dx2 = dx2 - com(dx1, dx2, dx3)
        dx3 = dx3 - com(dx1, dx2, dx3)
        
        dy1 = dy1 - com(dy1, dy2, dy3)
        dy2 = dy2 - com(dy1, dy2, dy3)
        dy3 = dy3 - com(dy1, dy2, dy3)
        
        #appending the data arrays
        xs.append(x1)
        ys.append(y1)
        t.append(t1)
        xj.append(x2)
        yj.append(y2)
        xsa.append(x3)
        ysa.append(y3)
        
        kdx1, kdy1, kx1, ky1, kdx2, kdy2, kx2, ky2, kdx3, kdy3, kx3, ky3 = rkstep(x1, x2, x3, y1, y2, y3, dx1, dx2, dx3, dy1, dy2, dy3, h)
        
        kdx4, kdy4, kx4, ky4, kdx5, kdy5, kx5, ky5, kdx6, kdy6, kx6, ky6 = rkstep(x1 + 0.5*kx1, x2 + 0.5*kx2, x3 + 0.5*kx3, y1 + 0.5*ky1, y2 + 0.5*ky2, y3 + 0.5*ky3, dx1 + 0.5*kdx1, dx2 + 0.5*kdx2, dx3 + 0.5*kdx3, dy1 + 0.5*kdy1, dy2 + 0.5*kdy2, dy3 + 0.5*kdy3, h)

        kdx7, kdy7, kx7, ky7, kdx8, kdy8, kx8, ky8, kdx9, kdy9, kx9, ky9 = rkstep(x1 + 0.5*kx4, x2 + 0.5*kx5, x3 + 0.5*kx6, y1 + 0.5*ky4, y2 + 0.5*ky5, y3 + 0.5*ky6, dx1 + 0.5*kdx4, dx2 + 0.5*kdx5, dx3 + 0.5*kdx6, dy1 + 0.5*kdy4, dy2 + 0.5*kdy5, dy3 + 0.5*kdy6, h)

        kdx10, kdy10, kx10, ky10, kdx11, kdy11, kx11, ky11, kdx12, kdy12, kx12, ky12 = rkstep(x1 + kx7, x2 + kx8, x3 + kx9, y1 + ky7, y2 + ky8, y3 + ky9, dx1 + kdx7, dx2 + kdx8, dx3 + kdx9, dy1 + kdy7, dy2 + kdy8, dy3 + kdy9, h)
        
        x1 = x1 + (1/6)*(kx1 + 2*kx4 + 2*kx7 + kx10)
        y1 = y1 + (1/6)*(ky1 + 2*ky4 + 2*ky7 + ky10)
        dx1 = dx1 + (1/6)*(kdx1 + 2*kdx4 + 2*kdx7 + kdx10)
        dy1 = dy1 + (1/6)*(kdy1 + 2*kdy4 + 2*kdy7 + kdy10)

        x2 = x2 + (1/6)*(kx2 + 2*kx5 + 2*kx8 + kx11)
        y2 = y2 + (1/6)*(ky2 + 2*ky5 + 2*ky8 + ky11)
        dx2 = dx2 + (1/6)*(kdx2 + 2*kdx5 + 2*kdx8 + kdx11)
        dy2 = dy2 + (1/6)*(kdy2 + 2*kdy5 + 2*kdy8 + kdy11)

        x3 = x3 + (1/6)*(kx3 + 2*kx6 + 2*kx9 + kx12)
        y3 = y3 + (1/6)*(ky3 + 2*ky6 + 2*ky9 + ky12)
        dx3 = dx3 + (1/6)*(kdx3 + 2*kdx6 + 2*kdx9 + kdx12)
        dy3 = dy3 + (1/6)*(kdy3 + 2*kdy6 + 2*kdy9 + kdy12)
        
        t1 = t1 + h

    #graphical output
    plt.plot(xs,ys)
    plt.plot(xj, yj)
    plt.plot(xsa, ysa)
    plt.xlabel("x/ au")
    plt.ylabel("y/ au")
    return plt.show

orbit(t1, x1, x2, x3, y1, y2, y3, dx1, dx2, dx3, dy1, dy2, dy3, h)
