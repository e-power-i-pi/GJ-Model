import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
import matplotlib.animation as anime
from matplotlib.animation import FuncAnimation
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors
import pandas as pd

# ------------------------------------ #

#constants
e = 1.6 * 1e-19
c = 3 * 1e8


def MagF(time, coords, P, R, B0, alpha):
    #x, y, z = coords
    r, t, p = coords #r, theta, phi
    o = 2* np.pi / P   # omega
    psi = p - (o*time)
    
    o_v = 2* np.pi / P  * np.array([cos(t), -sin(t), 0])
    
    factor = B0 * (R/r) * (R/r) * (R/r)
    
    Br = factor * ( (cos(alpha)*cos(t)) + (sin(alpha)*sin(t)*cos(psi)) )
    Bt = (factor/2) * ( (cos(alpha)*sin(t)) - (sin(alpha)*cos(t)*cos(psi)) )
    Bp = (factor/2) * sin(alpha) * sin(psi)
    
    relFac = 1 / (1 - (o*o*(r*1000)*(r*1000)*sin(t)*sin(t)/(c*c)) )
    
    return np.array([Br, Bt, Bp]), o_v, relFac

theta = np.linspace(-np.pi, np.pi, 180)
phi = np.linspace(-2*np.pi, 2*np.pi, 180)
r = np.linspace(10000, 200000, 1000)

choice = 50




# --------------------------------------------- #

# Get a simple image of charge density

fig = plt.figure()
map = Basemap(projection = "ortho", lat_0 = 40, lon_0 = 30)

map.drawmeridians(np.arange(0, 360, 30))
map.drawparallels(np.arange(0, 180, 30))

nGJ = np.empty((theta.size, phi.size))

for i in range(theta.size):
    for j in range(phi.size):
        period = 1
        omega = 2*np.pi / period
        x1, x2, rf = MagF(0, [choice, theta[i], phi[j]], period, 10, 1e10, 45 * np.pi/180)
        nGJ[i][j] = (2 * np.dot(x1,x2) / e) * rf   # (2*np.pi*c))
        
long, lat = np.meshgrid(180 / np.pi * theta, 180 / np.pi * phi)
x,y = map(long, lat)

cs = map.contourf(x,y, np.abs(nGJ), 15, cmap = "winter")

cb = map.colorbar(cs)

plt.show()





# ------------------------------- #

# For |nGJ| vs time

fig = plt.figure()

map = Basemap(projection = "ortho", lat_0 = 40, lon_0 = 30)

def init():
    fig
    return fig 


theta = np.linspace(-np.pi, np.pi, 180)
phi = np.linspace(-2*np.pi, 2*np.pi, 180)
long, lat = np.meshgrid(180 / np.pi * theta, 180 / np.pi * phi)
x,y = map(long, lat)

def animation(clock):
    nGJ = np.empty((theta.size, phi.size))

    plt.title("Absolute GJ charge density vs time - (P = 1s, r = 50km, $B = 10^{14}$ G, "+r"${\alpha} = 45^{\circ}$)"+
        "\n T = %.3f s" % (clock/20),
    pad = 15, fontsize = 10)

    for i in range(theta.size):
        for j in range(phi.size):
            period = 1
            omega = 2*np.pi / period
            x1, x2, rf = MagF(clock/20, [choice, theta[i], phi[j]], period, 10, 1e10, 45 * np.pi/180)
            nGJ[i][j] = (2 * np.dot(x1,x2) / e) * rf

    map.drawmeridians(np.arange(0, 360, 30))
    map.drawparallels(np.arange(0, 180, 30))
    cs = map.contourf(x,y, np.abs(nGJ), 15, cmap = "winter")
    map.colorbar(cs)
        
ani = FuncAnimation(fig, animation, init_func = init, frames = 21, interval = 1)
ani.save("|nGJ| vs time.gif", writer='imagemagick', fps = 2)

plt.show()





# ----------------------------------------- #

# For |nGJ| vs alpha along with colourbar

fig = plt.figure(figsize = (12,6))
fig.suptitle("Effect of GJ charge density on the observed power.")

plt.subplot(1,2,1)

path = "/Users/user/Desktop/Results/Interpolation data 2.0/"

alpha = np.array([20, 30, 40, 60, 80])
theta = np.array([20, 40, 50, 70, 80])
avgP = np.zeros((5,5))
p,q = 0,0

for i in alpha:
    for j in theta:
        f = pd.read_csv(path+"Results_"+str(i)+"/B=14_theta="+str(j)+"_deg_plas_alpha_"+str(i)+"_deg_profile.dat",
            header = None)

        avgP[p][q] = np.mean(np.array(f))
        q = (q+1)%5
    p = (p+1)%5

power = plt.imshow(avgP, origin = "lower", norm = "log", interpolation = "bicubic", cmap = "viridis")
plt.colorbar(power)
plt.title(r"$Luminosity$ vs $\theta, \alpha$")
plt.xlabel(r"$\theta$ (in degrees)")
plt.ylabel(r"$\alpha$ (in degrees)")
plt.xticks([0,1,2,3,4], labels = theta)
plt.yticks([0,1,2,3,4], labels = alpha)



# ----------------------------------------- #

plt.subplot(1,2,2)

map = Basemap(projection = "ortho", lat_0 = 30, lon_0 = 30)

def init():
    fig
    return fig 


theta = np.linspace(-np.pi, np.pi, 180)
phi = np.linspace(-2*np.pi, 2*np.pi, 180)
long, lat = np.meshgrid(180 / np.pi * theta, 180 / np.pi * phi)
x,y = map(long, lat)

# ************ For B and P axis (this and the bottom labelled code can be removed to get rid of these axes)

x10, y10 = map(15, 90)
x20, y20 = map(30, 15)
x30, y30 = map(0, 88)
line1, = plt.plot([x10, x20], [y10, y20], color = "black", linewidth = 2) # Rotation Axis

line2, =plt.plot([x10, x20], [y10, y20], color = "black", linewidth = 2) # Magnetic Axis

plt.annotate("B", xy = (x10, y10), color = "black")
a = plt.annotate("P", xy = (x30, y30), color = "white")

# ************

def animation(clock):
    nGJ = np.empty((theta.size, phi.size))

    plt.title("Absolute GJ charge density vs misalignment angle\n(R=10km, r = 50km, $B = 10^{14}$ G)\n"+
        r"${\alpha}$"+" = %d°" % (5 * clock), 
    pad = 15, fontsize = 9)

    for i in range(theta.size):
        for j in range(phi.size):
            period = 1
            omega = 2*np.pi / period
            x1, x2, rf = MagF(0, [choice, theta[i], phi[j]], period, 10, 1e10, 5 * clock * np.pi/180)
            nGJ[i][j] = (2 * np.dot(x1,x2) / e) * rf   # (2*np.pi*c))

    map.drawmeridians(np.arange(0, 360, 30))
    map.drawparallels(np.arange(0, 180, 30))
    cs = map.contourf(x,y, np.abs(nGJ), 15,extend = "both", levels = np.linspace(10**(25.5), 10**(27.8), 40),
        cmap = "viridis_r")
                     #norm=colors.LogNorm(vmin = 1e25, vmax = 1e28))
    cb = map.colorbar(cs, ticks = [1e26, 1e27, 10**27.5, 10**27.8])
    cb.set_ticklabels([r"$10^{26}$", r"$10^{27}$", r"$10^{27.5}$", r"$10^{27.8}$"])

    # ************ For B and P axis 

    # Rotation Axis
    x1, y1 = map(0 , 90 - (90/18 * clock))
    x2, y2 = map(30, 15/18 * clock)
    a.set_position((x1, y1))

    line1.set_data([x1, x20], [y1, y20])

    return line1,  

    # ************ 

        
ani = FuncAnimation(fig, animation, init_func = init, frames = 19, interval = 7) 
ani.save("|nGJ| vs alpha (with B&P axis).gif", writer='imagemagick', fps = 2)

plt.show()