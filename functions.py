import math as m
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import csv
import time as t
from parameters import *

###########################################################################
# Found Online Functions
###########################################################################

# https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

###########################################################################

# My Own Functions
def bodyForceGrav(dist): # RETURNS FORCE IN NEWTONS
    if dist >= CONST_R:
        return CONST_G*((CONST_MB*CONST_M)/dist**2)
    else:
        return (4/3)*CONST_PI*CONST_G*((CONST_MB*CONST_M)/((4/3)*CONST_PI*CONST_R**3))*dist
    
def bodyForceGravComponent(distAbsolute, distComponent):
    return CONST_G*((CONST_MB*CONST_M)/distAbsolute**3)*distComponent

def coeffecientDragCurve(mach): # Mach to CD based on "STABILITY CHARACTERISTICS OF THE APOLLO COMMAND MODULE" (NASA TN D-3890), Figures 6 & 7 (90 degree AoA)
    if mach >= 0 and mach < 0.2:
        return ((0.225-0.21)/0.2)*(mach-0.2)+0.225
    elif mach >= 0.2 and mach < 0.4:
        return ((0.25-0.225)/0.2)*(mach-0.4)+0.25
    elif mach >= 0.4 and mach < 0.7:
        return ((0.35-0.25)/0.3)*(mach-0.7)+0.35
    elif mach >= 0.7 and mach < 0.9:
        return ((0.45-0.35)/0.2)*(mach-0.9)+0.45
    elif mach >= 0.9 and mach < 1.2:
        return ((0.65-0.45)/0.3)*(mach-1.2)+0.65
    elif mach >= 1.2 and mach < 3:
        return ((0.575-0.65)/1.8)*(mach-3)+0.575
    elif mach >= 3 and mach < 4:
        return ((0.55-0.575)/1)*(mach-4)+0.55
    elif mach >= 4 and mach < 5:
        return 0.55
    elif mach >= 5 and mach < 6:
        return ((0.5625-0.55)/1)*(mach-6)+0.5625
    elif mach >= 6 and mach < 8:
        return ((0.5-0.5625)/2)*(mach-8)+0.5
    elif mach >= 8 and mach < 10:
        return ((0.55-0.5)/2)*(mach-10)+0.55
    elif mach >= 10:
        return 0.55
    else:
        return 0

def bodyTCurve(alt): # in Kelvin, 1976 US Standard Atmosphere 
    if alt >= 0 and alt < 11000:
        return ((216.65-288.15)/11000)*(alt-11000)+216.65
    elif alt >= 11000 and alt < 20000:
        return 216.65
    elif alt >= 20000 and alt < 32000:
        return ((228.65-216.65)/12000)*(alt-20000)+216.65
    elif alt >= 32000 and alt < 47000:
        return ((270.65-228.65)/15000)*(alt-32000)+228.65
    elif alt >= 47000 and alt < 51000:
        return 270.65
    elif alt >= 51000 and alt < 71000:
        return ((214.65-270.65)/20000)*(alt-51000)+270.65
    elif alt >= 71000 and alt < 84852:
        return ((187.15-214.65)/13852)*(alt-71000)+214.65
    elif alt >= 84852 and alt < 91000:
        return 187.15
    elif alt >= 91000 and alt < 107000:
        return ((202-187.15)/16000)*(alt-91000)+187.15
    elif alt >= 107000 and alt < 150000:
        return ((650-202)/43000)*(alt-107000)+202
    elif alt >= 150000 and alt < 300000:
        return ((1000-650)/150000)*(alt-150000)+650
    elif alt >= 300000:
        return 1000
    else:
        return 288

def molMass(alt): # Molecular Mass of the Air at an altitude, https://www.ngdc.noaa.gov/stp/space-weather/online-publications/miscellaneous/us-standard-atmosphere-1976/us-standard-atmosphere_st76-1562_noaa.pdf (page 15/figure 6)
    if alt < 85000:
        return 0.0289645
    elif alt >= 85000 and alt < 100000:
        return ((0.028695-0.0289645)/15000)*(alt-85000)+0.0289645
    elif alt >= 100000 and alt < 200000:
        return ((0.021368-0.028695)/100000)*(alt-100000)+0.028695
    elif alt >= 200000 and alt < 300000:
        return ((0.017750-0.021368)/100000)*(alt-200000)+0.021368
    elif alt >= 300000 and alt < 400000:
        return ((0.016043-0.017750)/100000)*(alt-300000)+0.017750
    elif alt >= 400000 and alt < 500000:
        return ((0.014379-0.016043)/100000)*(alt-400000)+0.016043
    elif alt >= 500000 and alt < 600000:
        return ((0.011608-0.014379)/100000)*(alt-500000)+0.014379
    elif alt >= 600000 and alt < 700000:
        return ((0.008083-0.011608)/100000)*(alt-600000)+0.011608
    elif alt >= 700000 and alt < 800000:
        return ((0.005702-0.008083)/100000)*(alt-700000)+0.008083
    elif alt >= 800000 and alt < 900000:
        return ((0.004526-0.005702)/100000)*(alt-800000)+0.005702
    elif alt >= 900000 and alt < 1000000:
        return ((0.004105-0.004526)/100000)*(alt-900000)+0.004526
    else:
        return 0.004105

def bodyScaleHeight(alt): # H = RT/Mg
    return (CONST_GAS*bodyTCurve(alt))/(molMass(alt)*(bodyForceGrav(alt+CONST_R)/CONST_M))

def bodyPCurve(alt): # in Pascals, huge credit to Dylan Hernquist for helping me come up with the function to properly find the pressure at an altitude
    x = np.linspace(alt, 0, 1000)
    y=[]
    for i in x:
        y.append(1/bodyScaleHeight(i))
    return 101325*m.exp(sp.integrate.trapezoid(y, x))

def density(alt):
    return (molMass(alt)*bodyPCurve(alt))/(CONST_GAS*bodyTCurve(alt))

def speedOfSound(alt):
    return m.sqrt((1.4*CONST_GAS*bodyTCurve(alt))/(molMass(alt)))

def machNumber(vel, alt):
    return vel/speedOfSound(alt)

def drag(vel, AbsVel, AbsAlt):
    return (0.5)*density(AbsAlt-CONST_R)*(vel*vel)*coeffecientDragCurve(machNumber(AbsVel, AbsAlt-CONST_R))*CONST_A

def calcInitOrbVel(ap, pe):
        return m.sqrt(CONST_G*CONST_MB*((2/(ap+CONST_R))-(1/(((2*CONST_R) + ap + pe)/2))))

def absolute(x, y):
    return np.sqrt((x**2)+(y**2))

class Projectile():
    def __init__(self, ap, pe):
        self.pos = [0, ap + CONST_R]
        self.vel = [calcInitOrbVel(ap, pe), 0]
        self.acc = [0,(-1*bodyForceGrav(ap + CONST_R))/CONST_M]
        self.drag = [drag(self.vel[0], absolute(self.vel[0], self.vel[1]),absolute(self.pos[0], self.pos[1])), 0]
        open("trajectory.csv", "w").close()
        self.write(["Time", "Positon X", "Position Y", "Absolute Altitude", "Velocity X", "Velocity Y", "Absolute Velocity", "Acceleration X", "Acceleration Y", "Absolute Acceleration", "Absolute Drag Force"])
    
    def write(self, x):
        csv.writer(open("trajectory.csv", "a", newline="")).writerow(x)


class Simulation():
    def __init__(self, projectile, dT, tLMax):
        self.projectile = projectile
        self.time = 0
        self.dT = dT
        self.tLMax = tLMax

    def run(self):
        start = t.perf_counter()
        while self.time <= self.tLMax:
            printProgressBar(self.time, self.tLMax, prefix = "Simulation:", length = 50)

            # record position
            line = []
            line.extend((self.time, self.projectile.pos[0], self.projectile.pos[1], absolute(self.projectile.pos[0], self.projectile.pos[1]), self.projectile.vel[0], self.projectile.vel[1], absolute(self.projectile.vel[0], self.projectile.vel[1]), self.projectile.acc[0], self.projectile.acc[1], absolute(self.projectile.acc[0], self.projectile.acc[1]), absolute(self.projectile.drag[0], self.projectile.drag[1])))
            self.projectile.write(line)
        
            # update position
            self.projectile.pos[0] = self.projectile.pos[0] + (self.projectile.vel[0]/(1/CONST_DT))
            self.projectile.pos[1] = self.projectile.pos[1] + (self.projectile.vel[1]/(1/CONST_DT))

            # update velocity
            self.projectile.vel[0] = self.projectile.vel[0] + (self.projectile.acc[0]/(1/CONST_DT))
            self.projectile.vel[1] = self.projectile.vel[1] + (self.projectile.acc[1]/(1/CONST_DT))

            # find drag force
            if self.projectile.vel[0] >= 0:
                a = -1
            else:
                a = 1
            
            if self.projectile.vel[1] >= 0:
                b = -1
            else:
                b = 1

            self.projectile.drag[0] = a*drag(self.projectile.vel[0], absolute(self.projectile.vel[0], self.projectile.vel[1]), absolute(self.projectile.pos[0], self.projectile.pos[1]))
            self.projectile.drag[1] = b*drag(self.projectile.vel[1], absolute(self.projectile.vel[0], self.projectile.vel[1]), absolute(self.projectile.pos[0], self.projectile.pos[1]))

            # find gravity force
            gravForceX = -1*bodyForceGravComponent(absolute(self.projectile.pos[0], self.projectile.pos[1]), self.projectile.pos[0])
            gravForceY = -1*bodyForceGravComponent(absolute(self.projectile.pos[0], self.projectile.pos[1]), self.projectile.pos[1])

            # compute new acceleration
            self.projectile.acc[0] = ((gravForceX)/CONST_M) + ((self.projectile.drag[0])/CONST_M)
            self.projectile.acc[1] = ((gravForceY)/CONST_M) + ((self.projectile.drag[1])/CONST_M)

            # add dT to time
            self.time = self.time + self.dT

            if absolute(self.projectile.pos[0], self.projectile.pos[1]) <= CONST_R:
                self.time = CONST_TL + self.dT
                print("\nProjectile Reached Surface")
        stop = t.perf_counter()
        print("\nSimulation Finished Running, Time Elapsed:", "{:.3f}".format(stop-start), "seconds")