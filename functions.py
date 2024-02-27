# ========================================
# FileName: functions.py
# Author: Kade Kirsch
# Email: kirschk@my.erau.edu
# GitHub: https://github.com/Kadulous0
#
# Description: contains vital functions to
# run a orbital decay/reentry simulation.
# ========================================

# disclaimer: in no way is this code effecient

import matplotlib.pyplot as plt
import matplotlib.colors as col
import colorama
colorama.init()
import pandas as pd
import numpy as np
import scipy as sp
import time as t
import math as m
import datetime as dt
import copy
import csv

from parameters import *

# ========================================
# Found Online Functions
# ========================================

# Print iterations progress
# https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
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

# ========================================
# End of Found Online Functions
# ========================================
        
# My Own Functions
def bodyAccGrav(dist):
    if dist >= CONST_R:
        return CONST_G*((CONST_MB)/dist**2)
    else:
        return (4/3)*m.pi*CONST_G*((CONST_MB)/((4/3)*m.pi*CONST_R**3))*dist
    
def bodyForceGrav(obj, dist):
    if dist >= CONST_R:
        return CONST_G*((CONST_MB*obj.mass)/dist**2)
    else:
        return (4/3)*m.pi*CONST_G*((CONST_MB*obj.mass)/((4/3)*m.pi*CONST_R**3))*dist

def geopotentialToMetric(alt_potential): # converts geopotential altitude to geometric altitude
    return (CONST_R*alt_potential)/(CONST_R-alt_potential)

def bodyTCurve(alt): # in Kelvin, 1976 US Standard Atmosphere 
    alt = geopotentialToMetric(alt)
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
    alt = geopotentialToMetric(alt)
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
    return (CONST_GAS*bodyTCurve(alt))/(molMass(alt)*(bodyAccGrav(alt+CONST_R)))

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

def drag(object):
    velocity = absolute(object.vel[0], object.vel[1])
    altitude = absolute(object.pos[0], object.pos[1])-CONST_R
    aoa = 180
    return (0.5)*density(altitude)*(velocity**2)*object.data.coeffecientDrag(aoa, machNumber(velocity, altitude))*object.area

def lift(object):
    velocity = absolute(object.vel[0], object.vel[1])
    altitude = absolute(object.acc[0], object.acc[1])-CONST_R
    aoa = 180
    return (0.5)*density(altitude)*(velocity**2)*object.data.coeffecientLift(aoa, machNumber(velocity, altitude))*object.area

def attitude(object):
    thetaPosition = m.atan2(object.pos[1], object.pos[0])
    thetaVelocity = m.atan2(object.vel[1], object.vel[0])
    return 180-((((thetaPosition-thetaVelocity)*(180/m.pi))+90)%360)

def calcInitOrbVel(ap, pe):
    return m.sqrt(CONST_MU*((2/(ap+CONST_R))-(1/(((2*CONST_R) + ap + pe)/2))))

def absolute(x, y):
    return np.sqrt((x**2)+(y**2))

def readData(filename, column):
    dataFrame = pd.read_csv(filename)
    return dataFrame[column]

def drawPlanet(radius, res, color, ax):
    planetX=[]
    planetY=[]
    rad = np.linspace(0, 2*m.pi, res)
    for i in rad:
        planetX.append(m.cos(i)*radius)
        planetY.append(m.sin(i)*radius)
    ax.plot(planetX, planetY, color="k")
    ax.fill(planetX, planetY, color=col.to_rgba(color))

def drawTrajectoryPlot(fileName, res, colorTraj, colorBody, colorFace, colorVar, margin):
    x=readData("trajectory.csv", "Position X").array
    y=readData("trajectory.csv", "Position Y").array
    c=readData("trajectory.csv", colorVar).array

    fig, ax = plt.subplots(figsize=(res/100,res/100))
    drawPlanet(CONST_R, 4096, colorBody, ax)
    im = ax.scatter(x, y, c=c, cmap=colorTraj, s=0.125, vmin=0)

    ax.set_xlim([-1.05*(CONST_R+margin), 1.05*(CONST_R+margin)])
    ax.set_ylim([-1.05*(CONST_R+margin), 1.05*(CONST_R+margin)])
    ax.set_box_aspect(1)
    ax.set_facecolor(colorFace)
    fig.tight_layout()
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1)
    #plt.colorbar(im, ax=ax, label="Velocity (m/s)", aspect=40, shrink=0.853, pad=0.01)
    plt.savefig(fileName)

def drawPlot(fileName, label, units, color):
    time = readData("trajectory.csv", "Time").array
    data = readData("trajectory.csv", label).array

    fig, ax = plt.subplots()
    ax.plot(time, data, color=color)
    ax.set_ylabel(label + " (" + units + ")")
    ax.set_xlabel("Time (seconds)")
    ax.set_title("Time vs. " + label)
    fig.tight_layout()
    plt.savefig(fileName)
    
def estTime(startTime, currentStep, maxSteps, dT):
    currentTime = t.perf_counter()
    elapsedTime = currentTime - startTime
    if currentStep==0:
        currentStep = 1
    estRemaining = (elapsedTime * (maxSteps / currentStep)) - elapsedTime
    print("ETA:", dt.timedelta(seconds=int(estRemaining)), " Delta T:", dT, end="    \r")

def dcopy(x):
    return copy.deepcopy(x)

class Apollo():
    def __init__(self):
        self.mass = 5560 #kg
        self.diameter = 3.91 #meters
        self.area = m.pi*(self.diameter/2)**2 #meters^2

    ## SHOULD BE USING 180 DEGREE AOA (UPDATE IT)
    # Mach to CD based on "STABILITY CHARACTERISTICS OF THE APOLLO COMMAND MODULE" (NASA TN D-3890), Figures 6 & 7 (90 degree AoA)
    def coeffecientDrag(self, aoa, mach): 
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

class Orion():
    def __init__(self):
        self.mass = 9300 #kg
        self.diameter = 5.02 #meters
        self.area = m.pi*(self.diameter/2)**2 #meters^2
    
    def coeffecientDrag(self, aoa, mach): 
        pass

    def coeffecientLift(self, aoa, mach):
        pass

class Capsule():
    def __init__(self, data, ap, pe):
        self.data = data
        self.mass = self.data.mass
        self.diameter = self.data.diameter
        self.area = self.data.area
        self.pos = [0, ap+CONST_R]
        self.vel = [calcInitOrbVel(ap, pe), 0]
        self.drag = [-drag(self), 0]
        self.grav = [0, -bodyForceGrav(self, self.pos[1])]
        self.acc = [self.drag[0]/self.mass, self.grav[1]/self.mass]

        open("trajectory.csv", "w").close()
        self.write(["Time", "Timestep", "Position Error", "Attitude", "Position X", "Position Y", "Absolute Altitude", "Absolute Velocity", "Absolute Acceleration", "Absolute Drag Force"])

    def write(self, x):
        csv.writer(open("trajectory.csv", "a", newline="")).writerow(x)


class Simulation():
    def __init__(self, object, dT_min, dT_max, tLMax, isDrag):
        # object
        self.object = object
        self.dragEnabled = isDrag

        # time
        self.time = 0
        self.dT_min = dT_min
        self.dT_max = dT_max
        self.dT = dT_min
        self.tLMax = tLMax

    def tickEuler(self, obj, timestepSize):
        # update forces
        theta_pos = m.atan2(obj.pos[1],obj.pos[0])
        obj.grav[0] = -bodyForceGrav(obj, absolute(obj.pos[0], obj.pos[1]))*m.cos(theta_pos)
        obj.grav[1] = -bodyForceGrav(obj, absolute(obj.pos[0], obj.pos[1]))*m.sin(theta_pos)
        if self.dragEnabled:
            # find drag force (newtons)
            theta_vel = m.atan2(obj.vel[1],obj.vel[0])
            obj.drag[0] = -drag(obj)*m.cos(theta_vel)
            obj.drag[1] = -drag(obj)*m.sin(theta_vel)

        # compute new acceleration
        obj.acc[0] = (obj.grav[0]/obj.mass) + (obj.drag[0]/obj.mass)
        obj.acc[1] = (obj.grav[1]/obj.mass) + (obj.drag[1]/obj.mass)
            
        # update velocity
        obj.vel[0] = obj.vel[0] + (obj.acc[0]*timestepSize)
        obj.vel[1] = obj.vel[1] + (obj.acc[1]*timestepSize)
            
        # update position
        obj.pos[0] = obj.pos[0] + (obj.vel[0]*timestepSize)
        obj.pos[1] = obj.pos[1] + (obj.vel[1]*timestepSize)

    def tickHeun(self, obj, timestepSize):
        initial = dcopy(obj)

        # update forces
        theta_pos = m.atan2(obj.pos[1],obj.pos[0])
        obj.grav[0] = -bodyForceGrav(obj, absolute(obj.pos[0], obj.pos[1]))*m.cos(theta_pos)
        obj.grav[1] = -bodyForceGrav(obj, absolute(obj.pos[0], obj.pos[1]))*m.sin(theta_pos)

        if self.dragEnabled and absolute(obj.pos[0], obj.pos[1])<150000000:
            # find drag force (newtons)
            theta_vel = m.atan2(obj.vel[1],obj.vel[0])
            obj.drag[0] = -drag(obj)*m.cos(theta_vel)
            obj.drag[1] = -drag(obj)*m.sin(theta_vel)

        # compute new acceleration
        obj.acc[0] = (obj.grav[0]+obj.drag[0])/obj.mass
        obj.acc[1] = (obj.grav[1]+obj.drag[1])/obj.mass
                
        # update velocity
        obj.vel[0] = obj.vel[0] + (timestepSize*(initial.acc[0]+obj.acc[0]))/2
        obj.vel[1] = obj.vel[1] + (timestepSize*(initial.acc[1]+obj.acc[1]))/2

        # update position
        obj.pos[0] = obj.pos[0] + (timestepSize*(initial.vel[0]+obj.vel[0]))/2
        obj.pos[1] = obj.pos[1] + (timestepSize*(initial.vel[1]+obj.vel[1]))/2

    def tickVerlet(self, obj, timestepSize):
        # position
        obj.pos[0] = obj.pos[0] + (obj.vel[0]*timestepSize) + (0.5*obj.acc[0]*timestepSize*timestepSize)
        obj.pos[1] = obj.pos[1] + (obj.vel[1]*timestepSize) + (0.5*obj.acc[1]*timestepSize*timestepSize)

        # forces
        theta_pos = m.atan2(obj.pos[1],obj.pos[0])
        obj.grav[0] = -bodyForceGrav(obj, absolute(obj.pos[0], obj.pos[1]))*m.cos(theta_pos)
        obj.grav[1] = -bodyForceGrav(obj, absolute(obj.pos[0], obj.pos[1]))*m.sin(theta_pos)
        if self.dragEnabled and absolute(obj.pos[0], obj.pos[1])<150000000:
            # find drag force (newtons)
            theta_vel = m.atan2(obj.vel[1],obj.vel[0])
            obj.drag[0] = -drag(obj)*m.cos(theta_vel)
            obj.drag[1] = -drag(obj)*m.sin(theta_vel)

        # acceleration
        old_acc = dcopy(obj.acc)
        obj.acc[0] = (obj.grav[0]+obj.drag[0])/obj.mass
        obj.acc[1] = (obj.grav[1]+obj.drag[1])/obj.mass

        # velocity
        obj.vel[0] = obj.vel[0] + (0.5*(obj.acc[0]+old_acc[0])*timestepSize)
        obj.vel[1] = obj.vel[1] + (0.5*(obj.acc[1]+old_acc[1])*timestepSize)

    def run(self):
        start = t.perf_counter() # performace calculation start

        while self.time <= self.tLMax:
            printProgressBar(self.time, self.tLMax, prefix = "Simulation:", length = 50, printEnd=" ")
            estTime(start, self.time, self.tLMax, self.dT)

            # do euler and verlet method
            initial = dcopy(self.object)
            self.tickEuler(self.object, self.dT)
            euler = dcopy(self.object)
            self.object = dcopy(initial)
            self.tickVerlet(self.object, self.dT)
            verlet = dcopy(self.object)

            # compare difference in positions
            positionError = absolute((euler.pos[0]-verlet.pos[0]),(euler.pos[1]-verlet.pos[1]))

            # determine new step size
            if positionError >= CONST_TOL_HIGH:
                # decrease step size
                if self.dT > self.dT_min:
                    self.dT = 0.5*self.dT
            elif positionError <= CONST_TOL_LOW:
                # increase step size
                if self.dT < self.dT_max:
                    self.dT = 2*self.dT

            # record position
            self.object = verlet
            line = []
            line.extend((self.time, self.dT, positionError, attitude(self.object), self.object.pos[0], self.object.pos[1], absolute(self.object.pos[0], self.object.pos[1])-CONST_R, absolute(self.object.vel[0], self.object.vel[1]), absolute(self.object.acc[0], self.object.acc[1]), absolute(self.object.drag[0], self.object.drag[1])))
            self.object.write(line)

            # print information to console
            print("\n\nAltitude:", round(absolute(self.object.pos[0], self.object.pos[1])-CONST_R,2), "meters", end="   \n")
            print("Velocity:", round(absolute(self.object.vel[0], self.object.vel[1]),2), "m/s", end="   \n")
            print("G-Forces:", round(((absolute(self.object.drag[0], self.object.drag[1])/self.object.mass)/bodyAccGrav(absolute(self.object.pos[0], self.object.pos[1]))),3), "g's", end="   ")
            print("\033[F\033[F\033[F\033[F\033[F")

            # add dT to time
            self.time = self.time + self.dT

            # collision detection of object and planet
            if absolute(self.object.pos[0], self.object.pos[1]) <= CONST_R:
                self.time = CONST_TL + m.e # sets time to timelimit to force simulation loop to end
                print("\n\n\n\n\n\nobject Reached Surface at Theta of", m.atan2(self.object.pos[1],self.object.pos[0])*(180/m.pi))

        stop = t.perf_counter() # performace calculation end
        if self.time == CONST_TL + m.e:
            print("Simulation Finished Running, Time Elapsed:", dt.timedelta(seconds=int((stop-start))), "(", (stop-start),"seconds )")
        else:
            print("\n\n\n\n\n\nSimulation Finished Running, Time Elapsed:", dt.timedelta(seconds=int((stop-start))), "(", (stop-start),"seconds )")
