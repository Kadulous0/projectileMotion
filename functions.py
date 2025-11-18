# ========================================
# FileName: functions.py
# Author: Kade Kirsch
# Email: kirschk@my.erau.edu
# GitHub: https://github.com/Kadulous0
#
# Description: contains vital functions to
# run an orbital decay/reentry simulation.
# ========================================

import matplotlib.pyplot as plt
import matplotlib.colors as col
import colorama
import pandas as pd
import numpy as np
import scipy as sp
import time as t
import math as m
import datetime as dt
import copy
import csv

from parameters import *

# initialize
colorama.init()


def gaussianIntegral(x, a):  # integral from -inf to inf of this is always == 1, a = max value
    return a * m.exp(-(x ** 2) / (2 * ((1 / (a * m.sqrt(2 * m.pi))) ** 2)))


def gaussianFunction(x, a, b, c):  # a = max value, b = offset, c = standard deviation
    return a * m.exp(-((x - b) ** 2) / (2 * (c ** 2)))


def convolve(function, value, peak, bounds, quality):
    x = np.linspace(-bounds, bounds, quality)
    y = []
    for j in x:
        y.append(function(value - j) * gaussianIntegral(j, peak))
    return sp.integrate.trapezoid(y, x)


def bodyAccGrav(dist):
    if dist >= CONST_R:
        return CONST_MU / dist ** 2
    else:
        return (CONST_MU / (CONST_R ** 3)) * dist


def bodyForceGrav(obj, dist):
    if dist >= CONST_R:
        return (CONST_MU * obj.mass) / dist ** 2
    else:
        return ((CONST_MU * obj.mass) / (CONST_R ** 3)) * dist


def geopotentialToMetric(alt_potential):  # converts geopotential altitude to geometric altitude
    return (CONST_R * alt_potential) / (CONST_R - alt_potential)


def bodyTCurve(alt):  # in Kelvin, 1976 US Standard Atmosphere
    alt = geopotentialToMetric(alt)
    if 0 <= alt < 11000:
        return ((216.65 - 288.15) / 11000) * (alt - 11000) + 216.65
    elif 11000 <= alt < 20000:
        return 216.65
    elif 20000 <= alt < 32000:
        return ((228.65 - 216.65) / 12000) * (alt - 20000) + 216.65
    elif 32000 <= alt < 47000:
        return ((270.65 - 228.65) / 15000) * (alt - 32000) + 228.65
    elif 47000 <= alt < 51000:
        return 270.65
    elif 51000 <= alt < 71000:
        return ((214.65 - 270.65) / 20000) * (alt - 51000) + 270.65
    elif 71000 <= alt < 84852:
        return ((187.15 - 214.65) / 13852) * (alt - 71000) + 214.65
    elif 84852 <= alt < 91000:
        return 187.15
    elif 91000 <= alt < 107000:
        return ((202 - 187.15) / 16000) * (alt - 91000) + 187.15
    elif 107000 <= alt < 150000:
        return ((650 - 202) / 43000) * (alt - 107000) + 202
    elif 150000 <= alt < 300000:
        return ((1000 - 650) / 150000) * (alt - 150000) + 650
    elif alt >= 300000:
        return 1000
    else:
        return ((216.65 - 288.15) / 11000) * (alt - 11000) + 216.65


def molMass(alt):  # Molecular Mass of the Air at an altitude,
    # https://www.ngdc.noaa.gov/stp/space-weather/online-publications/miscellaneous/us-standard-atmosphere-1976/us-standard-atmosphere_st76-1562_noaa.pdf (page 15/figure 6)
    alt = geopotentialToMetric(alt)
    if alt < 85000:
        return 0.0289645
    elif 85000 <= alt < 100000:
        return ((0.028695 - 0.0289645) / 15000) * (alt - 85000) + 0.0289645
    elif 100000 <= alt < 200000:
        return ((0.021368 - 0.028695) / 100000) * (alt - 100000) + 0.028695
    elif 200000 <= alt < 300000:
        return ((0.017750 - 0.021368) / 100000) * (alt - 200000) + 0.021368
    elif 300000 <= alt < 400000:
        return ((0.016043 - 0.017750) / 100000) * (alt - 300000) + 0.017750
    elif 400000 <= alt < 500000:
        return ((0.014379 - 0.016043) / 100000) * (alt - 400000) + 0.016043
    elif 500000 <= alt < 600000:
        return ((0.011608 - 0.014379) / 100000) * (alt - 500000) + 0.014379
    elif 600000 <= alt < 700000:
        return ((0.008083 - 0.011608) / 100000) * (alt - 600000) + 0.011608
    elif 700000 <= alt < 800000:
        return ((0.005702 - 0.008083) / 100000) * (alt - 700000) + 0.008083
    elif 800000 <= alt < 900000:
        return ((0.004526 - 0.005702) / 100000) * (alt - 800000) + 0.005702
    elif 900000 <= alt < 1000000:
        return ((0.004105 - 0.004526) / 100000) * (alt - 900000) + 0.004526
    else:
        return 0.004105


def bodyScaleHeight(alt):  # H = RT/Mg
    return (CONST_GAS * bodyTCurve(alt)) / (molMass(alt) * (bodyAccGrav(alt + CONST_R)))


def bodyPCurve(alt):  # in Pascals, huge credit to Dylan Hernquist for helping me come up with the function to properly
    # find the pressure at an altitude
    x = np.linspace(alt, 0, 1000)
    y = []
    for i in x:
        y.append(1 / bodyScaleHeight(i))
    return 101325 * m.exp(sp.integrate.trapezoid(y, x))


def density(alt):
    return (molMass(alt) * bodyPCurve(alt)) / (CONST_GAS * bodyTCurve(alt))


def speedOfSound(alt):
    return m.sqrt((1.4 * CONST_GAS * bodyTCurve(alt)) / (molMass(alt)))


def machNumber(vel, alt):
    return vel / speedOfSound(alt)


def drag(object):
    velocity = absolute(object.vel[0], object.vel[1])
    altitude = absolute(object.pos[0], object.pos[1]) - CONST_R
    aoa = 0
    return 0.5 * density(altitude) * (velocity ** 2) * object.data.coefficientDrag(aoa, machNumber(velocity, altitude)) * object.area


def lift(object):
    velocity = absolute(object.vel[0], object.vel[1])
    altitude = absolute(object.pos[0], object.pos[1]) - CONST_R
    aoa = 0
    return 0.5 * density(altitude) * (velocity ** 2) * object.data.coefficientLift(aoa, machNumber(velocity, altitude)) * object.area


def attitude(object):
    thetaPosition = m.atan2(object.pos[1], object.pos[0])
    thetaVelocity = m.atan2(object.vel[1], object.vel[0])
    return 180 - ((((thetaPosition - thetaVelocity) * (180 / m.pi)) + 90) % 360)


def calcInitOrbVel(ap, pe):
    return m.sqrt(CONST_MU * ((2 / (ap + CONST_R)) - (1 / (((2 * CONST_R) + ap + pe) / 2))))


def absolute(x, y):
    return np.sqrt((x ** 2) + (y ** 2))


def readData(filename, column):
    dataFrame = pd.read_csv(filename)
    return dataFrame[column]


def drawPlanet(radius, res, color, ax):
    planetX = []
    planetY = []
    rad = np.linspace(0, 2 * m.pi, res)
    for i in rad:
        planetX.append(m.cos(i) * radius)
        planetY.append(m.sin(i) * radius)
    ax.plot(planetX, planetY, color="k")
    ax.fill(planetX, planetY, color=col.to_rgba(color))


def drawTrajectoryPlot(fileName, res, colorTraj, colorBody, colorFace, colorVar, margin):
    x = readData("csvfiles/trajectory.csv", "Position X").array
    y = readData("csvfiles/trajectory.csv", "Position Y").array
    c = readData("csvfiles/trajectory.csv", colorVar).array

    fig, ax = plt.subplots(figsize=(res / 100, res / 100))
    drawPlanet(CONST_R, 4096, colorBody, ax)
    im = ax.scatter(x, y, c=c, cmap=colorTraj, s=0.125, vmin=0)

    ax.set_xlim([-1.05 * (CONST_R + margin), 1.05 * (CONST_R + margin)])
    ax.set_ylim([-1.05 * (CONST_R + margin), 1.05 * (CONST_R + margin)])
    ax.set_box_aspect(1)
    ax.set_facecolor(colorFace)
    fig.tight_layout()
    fig.subplots_adjust(left=0, bottom=0, right=1, top=1)
    # plt.colorbar(im, ax=ax, label="Velocity (m/s)", aspect=40, shrink=0.853, pad=0.01)
    plt.savefig("graphs/" + fileName)
    print(f"Trajectory plot saved to /graphs/{fileName}.png")


def drawPlot(fileName, label, units, color):
    time = readData("csvfiles/trajectory.csv", "Time").array
    data = readData("csvfiles/trajectory.csv", label).array

    fig, ax = plt.subplots()
    ax.plot(time, data, color=color)
    ax.set_ylabel(label + " (" + units + ")")
    ax.set_xlabel("Time (seconds)")
    ax.set_title("Time vs. " + label)
    fig.tight_layout()
    plt.savefig("graphs/" + fileName)
    print(f"Plot '{fileName}' saved to /graphs/{fileName}.png")


def printPressureAccuracy():
    print("Altitude : Percent Error")
    print(f"11 km : {(bodyPCurve(11000) - 22632.1) / 22632.1 * 100:.2f} %")
    print(f"20 km : {(bodyPCurve(20000) - 5474.89) / 5474.89 * 100:.2f} %")
    print(f"32 km : {(bodyPCurve(32000) - 868.019) / 868.019 * 100:.2f} %")
    print(f"47 km : {(bodyPCurve(47000) - 110.906) / 110.906 * 100:.2f} %")
    print(f"51 km : {(bodyPCurve(51000) - 66.9389) / 66.9389 * 100:.2f} %")
    print(f"71 km : {(bodyPCurve(71000) - 3.95642) / 3.95642 * 100:.2f} %\n")

def estTime(startTime, currentStep, maxSteps, dT):
    currentTime = t.perf_counter()
    elapsedTime = currentTime - startTime
    if currentStep == 0:
        currentStep = 1
    estRemaining = (elapsedTime * (maxSteps / currentStep)) - elapsedTime
    print("ETA:", dt.timedelta(seconds=int(estRemaining)), " Delta T:", dT, end="    \r")


def dcopy(x):
    return copy.deepcopy(x)


class Apollo:
    def __init__(self):
        self.mass = 5560  # kg
        self.diameter = 3.91  # meters
        self.area = m.pi * (self.diameter / 2) ** 2  # meters^2


    # Mach to CD based on "STABILITY CHARACTERISTICS OF THE APOLLO COMMAND MODULE" (NASA TN D-3890), Figures 6 & 7 (90 degree AoA)
    def coefficientDrag(self, aoa, mach):
        # TODO Reimplement using AOA=180 not AOA=90, drag data is currently incorrect
        # TODO Implement AOA dependant drag for APOLLO
        # TODO Parachute simulation
        if 0 <= mach < 0.2:
            return ((0.225 - 0.21) / 0.2) * (mach - 0.2) + 0.225
        elif 0.2 <= mach < 0.4:
            return ((0.25 - 0.225) / 0.2) * (mach - 0.4) + 0.25
        elif 0.4 <= mach < 0.7:
            return ((0.35 - 0.25) / 0.3) * (mach - 0.7) + 0.35
        elif 0.7 <= mach < 0.9:
            return ((0.45 - 0.35) / 0.2) * (mach - 0.9) + 0.45
        elif 0.9 <= mach < 1.2:
            return ((0.65 - 0.45) / 0.3) * (mach - 1.2) + 0.65
        elif 1.2 <= mach < 3:
            return ((0.575 - 0.65) / 1.8) * (mach - 3) + 0.575
        elif 3 <= mach < 4:
            return ((0.55 - 0.575) / 1) * (mach - 4) + 0.55
        elif 4 <= mach < 5:
            return 0.55
        elif 5 <= mach < 6:
            return ((0.5625 - 0.55) / 1) * (mach - 6) + 0.5625
        elif 6 <= mach < 8:
            return ((0.5 - 0.5625) / 2) * (mach - 8) + 0.5
        elif 8 <= mach < 10:
            return ((0.55 - 0.5) / 2) * (mach - 10) + 0.55
        elif mach >= 10:
            return 0.55
        else:
            return 0.225

    def coefficientLift(self, aoa, mach):
        # TODO Implement Lift for APOLLO
        return 0


class Orion:
    def __init__(self):
        self.mass = 9300  # kg
        self.diameter = 5.0292  # meters
        self.area = m.pi * (self.diameter / 2) ** 2  # meters^2

    def coefficientDrag(self, aoa, mach):
        # TODO Implement drag for ORION
        return 1.5

    def coefficientLift(self, aoa, mach):
        # TODO Implement lift for ORION
        return 0


class Capsule:
    def __init__(self, data, ap, pe):
        self.data = data
        self.mass = self.data.mass
        self.diameter = self.data.diameter
        self.area = self.data.area
        self.pos = [0, ap + CONST_R]
        self.vel = [calcInitOrbVel(ap, pe), 0]
        self.drag = [-drag(self), 0]
        self.lift = [0, lift(self)]
        self.grav = [0, -bodyForceGrav(self, self.pos[1])]
        self.acc = [self.drag[0] / self.mass, (self.grav[1] + self.lift[1]) / self.mass]

        open("csvfiles/trajectory.csv", "w").close()
        self.write(["Time", "Timestep", "Position Error", "Attitude", "Position X", "Position Y", "Absolute Altitude",
                    "Absolute Velocity", "Absolute Acceleration", "Absolute Drag Force"])

    def write(self, x):
        csv.writer(open("csvfiles/trajectory.csv", "a", newline="")).writerow(x)


class Simulation:
    def __init__(self, object, dT_min, dT_max, tLMax, dragEnabled):
        # object
        self.object = object
        self.dragEnabled = dragEnabled

        # time
        self.time = 0
        self.dT_min = dT_min
        self.dT_max = dT_max
        self.dT = dT_min
        self.tLMax = tLMax

    def tickEuler(self, obj, timestepSize):
        # update forces
        theta_pos = m.atan2(obj.pos[1], obj.pos[0])
        obj.grav[0] = -bodyForceGrav(obj, absolute(obj.pos[0], obj.pos[1])) * m.cos(theta_pos)
        obj.grav[1] = -bodyForceGrav(obj, absolute(obj.pos[0], obj.pos[1])) * m.sin(theta_pos)
        if self.dragEnabled:
            # find drag force (newtons)
            theta_vel = m.atan2(obj.vel[1], obj.vel[0])
            obj.drag[0] = -drag(obj) * m.cos(theta_vel)
            obj.drag[1] = -drag(obj) * m.sin(theta_vel)
            # obj.lift[0] = lift(obj) * m.cos(theta_vel + 90)
            # obj.lift[1] = lift(obj) * m.sin(theta_vel + 90)

        # compute new acceleration
        obj.acc[0] = (obj.grav[0] + obj.drag[0] + obj.lift[0]) / obj.mass
        obj.acc[1] = (obj.grav[1] + obj.drag[1] + obj.lift[1]) / obj.mass

        # update velocity
        obj.vel[0] = obj.vel[0] + (obj.acc[0] * timestepSize)
        obj.vel[1] = obj.vel[1] + (obj.acc[1] * timestepSize)

        # update position
        obj.pos[0] = obj.pos[0] + (obj.vel[0] * timestepSize)
        obj.pos[1] = obj.pos[1] + (obj.vel[1] * timestepSize)

    def tickVerlet(self, obj, timestepSize):
        # position
        obj.pos[0] = obj.pos[0] + (obj.vel[0] * timestepSize) + (0.5 * obj.acc[0] * timestepSize * timestepSize)
        obj.pos[1] = obj.pos[1] + (obj.vel[1] * timestepSize) + (0.5 * obj.acc[1] * timestepSize * timestepSize)

        # forces
        theta_pos = m.atan2(obj.pos[1], obj.pos[0])
        obj.grav[0] = -bodyForceGrav(obj, absolute(obj.pos[0], obj.pos[1])) * m.cos(theta_pos)
        obj.grav[1] = -bodyForceGrav(obj, absolute(obj.pos[0], obj.pos[1])) * m.sin(theta_pos)
        if self.dragEnabled and absolute(obj.pos[0], obj.pos[1]) < 150000000:
            # find drag force (newtons)
            theta_vel = m.atan2(obj.vel[1], obj.vel[0])
            obj.drag[0] = -drag(obj) * m.cos(theta_vel)
            obj.drag[1] = -drag(obj) * m.sin(theta_vel)
            obj.lift[0] = lift(obj) * -m.sin(theta_vel)
            obj.lift[1] = lift(obj) * m.cos(theta_vel)

        # acceleration
        old_acc = dcopy(obj.acc)
        obj.acc[0] = (obj.grav[0] + obj.drag[0] + obj.lift[0]) / obj.mass
        obj.acc[1] = (obj.grav[1] + obj.drag[1] + obj.lift[1]) / obj.mass

        # velocity
        obj.vel[0] = obj.vel[0] + (0.5 * (obj.acc[0] + old_acc[0]) * timestepSize)
        obj.vel[1] = obj.vel[1] + (0.5 * (obj.acc[1] + old_acc[1]) * timestepSize)

    def run(self):
        start = t.perf_counter()  # performance calculation start

        while self.time < self.tLMax:
            # adaptive timestepping
            # do euler and verlet method
            initial = dcopy(self.object)
            self.tickEuler(self.object, self.dT)
            euler = dcopy(self.object)
            self.object = dcopy(initial)
            self.tickVerlet(self.object, self.dT)
            verlet = dcopy(self.object)
    
            # compare difference in positions between euler and verlet
            positionError = absolute((euler.pos[0] - verlet.pos[0]), (euler.pos[1] - verlet.pos[1]))
    
            # determine new step size based on 'error'
            if positionError >= CONST_TOL_HIGH:
                # decrease step size
                if self.dT > self.dT_min:
                    self.dT = 0.5 * self.dT
            elif positionError <= CONST_TOL_LOW:
                # increase step size
                if self.dT < self.dT_max:
                    self.dT = 2 * self.dT
            
            # record position from verlet
            self.object = verlet

            line = []
            line.extend((self.time, self.dT, 0, attitude(self.object), self.object.pos[0], self.object.pos[1],
                         absolute(self.object.pos[0], self.object.pos[1]) - CONST_R,
                         absolute(self.object.vel[0], self.object.vel[1]),
                         absolute(self.object.acc[0], self.object.acc[1]),
                         absolute(self.object.drag[0], self.object.drag[1])))
            self.object.write(line)

            # collision detection of object and planet
            if absolute(self.object.pos[0], self.object.pos[1]) <= CONST_R:
                print(f"Object Reached Surface w/ Velocity of {absolute(self.object.vel[0], self.object.vel[1]):.2f} m/s")
                stop = t.perf_counter()  # performance calculation end
                print(f"Simulation Finished Running, Time Elapsed: {dt.timedelta(seconds=int((stop - start)))}\n")
                break

            # add dT to time
            self.time = self.time + self.dT
