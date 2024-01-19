# ========================================
# FileName: main.py
# Author: Kade Kirsch
# Email: kade2kirsch@gmail.com
# GitHub: https://github.com/Kadulous0
#
# Description: executes functions from
# functions.py to run an orbital simulation
# or to perform calculations based on
# various functions avaible.
# ========================================

from functions import *

#'''
capsule = Projectile(ap=160000, pe=160000)
sim = Simulation(projectile=capsule, dT=CONST_DT, tLMax=CONST_TL, isDrag=True)
sim.run()
#'''
drawTrajectoryPlot("figure", 4096, "inferno", (57/255, 116/255, 212/255), 160000)
drawPlot("altitude", "Absolute Altitude", "meters", "black")
drawPlot("velocity", "Absolute Velocity", "m/s", "black")
drawPlot("acceleration", "Absolute Acceleration", "m/s^2", "black")
drawPlot("drag", "Absolute Drag Force", "Newtons", "black")
#'''