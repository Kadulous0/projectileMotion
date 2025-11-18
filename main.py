# ========================================
# FileName: main.py
# Author: Kade Kirsch
# Email: kirschk@my.erau.edu
# GitHub: https://github.com/Kadulous0
#
# Description: executes functions from
# functions.py to run an orbital simulation
# or to perform calculations based on
# various functions avaible.
# ========================================

from functions import *

'''
Example Initial Condition's
Name: (apogee,perigee)

LEO Decay: (160000, 160000)
ISS Return: (420000,50000)
Moon Return: (385000000,42500)
'''

#'''
# Validate atmosphere accuracy.
printPressureAccuracy()

# Simulation Setup
capsuleType = Apollo()  # Apollo Spacecraft Capsule
object = Capsule(data=capsuleType, ap=150000, pe=32000)  # Initial Conditions
sim = Simulation(object=object, dT_min=CONST_DT_MIN, dT_max=CONST_DT_MAX, tLMax=CONST_TL, dragEnabled=True)  # Simulation Configuration

# Run Simulation
sim.run()

# Generate various graphs from simulation data
drawTrajectoryPlot("trajectory", 4096, "plasma", (57/255, 116/255, 212/255), "#111111", "Absolute Drag Force", 450000)
drawPlot("altitude", "Absolute Altitude", "meters", "black")
drawPlot("velocity", "Absolute Velocity", "m/s", "black")
drawPlot("acceleration", "Absolute Acceleration", "m/s^2", "black")
drawPlot("drag", "Absolute Drag Force", "Newtons", "black")
#'''

'''
# graphing function
x = np.arange(-10000, 130001, 100)  # bounds to plot
y = []
z = []
for i in x:
    y.append(convolve(bodyTCurve, i, 0.0003, 10000, 1000))
    z.append(bodyTCurve(i))

fig, ax = plt.subplots()
ax.plot(z, x, color="black")
ax.plot(y, x, color="tab:red", alpha=0.8)
ax.set_xlim([175, 300])
ax.set_ylim([0, 120000])
fig.tight_layout()
plt.savefig("graphs/test")
#'''