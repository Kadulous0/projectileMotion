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

''''''
# simulation setup
capsuleType = Apollo()
object = Capsule(data=capsuleType,ap=150000, pe=150000)
sim = Simulation(object=object, dT_min=CONST_DT_MIN, dT_max=CONST_DT_MAX, tLMax=CONST_TL, isDrag=True)

# run simulation
sim.run()

# generate graphs from simulation data
drawTrajectoryPlot("figure", 4096, "plasma", (57/255, 116/255, 212/255), "#111111", "Absolute Drag Force", 450000)
drawPlot("altitude", "Absolute Altitude", "meters", "black")
drawPlot("velocity", "Absolute Velocity", "m/s", "black")
drawPlot("acceleration", "Absolute Acceleration", "m/s^2", "black")
drawPlot("drag", "Absolute Drag Force", "Newtons", "black")
'''
# graph of function
x = np.linspace(0, 125001, 250)
y = []
for i in x:
    y.append(bodyScaleHeight(i))

fig, ax = plt.subplots()
ax.plot(y, x, color="black")
ax.set_xlabel("Body Scale Height (m)")
ax.set_ylabel("Altitude (m)")
fig.tight_layout()
plt.savefig("scaleheight")
'''