from functions import *

capsule = Projectile(CONST_AP, CONST_PE)
sim = Simulation(capsule, CONST_DT, CONST_TL)
sim.run()

drawTrajectoryPlot("figure", 4096, "red", "blue")
drawPlot("altitude", "Absolute Altitude", "meters", "black")
drawPlot("velocity", "Absolute Velocity", "m/s", "black")
drawPlot("acceleration", "Absolute Acceleration", "m/s^2", "black")
drawPlot("drag", "Absolute Drag Force", "Newtons", "black")