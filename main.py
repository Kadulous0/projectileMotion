from functions import *

capsule = Projectile(CONST_AP, CONST_PE)
sim = Simulation(capsule, CONST_DT, CONST_TL)
sim.run()

'''
dx = 1000
alt = np.arange(0,300001,dx)
v1=[]
v2=[]

for i in alt:
    v1.append(bodyTCurve(i))
    v2.append(np.log(bodyPCurve(i)))

fig, ax1 = plt.subplots()
color = 'blue'
ax1.set_xlabel('Temperature (K)', color=color)
ax1.set_ylabel('Altitude (m)')
ax1.plot(v1, alt, color=color)
ax1.tick_params(axis='x', labelcolor=color)

ax2 = ax1.twiny()
color = 'red'
ax2.set_xlabel('Natural Log Pressure (Pa)', color=color)
ax2.plot(v2, alt, color=color)
ax2.tick_params(axis='x', labelcolor=color)

fig.tight_layout()
plt.savefig("figure")
'''