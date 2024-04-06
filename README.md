# projectileMotion
2-D customizable simulation of atmospheric re-entry into Earths (maybe others too in the future) atmosphere. Using a custom implementation of the 1976 Standard Atmosphere model along with capsule data from the Apollo Command Module (and more in the future).

## Installation & Usage
1. Download the .zip and extract to a directory of your choosing.
2. Edit the "main.py" file
3. Customize the simulation by changing the apoapsis (ap) and periapsis (pe) located on line 25
4. (optional:) On line 26, change "isDrag=" to False to disable drag
5. (optional:) Edit the "parameters.py" file and change the "CONST_TL=" to a new preferred time limit in seconds.
6. Run "main.py" and wait til the capsule hits the Earth or the time limit is reached.
7. Check the directory for the graphs of the simulation.

## Work In Progress Features
- Multiple Types of Capsules
- Body Lift
- Improvements to Capsule Data
- Improvements to Atmospheric Data
- Increasing Simulation Speed

## Possible Future Plans
- Parachute Deployments
- Other Planets
- Sun/Moon Dynamics
- Weather
- 3-D Overhaul
