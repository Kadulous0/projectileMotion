# ========================================
# FileName: parameters.py
# Author: Kade Kirsch
# Email: kade2kirsch@gmail.com
# GitHub: https://github.com/Kadulous0
#
# Description: contains constants that
# functions.py uses for calculations.
# ========================================

import math as m

# GENERAL CONSTANTS (Rules of the Universe)
CONST_G = 6.674315e-11
CONST_GAS = 8.31446261815324

# BODY PARAMETERS (Earth)
CONST_MB = 5.972168*10**24 #kg
CONST_R = 6.371*10**6 #m
CONST_MU = CONST_G * CONST_MB

# PROJECTILE PARAMETERS (based off of Saturn V command module)
CONST_M = 5560 #kg
CONST_L = 3.91 #m
CONST_A = m.pi*CONST_L**2 #m**2

# SIMULATION PARAMETERS
CONST_DT = 1/60 #seconds (smaller number = more accurate + slower)
CONST_TL = 86400 #seconds (avoids infinite runtimes, 90 minutes suffices for a low perigee (25km))