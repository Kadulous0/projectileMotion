# ========================================
# FileName: parameters.py
# Author: Kade Kirsch
# Email: kirschk@my.erau.edu
# GitHub: https://github.com/Kadulous0
#
# Description: contains constants that
# functions.py uses for calculations.
# ========================================

# GENERAL CONSTANTS (Rules of the Universe)
CONST_G = 6.674315e-11
CONST_GAS = 8.31446261815324

# BODY PARAMETERS (Earth)
CONST_MB = 5.972168*10**24 #kg
CONST_R = 6.371*10**6 #m
CONST_MU = CONST_G * CONST_MB

# SIMULATION PARAMETERS
CONST_DT_MIN = 1/64 #seconds, smallest step size possible (smaller number = more accurate + slower)
CONST_DT_MAX = 64 #seconds, defines max step size
CONST_TL = 172800 #seconds (avoids infinite runtimes, 90 minutes usually suffices for a low perigee LEO return (< ~50km))
CONST_TOL_LOW = 5e-1
CONST_TOL_HIGH = 2e-0

'''
Time Quick Reference
30m:00s = 1,800s
01h:00m:00s = 3,600s
01h:30m:00s = 5,400s
03h:00m:00s = 10,800s
06h:00m:00s = 21,600s
12h:00m:00s = 43,200s
01d:00h:00m:00s = 86,400s
02d:00h:00m:00s = 172,800s
03d:00h:00m:00s = 259,200s
07d:00h:00m:00s = 604,800s
14d:00h:00m:00s = 1,209,600s
'''