"""
Utilities
"""
from math import *
from system_of_units import *

def wait():
    raw_input("Press a key...")


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step
def lin(x,x0,y0,x1,y1):
	return y0 + (x-x0)*(y1-y0)/(x1-x0)
