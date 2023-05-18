# Exp020_CSevolution.py
# Example of chemical shift evolution

from PO import *
PO.create(['I'])
from PO import *

rho = Iz
print(rho)
rho = rho.pulse(['I'], ['y'], [pi/2])# I90y-pulse
rho = rho.cs(['I'], [q])# CS evolution
print('Done!')