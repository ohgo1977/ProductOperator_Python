# Exp030_JCevolution.py
# Example of J-coupling evolution

from PO import *
PO.create(['I', 'S'])
from PO import *

rho = Ix
print(rho)
rho = rho.jc(['IS'], [pi*J12*t])# J-coupling eovlution
print('Done!')