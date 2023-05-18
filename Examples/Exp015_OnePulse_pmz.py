# Exp015_OnePulse_pmz.py
# One pulse experiment in the lowering/raising operator basis.

from PO import *
PO.create(['I'])
from PO import *

rho = Iz
rho = rho.xyz2pmz()
print(rho)
rho = rho.pulse(['I'], ['y'], [pi/2])# I90y-pulse
