# Exp100_INADEQUATE.py
# Levitt, M. H., Spin Dynamics(2nd Ed.), p.433.
# 2D-INADEQUATE using -45 deg phase shift

from PO import *
PO.create(['I', 'S'])
from PO import *

rho = Iz + Sz
print(rho)

States = 'cos'
if States == 'cos':
    phi = 0
elif States == 'sin':
    phi = -1/4*pi

rho = rho.pulse_phshift(['I', 'S'], [phi, phi], [3*pi/2, 3*pi/2])
rho = rho.jc(['IS'], [pi/2])
rho = rho.pulse_phshift(['I', 'S'], [phi, phi], [1*pi/2, 1*pi/2])
rho = rho.cs(['I', 'S'], [oI*t, oS*t])
rho = rho.pulse(['I', 'S'], ['y', 'y'], [1*pi/2, 1*pi/2])


phR = 0
a0_V,rho_V = rho.SigAmp(['I', 'S'], phR)

print('')
print(a0_V)
print(rho_V)
print('Done!')