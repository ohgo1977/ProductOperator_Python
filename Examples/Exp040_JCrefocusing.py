# Exp040_JCrefocusing.m
# Keeler, J., Understanding NMR Spectroscopy (1st Ed.), Wiley, 2005.
# pp. 168, Fig. 7.14
# I: t/2-   -t/2 => cs is not refocused
# S: t/2-180-t/2 => cs is refocused
#                   jc is refocused

from PO import *
PO.create(['I', 'S'])
from PO import *

rho = Ix + Sx
print(rho)

rho = rho.cs(['I'], [oI*t/2])
rho = rho.cs(['S'], [oS*t/2])
rho = rho.jc(['IS'], [pi*J12*t/2])

rho = rho.pulse(['S'], ['x'], [pi]) # Refocusing pulse on S

# What if refocusing pulse is also applied to I.
# rho = rho.pulse(['I'], ['x'], [pi]) # Refocusing pulse on I

rho = rho.cs(['I'], [oI*t/2])
rho = rho.cs(['S'], [oS*t/2])
rho = rho.jc(['IS'], [pi*J12*t/2])

print('Done!')