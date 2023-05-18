from PO import *
PO.create(['I'])
from PO import * # Need this line to access Ix, Iy, ... created by PO.create().

rho = Iz
print(rho) # print rho.txt automatically

# All four cases below are equivalent.
pulse_writing = 'exp4'

if pulse_writing == 'exp1':
    rho = rho.pulse(['I'], ['x'], [pi/2])# I90x-pulse
elif pulse_writing == 'exp2':
    rho = rho.pulse([1], ['x'], [pi/2])# I90x-pulse
elif pulse_writing == 'exp3':
    rho = rho.pulse(['I'], [0], [pi/2])# I90x-pulse
elif pulse_writing == 'exp4':
    rho = rho.pulse([1], [0], [pi/2])# I90x-pulse

print('')
print('pulse_writing is ',pulse_writing)
print('result is ', rho)