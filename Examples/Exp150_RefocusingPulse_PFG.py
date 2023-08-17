# Exp150_RefocusingPulse_PFG.m
# Keeler, J., Understanding NMR Spectroscopy, p. 406, 11.12.3
# Gradient G - 180+d pulse - Gradient G 
# The selection of p => -p pathway.
# "Cleaning up" the results of an imperfect 180 pulse.

from PO import *

pfg_switch = 1
init_status = 'TQ'

if init_status == 'SQ':
    spin_label_cell = ['I1']
    rho_info = 'I1p'

elif init_status == 'DQ':
    spin_label_cell = ['I1', 'I2']
    # It is relatively slow.
    rho_info = 'I1p*I2p'

elif init_status == 'TQ':
    spin_label_cell = ['I1', 'I2', 'I3']
    # It is slow.
    rho_info = 'I1p*I2p*I3p'


PO.create(spin_label_cell)
from PO import *

exec('rho = ' + rho_info)
print(rho)

gH_cell = [gH]*len(spin_label_cell)

if pfg_switch == 1:
    rho = rho.pfg(1, gH_cell)

# Imperfect 180 pulse (pi + d)
rho = rho.pulse(['I*'], ['x'], [pi + d])

if pfg_switch == 1:
    rho = rho.pfg(1, gH_cell)

rho.dispPO()
rho_dephase = rho.dephase()# Dephase has an issue
rho_dephase.dispPO()
print('Done!')