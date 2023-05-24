# Exp080_refocusedINEPT_InS.py
# Intensity calculation of refocused INEPT in InS system (n = 1,2 or 3)
# Levitt, M. H., Spin Dynamics (2nd Ed.), pp. 440 - 442, pp.488 - 491.

from PO import *
from sympy.plotting import plot

InS = 'IS'
if InS == 'IS':
    spin_label_cell = ['I1', 'S2']
    rho_info = 'I1z*B + S2z'
    jc_cell = ['I1S2']
    obs_cell = ['S2y']

elif InS == 'I2S':
    spin_label_cell = ['I1', 'I2', 'S3']
    rho_info = 'I1z*B + I2z*B + S3z'
    jc_cell = ['I1S3','I2S3']
    obs_cell = ['S3y']

elif InS == 'I3S':
    spin_label_cell = ['I1', 'I2', 'I3', 'S4']
    rho_info = 'I1z*B + I2z*B + I3z*B + S4z'
    jc_cell = ['I1S4','I2S4','I3S4']
    obs_cell = ['S4y']

PO.create(spin_label_cell)
from PO import *

q1 = 1/2*pi
q1_cell = [q1]*len(jc_cell)

q2 = pi*J*t
q2_cell = [q2]*len(jc_cell)

exec('rho = ' + rho_info)
print(rho)
rho = rho.pulse(['I*', 'S*'], ['x', 'x'], [3/2*pi, pi])
rho = rho.jc(jc_cell, q1_cell)

rho = rho.pulse(['I*', 'S*'], ['y', 'y'], [1/2*pi, 1/2*pi])
rho = rho.pulse(['I*', 'S*'], ['x', 'x'], [pi, pi])
rho = rho.jc(jc_cell, q2_cell)

rho_detect = rho.receiver('x')
rho_final = rho_detect.observable(['S*'])
rho_final.dispPO()
a0_V, rho_V = rho.SigAmp(['S*'], 'x')

# Plotting
J_Hz = 125
obs_id = rho_final.findterm(obs_cell) # obs_id: list of 0-based index
coef = rho_final.coef[obs_id[0]]
coef = coef.subs(J, J_Hz).subs(B, 1)
p1 = plot(coef, (t, 0, 1/J_Hz))
print('Done!')
