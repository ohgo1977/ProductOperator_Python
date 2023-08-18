# Exp130_PFG_Phase.py
# Keeler, J., Understanding NMR Spectroscopy, p. 399, 11.11.1
# Calculation of spatially dependent phase created by PFG for coherence order P.

from PO import *

# By implementing TR8() for simplification, the problem described here was solved.
# print('\nNote that there are cases that the coeffcient is not simplified nicely.')
# print('In the case of DQ, for example,')
# print('I1pI2m:  I*exp(I*GZ*gH)*sin(GZ*gH)/2 - exp(I*GZ*gH)*cos(GZ*gH)/2 + 1/2 that should be 0')
# print('I1mI2p: -I*exp(I*GZ*gH)*sin(GZ*gH)/2 + exp(I*GZ*gH)*cos(GZ*gH)/2 + 1/2 that should be 1')
# print('This is due to a weak capability of simplify() used in PO.CombPO()')
# print('In such case, the combination with rewrite(cos) will improve the result.')
# print('However, it is time consuming if rewrite(cos) is used in PO.CombPO().')
# print('Instead, the result can be simplified after finishing the calculation.')
# print('In this example, simplify_cos() can be applied by setting simplify_switch as 1\n')


init_status = 'DQ'
gH = symbols('gH')
simplify_switch = 1
print('\nsimplify_switch:', simplify_switch, '\n')

PO.simp = 'TR8'

if init_status == 'SQ':
    spin_label_cell = ['I1']

elif init_status == 'DQ': # It is slow.
    spin_label_cell = ['I1', 'I2']

elif init_status == 'TQ': # It is slow.
    spin_label_cell = ['I1', 'I2', 'I3']


ns = len(spin_label_cell)
gH_cell = [gH]*ns

for ii in range(2**ns):
    for jj in range(2**ns):
        M_in = np.matrix(np.zeros((2**ns, 2**ns)))
        if ii != jj:
            M_in[ii,jj] = 1
            rho = PO.M2pol(M_in, spin_label_cell)
            rho.disp = 0
            rho = rho.pfg(1, gH_cell)
            if simplify_switch == 1:
                rho = rho.simplify_exp()
            rho.dispPO()

print('Done!')