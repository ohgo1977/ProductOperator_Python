# Exp110_3QF_COSY_PS.py
# Guntert, P. et al., J. Magn. Reson. Ser. A, 101, 103-105, 1993.
# Guntert, P. Int. J. Quant. Chem., 106, 344-350, 2006.

# Para begin #
spin_label_cell = ['I1', 'I2', 'I3']
rho_ini = I1z
obs_cell = ['I*']
no_ph = 1 # The number of pulse phases used in the pulse sequence
ph_cell[1] = [0*pi/3,1*pi/3,2*pi/3,3*pi/3,4*pi/3,5*pi/3] # ph1
phRtab = [0,2]
phid = list(range(6)) # phid: 0-based index
# ph_cell[1] = [0*pi/3] # ph1
# phRtab = [0]
# phid = list(range(1)) # phid: 0-based index
coef_cell = []
disp_bin = 1
# Para end #

# PS begin #
PO.simp = 'fu'
rho = rho.pulse_phshift(['I*'], [ph1], [1/2*pi])
rho = rho.cs(['I*'], [o1*t])
rho = rho.jc(['I1I2', 'I1I3'], [pi*J12*t1, pi*J13*t1])
rho = rho.pulse_phshift(['I*'], [ph1], [1/2*pi])
rho = rho.pulse(['I*'], [0], [1/2*pi])
# PS end #