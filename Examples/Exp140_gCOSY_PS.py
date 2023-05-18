# Exp140_gCOSY_PS.py
# Berger, S.; Braun, S. 200 and More NMR Experiments A Practical Course, p. 526

# Para begin #
spin_label_cell = ['I1', 'I2']
rho_ini = I1z
obs_cell = ['I*']
no_ph = 2 # The number of pulse phases used in the pulse sequence
ph_cell[1] = [0, 2] # ph1
ph_cell[2] = [0, 0, 2, 2] # ph1
phRtab = [0,2]
phid = list(range(4)) # phid: 0-based index
coef_cell = []
disp_bin = 1
# Para end #

# PS begin #
rho = rho.pulse(['I*'], [ph1], [1/2*pi])
rho = rho.cs(['I1','I2'], [o1*t1, o2*t1])
rho = rho.jc(['I1I2'], [pi*J12*t1])
rho = rho.pfg(1, [gH, gH])
rho = rho.pulse(['I*'], [ph2], [1/2*pi])
rho = rho.pfg(1, [gH, gH])
rho = rho.dephase()
# PS end #