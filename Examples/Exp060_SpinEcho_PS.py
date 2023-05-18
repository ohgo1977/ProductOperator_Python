# Exp060_SpinEcho_PS.py
# Spin-echo (Hahn-echo) experiment with phase cycling.
# Effect of the miscalibration of 180 pulse can be checked.

# Para begin #
spin_label_cell = ['I']
rho_ini = Iz
obs_cell = ['I']
no_ph = 2 # The number of pulse phases used in the pulse sequence
ph_cell[1] = [1,2,3,0] # ph1
ph_cell[2] = [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3]  # ph2
phRtab = [0,3,2,1,2,1,0,3]
phid = list(range(16)) # phid: 0-based index
coef_cell = []
disp_bin = 1
# Para end #

# PS begin #
rho = rho.pulse(['I'], [ph1], [1/2*pi]) # 90-pulse
rho = rho.cs(['I'], [oI*t]) # Chemical shift evolution
# rho = rho.pulse(['I'],[ph2],[pi]) # 180-pulse,
rho = rho.pulse(['I'],[ph2],[pi + d]) # 180+d-pulse, where d indicates the miscalibration of 180-pulse
rho = rho.cs(['I'], [oI*t]) # Chemical shift evolution
# PS end #