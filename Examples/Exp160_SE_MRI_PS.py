# Exp160_SE_MRI_PS.py
# Spin-Echo MRI with phase-encoding.

# Para begin #
spin_label_cell = ['I']
rho_ini = Iz
obs_cell = ['I']
no_ph = 2 # The number of pulse phases used in the pulse sequence
ph_cell[1] = [1, 3, 3, 1] # ph1
ph_cell[2] = [0, 0, 1, 1] # ph2
phRtab = [2,0,2,0]
phid = list(range(4)) # phid: 0-based index
coef_cell = []
disp_bin = 1
# Para end #

# PS begin #
rho = rho.pulse(['I'], [ph1], [1/2*pi])
rho = rho.pfg(t, [gH], 'fe') # t is a half of echo time
rho = rho.pfg(a*t, [gH], 'pe') # Phase Encoding during t with an intensity a
rho = rho.pulse(['I'], [ph2], [pi]) # 180 pulse
rho = rho.pfg(t, [gH], 'fe')  # t is a half of echo time
# PS end #