# Exp050_OnePulse_PhaseCycling.py
# Example of writing a pulse sequence with phase cycling

# Para begin #
spin_label_cell = ['I']
rho_ini = Iz
obs_cell = ['I']
no_ph = 1 # The number of pulse phases used in the pulse sequence
ph_cell[1] = [1,2,3,0] # ph1
phRtab = [0,1,2,3]
phid = list(range(4)) # phid: 0-based index
coef_cell = []
disp_bin = 1
# Para end #

# PS begin #
rho = rho.pulse([1], [ph1], [1/2*pi]) # 90y-pulse
# PS end #