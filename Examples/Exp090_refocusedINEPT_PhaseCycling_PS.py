# Exp090_refocusedINEPT_PhaseCycling_PS.py
# refocused INEPT I => S
# Example to check phase cycling.
# Keeler, J., Understanding NMR Spectroscopy (1st Ed.), Wiley, 2005.
# pp. 174 - 175.

# Para begin #
spin_label_cell = ['I', 'S']
rho_ini = Iz*B + Sz
obs_cell = ['S']
no_ph = 7 # The number of pulse phases used in the pulse sequence
ph_cell[1] = [0,2] # ph1
ph_cell[2] = [0]   # ph2
ph_cell[3] = [0]   # ph3
ph_cell[4] = [0]   # ph4
ph_cell[5] = [1]   # ph5
ph_cell[6] = [0]   # ph6
ph_cell[7] = [0]   # ph7
phRtab = [0,2]
phid = list(range(2)) # phid: 0-based index
coef_cell = []
disp_bin = 1
# Para end #

# PS begin #
rho = rho.pulse(['I'], [ph1], [1/2*pi])
rho = rho.pulse(['I', 'S'], [ph3, ph2], [pi, pi])
rho = rho.jc(['IS'], [pi*J*2*t1])
rho = rho.pulse(['I', 'S'],[ph5, ph4],[1/2*pi, 1/2*pi])
rho = rho.pulse(['I', 'S'],[ph7, ph6],[pi, pi])     
rho = rho.jc(['IS'], [pi*J*2*t2])                       
# PS end #