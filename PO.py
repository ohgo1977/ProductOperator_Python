#  ------------------------------------------------------------------------
#  Class       : PO
#  Description : Functions for Product Operator Formalism of spin-1/2
#  Tested      : Python 3.8.5, SymPy 1.11.2, NumPy 1.23.3
#  Developer   : Dr. Kosuke Ohgo
#  ULR         : https://github.com/ohgo1977/PO_Python
#  Version     : 1.0.0
# 
#  Please read the manual (PO_Python_Manual.pdf) for details.
# 
#  ------------------------------------------------------------------------
# 
# MIT License
#
# Copyright (c) 2023 Kosuke Ohgo
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

print("Hello from PO.py!\n")
from sympy import exp, cos, sin, pi, symbols, I
from sympy import init_printing, Wild
import sympy as sym
import copy
from math import log2
import numpy as np
import sys
# from ttictoc import tic,toc # For time measurements using tic() and print(toc()).

# init_printing(use_unicode = False, wrap_line = False)
init_printing()

class PO:
    spin_label_cell_default = ['I', 'S', 'K', 'L', 'M']
    # Default labels for spin_label

    __index_switch = 1
    if __index_switch != 1 and __index_switch != 0:
         __index_switch = 1

    if __index_switch == 0:
        print('! CAUTION !')
        print('__index_switch is set as ', __index_switch)
        print('The input index for spin types becomes 0-based') 
        print('Example:')
        print('I-S spin system:  0 for I and 1 for S with 0-based index')
        print('sp2id(), jc1(), symcoef(), findterm(), observable(), dispPO() will include actual codes to convert from 1-based to 0-based.')
        # In pulse() and cs(), sp2id() converts 1-based index to 0-based index if PO.__index_switch is 1.
        # In selPO() and delPO(), findterm() converts 1-based index to 0-based index if PO.__index_switch is 1.

    asterisk_bin = 0
    if asterisk_bin != 1 and asterisk_bin != 0:
        asterisk_bin = 0
        print(asterisk_bin)
    # Control to the asterisk '*' between spin operators in the txt property.
    # 0 :w/o, 1: w/ '*'.

    def __init__(self, spin_no, spin_label, axis, coef, basis):
        self.axis = axis # np.matrix
        self.coef = sym.Matrix([coef]) # sym.Matrix
        self.spin_no = spin_no
        self.spin_label = spin_label # List
        self.basis = basis
        self.disp = 1
        self.PFGq = sym.Matrix([0]) # Initial value, is there a better initial value?
        self.Ncoef = PO.getNcoef(self) # np.array
        self.txt = PO.getTxt(self)
        self.logs = self.txt
        self.M = PO.getM(self) # sym.Matrix
        self.coherence = PO.getCoherence(self) # sym.Matrix

    def __str__(self):

        return self.txt

    def getTxt(self):
        # Generating a text describing the current spin operators.
        # If PO.asterisk_bin is 1, '*' will be added among spin operators and coeffcients, i.e.,
        # 2*Ix*Sx*cos(q).

        spin_label = self.spin_label
        axis = self.axis
        Ncoef = self.Ncoef
        coef = self.coef

        if coef == sym.Matrix([0]):
            txt = '0'

        else:
            for ii in range(axis.shape[0]):
                axis_vec = np.array(axis)[ii]
                pt = PO.axis2pt(spin_label, axis_vec)

                if Ncoef[ii] == 0.5:
                    Ncoef_t = '1/2'
                    if PO.asterisk_bin == 1:
                        Ncoef_t = Ncoef_t + '*' 
                elif Ncoef[ii] == 1.0:
                    Ncoef_t = ''
                else:
                    Ncoef_t = str(Ncoef[ii].astype(np.int64))
                    if PO.asterisk_bin == 1:
                        Ncoef_t = Ncoef_t + '*'

                # coef_t = '*' + '(' + str(sym.simplify(coef[ii])) + ')' # Makes program slow
                coef_t = '*' + '(' + str(coef[ii]) + ')'

                if coef_t == '*(1)':
                    coef_t = ''

                if ii == 0:
                    if coef_t != '*(0)':
                        txt = Ncoef_t + pt + coef_t
                    else:
                        txt = ''
                else:
                    if coef_t != '*(0)':
                        txt = txt + ' + ' + Ncoef_t + pt + coef_t
        return txt
    
    def axis2pt(spin_label, axis_vec):
        # Generating a text of a spin operator corresponding to axis_vec.
        # axis_vec should be a np.array.
        # Example. 
        # spin_label = ['I', 'S'] and axis_vec = [3, 3] => pt = 'IzSz'
        #
        # If PO.asterisk_bin is 1, '*' will be added between spin operators, i.e.,
        # pt = 'Iz*Sz'

        pt = ''
        if np.count_nonzero(axis_vec, axis=0) == 0:
            pt = 'E'
        else:
            jj_int = 0
            for jj in range(len(np.array(axis_vec))):
                axis_v = axis_vec[jj]
                st = spin_label[jj]

                if axis_v != 0:
                    jj_int += 1
                    if axis_v == 1:
                        at = 'x'
                    elif axis_v == 2:
                        at = 'y'
                    elif axis_v == 3:
                        at = 'z'
                    elif axis_v == 4:
                        at = 'p'
                    elif axis_v == 5:
                        at = 'm'
                    elif axis_v == 6:
                        at = 'a'
                    elif axis_v == 7:
                        at = 'b'
                    # end

                    if PO.asterisk_bin == 0:
                        pt = pt + st + at

                    elif PO.asterisk_bin == 1: # Iz*Sz style
                        if jj_int == 1:
                            pt = pt + st + at
                        else:
                            pt = pt + '*' + st + at 
                                  
                # end
            # end
        return pt
    
    def getNcoef(self):
        # Calculating coeffcients corresponding to 2^(N-1) where N is the number of the active spins.
        axis = self.axis
        basis = self.basis

        if basis == 'xyz':
            Ncoef = np.float_power(2, (np.count_nonzero(np.array(axis), axis=1)-1))
            # 2^(N-1) 
            # Examples 
            # IxSx   => N = 2 => Ncoef = 2
            # IxSxKx => N = 3 => Ncoef = 4
        elif basis == 'pmz' or basis == 'pol':
            Ncoef = np.array([1]*axis.shape[0])

        # Ncoef = np.transpose(np.matrix(Ncoef)) # This line will affet getM().
        return Ncoef

    def getM(self):
        # Calculating a matrix form of the current density operator.
        # The output matrix M_out is a sym class.

        axis = self.axis
        Ncoef = self.Ncoef
        coef = self.coef

        for ii in range(axis.shape[0]):
            axis_tmp = axis[ii]
            for jj in range(axis_tmp.shape[1]):
                v_tmp = axis_tmp[0,jj]
                Ma = PO.axis2M(v_tmp)

                if jj == 0:
                    M_tmp = Ma
                else:
                    M_tmp = np.kron(M_tmp, Ma)


            M_tmp = sym.Matrix(M_tmp) # Convert to Symbolic.
            Mo = M_tmp * 2 ** (axis_tmp.shape[1] - np.count_nonzero(axis_tmp)) * Ncoef[ii] * coef[ii]

            if ii == 0:
                M_out = Mo
            else:
                M_out = M_out + Mo
            
        # M_out = sym.simplify(M_out) # It makes the program slow.
        # M_out = sym.nsimplify(M_out, rational=True) # It makes the program slow.

        return M_out

    def getCoherence(self):
        # Calculating a matrix swhoing coherences included in the current density operator.
        # Note that 'u' (up) and 'd' (down) are used instead of 'p' (+) and 'm' (-) in this display.
        # Details are in PO.rho_box().

        spin_no = self.spin_no
        coherence_out = PO.rho_box(spin_no)[0]
        for ii in range(2**spin_no):
            for jj in range(2**spin_no):
                if self.M[ii,jj] == 0:
                    coherence_out[ii,jj] = 0
        return coherence_out

    def CombPO(self):

        # Combining coeffcieints of same type of terms in a PO-class object.
        axis_in = self.axis
        coef_in = self.coef
        # https://stackoverflow.com/questions/36068437/is-there-python-equivalent-to-matlab-unique
        c, IA, IC = np.unique(axis_in,return_index=True,return_inverse=True, axis=0)

        coef_out = sym.Matrix(np.zeros((len(IA), 1)))
        axis_out = np.matrix(np.zeros((len(IA), axis_in.shape[1]))).astype(int)

        for ii in range(len(IA)):
            IA_tmp = IA[ii]
            IC_tmp = np.where(IC == IC[IA_tmp])[0]# There are two output
            axis_out[ii,:] = axis_in[IA_tmp,:]
            coef_out[ii,0] = sum(coef_in[tuple(IC_tmp),0])# Need to convert to tuple

        coef_out  = sym.simplify(sym.nsimplify(sym.expand(coef_out), rational=True)) # This line will be good for dephase().
        len_coef_out = len(coef_out)

        ii_int = 0
        # Remove terms with coef = 0
        for ii in range(len_coef_out - 1, -1, -1):
            if isinstance(coef_out[ii], sym.core.numbers.Zero):
                axis_out = np.delete(axis_out, ii, axis=0)
                coef_out = coef_out.row_del(ii)
                ii_int += 1

        # All zero condition.
        if ii_int == len_coef_out:
            axis_out = np.matrix([0]*self.spin_no)
            coef_out = sym.Matrix([0])

        return PO(self.spin_no, self.spin_label, axis_out, coef_out, self.basis)

    def pulse(self, sp_cell, ph_cell, q_cell):
        # Calculation of the change of rho under simultaneous pulses.
        # sp_cell: type of spins in a list.
        # Character of spin ('I' or 'S' etc.) or
        # the order number of spin (0 for 'I', 1 for 'S' etc.). # 0-based index
        #                          (1 for 'I', 2 for 'S' etc.). # 1-based index
        # ph_cell: phases of the pulses in a list. 
        # Character ('x','X','-x' etc.) or number (0,1,2,3)
        # quadrature phases are only accepted in this function.
        # q_cell: flip angles (rad) in a list (symbolic or double)
        # Exmaple:
        # rho_new = pulse(rho,['I', 'S'],['x', 'y'],[pi/2, pi/2])
        spin_label_cell = self.spin_label
        id_vec, ii_vec = PO.sp2id(sp_cell, spin_label_cell)
        for ii in range(len(id_vec)):
            sp = id_vec[ii]
            ph = ph_cell[ii_vec[ii]]
            q = q_cell[ii_vec[ii]]
            self = PO.pulse1(self, sp, ph, q)[0]

        return self

    def pulse1(self, sp, ph, q):
        #  Calculation of the change of rho under a pulse.
        #  sp: type of spin, character of spin ('I' or 'S' etc.) or
        #  the order number of spin (0 for 'I', 1 for 'S' etc.). # 0-based index
        #  ph: phase of the pulse, character ('x','X','-x' etc.) or number (0,1,2,3)
        #      quadrature phase only.
        #  q: flip angle in radian (symbolic or double)
        basis_org = self.basis
        spin_no = self.spin_no
        spin_label_cell = self.spin_label
        s0 = self.logs
        disp0 = self.disp
        self = PO.set_basis(self,'xyz')
        axis_tmp = np.matrix([0]*spin_no)

        if isinstance(sp, str):
            for ii in range(len(spin_label_cell)):
                if sp == spin_label_cell[ii]:
                    id_sp = ii # 0-based index
        elif isinstance(sp, int):
            id_sp = sp # 0-based index

        if ph == 'x' or ph == 'X' or (isinstance(ph, int) and ph == 0):
            phase_id = 1
            coef_tmp = sym.Matrix([1])
            ph_tmp = 'x'
        elif ph == 'y' or ph == 'Y' or (isinstance(ph, int) and ph == 1):
            phase_id = 2
            coef_tmp = sym.Matrix([1])
            ph_tmp = 'y'
        elif ph == '-x' or ph == '-X' or (isinstance(ph, int) and ph == 2):
            phase_id = 1
            coef_tmp = sym.Matrix([-1])
            ph_tmp = '-x'
        elif ph == '-y' or ph == '-Y' or (isinstance(ph, int) and ph == 3):
            phase_id = 2
            coef_tmp = sym.Matrix([-1])
            ph_tmp = '-y'


        axis_tmp[0, id_sp] = phase_id

        H = PO(self.spin_no, self.spin_label, axis_tmp, coef_tmp, 'xyz')
        obj = PO.UrhoUinv_mt(self, H, q)
        obj = PO.set_basis(obj, basis_org)
        obj.disp = disp0

        q_s = str(q)
        flag_tmp = q_s.find('pi')
        if flag_tmp == -1: # Not including 'pi' in q_s
            s_out = 'Pulse: %s %s %s' % (spin_label_cell[id_sp],str(q),ph_tmp)
        else:
            try:
                s_out = 'Pulse: %s%4d%s' % (spin_label_cell[id_sp],round(q/pi*180),ph_tmp)
            except:
                s_out = 'Pulse: %s %s %s' % (spin_label_cell[id_sp],str(q),ph_tmp)

        s1 = '%s' % (s_out)
        s2 = '    %s' % (obj.txt)
        obj.logs = '%s\n%s\n%s' %(s0, s1, s2)

        if obj.disp == 1:
            print('%s' % (s1))
            print('%s' % (s2))

        return obj, id_sp


    def pulse_phshift(self, sp_cell, ph_cell, q_cell):
        # Calculation of the change of rho under simultaneous pulses with arbitrary phases.
        # ph_cell: phases of the pulses in a list. 
        # Arbitrary phases in radian are allowed.
        # rho_new = pulse(rho, ['I', 'S'], [pi/4, pi/4], [pi/2, pi/2])
        # rho_new = pulse(rho, ['I', 'S'], [q, q], [pi/2, pi/2])

        spin_label_cell = self.spin_label
        id_vec, ii_vec = PO.sp2id(sp_cell, spin_label_cell)
        for ii in range(len(id_vec)):
            sp = id_vec[ii]
            ph = ph_cell[ii_vec[ii]]
            q = q_cell[ii_vec[ii]]
            self = PO.pulse_phshift1(self, sp, ph, q)[0]

        return self

    def pulse_phshift1(self, sp, ph, q):
        #  Calculation of the change of rho under a pulse with an arbitrary phase.
        #  sp: type of spin, character of spin ('I' or 'S' etc.) or
        # the order number of spin (0 for 'I', 1 for 'S' etc.). # 0-based index
        #  ph: Arbitrary phase in radian.
        #  q: flip angle in radian (symbolic or double)

        s0 = self.logs
        disp0 = self.disp
        spin_label_cell = self.spin_label

        self.disp = 0        
        self = PO.cs1(self, sp, -ph)[0]
        self = PO.pulse1(self, sp, 'x', q)[0]
        self, id_sp = PO.cs1(self, sp, ph) # id_sp: 0-based
        self.disp = disp0

        q_s = str(q)
        flag_tmp = q_s.find('pi')
        if flag_tmp == -1: # Not including 'pi' in q_s
            s_out = 'Pulse: %s %s %s' % (spin_label_cell[id_sp], str(q), str(ph))
        else:
            try:
                s_out = 'Pulse: %s%4d %s' % (spin_label_cell[id_sp], round(q/pi*180), str(ph))
            except:
                s_out = 'Pulse: %s %s %s' % (spin_label_cell[id_sp], str(q), str(ph))

        s1 = '%s' % (s_out)
        s2 = '    %s' % (self.txt)
        self.logs = '%s\n%s\n%s' %(s0, s1, s2)

        if self.disp == 1:
            print('%s' % (s1))
            print('%s' % (s2))

        return self, id_sp

    def cs(self, sp_cell, q_cell):

        spin_label_cell = self.spin_label
        id_vec, ii_vec = PO.sp2id(sp_cell, spin_label_cell)
        for ii in range(len(id_vec)):
            sp = id_vec[ii]# int
            q = q_cell[ii_vec[ii]]
            self = PO.cs1(self, sp, q)[0] #cs1 outputs two arugements as a tuple

        return self


    def cs1(self, sp, q):
        # obj = cs1(obj,sp,q)
        # Calculation of the chemical shift evolution of rho.
        # obj: PO class object
        # sp: type of spin, character of spin ('I' or 'S' etc.) or
        # the order number of spin (0 for 'I', 1 for 'S' etc.).
        # q: flip angle (symbolic or double)

        basis_org = self.basis
        spin_no = self.spin_no
        s0 = self.logs
        disp0 = self.disp
        spin_label_cell = self.spin_label
        self = PO.set_basis(self,'xyz')
        axis_tmp = np.matrix([0]*spin_no)

        if isinstance(sp, str):
            for ii in range(len(spin_label_cell)):
                if sp == spin_label_cell[ii]:
                    id_sp = ii
        elif isinstance(sp, int):
            id_sp = sp # 0-based index

        coef_tmp = sym.Matrix([1])
        axis_tmp[0, id_sp] = 3

        H = PO(self.spin_no, self.spin_label, axis_tmp, coef_tmp, 'xyz')
        obj = PO.UrhoUinv_mt(self, H, q)
        obj = PO.set_basis(obj, basis_org)
        obj.disp = disp0

        q_s = str(q)
        flag_tmp = q_s.find('pi')
        if flag_tmp == -1: # Not including 'pi' in q_s
            s_out = 'CS: %s %s' % (spin_label_cell[id_sp],str(q))
        else:
            try:
                s_out = 'CS: %s%4d' % (spin_label_cell[id_sp], round(q/pi*180))
            except:
                s_out = 'CS: %s %s' % (spin_label_cell[id_sp], str(q))

        s1 = '%s' % (s_out)
        s2 = '    %s' % (obj.txt)
        obj.logs = '%s\n%s\n%s' %(s0, s1, s2)

        if obj.disp == 1:
            print('%s' % (s1))
            print('%s' % (s2))

        return obj, id_sp

    def jc(self, sp_cell, q_cell):

        for ii in range(len(sp_cell)):
            sp = sp_cell[ii]
            q = q_cell[ii]
            self = PO.jc1(self, sp, q)[0] #jc1 outputs two arugements as a tuple

        return self

    def jc1(self, sp, q):
        basis_org = self.basis
        spin_no = self.spin_no
        spin_label_cell = self.spin_label
        s0 = self.logs
        disp0 = self.disp
        self = PO.set_basis(self,'xyz')
        axis_tmp = np.matrix([0]*spin_no)


        if isinstance(sp, list): # This condiion means that sp is a list of Int, i.e., sp = [1,2]
            sp_tmp = ''
            for ii in range(2):
                if PO.__index_switch == 1: # 1-based index
                    id_tmp = sp[ii] - 1 # 1-based index => 0-based index
                else:
                    id_tmp = sp[ii] # 0-based index

                axis_tmp[0, id_tmp] = 3
                sp_tmp = sp_tmp + spin_label_cell[id_tmp]

            id_sp = sp

        elif isinstance(sp, str):
            id_sp = [0]*2
            ii_int = -1
            for ii in range(len(spin_label_cell)):
                spin_label_tmp = spin_label_cell[ii]
                if sp.find(spin_label_tmp) >= 0:
                    id_tmp = ii
                    axis_tmp[0, id_tmp] = 3
                    ii_int += 1
                    id_sp[ii_int] = id_tmp
            sp_tmp = sp

        coef_tmp = sym.Matrix([1])

        H = PO(self.spin_no, self.spin_label, axis_tmp, coef_tmp, 'xyz')
        obj = PO.UrhoUinv_mt(self, H, q)
        obj = PO.set_basis(obj, basis_org)
        obj.disp = disp0

        q_s = str(q)
        flag_tmp = q_s.find('pi')
        if flag_tmp == -1: # Not including 'pi' in q_s
            s_out = 'JC: %s %s' % (sp_tmp,str(q))
        else:
            try:
                s_out = 'JC: %s%4d' % (sp_tmp,round(q/pi*180))
            except:
                s_out = 'JC: %s %s' % (sp_tmp,str(q))

        s1 = '%s' % (s_out)
        s2 = '    %s' % (obj.txt)
        obj.logs = '%s\n%s\n%s' %(s0, s1, s2)

        if obj.disp == 1:
            print('%s' % (s1))
            print('%s' % (s2))

        return obj, id_sp

    def sp2id(sp_cell, spin_label_cell):
        
        for ii in range(len(sp_cell)):
            sp = sp_cell[ii]
            if isinstance(sp, int):
                if PO.__index_switch == 1:
                    id_vec_tmp = np.array([sp - 1])# 1-based index => 0-based index
                else:
                    id_vec_tmp = np.array([sp]) # 0-based index  

            elif isinstance(sp, str):
                if sp.find('*') == 0 or sp.find('*') == 1:# Wildcard 'I*' or '*'
                    if len(sp) == 1: # '*'
                        id_vec_tmp = np.arange(0,len(spin_label_cell))
                    elif len(sp) == 2: # 'I*'
                        # https://stackoverflow.com/questions/14849293/find-all-index-position-in-list-based-on-partial-string-inside-item-in-list
                        id_vec_tmp = [i for i, s in enumerate(spin_label_cell) if sp[0] in s]
                else:
                    id_vec_tmp = [i for i, s in enumerate(spin_label_cell) if sp in s]

            ii_vec_tmp = [ii]*len(id_vec_tmp)

            if ii == 0:
                id_vec = id_vec_tmp
                ii_vec = ii_vec_tmp
            else:
                id_vec = np.concatenate((id_vec, id_vec_tmp),0)
                ii_vec = np.concatenate((ii_vec, ii_vec_tmp),0)
            
        id_vec, id_tmp, tmp = np.unique(id_vec, return_index=True, return_inverse=True, axis=0)
        id_vec = id_vec.astype(int).tolist() # Convert elements to int => list, 0-based index

        ii_vec = np.array(ii_vec)
        ii_vec = ii_vec[id_tmp]
        return id_vec, ii_vec

    def pfg(self, G_coef, gamma_cell):
        # obj = pfg(obj,G_coef,gamma_cell)
        # applys pulse field gradient to all spins.
        # G_coef is a relative strengh of the gradient field and 
        # gamma_cell a cell array including gyromagnetic ratio of the spins.
        # Symbolic constant GZ is used as a stamp to show terms affected by pfg().
        # This information is used in dephase().
        # This method was obitaned from POMA by Gunter (2006).

        GZ = symbols('GZ')
        obj_tmp = self
        spin_label_cell = obj_tmp.spin_label
        disp0 = obj_tmp.disp

        q_tmp = sym.Matrix([0])
        for ii in range(self.spin_no):
            s0 = obj_tmp.logs
            sp_tmp = spin_label_cell[ii]
            q = G_coef*GZ*gamma_cell[ii]
            q0 = obj_tmp.PFGq # Initial PFGq
            obj_tmp.disp = 0
            obj_tmp = PO.cs1(obj_tmp, sp_tmp, q)[0] # PFGq of the output obj_tmp is reset to sym.Matrix([0]).
            obj_tmp.disp = disp0

            q_s = str(q)
            flag_tmp = q_s.find('pi')
            if flag_tmp == -1: # Not including 'pi' in q_s
                s_out = 'PFG: %s %s' % (spin_label_cell[ii],str(q))
            else:
                s_out = 'PFG: %s%4d' % (spin_label_cell[ii],round(q/pi*180))

            s1 = '%s' % (s_out)
            s2 = '    %s' % (obj_tmp.txt)
            obj_tmp.logs = '%s\n%s\n%s' %(s0, s1, s2)

            if obj_tmp.disp == 1:
                print('%s' % (s1))
                print('%s' % (s2))

            if q0 == sym.Matrix([0]): # Initial condition
                # obj_tmp.PFGq[0] = q
                obj_tmp.PFGq[0] = GZ*gamma_cell[ii] # Ignore G_coef
            else:
                # q_tmp[0] = q
                # obj_tmp.PFGq = q0.row_join(q_tmp)
                q_tmp[0] = GZ*gamma_cell[ii] # Ignore G_coef
                obj_tmp.PFGq = q0.row_join(q_tmp)

        return obj_tmp

    def dephase(self):
        # obj = dephase(obj)
        # delete terms affected by pfg().
        # This method was obitaned from POMA by Gunter (2006).

        obj = copy.deepcopy(self)# Give a differnt ID to obj from self
        PFGq_in = obj.PFGq
        coef_in = obj.coef
        # coef_in  = sym.nsimplify(sym.expand(coef_in), rational=True) # This line will be good for dephase().
        s0 = obj.logs
        Zpfg = symbols('Zpfg')

        if PFGq_in != sym.Matrix([0]):
            for ii in range(len(PFGq_in)):
                PFGq_tmp = PFGq_in[ii]
                for jj in range(len(coef_in)):
                    # Replace PFGq_tmp to Zpfg
                    coef_in[jj] = coef_in[jj].simplify()
                    coef_in[jj] = coef_in[jj].subs(PFGq_tmp, Zpfg)

        # https://docs.sympy.org/latest/modules/core.html
        # 2.1. Pattern -> expr
        a, b = map(Wild, 'ab')
        for ii in range(len(coef_in)):
            coef_in[ii] = coef_in[ii].rewrite(exp).expand() # This line is important to choose exp(I*a*Zpfg) correctly
            coef_in[ii] = coef_in[ii].replace(exp(a*Zpfg),0) # exp(I*a*Zpfg) => 0
            coef_in[ii] = coef_in[ii].replace(exp(a*Zpfg + b),0) # exp(I*a*Zpfg + b) => 0
            # coef_in[ii] = sym.simplify(coef_in[ii]) # Conversion from exp to cos and sin
            coef_in[ii] = coef_in[ii].rewrite(cos).simplify() # Conversion from exp to cos and sin

        obj.coef = coef_in
        obj_out = PO.CombPO(obj) #PFGq is reset to 0
        # obj_out.PFGq = PFGq_in # No need to keep the record of PFGq_in

        s_out = 'Dephasing of terms including %s' % (str(PFGq_in)[8:-2]) # Matrix([[GZ*gH, GZ*gC]]) => [GZ*gH, GZ*gC]
        s1 = '%s' % (s_out)
        s2 = '    %s' % (obj_out.txt)
        obj_out.logs = '%s\n%s\n%s' %(s0, s1, s2)

        if obj_out.disp == 1:
            print('%s' % (s1))
            print('%s' % (s2))

        return obj_out

    def UrhoUinv_mt(self, H, q):

        if not(isinstance(self, PO)) or not(isinstance(H, PO)):
            sys.exit('Both obj and H should be the PO object!')

        if H.basis != 'xyz':
            sys.exit('The basis of H should be xyz!')

        if H.axis.shape[0] > 1:
            sys.exit('H must be a single term!')

        basis_org = self.basis
        self = PO.set_basis(self,'xyz')

        mt = np.matrix([[0, 3, -2],[-3, 0, 1],[2, -1, 0]])
        mt_large = np.matrix([[0, 3, -2, 1],[-3, 0, 1, 2],[2, -1, 0, 3],[1, 2, 3, 0]])
        q = q*H.coef # Matrix
        q = q[0,0]# Matrix to a single symbolic value

        H_axis = np.array(H.axis)

        type_mask_mat = (np.multiply(self.axis, H_axis) != 0)*1 # [True, False] => [0, 1] 
        axis_diff_mat = (self.axis != H_axis)*1 # [True, False] => [0, 1] 
        axis_mask_mat = np.multiply(type_mask_mat, axis_diff_mat)
        axis_mask_vec = axis_mask_mat.sum(axis=1)

        # Python Assigment
        # a = [1, 2, 3, 4]
        # b = a
        # b[0] = 'X'
        # Then not only b but also a become
        # ['X', 2, 3, 4]
        axis_mask_vec2 = axis_mask_vec
        axis_mask_vec2[axis_mask_vec2 != 1] = 0
        axis_new_tmp = np.zeros((self.axis.shape[0] + axis_mask_vec2.sum(axis=0)[0,0], self.axis.shape[1]))
        coef_new_tmp = sym.Matrix(np.zeros((self.axis.shape[0] + axis_mask_vec2.sum(axis=0)[0,0], 1)))

        axis_mask_vec3 = axis_mask_vec2 # n x 1 matrix
        id_tmp_vec = np.where(axis_mask_vec3 == 1)[0] # Info of row, array

        # Python Index
        # MATLAB      => Python
        # a(1:3,5:9)  => a[0:3, 4:9]
        # a(1:5,:)    => a[0:5] or a[:5] or a[0:5, :]
        # a(3:2:21,:) => a[2:21:2,:]
        for ii in range(len(id_tmp_vec) - 1, -1, -1):
            id_tmp = id_tmp_vec[ii]

            if id_tmp != len(axis_mask_vec3) - 1:
                axis_mask_vec3 = np.concatenate((np.concatenate((axis_mask_vec3[0:id_tmp + 1], np.matrix([2])),axis=0), axis_mask_vec3[id_tmp + 1:len(axis_mask_vec3)]),axis=0)

            elif id_tmp == len(axis_mask_vec3) - 1:
                axis_mask_vec3 = np.concatenate((axis_mask_vec3[0:id_tmp + 1], np.matrix([2])),axis=0)

        # Non-sin terms
        axis_new_tmp[np.where(axis_mask_vec3 != 2)[0],:] = self.axis

        arr = np.where(axis_mask_vec3 != 2)[0]
        id_tmp = arr.tolist()
        ii_int = 0
        for ii in id_tmp:
            coef_new_tmp[ii,0] = self.coef[ii_int]
            # ii_int = ii_int + 1
            ii_int += 1

        # sin terms
        arr1 = np.where(axis_mask_vec3 == 1)[0]
        id_tmp1 = arr1.tolist()
        arr2 = np.where(axis_mask_vec3 == 2)[0]
        id_tmp2 = arr2.tolist()
        # print(coef_new_tmp[id_tmp, 0]) # This works.
        # coef_new_tmp[id_tmp2,0] = coef_new_tmp[id_tmp1,0] # This doesn't work
        for ii in range(len(id_tmp1)):
            coef_new_tmp[id_tmp2[ii],0] = coef_new_tmp[id_tmp1[ii],0]

        # cos terms
        axis_cos = axis_new_tmp[np.where(axis_mask_vec3 == 1)[0],:]
        for ii in range(len(id_tmp1)):
            coef_new_tmp[id_tmp1[ii],0] = cos(q)*coef_new_tmp[id_tmp1[ii],0]

        H_axis_mat = np.repeat(H_axis, axis_cos.shape[0], axis=0)

        axis_cos4 = axis_cos
        H_axis_mat4 = H_axis_mat
        axis_cos4[axis_cos4 == 0] = 4
        H_axis_mat4[H_axis_mat4 == 0] = 4

        axis_sin = np.zeros((axis_cos4.shape[0], axis_cos4.shape[1]))
        for ii in range(axis_sin.shape[0]):
            id1 = H_axis_mat4[ii,:]
            id1 = id1.astype(np.int64)
            id2 = axis_cos4[ii,:]
            id2 = id2.astype(np.int64)

            axis_sin_tmp = mt_large[np.ix_(id1 - 1, id2 - 1)]
            axis_sin[ii,:] = np.diag(axis_sin_tmp)

        axis_new_tmp[np.where(axis_mask_vec3 == 2)[0],:] = axis_sin
        axis_new_tmp = abs(axis_new_tmp)

        # It looks that it is not necessary to add ~isempty(axis_sin) necessary.
        axis_sin_sign = axis_sin
        axis_sin_sign[axis_sin_sign == 0] = 1
        axis_sin_sign = np.prod(np.sign(axis_sin_sign),axis=1)

        arr2 = np.where(axis_mask_vec3 == 2)[0]
        id_tmp2 = arr2.tolist()
        for ii in range(len(id_tmp2)):
            coef_new_tmp[id_tmp2[ii],0] = sin(q)*coef_new_tmp[id_tmp2[ii],0]*axis_sin_sign[ii]

        axis_new_tmp = np.matrix(axis_new_tmp) # if axis_new_tmp is array, it causes an error in getM
        obj_out = PO(self.spin_no, self.spin_label, axis_new_tmp, coef_new_tmp, 'xyz') # basis should be 'xyz'
        obj_out = PO.CombPO(obj_out)
        obj_out = PO.set_basis(obj_out, basis_org)
        return obj_out

    def xyz2pmz(self):
        # conversion from Cartesian operator basis to lowring/raising operator basis
        if self.basis != 'xyz':
            sys.exit('The basis of the object should be xyz')

        axis_in = self.axis
        coef_in = self.coef
        Ncoef_in = self.Ncoef
        spin_no = self.spin_no
        disp0 = self.disp

        for ii in range(axis_in.shape[0]):
            axis_tmp = axis_in[ii,:]# np.matrix

            if len(np.where(axis_tmp == 1)[0]) == 0 and len(np.where(axis_tmp == 2)[0]) == 0:
                # Only *z operators Iz, 2IzSz, 4IzSzKz,...
                axis_out_tmp = axis_tmp
                coef_out_tmp = sym.Matrix([1])
            else:
                # Conversion from Ix and Iy to Ip and Im
                xn = len(np.where(axis_tmp == 1)[0]) # Example: IxSyKxMz =>(Ip + Im)(Sp - Sm)(Kp + Km)Mz
                yn = len(np.where(axis_tmp == 2)[0]) # xn =2, yn = 1 => 2^(2+1) = 8 terms.
                xyn = 2 ** (xn + yn)
                axis_out_tmp = np.tile(axis_tmp, (xyn, 1)) # repmat([1 2 1 3],8,1)
                axis_out_tmp[axis_out_tmp == 1] = 0        # Remove x operators
                axis_out_tmp[axis_out_tmp == 2] = 0        # Remove y operators => [0 0 0 3; 0 0 0 3;... 0 0 0 3]
                coef_out_tmp = sym.Matrix([1]*xyn)

                dec = np.r_[0:xyn]# 0, 1, ..., xyn
                bin_mat = PO.de2bi(dec, (xn + yn)) # np.array, 2D
                bin_mat[bin_mat == 0] = 4
                bin_mat[bin_mat == 1] = 5
                bin_mat = np.matrix(bin_mat) # np.matrix, 2D

                int_count = -1
                for jj in range(spin_no):
                    axis_v = axis_tmp[0,jj]
                    if axis_v == 1 or axis_v == 2:
                        # int_count = int_count + 1
                        int_count += 1
                        bin_vec = bin_mat[:,int_count] # np.matrix, 2D
                        axis_out_tmp[:,jj] = bin_vec
                        bin_vec = np.array(bin_vec).transpose()[0] # np.array, 1 x n, 1D

                        c_tmp = np.array([complex(1 + 0*1j)]*bin_vec.shape[0])
                    
                        if axis_v == 1: # Ix = 1/2*Ip + 1/2*Im
                            c_tmp[bin_vec == 4] = 1/2
                            c_tmp[bin_vec == 5] = 1/2

                        elif axis_v == 2: # Iy = 1/(2i)*Ip - 1/(2i)*Im
                            c_tmp[bin_vec == 4] = 1/(2*1j)
                            c_tmp[bin_vec == 5] = -1/(2*1j)

                        for kk in range(len(c_tmp)):
                            coef_out_tmp[kk] = coef_out_tmp[kk]*c_tmp[kk]
                        # coef_out_tmp = np.multiply(coef_out_tmp,c_tmp[0,:])

            coef_out_tmp = coef_out_tmp*coef_in[ii]*Ncoef_in[ii] # Ncoef => coef

            if ii == 0:
                axis_out = axis_out_tmp
                coef_out = coef_out_tmp
            else:
                axis_out = np.concatenate((axis_out, axis_out_tmp),axis=0)
                coef_out = np.concatenate((coef_out, coef_out_tmp),axis=0) # Change to col_join?

        axis_out = np.matrix(axis_out)
        coef_out = sym.Matrix(coef_out)
        obj_out = PO(self.spin_no, self.spin_label, axis_out, coef_out, 'pmz')
        obj_out = PO.CombPO(obj_out)

        obj_out.disp = disp0
        obj_out.logs = '%s\n%s' % (self.logs, obj_out.txt) # Original use char.
        return obj_out

    def pmz2xyz(self):
        # conversion from Cartesian operator basis to lowring/raising operator basis
        if self.basis != 'pmz':
            sys.exit('The basis of the object should be pmz')

        axis_in = self.axis
        coef_in = self.coef
        Ncoef_in = self.Ncoef
        spin_no = self.spin_no
        disp0 = self.disp

        for ii in range(axis_in.shape[0]):
            axis_tmp = axis_in[ii,:]# np.matrix

            if len(np.where(axis_tmp == 4)[0]) == 0 and len(np.where(axis_tmp == 5)[0]) == 0:
                # Only *z operators Iz, 2IzSz, 4IzSzKz,...
                axis_out_tmp = axis_tmp
                coef_out_tmp = sym.Matrix([1])
                xn = 0
                yn = 0
                zn = len(np.where(axis_tmp == 3)[0])

            else:
                # Conversion from Ip and Im to Ix and Iy
                xn = len(np.where(axis_tmp == 4)[0]) # p
                yn = len(np.where(axis_tmp == 5)[0]) # m
                zn = len(np.where(axis_tmp == 3)[0]) # m
                xyn = 2 ** (xn + yn)
                axis_out_tmp = np.tile(axis_tmp, (xyn, 1))
                axis_out_tmp[axis_out_tmp == 4] = 0        # Remove Ip
                axis_out_tmp[axis_out_tmp == 5] = 0        # Remove Im
                coef_out_tmp = sym.Matrix([1]*xyn)

                dec = np.r_[0:xyn]# 0, 1, ..., xyn
                bin_mat = PO.de2bi(dec, (xn + yn)) # np.array, 2D
                bin_mat[bin_mat == 0] = 2
                bin_mat[bin_mat == 1] = 1
                bin_mat = np.matrix(bin_mat) # np.matrix, 2D

                int_count = -1
                for jj in range(spin_no):
                    axis_v = axis_tmp[0,jj]
                    if axis_v == 4 or axis_v == 5:
                        # int_count = int_count + 1
                        int_count += 1
                        bin_vec = bin_mat[:,int_count] # np.matrix, 2D
                        axis_out_tmp[:,jj] = bin_vec
                        bin_vec = np.array(bin_vec).transpose()[0] # np.array, 1 x n, 1D

                        c_tmp = np.array([complex(1 + 0*1j)]*bin_vec.shape[0])
                    
                        if axis_v == 4: # Ip = Ix + 1i*Iy
                            c_tmp[bin_vec == 1] = 1
                            c_tmp[bin_vec == 2] = 1j

                        elif axis_v == 5: # Im = Ix - 1i*Iy
                            c_tmp[bin_vec == 1] = 1
                            c_tmp[bin_vec == 2] = -1j

                        for kk in range(len(c_tmp)):
                            coef_out_tmp[kk] = coef_out_tmp[kk]*c_tmp[kk]
                        # coef_out_tmp = np.multiply(coef_out_tmp,c_tmp[0,:])

            coef_out_tmp = coef_out_tmp*coef_in[ii]*Ncoef_in[ii]*(1/2)**(xn + yn + zn - 1) # Ncoef => coef
            # From 1*IpSpIz (Ncoef =1), IxSxIz IxSyIz etc. are created.
            # Then Ncoef for IxSxIz in xyz-basis is calculated as 4 automatically. 
            # To compensate 4 in the new Ncoef, it is necessary to apply (1/2)*(xn + yn + zn -1) to the new coef.
            # In this case (1/2)*(xn + yn + zn -1) = (1/2)*(2 + 0 + 1 -1) = 1/4.
            # This line is a main difference from the one in xyz2pmz().           
            if ii == 0:
                axis_out = axis_out_tmp
                coef_out = coef_out_tmp
            else:
                axis_out = np.concatenate((axis_out, axis_out_tmp),axis=0)
                coef_out = np.concatenate((coef_out, coef_out_tmp),axis=0) # Change to col_join?

        axis_out = np.matrix(axis_out)
        coef_out = sym.Matrix(coef_out)

        obj_out = PO(self.spin_no, self.spin_label, axis_out, coef_out, 'xyz')
        obj_out = PO.CombPO(obj_out)

        obj_out.logs = '%s\n%s' % (self.logs, obj_out.txt) # Original use char.
        obj_out.disp = disp0
        return obj_out

    def xyz2pol(self):
        # conversion from Cartesian operator basis to Polarization operator basis
        if self.basis != 'xyz':
            sys.exit('The basis of the object should be xyz')

        axis_in = self.axis
        coef_in = self.coef
        spin_no = self.spin_no
        disp0 = self.disp

        for ii in range(axis_in.shape[0]):
            axis_tmp = axis_in[ii,:]# np.matrix
            xyzen = 2 ** spin_no
            coef_out_tmp = sym.Matrix([1]*xyzen)

            dec = np.r_[0:xyzen]# 0, 1, ..., xyzen -1
            bin_mat = PO.de2bi(dec, spin_no) # np.array

            for jj in range(spin_no):
                axis_v = axis_tmp[0,jj]
                if axis_v == 1 or axis_v == 2:
                    bin_mat[bin_mat[:,jj] == 0, jj] = 4 # p
                    bin_mat[bin_mat[:,jj] == 1, jj] = 5 # m
                elif axis_v == 3 or axis_v == 0:
                    bin_mat[bin_mat[:,jj] == 0, jj] = 6 # a
                    bin_mat[bin_mat[:,jj] == 1, jj] = 7 # b

            axis_out_tmp = bin_mat

            for jj in range(spin_no):
                axis_v = axis_tmp[0,jj]
                bin_vec = bin_mat[:,jj]

                c_tmp = np.array([complex(1 + 0*1j)]*bin_vec.shape[0])
                
                if axis_v == 1: # Ix = 1/2*Ip + 1/2*Im
                    c_tmp[bin_vec == 4] = 1/2
                    c_tmp[bin_vec == 5] = 1/2

                elif axis_v == 2: # Iy = 1/(2i)*Ip - 1/(2i)*Im
                    c_tmp[bin_vec == 4] = 1/(2*1j)
                    c_tmp[bin_vec == 5] = -1/(2*1j)

                elif axis_v == 3: # Iz = 1/2*Ia - 1/2*Ib
                    c_tmp[bin_vec == 6] = 1/2
                    c_tmp[bin_vec == 7] = -1/2

                elif axis_v == 0: # hE = 1/2*Ia + 1/2*Ib
                    c_tmp[bin_vec == 6] = 1/2
                    c_tmp[bin_vec == 7] = 1/2

                for kk in range(len(c_tmp)):
                    coef_out_tmp[kk] = coef_out_tmp[kk]*c_tmp[kk]
                # coef_out_tmp = np.multiply(coef_out_tmp,c_tmp[0,:])

            coef_out_tmp = coef_out_tmp*coef_in[ii]*2**(spin_no - 1)

            if ii == 0:
                axis_out = axis_out_tmp
                coef_out = coef_out_tmp
            else:
                axis_out = np.concatenate((axis_out, axis_out_tmp),axis=0)
                coef_out = np.concatenate((coef_out, coef_out_tmp),axis=0) # Change to col_join?

        axis_out = np.matrix(axis_out)
        coef_out = sym.Matrix(coef_out)
        obj_out = PO(self.spin_no, self.spin_label, axis_out, coef_out, 'pol')
        obj_out = PO.CombPO(obj_out)
        obj_out.logs = '%s\n%s' % (self.logs, obj_out.txt) # Original use char.
        obj_out.disp = disp0
        return obj_out

    def pol2xyz(self):
        # conversion from Polarization operator basis to Cartesian operator basis.
        if self.basis != 'pol':
            sys.exit('The basis of the object should be pol')

        axis_in = self.axis
        coef_in = self.coef
        spin_no = self.spin_no
        disp0 = self.disp

        for ii in range(axis_in.shape[0]):
            axis_tmp = axis_in[ii,:]# np.matrix
            xyzen = 2 ** spin_no
            coef_out_tmp = sym.Matrix([1]*xyzen) # sym.ones()

            dec = np.r_[0:xyzen]# 0, 1, ..., xyzen -1
            bin_mat = PO.de2bi(dec, spin_no) # np.array

            # Ia =>  Iz + hE
            # Ib => -Iz + hE
            # Ip =>  Ix + Iy*1i
            # Im =>  Ix - Iy*1i
            for jj in range(spin_no):
                axis_v = axis_tmp[0,jj]
                if axis_v == 4 or axis_v == 5: # p or m
                    bin_mat[bin_mat[:,jj] == 1, jj] = 2 # y
                    bin_mat[bin_mat[:,jj] == 0, jj] = 1 # x
                elif axis_v == 6 or axis_v == 7: # a or b
                    bin_mat[bin_mat[:,jj] == 0, jj] = 3 # z
                    bin_mat[bin_mat[:,jj] == 1, jj] = 0 # hE
                else:
                    bin_mat[:, jj] = 0 # Convert any remained 1 to 0.


            axis_out_tmp = bin_mat

            for jj in range(spin_no):
                axis_v = axis_tmp[0,jj]
                bin_vec = bin_mat[:,jj]

                c_tmp = np.array([complex(1 + 0*1j)]*bin_vec.shape[0])
                
                if axis_v == 4: # Ip = Ix + 1i*Iy
                    c_tmp[bin_vec == 1] = 1
                    c_tmp[bin_vec == 2] = 1j

                elif axis_v == 5: # Im = Ix - 1i*Iy
                    c_tmp[bin_vec == 1] = 1
                    c_tmp[bin_vec == 2] = -1*1j

                elif axis_v == 6: # Ia = hE + Iz
                    c_tmp[bin_vec == 0] = 1
                    c_tmp[bin_vec == 3] = 1

                elif axis_v == 7: # Ib = hE - Iz
                    c_tmp[bin_vec == 0] = 1
                    c_tmp[bin_vec == 3] = -1

                for kk in range(len(c_tmp)):
                    coef_out_tmp[kk] = coef_out_tmp[kk]*c_tmp[kk]
                # coef_out_tmp = np.multiply(coef_out_tmp,c_tmp[0,:])

            coef_out_tmp = coef_out_tmp*coef_in[ii]*(1/2)**(spin_no - 1)

            if ii == 0:
                axis_out = axis_out_tmp
                coef_out = coef_out_tmp
            else:
                axis_out = np.concatenate((axis_out, axis_out_tmp),axis=0)
                coef_out = np.concatenate((coef_out, coef_out_tmp),axis=0) # Change to col_join?

        axis_out = np.matrix(axis_out)
        coef_out = sym.Matrix(coef_out)
        obj_out = PO(self.spin_no, self.spin_label, axis_out, coef_out, 'xyz')
        obj_out = PO.CombPO(obj_out)
        obj_out.logs = '%s\n%s' % (self.logs, obj_out.txt) # Original use char.
        obj_out.disp = disp0
        return obj_out

    def pmz2pol(self):
        # conversion from lowring/raising operator basis to Polarization operator basis
        if self.basis != 'pmz':
            sys.exit('The basis of the object should be pmz')

        axis_in = self.axis
        coef_in = self.coef
        Ncoef_in = self.Ncoef
        spin_no = self.spin_no
        disp0 = self.disp

        for ii in range(axis_in.shape[0]):
            axis_tmp = axis_in[ii,:]# np.matrix

            if len(np.where(axis_tmp == 3)[0]) == 0 and len(np.where(axis_tmp == 0)[0]) == 0:
                # Not including Iz (3) and E (0)
                axis_out_tmp = axis_tmp
                coef_out_tmp = sym.Matrix([1])
            else:
                # Conversion from Iz and E to Ia and Ib
                xn = len(np.where(axis_tmp == 0)[0])
                yn = len(np.where(axis_tmp == 3)[0])
                xyn = 2 ** (xn + yn)
                axis_out_tmp = np.tile(axis_tmp, (xyn, 1))
                axis_out_tmp[axis_out_tmp == 0] = 0    
                axis_out_tmp[axis_out_tmp == 3] = 0  
                coef_out_tmp = sym.Matrix([1]*xyn)

                dec = np.r_[0:xyn]# 0, 1, ..., xyn
                bin_mat = PO.de2bi(dec, (xn + yn)) # np.array, 2D
                bin_mat[bin_mat == 0] = 6
                bin_mat[bin_mat == 1] = 7
                bin_mat = np.matrix(bin_mat) # np.matrix, 2D

                int_count = -1
                for jj in range(spin_no):
                    axis_v = axis_tmp[0,jj]
                    if axis_v == 3 or axis_v == 0:
                        # int_count = int_count + 1
                        int_count += 1
                        bin_vec = bin_mat[:,int_count] # np.matrix, 2D
                        axis_out_tmp[:,jj] = bin_vec
                        bin_vec = np.array(bin_vec).transpose()[0] # np.array, 1 x n, 1D

                        c_tmp = np.array([complex(1 + 0*1j)]*bin_vec.shape[0])
                    
                        if axis_v == 3: # Iz = 1/2*Ia - 1/2*Ib
                            c_tmp[bin_vec == 6] = 1/2 # No Nceof in pmz, thus need to define
                            c_tmp[bin_vec == 7] = -1/2

                        elif axis_v == 2: # E = Ia + Ib
                            c_tmp[bin_vec == 4] = 1 # No Nceof in pmz, thus need to define
                            c_tmp[bin_vec == 5] = 1

                        for kk in range(len(c_tmp)):
                            coef_out_tmp[kk] = coef_out_tmp[kk]*c_tmp[kk]
                        # coef_out_tmp = np.multiply(coef_out_tmp,c_tmp[0,:])

            coef_out_tmp = coef_out_tmp*coef_in[ii]*Ncoef_in[ii] # Ncoef => coef

            if ii == 0:
                axis_out = axis_out_tmp
                coef_out = coef_out_tmp
            else:
                axis_out = np.concatenate((axis_out, axis_out_tmp),axis=0)
                coef_out = np.concatenate((coef_out, coef_out_tmp),axis=0) # Change to col_join?

        axis_out = np.matrix(axis_out)
        coef_out = sym.Matrix(coef_out)
        obj_out = PO(self.spin_no, self.spin_label, axis_out, coef_out, 'pol')
        obj_out = PO.CombPO(obj_out)
        obj_out.logs = '%s\n%s' % (self.logs, obj_out.txt) # Original use char.
        obj_out.disp = disp0
        return obj_out

    def pol2pmz(self):
        # conversion from Polarization operator basis to lowring/raising operator basis
        if self.basis != 'pol':
            sys.exit('The basis of the object should be pol')

        axis_in = self.axis
        coef_in = self.coef
        Ncoef_in = self.Ncoef
        spin_no = self.spin_no
        disp0 = self.disp

        for ii in range(axis_in.shape[0]):
            axis_tmp = axis_in[ii,:]# np.matrix

            if len(np.where(axis_tmp == 6)[0]) == 0 and len(np.where(axis_tmp == 7)[0]) == 0:
                # Not including Ia (6) and Ib (7)
                axis_out_tmp = axis_tmp
                coef_out_tmp = sym.Matrix([1])
            else:
                # Conversion from Ia and Ib to Iz and E
                xn = len(np.where(axis_tmp == 6)[0])
                yn = len(np.where(axis_tmp == 7)[0])
                xyn = 2 ** (xn + yn)
                axis_out_tmp = np.tile(axis_tmp, (xyn, 1))
                axis_out_tmp[axis_out_tmp == 6] = 0    
                axis_out_tmp[axis_out_tmp == 7] = 0  
                coef_out_tmp = sym.Matrix([1]*xyn)

                dec = np.r_[0:xyn]# 0, 1, ..., xyn
                bin_mat = PO.de2bi(dec, (xn + yn)) # np.array, 2D
                bin_mat[bin_mat == 0] = 0
                bin_mat[bin_mat == 1] = 3
                bin_mat = np.matrix(bin_mat) # np.matrix, 2D

                int_count = -1
                for jj in range(spin_no):
                    axis_v = axis_tmp[0,jj]
                    if axis_v == 6 or axis_v == 7:
                        # int_count = int_count + 1
                        int_count += 1
                        bin_vec = bin_mat[:,int_count] # np.matrix, 2D
                        axis_out_tmp[:,jj] = bin_vec
                        bin_vec = np.array(bin_vec).transpose()[0] # np.array, 1 x n, 1D

                        c_tmp = np.array([complex(1 + 0*1j)]*bin_vec.shape[0])
                    
                        if axis_v == 6: # Ia = 1/2*E + Iz
                            c_tmp[bin_vec == 0] = 1/2 # No Nceof with pmz basis, thus need to provide 1/2
                            c_tmp[bin_vec == 3] = 1

                        elif axis_v == 7: # Ib = 1/2*E + Iz
                            c_tmp[bin_vec == 0] = 1/2 # No Nceof with pmz basis, thus need to provide 1/2
                            c_tmp[bin_vec == 3] = -1

                        for kk in range(len(c_tmp)):
                            coef_out_tmp[kk] = coef_out_tmp[kk]*c_tmp[kk]
                        # coef_out_tmp = np.multiply(coef_out_tmp,c_tmp[0,:])

            coef_out_tmp = coef_out_tmp*coef_in[ii]*Ncoef_in[ii] # Ncoef => coef

            if ii == 0:
                axis_out = axis_out_tmp
                coef_out = coef_out_tmp
            else:
                axis_out = np.concatenate((axis_out, axis_out_tmp),axis=0)
                coef_out = np.concatenate((coef_out, coef_out_tmp),axis=0) # Change to col_join?

        axis_out = np.matrix(axis_out)
        coef_out = sym.Matrix(coef_out)
        obj_out = PO(self.spin_no, self.spin_label, axis_out, coef_out, 'pmz')
        obj_out = PO.CombPO(obj_out)
        obj_out.logs = '%s\n%s' % (self.logs, obj_out.txt) # Original use char.
        obj_out.disp = disp0
        return obj_out

    def de2bi(d, n):
        # https://stackoverflow.com/questions/44449106/using-numpy-arrays-for-integer-and-array-inputs
        d = np.array(d)
        d = np.reshape(d, (1, -1))
        power = np.flipud(2**np.arange(n))
        g = np.zeros((np.shape(d)[1], n))

        for i, num in enumerate(d[0]):
            g[i] = num * np.ones((1,n))
        b = np.floor((g%(2*power))/power)
        return b

    def set_basis(self, basis_in):

        # change the expression of obj using basis_in
        if self.basis != basis_in:
            if self.basis == 'xyz' and basis_in == 'pmz':
                obj = PO.xyz2pmz(self)
                branch_v = 1
            elif self.basis == 'xyz' and basis_in == 'pol':
                obj = PO.xyz2pol(self)
                branch_v = 2
            elif self.basis == 'pmz' and basis_in == 'xyz':
                obj = PO.pmz2xyz(self)
                branch_v = 3
            elif self.basis == 'pmz' and basis_in == 'pol':
                obj = PO.pmz2pol(self)
                branch_v = 4
            elif self.basis == 'pol' and basis_in == 'xyz':
                obj = PO.pol2xyz(self)
                branch_v = 5
            elif self.basis == 'pol' and basis_in == 'pmz':
                obj = PO.pol2pmz(self)
                branch_v = 6
        else:
            obj = self
            branch_v = 7

        # print('branch_v',branch_v)
        return obj


    def add_sub_fun(self, other, add_sub_coef):
        spin_no = self.spin_no
        spin_label = self.spin_label
        basis = self.basis
        self_fc = add_sub_coef[0]
        other_fc = add_sub_coef[1]

        if isinstance(other, PO):
            if self.axis.shape[1] != other.axis.shape[1]:
                sys.exit('The number of spin types for obj1 and obj2 must be same!')

            # xyz + pmz => pmz
            if (self.basis == 'xyz' and other.basis == 'pmz') or \
                (self.basis == 'pmz' and other.basis == 'xyz'):

                self = PO.set_basis(self,'pmz')
                other = PO.set_basis(other,'pmz')

            # xyz + pol => pol
            # pmz + pol => pol
            elif (self.basis == 'xyz' and other.basis == 'pol') or \
                (self.basis == 'pol' and other.basis == 'xyz') or \
                (self.basis == 'pmz' and other.basis == 'pol') or \
                (self.basis == 'pol' and other.basis == 'pmz'):

                self = PO.set_basis(self,'pol')
                other = PO.set_basis(other,'pol')
            
            axis = np.concatenate((self.axis, other.axis),axis=0)
            coef = sym.Matrix([self_fc*self.coef, other_fc*other.coef])
        else:
            if self.basis == 'xyz':
                coef_fc = sym.Integer(2)

            else:
                coef_fc = sym.Integer(1)

            axis = np.concatenate((self.axis, np.matrix([0]*spin_no)),axis=0)
            coef = sym.Matrix([self_fc*self.coef, other_fc*coef_fc*other]) # Case of 'xyz', 1 = hE*2

        obj_out = PO(spin_no, spin_label, axis, coef, self.basis)
        obj_out = PO.CombPO(obj_out)
        return obj_out

    def __add__(self, other):
    # __add__( self , other )	self + other
        return PO.add_sub_fun(self, other, [1,1])
    
    def __radd__(self, other):
    # __radd__( self , other )	other + self
        return PO.add_sub_fun(self, other, [1,1])

    def __sub__(self, other):
        # __sub__( self , other )	 self - other
        return PO.add_sub_fun(self, other, [1,-1])
    
    def __rsub__(self, other):
        # __rsub__( self , other )	 other - self
        return PO.add_sub_fun(self, other, [-1,1])

    def __pos__(self):
        return self   

    def __neg__(self):
        spin_no = self.spin_no
        spin_label = self.spin_label
        basis = self.basis
        axis = self.axis
        coef = -1*self.coef

        obj_out = PO(spin_no, spin_label, axis, coef, basis)
        return PO.CombPO(obj_out)      

    def __mul__(self, other):
        #  __mul__( self , other )	 self * other
        if isinstance(other, PO) == 0:
            spin_no = self.spin_no
            spin_label = self.spin_label
            basis = self.basis
            axis = self.axis
            coef = other*self.coef

            obj_out = PO(spin_no, spin_label, axis, coef, basis)
            return PO.CombPO(obj_out)     

        elif isinstance(other, PO) == 1:
            if self.axis.shape[1] != other.axis.shape[1]:
                sys.exit('The number of spin types for obj1 and obj2 must be same!')

            # xyz * xyz
            if self.basis == 'xyz' and other.basis == 'xyz':
                basis_org = 'xyz'

            # xyz * pmz or pmz * xyz
            elif (self.basis == 'xyz' and other.basis == 'pmz') or \
                (self.basis == 'pmz' and other.basis == 'xyz'):
                basis_org = 'pmz'

            # xyz * pol or pol * xyz 
            elif (self.basis == 'xyz' and other.basis == 'pol') or \
                (self.basis == 'pol' and other.basis == 'xyz'):
                basis_org = 'pol'

            # pmz * pmz
            elif self.basis == 'pmz' and other.basis == 'pmz':
                basis_org = 'pmz'
            
            # pmz * pol or pol * pmz
            elif (self.basis == 'pmz' and other.basis == 'pol') or \
                (self.basis == 'pol' and other.basis == 'pmz'):
                basis_org = 'pol'

            # pol * pol
            elif self.basis == 'pol' and other.basis == 'pol':
                basis_org = 'pol'

            self = PO.set_basis(self,'xyz')
            other = PO.set_basis(other,'xyz')

            axis1 = self.axis
            coef1 = self.coef
            Ncoef1 = self.Ncoef

            axis2 = other.axis
            coef2 = other.coef
            Ncoef2 = other.Ncoef

            row1 = axis1.shape[0]
            row2 = axis2.shape[0]

            axis1M = np.tile(axis1, (row2, 1))
            coef1V = np.tile(coef1, (row2, 1))# tile: numpy, coef: sympy, coef1V: numpy
            Ncoef1V = np.tile(Ncoef1, row2) # Since Ncoef1V is array, order is (1, row2)

            axis2M = np.repeat(axis2, row1, axis=0)
            coef2V = np.repeat(coef2, row1, axis=0)# repeat: numpy, coef: sympy, coef2V: numpy
            Ncoef2V = np.repeat(Ncoef2, row1, axis=0)

            comp_M = np.multiply(axis1M, axis2M) > 0
            comp_M = comp_M * 1
            comp_V = comp_M.sum(axis=1)

            axis_new = axis1M + axis2M
            coef_new = np.multiply(coef1V, coef2V)# multiply: numpy
            Ncoef_new = np.multiply(Ncoef1V, Ncoef2V)

            at = np.matrix([[0, 3, 2],[3, 0, 1],[2, 1, 0]])
            ct = 1/2*np.matrix([[1/2, 1*I, -1*I],[-1*I, 1/2, 1*I],[1*I, -1*I, 1/2]])

            for ii in range(comp_V.shape[0]):
                if comp_V[ii] == 1:
                    comp_M_row = np.array(comp_M[ii,:])
                    jj_vec = np.where(comp_M_row == 1)[1]# There are two output

                    for jj in jj_vec:
                        id1 = axis1M[ii,jj]
                        id2 = axis2M[ii,jj]
                        id1 = id1.astype(np.int64)
                        id2 = id2.astype(np.int64)
                        axis_new[ii,jj] = at[id1 - 1, id2 - 1]
                        coef_new[ii] = coef_new[ii]*ct[id1 - 1, id2 - 1]

            coef_new = sym.Matrix(coef_new) # The conversion to symbolic matrix is necessary for getM().
            obj_base = PO(self.spin_no, self.spin_label, axis_new, coef_new, self.basis)  

            Ncoef_in = Ncoef_new
            Ncoef_out = obj_base.Ncoef
            Ncoef_rate = np.divide(Ncoef_in,Ncoef_out)

            id_coef = np.where(Ncoef_rate != 1)[0]# There is one output

            for ii in id_coef:
                coef_new[ii] = coef_new[ii] * Ncoef_rate[ii]

            obj_out = PO(self.spin_no, self.spin_label, axis_new, coef_new, self.basis)
            obj_out = PO.CombPO(obj_out)
            return PO.set_basis(obj_out,basis_org)
          
    def __rmul__(self, other):
        # __rmul__( self , other )	other * self
        if isinstance(other, PO) == 0:#  this is for a*Ix style.
            return self * other
        
    def __truediv__(self, other):
        if isinstance(self, PO) == 1 and isinstance(other, PO) != 1:
            obj = 1/other*self
        else:
            sys.exit('PO-class object cannot the be the divisor!') 
            # This branch works when both self and other are PO.
        return obj

    def __pow__(self, other):
        if (other % 1) == 0:
            if other == 0:
                    obj = 2*hE
            elif other > 0:
                for ii in range(int(other)):
                    if ii == 0:
                        obj = self
                    else:
                        obj = obj*self
            else:
                sys.exit('Can''t calculate obj^-n!')
        else:
            sys.exit('n in obj^n should be 0 or a positive integer!')
        return obj



    def create(spin_label, *args):
        # Function to create basic product operators.

        if len(args) == 0: # add_cell and symcoef_switch are not given
            add_cell = []
            symcoef_switch = 'on'
        elif len(args) == 1: # add_cell is given as args
            add_cell = args[0]
            symcoef_switch = 'on'

        spin_no = len(spin_label)
        for ii in range(spin_no):
            for jj in range(7):
                if jj == 0:
                    phase_s = 'x'
                    basis = 'xyz'
                elif jj == 1:
                    phase_s = 'y'
                    basis = 'xyz'
                elif jj == 2:
                    phase_s = 'z'
                    basis = 'xyz' 
                elif jj == 3:
                    phase_s = 'p'
                    basis = 'pmz' 
                elif jj == 4:
                    phase_s = 'm'
                    basis = 'pmz'
                elif jj == 5:
                    phase_s = 'a'
                    basis = 'pol' 
                elif jj == 6:
                    phase_s = 'b'
                    basis = 'pol'
                # end 

                axis = [0]*spin_no
                axis[ii] = jj + 1
                axis = np.matrix(axis)
                sp = spin_label[ii] + phase_s
                coef = 1

                cmd_s = sp + ' = PO(spin_no, spin_label, axis, coef, basis)'
                # https://qiita.com/Masahiro_T/items/1e5ff567fa5a3c5aebf9
                exec(cmd_s, locals(), globals())
            # end
        # end

        axis = [0]*spin_no
        axis = np.matrix(axis)
        sp = 'hE'
        basis = 'xyz'
        coef = 1
        cmd_s = sp + ' = PO(spin_no, spin_label, axis, coef, basis)'
        exec(cmd_s, locals(), globals())        

        if symcoef_switch == 'on':
            PO.symcoef(spin_label, add_cell)


    def symcoef(spin_label_cell, *args):

        spin_no = len(spin_label_cell)

        # Frequency o & gyromagnetic ratio g
        header_list = ['o', 'g']
        for jj in range(len(header_list)):
            header_str = header_list[jj]
            for ii in range(spin_no):
                sp = spin_label_cell[ii]
                varname = header_str + sp
                PO.search2create(varname)

                varname = header_str + sp[0]
                PO.search2create(varname)

                varname = header_str + sp[-1]
                PO.search2create(varname)

                if PO.__index_switch == 1: # 1-based index
                    varname = header_str + str(ii + 1) # MATLAB ahd Human: 1-based index
                else:
                    varname = header_str + str(ii) # Python: 0-based index
                PO.search2create(varname)

        # coupling J
        for ii in range(spin_no):
            for jj in range(spin_no):
                if ii < jj:
                    sp_ii = spin_label_cell[ii]
                    sp_jj = spin_label_cell[jj]

                    varname = 'J' + sp_ii + sp_jj
                    PO.search2create(varname)

                    varname = 'J' + sp_ii[0] + sp_jj[0]
                    PO.search2create(varname)

                    varname = 'J' + sp_ii[-1] + sp_jj[-1]
                    PO.search2create(varname)

                    if PO.__index_switch == 1: # 1-based index
                        varname = 'J' + str(ii + 1) + str(jj + 1) # MATLAB ahd Human: 1-based index
                    else:
                        varname = 'J' + str(ii) + str(jj) # Python: 0-based index

                    PO.search2create(varname)
                    
        symcoef_list = ['a', 'b', 'c', 'd', 'f', 'gH', 'gC', 'gN', 'q', 't1', 't2', 't3', 't4', 't', 'w', 'B', 'J', 'G']

        if len(args) == 0:
            add_cell = []
        elif len(args) == 1: # add_cell was given
            add_cell = args[0]

        symcoef_list = symcoef_list + add_cell

        for ii in range(len(symcoef_list)):
            varname = symcoef_list[ii]
            PO.search2create(varname)


    def search2create(varname):
        if varname not in globals():
            cmd_s = varname + " = symbols('" + varname + "')"
            exec(cmd_s, locals(), globals())

    def load_PS(fname):
        f = open(fname, 'r')
        mlines = f.readlines()
        para_bin = 0
        ps_bin = 0
        rho_ini_bin = 0
        para_lines = []
        ps_lines = []

        for ii in range(len(mlines)):
            mline = mlines[ii]
            if mline.find('# Para begin #') != -1:
                para_bin = 1
                continue
            elif mline.find('# Para end #') != -1:
                para_bin = 0
            elif mline.find('# PS begin #') != -1:
                ps_bin = 1
                continue
            elif mline.find('# PS end #') != -1:
                ps_bin = 0
            elif mline.find('rho_ini') == 0:
                rho_ini_bin = 1

            if para_bin == 1:
                if rho_ini_bin == 0:
                    para_lines = para_lines + [mline[0:-1]]
                elif rho_ini_bin == 1:
                    rho_ini_line = mline[0:-1]
                    rho_ini_bin = 0
            elif ps_bin == 1:
                ps_lines = ps_lines + [mline[0:-1]]
        
        return para_lines, ps_lines, rho_ini_line

    def run_PS(fname):
        # MATLAB version
        # Variables are only effectvie in run_PS().
        #
        # Python version
        # Can exec() in Python handle variable-assignment to store variables in a local workspace only?
        # => It seems NO.
        # https://discuss.python.org/t/defining-local-variable-with-exec/8665/
        # Due to this reasons, any variables used for the actual spin dynamics calculation 
        # should be generated by exec() and they should be stored in the global workspace.

        para_lines, ps_lines, rho_ini_line = PO.load_PS(fname)

        for ii in range(len(para_lines)):
            exec(para_lines[ii], locals(), globals())
            if para_lines[ii].find('no_ph') >= 0:
                exec('ph_cell = [0]*(no_ph + 1)',locals(), globals())

        # Optional parameters
        varname = 'phid'
        if varname not in globals():
            ph_length = 0
            for ii in range(len(ph_cell)):
                ph_length = max(ph_length, len(ph_cell[ii]))

            h_length = max(ph_length, len(phRtab))
            exec('h_lenght = ' + str(h_length), locals(), globals()) 
            exec('phid = list(range(h_length))', locals(), globals())

        varname = 'coef_cell'
        if varname not in globals():
            exec('coef_cell = []', locals(), globals())

        varname = 'spin_label_cell'
        if varname not in globals():
            exec('spin_label_cell = PO.spin_label_cell_default', locals(), globals())

        varname = 'disp_bin'
        if varname not in globals():
            exec('disp_bin = 1', locals(), globals())

        varname = 'obs_cell'
        if varname not in globals():
            exec("obs_cell = ['*']", locals(), globals())
        
        PO.create(spin_label_cell, coef_cell) # At what scope variables are created?

        exec(rho_ini_line, locals(), globals())
        exec('rho_ini.disp = disp_bin', locals(), globals())

        int_ii = 0
        rho_cell = [0]*len(phid)
        rho_detect_cell = [0]*len(phid)
        rho_total = 0
        a0_M = sym.Matrix([])
        rho_M = sym.Matrix([])

        for idx, ii in enumerate(phid): # idx: index, ii: element of phid
            print('')
            print('ii:', ii)
            if 'ph_cell' in globals():
                # for jj in range(len(ph_cell)): # Let me leave this block although it won't be used.
                    # if PO.__index_switch == 1:
                    #     cmd_s = 'ph' + str(jj + 1) + " = PO.phmod(ph_cell[jj], ii)" # ph_cell[0] => ph1
                    # else:
                    #     cmd_s = 'ph' + str(jj) + " = PO.phmod(ph_cell[jj], ii)" # ph_cell[0] => ph0

                for jj in range(len(ph_cell) - 1):
                    cmd_s = 'ph' + str(jj + 1) + " = PO.phmod(ph_cell[jj + 1], ii)" # ph_cell[1] => ph1
                    exec(cmd_s, locals(), globals())

            exec('phR = PO.phmod(phRtab, ii)', locals(), globals())
            exec('rho = copy.deepcopy(rho_ini)', locals(), globals()) # rho should be deep-copied from rho_ini

            if rho_ini.disp == 1:
                print(rho)
            
            # Run pulse sequnce
            for jj in range(len(ps_lines)):
                try:
                    exec(ps_lines[jj], locals(), globals()) # At what scope this line is executed?
                except:
                     sys.exit('PS Line ' + str(jj))

            # Store rho
            rho_cell[int_ii] = rho # rho can be accessible from run_PS, global

            # Receiver
            rho_detect = PO.receiver(rho, phR) # Need deepcopy of rho in the receiver. Otherwise, the stored rho in rho_cell is changed by receiver.
            rho_detect_cell[int_ii] = rho_detect
            rho_total = rho_detect + rho_total

            # Create a0_M, rho_M
            a0_V, rho_V = rho.SigAmp(obs_cell, phR)
            a0_M = a0_M.col_join(a0_V)
            rho_M = rho_M.col_join(rho_V)
            
            int_ii += 1 # End of For loop

        # Observable
        rho_obs = PO.observable(rho_total, obs_cell)
        print('\nPulse Seuqnce: Done!')
        return rho_cell, rho_detect_cell, rho_total, rho_obs, a0_M, rho_M

    def observable(self, *args):

        spin_label_cell = self.spin_label
        if len(args) < 1:
            sp_cell = spin_label_cell
        elif len(args) == 1:
            sp_cell = args[0]
        
        s0 = self.logs
        disp0 = self.disp
        obj = copy.deepcopy(self)# Give a differnt ID to obj from self
        basis_org = obj.basis
        if basis_org != 'xyz':
            obj = PO.set_basis(obj, 'xyz')

        # Selection of observable terms
        axis_in = obj.axis
        id_in = ((axis_in == 1)*1 + (axis_in == 2)*1).sum(1) # terms with only one x or one y
        id_in = np.where(id_in == 1)[0].astype(int).tolist() # 0-based Index

        if PO.__index_switch == 1:
            id_in = [x + 1 for x in id_in] # 0-based index => 1-based index

        obj = PO.selPO(obj, id_in)        

        id_vec, _ = PO.sp2id(sp_cell, spin_label_cell) # id_vec: 0-based index

        id_in = np.array([])
        for ii in range(obj.axis.shape[0]):
            axis_tmp = obj.axis[ii,:]
            if ((axis_tmp[0,id_vec] == 1)*1 + (axis_tmp[0,id_vec] == 2)*1).sum(1) == 1:
                id_in = np.concatenate((id_in,[ii]),0)

        id_in = id_in.astype(int).tolist() # 0-based ndex

        if PO.__index_switch == 1:
            id_in = [x + 1 for x in id_in] # 0-based index => 1-based index

        obj = PO.selPO(obj, id_in)
        obj = PO.set_basis(obj, basis_org)
        obj.disp = disp0

        # https://docs.python.org/3/tutorial/datastructures.html#list-comprehensions
        # In cotnrast with Matlab, spin_label_cell[id_in] won't work in Python.
        # instead, use List Comprehensions
        s_out = 'Selecting Observable of %s' % ' '.join([spin_label_cell[x] for x in id_vec])
        s1 = '%s' % (s_out)
        s2 = '    %s' % (obj.txt)
        obj.logs = '%s\n%s\n%s' %(s0, s1, s2)

        if disp0 == 1:
            print('')
            print('%s' % (s1))
            print('%s' % (s2))

        return obj

    def selPO(self, id_in):
        # id_in: 1-based index if PO.__index_switch = 1.
        #        0-based index if PO.__index_switch = 0.
        # PO.findterm() will handle the conversion from 1-based index to 0-based index if PO.__index_switch is 1.

        s0 = self.logs
        id_in = self.findterm(id_in) # The outpted id_in is 0-based index.
        axis_tmp = self.axis[id_in,:]
        coef_tmp = self.coef[id_in,:]
        self = PO(self.spin_no, self.spin_label, axis_tmp, coef_tmp, self.basis)

        s_out = 'selPO'
        s1 = '%s' % (s_out)
        s2 = '    %s' % (self.txt)
        self.logs = '%s\n%s\n%s' %(s0, s1, s2)

        return self


    def delPO(self, id_in):
        # id_in: 1-based index if PO.__index_switch = 1.
        #        0-based index if PO.__index_switch = 0.
        # PO.findterm() will handle the conversion from 1-based index to 0-based index if PO.__index_switch is 1.

        s0 = self.logs
        id_in = self.findterm(id_in) # The outpted id_in is 0-based index.

        axis_tmp = self.axis
        axis_tmp = np.delete(axis_tmp, id_in, axis=0)

        if axis_tmp.shape[0] == 0: # Delete everything
            axis_tmp = np.matrix(np.zeros((1, self.spin_no))).astype(int)
            coef_tmp = sym.Matrix([0])
        else:
            coef_tmp = copy.deepcopy(self.coef) # Need to deepcopy
            for ii in range(len(id_in) - 1, -1, -1):
                coef_tmp.row_del(id_in[ii])

        self = PO(self.spin_no, self.spin_label, axis_tmp, coef_tmp, self.basis)

        s_out = 'delPO'
        s1 = '%s' % (s_out)
        s2 = '    %s' % (self.txt)
        self.logs = '%s\n%s\n%s' %(s0, s1, s2)
        
        return self

    def findterm(self, id_in):
        # id_out: 0-based index
        # id_in: 1-based index if PO.__index_switch = 1.
        #        0-based index if PO.__index_switch = 0.
        # id_out = findterm(rho, [1, 2]) # These indexes are shown in PO.dispPO().
        # id_out = findterm(rho, ['Ix'])
        # id_out = findterm(rho, ['IzSz'])
        # id_out = findterm(rho, ['Ix', 'IzSz'])
        # id_out = findterm(rho, ['IxS*'])
        # id_out = findterm(rho, ['I*S*', 'I*S*K*'])

        spin_label_cell = self.spin_label
        spin_no = self.spin_no
        sp_cell = id_in
        id_out = np.array([])

        for ii in range(len(sp_cell)):
            sp = sp_cell[ii]

            if type(sp) == str:
                axis_tmp = np.matrix([0]*spin_no)
                for jj in range(len(spin_label_cell)):

                    spin_label_tmp = spin_label_cell[jj]

                    if sp.find(spin_label_tmp) >= 0:
                        id_tmp = jj
                        phase_s = sp[sp.find(spin_label_tmp) + 1]
                        
                        if phase_s == 'x':
                            phase_id = 1
                        elif phase_s == 'y':
                            phase_id = 2
                        elif phase_s == 'z':
                            phase_id = 3
                        elif phase_s == 'p':
                            phase_id = 4
                        elif phase_s == 'm':
                            phase_id = 5
                        elif phase_s == 'a':
                            phase_id = 6
                        elif phase_s == 'b':
                            phase_id = 7
                        elif phase_s == '*':
                            phase_id = np.matrix(list(range(1,8))).transpose()
                        else:
                            phase_id = 0

                        if type(phase_id) == int:
                            axis_tmp[:, id_tmp] = phase_id

                        elif phase_id.shape[0] == 7:
                            rate_tmp = axis_tmp.shape[0]
                            axis_tmp = np.tile(axis_tmp, (7, 1)) # [1,2,3,..,6,7,1,2,3,...,6,7,]
                            phase_id_vec = np.repeat(phase_id, rate_tmp, axis=0) # [1,1,1,1,...,7,7,7]
                            axis_tmp[:, id_tmp] = phase_id_vec

                # Detect the IDs of obj.axis corresponding to axis_tmp
                id_tmp = np.where(PO.ismember_row(self.axis, axis_tmp))[0]

            elif type(sp) == int:
                if PO.__index_switch == 1:
                    id_tmp = [sp - 1] # Converting 1-based index to 0-based index
                else:
                    id_tmp = [sp] # Keeping 0-based index
            
            id_out = np.concatenate((id_out,id_tmp),0)

        # sort => int => list
        id_out = np.sort(id_out).astype(int).tolist() # 0-based index

        return id_out

    def ismember_row(a,b):
        # https://stackoverflow.com/questions/71392305/python-equivalent-of-matlab-ismember-on-rows-for-large-arrays
        # Get the unique row index
        _, rev = np.unique(np.concatenate((b,a)),axis=0,return_inverse=True)
        # Split the index
        a_rev = rev[len(b):]
        b_rev = rev[:len(b)]
        # Return the result:
        return np.isin(a_rev,b_rev)
    
    def receiver(self, phR):

        q = PO.ph2q(phR)
        disp0 = self.disp
        obj = copy.deepcopy(self)# Give a differnt ID to obj from self
        obj.disp = 0    
        s0 = obj.logs
        sp_cell = obj.spin_label
        q_cell = [-q]*len(sp_cell)
        obj = PO.cs(obj, sp_cell, q_cell)
        obj.disp = disp0

        s_out = 'Receiver with %s' % (PO.ph_num2str(phR))

        s1 = '%s' % (s_out)
        s2 = '    %s' % (obj.txt)
        self.logs = '%s\n%s\n%s' %(s0, s1, s2)

        if disp0 == 1:
            print('%s' % (s1))
            print('%s' % (s2))

        return obj

    def dispPO(self):
        print('')
        if PO.asterisk_bin == 1:
            s_header = '%-4s%7s%-20s %-s' % ('ID', '', 'Term', 'Coef')
        else:
            s_header = '%-4s%5s%-20s %-s' % ('ID', '',  'Term', 'Coef')

        print(s_header)
        Ncoef = self.Ncoef
        coef = self.coef
        axis = self.axis

        for ii in range(axis.shape[0]):
            axis_tmp = np.array(axis)[ii]
            pt = PO.axis2pt(self.spin_label, axis_tmp)

            asterisk_str = '*'
            if Ncoef[ii] == 0.5:
                Ncoef_t = '1/2'
            elif Ncoef[ii] == 1.0:
                Ncoef_t = ''
                asterisk_str = ' '
            else:
                Ncoef_t = str(Ncoef[ii].astype(np.int64))

            coef_t = str(sym.simplify(coef[ii]))

            if PO.__index_switch == 1:
                line_id = ii + 1
            else:
                line_id = ii

            if PO.asterisk_bin == 1:
                s_out = '%-4d%6s%s%-20s %s' % (line_id, Ncoef_t, asterisk_str, pt, coef_t)
            else:
                s_out = '%-4d%5s%-20s %s' % (line_id, Ncoef_t, pt, coef_t)

            print(s_out)
        print('')

    def SigAmp(self, sp_cell, phR):

        spin_label_cell = self.spin_label
        id_vec,_ = PO.sp2id(sp_cell, spin_label_cell)

        sp_cell_tmp = [spin_label_cell[x] for x in id_vec]

        for ii in range(len(sp_cell_tmp)):
            sp_tmp = sp_cell_tmp[ii]
            sp_m = sp_tmp + 'm'
            exec('obsPO = '+ sp_m, {}, globals())
            if ii == 0:
                obsPO_M = obsPO.M
            else:
                obsPO_M = obsPO_M + obsPO.M

        a0_V = sym.Matrix([])
        rho_V = sym.Matrix([])

        rho_M = self.coherence
        self_M = self.M

        for ii in range(obsPO_M.shape[0]):
            for jj in range(obsPO_M.shape[1]):
                v_tmp = self_M[ii, jj]*obsPO_M[ii,jj]
                if not isinstance(obsPO_M[ii,jj], sym.core.numbers.Zero):# This selection should be based on obsPO_M[ii,jj], not v_tmp
                    a0_V = a0_V.row_join(sym.Matrix([v_tmp]))
                    rho_V = rho_V.row_join(sym.Matrix([rho_M[ii,jj]]))

        a0_V = 2*I*PO.rec_coef(phR)*a0_V
        a0_V = sym.nsimplify(sym.simplify(a0_V)).as_mutable() # sym.simplify makes a0_V as immutable

        # if self.disp == 1:
        #     ph_s = PO.ph_num2str(phR)
        #     print('phRec: ', ph_s)

        return a0_V, rho_V

    def rec_coef(ph):

        if ph == 'x' or ph == 'X' or (isinstance(ph, int) and ph == 0):
            coef = 1
        elif ph == 'y' or ph == 'Y' or (isinstance(ph, int) and ph == 1):
            coef = -I
        elif ph == '-x' or ph == '-X' or (isinstance(ph, int) and ph == 2):
            coef = -1
        elif ph == '-y' or ph == '-Y' or (isinstance(ph, int) and ph == 3):
            coef = I
        else:
            coef = exp(-I*ph)

        return coef

    def ph_num2str(ph_n):
        if not isinstance(ph_n, int):
            ph_s = str(ph_n)
        else:
            if ph_n == 0:
                ph_s = 'x'
            elif ph_n == 1:
                ph_s = 'y'
            elif ph_n == 2:
                ph_s = '-x'
            elif ph_n == 3:
                ph_s = '-y'

        return ph_s

    def ph2q(ph):
        if ph == 'x' or ph == 'X' or (isinstance(ph, int) and ph == 0):
            q = 0
        elif ph == 'y' or ph == 'Y' or (isinstance(ph, int) and ph == 1):
            q = pi/2
        elif ph == '-x' or ph == '-X' or (isinstance(ph, int) and ph == 2):
            q = pi
        elif ph == '-y' or ph == '-Y' or (isinstance(ph, int) and ph == 3):
            q = 3/2*pi
        else:
            q = ph # Not used in MATLAB version

        return q

    def phmod(phx, ii):
        # ii is 0-based index
        # phx   = [3,2,1,0]
        # ii    = [0,1,2,3,4,5,6,7]
        # phid  = [0,1,2,3,0,1,2,3]
        # phout = [3,2,1,0,3,2,1,0]

        phid = ii % len(phx)
        phout = phx[phid]
        return phout

    def axis2M(axis_v):
        if axis_v == 0: # 1/2 E
            M_tmp = ([[0.5, 0],[0, 0.5]])
        elif axis_v == 1: # x
            M_tmp = ([[0, 0.5],[0.5, 0]])
        elif axis_v == 2: # y
            M_tmp = ([[0, 1/(2*I)],[-1/(2*I), 0]])
        elif axis_v == 3: # z
            M_tmp = ([[0.5, 0],[0, -0.5]])
        elif axis_v == 4: # p
            M_tmp = ([[0, 1],[0, 0]])
        elif axis_v == 5: # m
            M_tmp = ([[0, 0],[1, 0]])
        elif axis_v == 6: # a
            M_tmp = ([[1, 0],[0, 0]])
        elif axis_v == 7: # b
            M_tmp = ([[0, 0],[0, 1]])

        M_tmp = np.matrix(M_tmp)

        return M_tmp

    def M2pol(*args):
        M_in = args[0]
        if len(args) < 2:
            spin_label_cell = PO.spin_label_cell_default
        else:    
            spin_label_cell =  args[1]
        M_in = sym.Matrix(M_in)

        if M_in.shape[0] != M_in.shape[1]:
            sys.exit('M_in should be 2^n x 2^n !')
    
        spin_no = log2(M_in.shape[0])
        if spin_no%1 != 0:
            sys.exit('M_in should be 2^n x 2^n !')
        spin_no = int(spin_no)

        if len(spin_label_cell) < spin_no:
            sys.exit('the size of spin_label_cell must be same as or bigger than spin_no')

        spin_label = spin_label_cell[0:spin_no]
        rho_num_cell = PO.rho_box(spin_no)[1] # Use for axis
        int_ii = 0

        for ii in range(2**spin_no):
            for jj in range(2**spin_no):
                if M_in[ii,jj] != 0:
                    if int_ii == 0:
                        coef_tmp = sym.Matrix([M_in[ii,jj]]) # sym.Matrix([]), not sym.Matrix(value)
                    else:
                        coef_tmp = coef_tmp.col_join(sym.Matrix([M_in[ii,jj]]))
                    
                    for kk in range(spin_no):
                        axis_v_tmp = np.matrix([rho_num_cell[kk][ii][jj]])
                        if kk == 0:
                            axis_vec = axis_v_tmp
                        else:
                            axis_vec = np.concatenate((axis_vec, axis_v_tmp),1)

                    if int_ii == 0:
                        axis_tmp = axis_vec
                    else:
                        axis_tmp = np.concatenate((axis_tmp, axis_vec),0)
                    int_ii += 1

        obj = PO(spin_no, spin_label, axis_tmp, coef_tmp, 'pol')
        return obj

    def M2pmz(*args):
        obj = PO.M2pol(*args)
        obj = PO.pol2pmz(obj)

        return obj

    def M2xyz(*args):
        obj = PO.M2pol(*args)
        obj = PO.pol2xyz(obj)

        return obj

    def rho_box(n):        
        n_s = 2**n
        dec = np.r_[n_s-1:-1:-1] # n_s-1, n_s-2, ..., 1, 0
        bin_mat = PO.de2bi(dec, n)
        rho_cell = sym.zeros(n_s,n_s)
        rho_num_cell = np.zeros((n,n_s,n_s))

        for ii in range(n_s):
            for jj in range(n_s):
                r_vec = bin_mat[jj,:]
                c_vec = bin_mat[ii,:]
                rho_vec = c_vec - r_vec
                rho_tmp = list("".join([x*n for x in 'b'])) 
                # ['b', 'b', 'b'] instead of 'bbb'
                # if rho_tmp = 'bbb', rho_tmp[0] = 'a' doesn't work.
                # https://stackoverflow.com/questions/41752946/replacing-a-character-from-a-certain-index#:~:text=As%20strings%20are%20immutable%20in,value%20at%20the%20desired%20index.&text=You%20can%20quickly%20(and%20obviously,%22slices%22%20of%20the%20original.

                for kk in range(n):
                    rho_num_tmp = 7
                    if r_vec[kk] == 1:
                        rho_tmp[kk] = 'a'
                        rho_num_tmp = 6

                    if rho_vec[kk] == 1:
                        rho_tmp[kk] = 'u'
                        rho_num_tmp = 4
                                              
                    elif rho_vec[kk] == -1:
                        rho_tmp[kk] = 'd'
                        rho_num_tmp = 5

                    # Symbols abm, bbm, pbm, mbm will be displayed as bold a, b, p, m with Latex printer
                    # because bm is a Latx format for bold.
                    # Thus, 'u' and 'd' are used instead of 'p' and 'm', respectively.

                    rho_num_cell[kk][ii][jj] = rho_num_tmp

                rho_tmp = "".join(rho_tmp)
                rho_cell[ii, jj] = symbols(rho_tmp)


        return rho_cell, rho_num_cell



