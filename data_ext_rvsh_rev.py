##################################################################################
#                                                                                #
#         Data extract from Hall Measurement System output                       #
#                 - The City College of New York -                               #
#                                                                                #
#         Created by Jisoo Moon                                                  #
#         Version 1.0 (Apr. 02, 2025)                                            #
#                                                                                #
##################################################################################
'''
   0: Set Field [G]
   1: Resistivity (Bavg) [ohm cm]
   2: Sheet Resistivity (Bavg) [ohm/sqr]
   3: Resistivity (Gavg) [ohm cm]
   4: Sheet Resistivity (Gavg) [ohm/sqr]
   5: Resistivity A [ohm cm]
   6: Sheet Resistivity A [ohm/sqr]
   7: Resistance A1 [ohm]
   8: Resistance A2 [ohm]
   9: Factor A
  10: Resistivity B [ohm cm]
  11: Sheet Resistivity B [ohm/sqr]
  12: Resistance B1 [ohm]
  13: Resistance B2 [ohm]
  14: Factor B
'''

file_in = 'cdi1010_rvsh_11K.txt'
file_out= 'cdi1010_rvsh_11K_processed.txt'

tolerance = 1e-15
max_loop = int(500)

from numpy import genfromtxt,log,pi,exp
from math import isnan

data = genfromtxt(file_in,skip_header=1,delimiter=',',encoding='unicode_escape')

period = int(7)
Ash_raw = tuple(data[:,0])
A1_raw  = tuple(data[:,1])
A2_raw  = tuple(data[:,2])
Bsh_raw = tuple(data[:,3])
B1_raw  = tuple(data[:,4])
B2_raw  = tuple(data[:,5])
field_raw = tuple(data[:,6]/1e4)

ncheck = len(Ash_raw)
#print('ncheck =',ncheck)
if ncheck % period != 0:
    print('The data structure is different from extected.')
    exit()

ndata = ncheck/period

A1 = ()
A2 = ()
A_avg = ()
Ash_hms = ()
Ash_sqr = ()
Ash_cal = ()
B1 = ()
B2 = ()
B_avg = ()
Bsh_hms = ()
Bsh_sqr = ()
Bsh_cal = ()
field  = ()

for i in range(int(ndata)):
    j = i*period

    # Positive field
    A1_temp = A1_raw[j]
    A2_temp = A2_raw[j+1]
    A_temp = (A1_temp + A2_temp) / float(2)
    
    A1 += (A1_temp,)
    A2 += (A2_temp,)
    A_avg += (A_temp,)

    Ash_hms += (Ash_raw[j+2],)
    Ash_sqr += (A_temp * pi / log(float(2)),)

    ind = int(0)
    delta = 0
    Z_temp = float(2)*log(2)/(pi*A1_temp + pi*A2_temp)
    Y_ = 1 / exp(pi*Z_temp*A1_temp) + 1/exp(pi*Z_temp*A2_temp)
    Z_= Z_temp - ((1-Y_)/pi)/(A1_temp/exp (pi * Z_temp*A1_temp)  + A2_temp/exp(pi*Z_temp*A2_temp))
    while ((abs(Z_-Z_temp) > tolerance) and (ind <= max_loop)):
        ind = ind + int(1)
        Z_temp = Z_
        Y_ = 1 / exp (pi * Z_temp*A1_temp) + 1/exp(pi*Z_temp*A2_temp)
        Z_ = Z_temp - ((1-Y_)/pi)/(A1_temp/exp (pi * Z_temp*A1_temp)  + A2_temp/exp(pi*Z_temp*A2_temp))
        #print('[', ind, ']\t', 'Convergence (sheet conductance difference, /Ohms) = %.5e' %(Z_-Z_temp))
        if ind >= max_loop-int(1):
            print('WARNING: Sheet resistance does not convergence')
    sheet = 1/Z_
    Ash_cal += (sheet,)

    B1_temp = B1_raw[j+3]
    B2_temp = B2_raw[j+4]
    B_temp = (B1_temp + B2_temp) / float(2)
    
    B1 += (B1_temp,)
    B2 += (B2_temp,)
    B_avg += (B_temp,)

    Bsh_hms += (Bsh_raw[j+5],)
    Bsh_sqr += (B_temp * pi / log(float(2)),)

    ind = int(0)
    delta = 0
    Z_temp = float(2)*log(2)/(pi*B1_temp + pi*B2_temp)
    Y_ = 1 / exp(pi*Z_temp*B1_temp) + 1/exp(pi*Z_temp*B2_temp)
    Z_= Z_temp - ((1-Y_)/pi)/(B1_temp/exp (pi * Z_temp*B1_temp)  + B2_temp/exp(pi*Z_temp*B2_temp))
    while ((abs(Z_-Z_temp) > tolerance) and (ind <= max_loop)):
        ind = ind + int(1)
        Z_temp = Z_
        Y_ = 1 / exp (pi * Z_temp*B1_temp) + 1/exp(pi*Z_temp*B2_temp)
        Z_ = Z_temp - ((1-Y_)/pi)/(B1_temp/exp (pi * Z_temp*B1_temp)  + B2_temp/exp(pi*Z_temp*B2_temp))
        #print('[', ind, ']\t', 'Convergence (sheet conductance difference, /Ohms) = %.5e' %(Z_-Z_temp))
        if ind >= max_loop-int(1):
            print('WARNING: Sheet resistance does not convergence')
    sheet = 1/Z_
    Bsh_cal += (sheet,) 

    field += (field_raw[j+6],)

'''
res_gavg     = data[:,3]
shres_gavg   = data[:,4]
resA         = data[:,5]
shresA       = data[:,6]
resistanceA1 = data[:,7]
resistanceA2 = data[:,8]
factorA      = data[:,9]
resB         = data[:,10]
shresB       = data[:,11]
resistanceB1 = data[:,12]
resistanceB2 = data[:,13]
factorB      = data[:,14]
'''

fout = open(file_out, 'w')

# Header
fout.write(
    'B(T)\t'	# in Tesla
    'RA1(Ohms)\t'	# Temperature (K)
    'RA2(Ohms)\t'	# Temperature (K)
    'RxxA(Ohms)\t'	# Magnetic Field (Oe)
    'RshA_sq(O/sq)\t'	# Magnetic Field (Oe)
    'RshA_cal(O/sq)\t'	# Temperature (K)
    'RshA_hms(O/sq)\t'	# Temperature (K)
    'RxxB1(Ohms)\t'	# Temperature (K)
    'RxxB2(Ohms)\t'	# Temperature (K)
    'RxxB(Ohms)\t'	# Magnetic Field (Oe)
    'RshB_sq(O/sq)\t'	# Magnetic Field (Oe)
    'RshB_cal(O/sq)\t'	# Temperature (K)
    'RshB_hms(O/sq)\t'	# Temperature (K)
    'Rsh(Ohms/sq)'	# Magnetic Field (Oe)
    )

nlength = len(field)
for i in range(nlength):
    fout.write('\r\n')
    fout.write(str('{:.3f}'.format(field[i])))	# Moment (emu)
    fout.write('\t')
    fout.write(str('{:.4f}'.format(A1[i])))# Temperature (K)
    fout.write('\t')
    fout.write(str('{:.4f}'.format(A2[i])))	# Magnetic Field (Oe)
    fout.write('\t')
    fout.write(str('{:.4f}'.format(A_avg[i])))	# Magnetic Field (Oe)
    fout.write('\t')
    fout.write(str('{:.4f}'.format(Ash_sqr[i])))	# Magnetic Field (Oe)
    fout.write('\t')
    fout.write(str('{:.4f}'.format(Ash_cal[i])))# Temperature (K)
    fout.write('\t')
    fout.write(str('{:.4f}'.format(Ash_hms[i])))	# Magnetic Field (Oe)
    fout.write('\t')
    fout.write(str('{:.4f}'.format(B1[i])))# Temperature (K)
    fout.write('\t')
    fout.write(str('{:.4f}'.format(B2[i])))	# Magnetic Field (Oe)
    fout.write('\t')
    fout.write(str('{:.4f}'.format(B_avg[i])))	# Magnetic Field (Oe)
    fout.write('\t')
    fout.write(str('{:.4f}'.format(Bsh_sqr[i])))	# Magnetic Field (Oe)
    fout.write('\t')
    fout.write(str('{:.4f}'.format(Bsh_cal[i])))# Temperature (K)
    fout.write('\t')
    fout.write(str('{:.4f}'.format(Bsh_hms[i])))	# Magnetic Field (Oe)
    fout.write('\t')
    fout.write(str('{:.4f}'.format((Ash_cal[i]+Bsh_cal[i])/float(2))))	# Magnetic Field (Oe)

fout.close()

exit()
