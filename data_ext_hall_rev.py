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

file_in = 'cdi1010_hall_300K.txt'
file_out= 'cdi1010_hall_300K_processed.txt'

from numpy import genfromtxt
from math import isnan

data = genfromtxt(file_in,skip_header=1,delimiter=',',encoding='unicode_escape')

perio = int(3)
rh_A_temp  = tuple(data[:,0])
rh_B_temp  = tuple(data[:,1])
field_temp = tuple(data[:,2]/1e4)

ncheck = len(rh_A_temp)
#print('ncheck =',ncheck)
if ncheck % perio != 0:
    print('The data structure is different from extected.')
    exit()

ndata = ncheck/perio

rh_A = ()
rh_B = ()
field = ()

for i in range(int(ndata)):
    j = i*perio
    
    rh_A += (rh_A_temp[j],)
    rh_B += (rh_B_temp[j+1],)
    field += (field_temp[j+2],)


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
    'B(T)\t'	# Moment (emu)
    'R_A(Ohms)\t'	# Temperature (K)
    'R_B(Ohms)\t'	# Temperature (K)
    'R(Ohms)'	# Magnetic Field (Oe)
    )

nlength = len(rh_A)
for i in range(nlength):
    fout.write('\r\n')
    fout.write(str('{:.4f}'.format(field[i])))	# Moment (emu)
    fout.write('\t')
    fout.write(str('{:.5f}'.format(rh_A[i])))# Temperature (K)
    fout.write('\t')
    fout.write(str('{:.5f}'.format(rh_B[i])))	# Magnetic Field (Oe)
    fout.write('\t')
    fout.write(str('{:.5f}'.format(rh_A[i]+rh_B[i])))	# Magnetic Field (Oe)

fout.close()

exit()
