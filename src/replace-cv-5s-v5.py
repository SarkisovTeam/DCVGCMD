# This code is written by Tim Yuan
# Created on the 19th Sept, 2023
# Last modified on the 19th June, 2024

# This code take a configuration from the
# bulk simulation gro file, cut a box size
# according to the input volume. 
# This replaces the the control volume
# in the DCVGCMD and generate a new 
# configuration. 

# Lev Sarkisov Research Group

code="replace-cv-5s-v5.py"
modified="19th June, 2024"

import sys
import os.path
import numpy as np
import argparse
import subprocess
import math

np.set_printoptions(threshold=sys.maxsize)


###########################################################################
# Define input options
###########################################################################
def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                            description="Tim Yuan\n\nLev Sarkisov Research Group\n\n"+modified+"\n\n"
                            "This program is used to perform DCVGCMD"
                            )
    parser.add_argument("-i", "--input_i", help="Specify the input file DCVGCMD setup",
                       default="", metavar='')

    parser.add_argument("-b", "--input_b", help="Specify the input bulk configuration",
                       default="", metavar='')

    parser.add_argument("-p", "--input_p", help="Specify the input topology",
                       default="", metavar='')

    parser.add_argument("-o", "--output_o", help="specify the name for output grofile",
                        default="dcvgcmd.gro", metavar='')

    parser.add_argument("-op", "--oposres", help="specify the name for output grofile for posres",
                        default="posres.gro", metavar='')

    parser.add_argument("-ip", "--iposres", help="specify the name for original grofile for posres",
                        default="membrane.gro", metavar='')

    parser.add_argument("-t", "--output_t", help="specify the name for output topology",
                        default="dcvgcmd.top", metavar='')

    parser.add_argument("-s", "--output_s", help="specify the name for output info store file",
                        default="info.store", metavar='')

    parser.add_argument("-box", "--box", help="specify the size of the cv box, 0x0x0 format",
                        default="5x5x5", metavar='')

    parser.add_argument("-pos", "--pos", help="specify the starting position of the cv, 0x0x0 format",
                        default="0x0x1", metavar='')


    args = parser.parse_args()
    return args

args = parseargs()


###########################################################################
# Main function                                                           #
###########################################################################
def main():
  print("\n\nThe command python reads\n\n")
  print(subprocess.list2cmdline(sys.argv[1:]))
  pm.n2num = 3
  pm.co2num = 5
  pm.o2num = 3
  pm.solnum = 4
  pm.gasspecies = ['n2', 'N2', 'co2', 'CO2', 'o2', 'O2', 'SOL']
  fini = args.input_i
  finb = args.input_b
  finp = args.input_p
  fout = args.output_o
  ftop = args.output_t
  fopr = args.oposres
  fipr = args.iposres
  cbox = np.asarray(args.box.split('x'), dtype=np.float64)
  pos  = np.asarray(args.pos.split('x'), dtype=np.float64)

  ###########################################################################
  # Read the bulk and system gro file
  ###########################################################################
  bstore, bstore_true, bbox = readgro(finb)
  sstore, sstore_true, sbox = readgro(fini)

  pm.maxindex =max(len(bstore), len(sstore))

  ###########################################################################
  # GCMC move
  ###########################################################################
  other, co2, n2, o2, sol = gcmd(bstore, bstore_true, bbox, sstore, sstore_true, sbox, cbox, pos)

  totnum = len(other)+len(co2)+len(n2)+len(o2)+len(sol)
  ###########################################################################
  # Write an output 
  ###########################################################################
  print("\n\nWriting files\n\n")
  wout(fout, finp, ftop, other, co2, n2, o2, sol, sbox)
  posres(fipr, fopr, totnum, co2, n2, o2, sol)
  print("\n\nFile writen complete\n\n")



###########################################################################
# read original membrane gro file
###########################################################################
def posres(fipr, fopr, totnum, co2, n2, o2, sol):
  if not os.path.isfile(fipr):
    print("Missing pos res file, check input")
    exit()
  fiprhand = open(fipr, 'r').readlines()
  foprhand = open(fopr, 'w')
  foprhand.write("posres\n")
  foprhand.write(str(totnum)+"\n")

  for line in fiprhand[2:-1]:
    foprhand.write(line)
  for i in range (len(co2)):
    foprhand.write("{:>5.5s}{:5s}{:>5s}{:>5.5s}{:8.3f}{:8.3f}"\
                    "{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".\
            format(co2[i][0], co2[i][1], co2[i][2], co2[i][3],\
                   float(co2[i][4]), float(co2[i][5]), float(co2[i][6]),\
                   float(co2[i][7]), float(co2[i][8]), float(co2[i][9])))

  for i in range (len(n2)):
    foprhand.write("{:>5.5s}{:5s}{:>5s}{:>5.5s}{:8.3f}{:8.3f}"\
                    "{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".\
            format(n2[i][0], n2[i][1], n2[i][2], n2[i][3],\
                   float(n2[i][4]), float(n2[i][5]), float(n2[i][6]),\
                   float(n2[i][7]), float(n2[i][8]), float(n2[i][9])))

  for i in range (len(o2)):
    foprhand.write("{:>5.5s}{:5s}{:>5s}{:>5.5s}{:8.3f}{:8.3f}"\
                    "{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".\
            format(o2[i][0], o2[i][1], o2[i][2], o2[i][3],\
                   float(o2[i][4]), float(o2[i][5]), float(o2[i][6]),\
                   float(o2[i][7]), float(o2[i][8]), float(o2[i][9])))
 
  for i in range (len(sol)):
    foprhand.write("{:>5.5s}{:5s}{:>5s}{:>5.5s}{:8.3f}{:8.3f}"\
                    "{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".\
            format(sol[i][0], sol[i][1], sol[i][2], sol[i][3],\
                   float(sol[i][4]), float(sol[i][5]), float(sol[i][6]),\
                   float(sol[i][7]), float(sol[i][8]), float(sol[i][9])))

  foprhand.write(fiprhand[-1])
  
###########################################################################
# Write an output
###########################################################################
def wout(fout, finp, ftop, other, co2, n2, o2, sol, sbox):
  f_o_hand = open(fout, 'w')
  f_o_hand.write("After GCMC move\n")
  totnum = len(other)+len(co2)+len(n2)+len(o2)+len(sol)
  f_o_hand.write(str(totnum)+"\n")


  for i in range (len(other)):
    f_o_hand.write("{:>5.5s}{:5s}{:>5s}{:>5.5s}{:8.3f}{:8.3f}"\
                    "{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".\
            format(other[i][0], other[i][1], other[i][2], other[i][3],\
                   float(other[i][4]), float(other[i][5]), float(other[i][6]),\
                   float(other[i][7]), float(other[i][8]), float(other[i][9])))

  for i in range (len(co2)):
    f_o_hand.write("{:>5.5s}{:5s}{:>5s}{:>5.5s}{:8.3f}{:8.3f}"\
                    "{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".\
            format(co2[i][0], co2[i][1], co2[i][2], co2[i][3],\
                   float(co2[i][4]), float(co2[i][5]), float(co2[i][6]),\
                   float(co2[i][7]), float(co2[i][8]), float(co2[i][9])))

  for i in range (len(n2)):
    f_o_hand.write("{:>5.5s}{:5s}{:>5s}{:>5.5s}{:8.3f}{:8.3f}"\
                    "{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".\
            format(n2[i][0], n2[i][1], n2[i][2], n2[i][3],\
                   float(n2[i][4]), float(n2[i][5]), float(n2[i][6]),\
                   float(n2[i][7]), float(n2[i][8]), float(n2[i][9])))

  for i in range (len(o2)):
    f_o_hand.write("{:>5.5s}{:5s}{:>5s}{:>5.5s}{:8.3f}{:8.3f}"\
                    "{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".\
            format(o2[i][0], o2[i][1], o2[i][2], o2[i][3],\
                   float(o2[i][4]), float(o2[i][5]), float(o2[i][6]),\
                   float(o2[i][7]), float(o2[i][8]), float(o2[i][9])))

  for i in range (len(sol)):
    f_o_hand.write("{:>5.5s}{:5s}{:>5s}{:>5.5s}{:8.3f}{:8.3f}"\
                    "{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".\
            format(sol[i][0], sol[i][1], sol[i][2], sol[i][3],\
                   float(sol[i][4]), float(sol[i][5]), float(sol[i][6]),\
                   float(sol[i][7]), float(sol[i][8]), float(sol[i][9])))

  f_o_hand.write(str(sbox[0])+"  "+str(sbox[1])+"  "+str(sbox[2])+"\n")
  f_o_hand.close()


  ###########################################################################
  # read and write a new topology file
  ###########################################################################
  newtop = []
  if os.path.isfile(finp):
    handle_fin_s = open(finp, 'r').readlines()
    co2_check = 0
    n2_check = 0
    o2_check = 0
    sol_check = 0
    for i in range(len(handle_fin_s)):
      if '[ molecules ]' in handle_fin_s[i]:
        index_mol = i
    for i in range(index_mol+2):
        newtop.append(handle_fin_s[i])

    for line in handle_fin_s:
      if line.startswith('CO2' or 'co2'):
        co2_tbf = line.split()
        if len(co2_tbf) == 2 and co2_tbf[1].isdigit():
          if len(co2) != 0:
            co2_tbf[1] = len(co2)
            newtop.append(str(co2_tbf[0])+"   "+str(int(round(co2_tbf[1]/pm.co2num)))+"\n")
        co2_check = 1
    if co2_check == 0 and  len(co2) != 0:
      newtop.append("CO2   "+str(int(round(len(co2)/pm.co2num)))+"\n")



    for line in handle_fin_s:
      if line.startswith('N2' or 'n2'):
        n2_tbf = line.split()
        if len(n2_tbf) ==2 and n2_tbf[1].isdigit():
          if len(n2) != 0:
            n2_tbf[1] = len(n2)
            newtop.append(str(n2_tbf[0])+"    "+str(int(round(n2_tbf[1]/pm.n2num)))+"\n")
        else:
          newtop.append(line)
        n2_check = 1
    if n2_check == 0 and  len(n2) != 0:
      newtop.append("N2    "+str(int(round(len(n2)/pm.n2num)))+"\n")

    for line in handle_fin_s:
      if line.startswith('O2' or 'o2'):
        o2_tbf = line.split()
        if len(o2_tbf) == 2 and o2_tbf[1].isdigit():
          if len(o2) != 0:
            o2_tbf[1] = len(o2)
            newtop.append(str(o2_tbf[0])+"   "+str(int(round(o2_tbf[1]/pm.o2num)))+"\n")
        o2_check = 1
    if o2_check == 0 and  len(o2) != 0:
      newtop.append("O2   "+str(int(round(len(o2)/pm.o2num)))+"\n")

    for line in handle_fin_s:
      if line.startswith('SOL' or 'sol'):
        sol_tbf = line.split()
        if len(sol_tbf) == 2 and sol_tbf[1].isdigit():
          if len(sol) != 0:
            sol_tbf[1] = len(sol)
            newtop.append(str(sol_tbf[0])+"   "+str(int(round(sol_tbf[1]/pm.solnum)))+"\n")
        sol_check = 1
    if sol_check == 0 and  len(sol) != 0:
      newtop.append("SOL   "+str(int(round(len(sol)/pm.solnum)))+"\n")


    f_o_hand = open(ftop, 'w')
    for i in range(len(newtop)):
      f_o_hand.write(newtop[i])
    f_o_hand.close()
  else:
    print("No input topology file found")
 


###########################################################################
# Replace the control volume with bulk
###########################################################################
def gcmd(bstore, bstore_true, bbox, sstore, sstore_true, sbox, cbox, pos):
  assert(bbox[0]>=(cbox[0]-0.00001)), "The x dimension of the bulk should exceed the CV"
  assert(bbox[1]>=(cbox[1]-0.00001)), "The y dimension of the bulk should exceed the CV"
  assert(bbox[2]>=(cbox[2]-0.00001)), "The z dimension of the bulk should exceed the CV"

  assert(sbox[0]>=(cbox[0]-0.00001)), "The x dimension of the system should exceed the CV"
  assert(sbox[1]>=(cbox[1]-0.00001)), "The y dimension of the system should exceed the CV"
  assert(sbox[2]>=(cbox[2]-0.00001)), "The z dimension of the system should exceed the CV"

  mlenx = bbox[0] - cbox[0]
  mleny = bbox[1] - cbox[1]
  mlenz = bbox[2] - cbox[2]

  posx = np.random.uniform(0, mlenx, 1)[0]
  posy = np.random.uniform(0, mleny, 1)[0]
  posz = np.random.uniform(0, mlenz, 1)[0]

  ###########################################################################
  # Cut a box and get all the molecules in the bulk
  ###########################################################################
  bulk_molecules = []
  for i in range(len(bstore)):
    if float(bstore[i][4]) >= posx and float(bstore[i][4]) < (posx+cbox[0]):
      if float(bstore[i][5]) >= posy and float(bstore[i][5]) < (posy+cbox[1]):
        if float(bstore[i][6]) >= posz and float(bstore[i][6]) < (posz+cbox[2]):
          bstore[i][4] = float(bstore[i][4]) - posx + pos[0]
          bstore[i][5] = float(bstore[i][5]) - posy + pos[1]
          bstore[i][6] = float(bstore[i][6]) - posz + pos[2]
          bulk_molecules.append(bstore[i])


  ###########################################################################
  # Check the integrety of the molecules
  # Note both the CO2 and N2 has 3 atoms
  ###########################################################################
  #n2------------------------------------------------------------------------
  del_list = []
  int_check = np.zeros(pm.maxindex)

  for i in range (len(bulk_molecules)):
    if bulk_molecules[i][1].strip()  in ['n2', 'N2']:
      index_i = int(bulk_molecules[i][0]) 
      int_check[index_i] += 1

  for i in range(len(bulk_molecules)):
    if bulk_molecules[i][1].strip() in ['n2', 'N2']:
      index_i = int(bulk_molecules[i][0]) 
      if int_check[index_i] != pm.n2num :
        del_list.append(i)
  del_list.sort(reverse=True)
  for i in range(len(del_list)):
    del bulk_molecules[del_list[i]]

 
  #co2-----------------------------------------------------------------------
  del_list = []
  int_check = np.zeros(pm.maxindex)

  for i in range (len(bulk_molecules)):
    if bulk_molecules[i][1].strip() in ['co2','CO2']:
      index_i = int(bulk_molecules[i][0]) 
      int_check[index_i] += 1
  for i in range(len(bulk_molecules)):
    if bulk_molecules[i][1].strip() in ['co2','CO2']:
      index_i = int(bulk_molecules[i][0]) 
      if int_check[index_i] != pm.co2num :
        del_list.append(i)
  del_list.sort(reverse=True)
  for i in range(len(del_list)):
    del bulk_molecules[del_list[i]]


  #o2------------------------------------------------------------------------
  del_list = []
  int_check = np.zeros(pm.maxindex)

  for i in range (len(bulk_molecules)):
    if bulk_molecules[i][1].strip()  in ['o2', 'O2']:
      index_i = int(bulk_molecules[i][0]) 
      int_check[index_i] += 1

  for i in range(len(bulk_molecules)):
    if bulk_molecules[i][1].strip() in ['o2', 'O2']:
      index_i = int(bulk_molecules[i][0]) 
      if int_check[index_i] != pm.o2num :
        del_list.append(i)
  del_list.sort(reverse=True)
  for i in range(len(del_list)):
    del bulk_molecules[del_list[i]]


  #co2-----------------------------------------------------------------------
  del_list = []
  int_check = np.zeros(pm.maxindex)

  for i in range (len(bulk_molecules)):
    if bulk_molecules[i][1].strip()  in ['SOL', 'sol']:
      index_i = int(bulk_molecules[i][0]) 
      int_check[index_i] += 1

  for i in range(len(bulk_molecules)):
    if bulk_molecules[i][1].strip() in ['SOL', 'sol']:
      index_i = int(bulk_molecules[i][0]) 
      if int_check[index_i] != pm.solnum :
        del_list.append(i)
  del_list.sort(reverse=True)
  for i in range(len(del_list)):
    del bulk_molecules[del_list[i]]

 
  ###########################################################################
  # Cut a box and delete all the molecules in the cv
  ###########################################################################
  cv_molecules = []
  cd_molecules = []
  for i in range(len(sstore)):
    if sstore[i][1].strip() in pm.gasspecies:
      if float(sstore[i][4]) < pos[0] or float(sstore[i][4]) >= (pos[0]+cbox[0]) \
      or float(sstore[i][5]) < pos[1] or float(sstore[i][5]) >= (pos[1]+cbox[1]) \
      or float(sstore[i][6]) < pos[2] or float(sstore[i][6]) >= (pos[2]+cbox[2]):
            cv_molecules.append(sstore_true[i])
      else:
          cd_molecules.append(sstore_true[i])
    else:
      cv_molecules.append(sstore_true[i])

  ###########################################################################
  # Check the integrety of the molecules
  # Note both the CO2 and N2 has 3 atoms
  ###########################################################################

  del_list = []
  int_check = np.zeros(pm.maxindex)


  for i in range (len(cv_molecules)):
    if cv_molecules[i][1].strip()  in ['n2', 'N2']:
      index_i = int(cv_molecules[i][0]) 
      int_check[index_i] += 1
  for i in range(len(cv_molecules)):
    if cv_molecules[i][1].strip() in ['n2', 'N2']:
      index_i = int(cv_molecules[i][0]) 
      if (int_check[index_i]%pm.n2num) != 0 :
        del_list.append(i)
  del_list.sort(reverse=True)
  for i in range(len(del_list)):
    cd_molecules.append(cv_molecules[del_list[i]])
    del cv_molecules[del_list[i]]

 
  del_list = []
  int_check = np.zeros(pm.maxindex)

  for i in range (len(cv_molecules)):
    if cv_molecules[i][1].strip() in ['co2','CO2']:
      index_i = int(cv_molecules[i][0]) 
      int_check[index_i] += 1
  for i in range(len(cv_molecules)):
    if cv_molecules[i][1].strip() in ['co2','CO2']:
      index_i = int(cv_molecules[i][0]) 
      if (int_check[index_i]%pm.co2num) != 0:
        del_list.append(i)
  del_list.sort(reverse=True)
  for i in range(len(del_list)):
    cd_molecules.append(cv_molecules[del_list[i]])
    del cv_molecules[del_list[i]]


 
  del_list = []
  int_check = np.zeros(pm.maxindex)

  for i in range (len(cv_molecules)):
    if cv_molecules[i][1].strip() in ['o2','O2']:
      index_i = int(cv_molecules[i][0]) 
      int_check[index_i] += 1
  for i in range(len(cv_molecules)):
    if cv_molecules[i][1].strip() in ['o2','O2']:
      index_i = int(cv_molecules[i][0]) 
      if (int_check[index_i] % pm.o2num)!= 0:
        del_list.append(i)
  del_list.sort(reverse=True)
  for i in range(len(del_list)):
    cd_molecules.append(cv_molecules[del_list[i]])
    del cv_molecules[del_list[i]]


 
  del_list = []
  int_check = np.zeros(pm.maxindex)

  for i in range (len(cv_molecules)):
    if cv_molecules[i][1].strip() in ['SOL','sol']:
      index_i = int(cv_molecules[i][0]) 
      int_check[index_i] += 1
  for i in range(len(cv_molecules)):
    if cv_molecules[i][1].strip() in ['SOL','sol']:
      index_i = int(cv_molecules[i][0]) 
      if (int_check[index_i] % pm.solnum) != 0:
        del_list.append(i)
  del_list.sort(reverse=True)
  for i in range(len(del_list)):
    cd_molecules.append(cv_molecules[del_list[i]])
    del cv_molecules[del_list[i]]




  ###########################################################################
  # Write the atoms info into a file
  ###########################################################################
  f_o_hand = open(args.output_s, 'w')
  f_o_hand.write("Atoms added from the bulk:  "+str(len(bulk_molecules))+ " " + \
    str(posx)+ " " + str(posy)+ " " +str(posz)+ " " +"\n")
  for i in range (len(bulk_molecules)):
    f_o_hand.write("{:>5.5s}{:5s}{:>5s}{:>5.5s}{:8.3f}{:8.3f}"\
                    "{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".\
     format(bulk_molecules[i][0], bulk_molecules[i][1], bulk_molecules[i][2], bulk_molecules[i][3],\
      float(bulk_molecules[i][4]), float(bulk_molecules[i][5]), float(bulk_molecules[i][6]),\
      float(bulk_molecules[i][7]), float(bulk_molecules[i][8]), float(bulk_molecules[i][9])))


  f_o_hand.write("Atoms deleted from the CV:  "+str(len(cd_molecules))+"\n")
  for i in range(len(cd_molecules)):
    f_o_hand.write("{:>5.5s}{:5s}{:>5s}{:>5.5s}{:8.3f}{:8.3f}"\
                    "{:8.3f}{:8.4f}{:8.4f}{:8.4f}\n".\
     format(cd_molecules[i][0], cd_molecules[i][1], cd_molecules[i][2], cd_molecules[i][3],\
      float(cd_molecules[i][4]), float(cd_molecules[i][5]), float(cd_molecules[i][6]),\
      float(cd_molecules[i][7]), float(cd_molecules[i][8]), float(cd_molecules[i][9])))


  f_o_hand.close()

 

   
 
  ###########################################################################
  # Combine the bulk and cv 
  ###########################################################################
  other = []
  co2 = []
  n2 = []
  o2 = []
  sol = []

  for i in range(len(cv_molecules)):
    if   cv_molecules[i][1].strip() in ['n2' , 'N2']:
      n2.append(cv_molecules[i])
    elif cv_molecules[i][1].strip() in ['co2' , 'CO2']:
      co2.append(cv_molecules[i])
    elif cv_molecules[i][1].strip() in ['o2' , 'O2']:
      o2.append(cv_molecules[i])
    elif cv_molecules[i][1].strip() in ['SOL' , 'sol']:
      sol.append(cv_molecules[i])
    else:
      other.append(cv_molecules[i])
  for i in range(len(bulk_molecules)):
    if   bulk_molecules[i][1].strip() in ['n2' , 'N2']:
      n2.append(bulk_molecules[i])
    elif bulk_molecules[i][1].strip() in ['co2' , 'CO2']:
      co2.append(bulk_molecules[i])
    elif bulk_molecules[i][1].strip() in ['o2' , 'O2']:
      o2.append(bulk_molecules[i])
    elif bulk_molecules[i][1].strip() in ['SOL' , 'sol']:
      sol.append(bulk_molecules[i])
 
    else:
      other.append(bulk_molecules[i])

  return other, co2, n2, o2, sol


###########################################################################
# Read the gro file
###########################################################################
def readgro(fin):
    if os.path.isfile(fin):
        handle_fin_s = open(fin, 'r').readlines()
        sitenum = int(handle_fin_s[1])
        site = np.zeros((sitenum, 6 ), dtype = np.float64)
        count = 0
        store = []
        store_true = []
        dimen = handle_fin_s[-1].split()
        box = [float(dimen[0]),float(dimen[1]), float(dimen[2])]
        for i in handle_fin_s[2:len(handle_fin_s)-1]:
            count += 1
            x = [i[0:5].strip(), i[5:10].strip(),i[10:15].strip(),\
            i[15:20].strip(),i[20:28].strip(),i[28:36].strip(),i[36:44].strip(),\
            i[44:52].strip(),i[52:60].strip(),i[60:68].strip()]
            if x[7] == '':
              x[7] = 0
              x[8] = 0
              x[9] = 0
            xtrue = list(x)
            store_true.append(xtrue)
            #warp x y and z back to the box
            x[4] = float(x[4])
            x[5] = float(x[5])
            x[6] = float(x[6])
            x[4] = x[4] - math.floor(x[4]/box[0])*box[0]
            x[5] = x[5] - math.floor(x[5]/box[1])*box[1]
            x[6] = x[6] - math.floor(x[6]/box[2])*box[2]
            store.append(x)
    else:
        print("Missing configuration")
        exit()
    return store, store_true, box

###########################################################################
# Class for the parameters                                                #
###########################################################################
class pm:
    '''This class includes the parameters'''



###########################################################################
# Boilderplate code to call main() function                               #
###########################################################################
if __name__ == '__main__':
    main()

