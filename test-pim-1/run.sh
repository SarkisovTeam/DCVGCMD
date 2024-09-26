#!/bin/bash

# This is the script for performing DCVGCMD using GROMACS
#-------------------------------------------------------------------
# We need to load the necessary executables, in this case GROMACS
# This scrip assumes it runs with GPU and ntmpi, change the commands
# accordingly to fit your own machine. 
#-------------------------------------------------------------------
gmx2023
export GMX_MAXCONSTRWARN=-1
export OMP_NUM_THREADS=8

#-------------------------------------------------------------------
# The default prefix of GROMACS is gmx, please change it accordingly
#-------------------------------------------------------------------
gmx=gmx
#gmx=gmx_mpi
#-------------------------------------------------------------------
a=100   #we have 100 moves, and each simulation have 50 ps 


for (( i = 1 ; i<=$a; i++ ))
do
  echo $i
  j=$(($i-1))
  echo $j
  mkdir $i
  cd $i
  #-------------------------------------------------------------------
  # In the actual simulation, we can generate a random number to choose a
  # configuration for the control volume, but here, we only have 1 config
  # for the gas mixture which is named 1.gro
  #-------------------------------------------------------------------
  #rand=`shuf -i 1-5000 -n 1`
  rand=1

  #-------------------------------------------------------------------
  # Check if a simulation is finished
  # If there is a blown-up configuration
  # Reject the configuration and reconstruct it
  #-------------------------------------------------------------------
  FILE=$'nvt-1bar-308K-mix.gro'
  while true
  do
  if [ -f $FILE ]; then
   echo "File $FILE exists."
    break
  else
   echo "File $FILE does not exist."

  #-------------------------------------------------------------------
  # Inlet control volume exchange
  #-------------------------------------------------------------------
   python ../../src/replace-cv-5s-v5.py -b ../../data/${rand}.gro -i ../$j/membrane.gro -p ../0/mix.top  -t mix.top -s info_inlet -box 8.09715x5.64766x5 -pos 0x0x1 -o inlet.gro -ip ../0/membrane.gro -op posres.gro

  $gmx editconf  -f inlet.gro -o inlet.gro -resnr 1

  #-------------------------------------------------------------------
  # Outlet control volume exchange
  #-------------------------------------------------------------------
  python ../../src/replace-cv-5s-v5.py -b ../0/empty.gro -i inlet.gro -p ../0/mix.top  -t mix.top -s info_outlet -box 8.09715x5.64766x5  -pos 0x0x20 -o outlet.gro -ip ../0/membrane.gro -op posres.gro

  $gmx editconf  -f outlet.gro -o outlet.gro -resnr 1

$gmx make_ndx -f outlet.gro -o index-run.ndx <<EOF
! r PIM
a CO
a OC1
a OC2
a MN1
a NM1
a NM2
q
EOF


  $gmx grompp -f ../0/nvteq-t308.mdp -c  outlet.gro -p mix.top -o nvt-1bar-308K-mix.tpr  -r posres.gro -n index-run.ndx -maxwarn 100
  $gmx mdrun -s nvt-1bar-308K-mix.tpr -v -deffnm nvt-1bar-308K-mix -ntmpi 1  -notunepme

  fi
  done

  #-------------------------------------------------------------------
  # Analysis preparation
  #-------------------------------------------------------------------
  cp nvt-1bar-308K-mix.gro membrane.gro

  #-------------------------------------------------------------------
  # Remove unnecessary files
  #-------------------------------------------------------------------
  rm \#*    
cd ..
done

