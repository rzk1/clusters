import sys, os, subprocess
from sys import argv
from math import sqrt, floor, ceil

# import all user-defined functions for working with water files
import clusters
import clusterdata
import recursive

# =============================================================
# ------------- the main program ------------------------------
# =============================================================
# read command line arguments
script, infile = argv

clusterdata.Rcutoff=float(clusterdata.Rcutoff)
#filename = os.path.basename(infile)
clusterdata.snapshotdir = os.path.dirname(infile)

# read the box size from the cp2k input template
# the size of this box determines the shift in gas_phase calculations
clusterdata.abc_gasphase[:]=clusters.get_abc_from_template(clusterdata.cp2kfname+".inp")

# determine the largest cluster we have to do
isize=0
for icluster in clusterdata.docluster:
 isize += 1
 if (icluster==1):
  clusterdata.largest_cluster = isize
if (clusterdata.largest_cluster < 1):
 exit(3)

# check if mode_epsilon_file array is the same length as docluster
if (len(clusterdata.docluster) != len(clusterdata.mode_epsilon_file) ):
 print "Lenght of docluster and mode_epsilon_file must be equal"
 exit(3)
 
clusterdata.ncells = [0,0,0] # how many periodic cells in each direction we should consider
for idim in range(0,3):
 clusterdata.ncells[idim] = int(ceil( ((clusterdata.largest_cluster-1)*clusterdata.Rcutoff)/clusterdata.abc[idim] ))
 # it is important that Rcutoff is such that the largest clusters is less than half of a cell
 # otherwise spurious convergence efffect will be observed
 if ( (clusterdata.largest_cluster-1)*clusterdata.Rcutoff > clusterdata.abc[idim]/2.0 ):
  print "Rcutoff is too large for this system. Rcutoff must be smaller than %10.5f" % (clusterdata.abc[idim]/(2.0*(clusterdata.largest_cluster-1))) 
  sys.exit(1)  

# get atomic coordinates from the file (they will be shifted to the 0th cell)
array_of_lines_zero = clusters.get_0cell_coordinates(infile,clusterdata.atoms_per_molecule,clusterdata.abc)

natoms = len(array_of_lines_zero)
if (natoms%clusterdata.atoms_per_molecule != 0): "Something is wrong: natoms = %10d" % natoms
clusterdata.nmols_0cell = int(floor(natoms/clusterdata.atoms_per_molecule))

nmols_central=clusterdata.nmols_0cell
if (clusterdata.only_first_N_molecules>0):
 nmols_central=clusterdata.only_first_N_molecules

clusterdata.array_of_lines = clusters.expand_0cell_to_multiple_cells(array_of_lines_zero,clusterdata.ncells,clusterdata.abc)

natoms = len(clusterdata.array_of_lines)
if (natoms%clusterdata.atoms_per_molecule != 0): "Something is wrong: natoms = %10d" % natoms
nmols = int(floor(natoms/clusterdata.atoms_per_molecule))

print "Molecules: %6d, Zero-cell molecules: %6d, Central molecules: %6d" % (nmols, clusterdata.nmols_0cell, nmols_central)

connectivityFile=infile+".connectivity-"+"%.7f" % clusterdata.Rcutoff
if ( os.path.isfile(connectivityFile) ):
 clusterdata.connectivity = clusters.read_connectivity_matrix(connectivityFile,nmols)
else:
 clusterdata.connectivity = clusters.create_connectivity_matrix(clusterdata.array_of_lines,nmols,clusterdata.Rcutoff)
 clusters.write_connectivity_matrix(connectivityFile,clusterdata.connectivity,nmols)

# initialize all bookkeeping arrays
if (clusterdata.action==clusterdata.action_submit):
 clusters.init_bookkeeping_data()

if (clusterdata.action==clusterdata.action_readenergy):
 clusters.extract_DFT_energies_to_one_file(clusterdata.docluster,clusterdata.snapshotdir,clusterdata.tempdir,clusterdata.indivdir,clusterdata.energyfile)
 clusters.delete_epsilon_files(clusterdata.mode_epsilon_file,clusterdata.snapshotdir,clusterdata.epsilonfile)
 # get charges on atoms (to compute purely electrostatic energy)
 #clusters.unpack_arc(1,clusterdata.snapshotdir,clusterdata.tempdir)
 #clusterdata.array_of_charges = clusters.get_atomic_charges(nmols_unique,clusterdata.snapshotdir,clusterdata.tempdir,clusterdata.indivdir)
 #clusters.delete_temp_dir(1,clusterdata.snapshotdir,clusterdata.tempdir)

# =============================================================
# now go over all clusters
if (clusterdata.action==clusterdata.action_submit): 
 cluster_index=[]
 recursive.loop_over_all_clusters(cluster_index, 0, nmols_central, nmols, clusterdata.largest_cluster)
elif (clusterdata.action==clusterdata.action_readenergy):
 for imax in range(clusterdata.largest_cluster):
  if (clusterdata.docluster[imax]==1):
   cluster_index=[]
   recursive.loop_over_all_clusters(cluster_index, 0, nmols_central, nmols, imax+1)
 print 'Total energy: %20.10f' % clusterdata.total_energy
# =============================================================

############# wrap up ################
if (clusterdata.action==clusterdata.action_submit):
 clusters.close_submit_report()

