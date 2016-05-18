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
 
ncells = [0,0,0] # how many periodic cells in each direction we should consider
for idim in range(0,3):
 ncells[idim] = int(ceil( ((clusterdata.largest_cluster-1)*clusterdata.Rcutoff)/clusterdata.abc[idim] ))
 # it is important that Rcutoff is less than half of a cell
 # otherwise spurious convergence efffect will be observed
 if ( clusterdata.Rcutoff > clusterdata.abc[idim]/2.0 ):
  print "Rcutoff is too large for this system. Rcutoff must be smaller than %10.5f" % (clusterdata.abc[idim]/2.0) 
  sys.exit(1)  

# get atomic coordinates from the file (they will be shifted to the 0th cell)
array_of_lines_zero = clusters.get_0cell_coordinates(infile,clusterdata.atoms_per_molecule,clusterdata.abc)

natoms = len(array_of_lines_zero)
if (natoms%clusterdata.atoms_per_molecule != 0): "Something is wrong: natoms = %10d" % natoms
nmols_unique = int(floor(natoms/clusterdata.atoms_per_molecule))

if (clusterdata.only_first_N_molecules>0):
 nmols_unique=clusterdata.only_first_N_molecules

clusterdata.array_of_lines = clusters.expand_0cell_to_multiple_cells(array_of_lines_zero,ncells,clusterdata.abc)

natoms = len(clusterdata.array_of_lines)
if (natoms%clusterdata.atoms_per_molecule != 0): "Something is wrong: natoms = %10d" % natoms
nmols = int(floor(natoms/clusterdata.atoms_per_molecule))

#print "%10d%10d%10d" % (natoms, nmols, nmols_unique)

print "%10d%10d%10d" % (natoms, nmols, nmols_unique)
clusterdata.connectivity = clusters.create_connectivity_matrix(clusterdata.array_of_lines,nmols,clusterdata.Rcutoff)

# initialize all bookkeeping arrays
clusterdata.init_bookkeeping_data()
# store bookkeeping data
#clusterdata.record_number_of_molecules(nmols_unique,nmols)

# =============================================================
# now create clusters
cluster_index=[]
recursive.loop_over_all_clusters(cluster_index, 0, nmols_unique, nmols, clusterdata.largest_cluster)
# =============================================================

#TEMP_COMM
#TEMP_COMMfor molA in range(0,nmols_unique):
#TEMP_COMM  for molB in range(molA+1,nmols):
#TEMP_COMM   if (rAB < Rcutoff):
#TEMP_COMM   # do larger clusters only if the current A-B pair has a potential to form
#TEMP_COMM   # the cluster of the appropriate size (the largest is the linear cluster)
#TEMP_COMM   if ( rAB < ((largest_cluster-1.0+0.001)*Rcutoff) and do_this_cluster_size[2]==1):
#TEMP_COMM
#TEMP_COMM   #end (Max_Cluster_Size-1)*Rcut criterion
#TEMP_COMM  #end molB loop 
#TEMP_COMM#end molA loop
#TEMP_COMM

############# wrap up ################
clusterdata.close_submit_report()

