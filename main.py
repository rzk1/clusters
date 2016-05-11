import sys, os, subprocess
from sys import argv
from math import sqrt, floor, ceil

# import all user-defined functions for working with water files
import clusters
import recursive

# =============================================================
# ------------- the main program ------------------------------
# =============================================================
# set important variables manually
indivdir = "cluster"
cp2kfname = "standard"
abc = [15.5356853362,15.5356853362,15.5356853362] # periodic box lengths
Rcutoff = 4.0
largest_cluster=4
atoms_per_molecule=3 # how many atoms are in a molecule (not all subrout are generalized)

# read command line arguments
script, infile = argv

Rcutoff=float(Rcutoff)
#filename = os.path.basename(infile)
snapshotdir = os.path.dirname(infile)

# read the box size from the cp2k input template
# the size of this box determines the shift in gas_phase calculations
abc_gasphase = [0.0,0.0,0.0]
abc_gasphase[:]=clusters.get_abc_from_template(cp2kfname+".inp")
 
ncells = [0,0,0] # how many periodic cells in each direction we should consider
for idim in range(0,3):
 ncells[idim] = int(ceil( ((largest_cluster-1)*Rcutoff)/abc[idim] ))
 # it is important that Rcutoff is less than half of a cell
 # otherwise spurious convergence efffect will be observed
 if ( Rcutoff > abc[idim]/2.0 ):
  print "Rcutoff is too large for this system. Rcutoff must be smaller than %10.5f" % (abc[idim]/2.0) 
  sys.exit(1)  

# get atomic coordinates from the file (they will be shifted to the 0th cell)
array_of_lines_zero = clusters.get_0cell_coordinates(infile,atoms_per_molecule,abc)

natoms = len(array_of_lines_zero)
if (natoms%atoms_per_molecule != 0): "Something is wrong: natoms = %10d" % natoms
nmols_unique = int(floor(natoms/atoms_per_molecule))

array_of_lines = clusters.expand_0cell_to_multiple_cells(array_of_lines_zero,ncells,abc)

natoms = len(array_of_lines)
if (natoms%atoms_per_molecule != 0): "Something is wrong: natoms = %10d" % natoms
nmols = int(floor(natoms/atoms_per_molecule))

#print "%10d%10d%10d" % (natoms, nmols, nmols_unique)

print "%10d%10d%10d" % (natoms, nmols, nmols_unique)
print "Creating connectivity matrix..."
connectivity = clusters.create_connectivity_matrix(array_of_lines,nmols,Rcutoff)

#TEMP_COMM# bookkeeping
#TEMP_COMM# init farming files and counters for different cluster size
#TEMP_COMMfarming_file = []
#TEMP_COMMntuples = []
#TEMP_COMMntuples_kept = []
#TEMP_COMMntuples_not_submitted = []
#TEMP_COMMibatch = []
#TEMP_COMMdo_this_cluster_size=[0,0,0,0,0,0,0] 
#TEMP_COMMfor icluster in range(0,largest_cluster):
#TEMP_COMM ntuples.append(0)
#TEMP_COMM ntuples_kept.append(0)
#TEMP_COMM ntuples_not_submitted.append(0)
#TEMP_COMM ibatch.append(1)
#TEMP_COMM do_this_cluster_size[icluster] = 1
#TEMP_COMM farming_file.append( clusters.farming_file_start_writing( icluster+1, ibatch[icluster], snapshotdir ) )
#TEMP_COMM #print farming_file
#TEMP_COMM
#TEMP_COMM# use 0th element to keep the number of molecules
#TEMP_COMMntuples[0] = nmols
#TEMP_COMMntuples_kept[0] = nmols_unique
#TEMP_COMM

# =============================================================
# now create clusters
cluster_index=[]
recursive.loop_over_all_clusters(cluster_index, 0, nmols_unique, nmols, largest_cluster, connectivity)

#TEMP_COMM
#TEMP_COMMfor molA in range(0,nmols_unique):
#TEMP_COMM
#TEMP_COMM print "%05d" % (molA+1)
#TEMP_COMM 
#TEMP_COMM molecules = [molA]
#TEMP_COMM farming_file[0], ntuples_not_submitted[0], ibatch[0] = clusters.create_new_cluster_record( molecules, array_of_lines, abc_gasphase, indivdir, snapshotdir, cp2kfname, farming_file[0], ntuples_not_submitted[0], ibatch[0] )
#TEMP_COMM
#TEMP_COMM if (do_this_cluster_size[1]==1):
#TEMP_COMM
#TEMP_COMM  for molB in range(molA+1,nmols):
#TEMP_COMM
#TEMP_COMM   ntuples[1] += 1
#TEMP_COMM
#TEMP_COMM   rAB = clusters.get_OO_distance(array_of_lines,molA,molB) 
#TEMP_COMM   #print "Candidate: %05d - %05d, Distance: %10.5f" % (molA+1,molB+1,rAB)
#TEMP_COMM   
#TEMP_COMM   # if the molecules in the AB-pair are close write their file
#TEMP_COMM   if (rAB < Rcutoff):
#TEMP_COMM    molecules=[molA,molB]
#TEMP_COMM    farming_file[1], ntuples_not_submitted[1], ibatch[1] = clusters.create_new_cluster_record( molecules, array_of_lines, abc_gasphase, indivdir, snapshotdir, cp2kfname, farming_file[1], ntuples_not_submitted[1], ibatch[1] )
#TEMP_COMM    ntuples_kept[1] += 1
#TEMP_COMM
#TEMP_COMM   # do larger clusters only if the current A-B pair has a potential to form
#TEMP_COMM   # the cluster of the appropriate size (the largest is the linear cluster)
#TEMP_COMM   if ( rAB < ((largest_cluster-1.0+0.001)*Rcutoff) and do_this_cluster_size[2]==1):
#TEMP_COMM
#TEMP_COMM    # triples
#TEMP_COMM    for molC in range(molB+1,nmols):
#TEMP_COMM
#TEMP_COMM     ntuples[2] += 1
#TEMP_COMM
#TEMP_COMM     rAC = clusters.get_OO_distance(array_of_lines,molA,molC) 
#TEMP_COMM     rBC = clusters.get_OO_distance(array_of_lines,molB,molC) 
#TEMP_COMM    
#TEMP_COMM     # count the number of short distances in the triple
#TEMP_COMM     nshortd = 0
#TEMP_COMM     if (rAB < Rcutoff):
#TEMP_COMM      nshortd+=1
#TEMP_COMM     if (rAC < Rcutoff):
#TEMP_COMM      nshortd+=1
#TEMP_COMM     if (rBC < Rcutoff):
#TEMP_COMM      nshortd+=1
#TEMP_COMM     # include the triple only if there are at least 2 short distances 
#TEMP_COMM     if (nshortd >= 2):
#TEMP_COMM      molecules=[molA,molB,molC]
#TEMP_COMM      farming_file[2], ntuples_not_submitted[2], ibatch[2] = clusters.create_new_cluster_record( molecules, array_of_lines, abc_gasphase, indivdir, snapshotdir, cp2kfname, farming_file[2], ntuples_not_submitted[2], ibatch[2] )
#TEMP_COMM      ntuples_kept[2] += 1
#TEMP_COMM    
#TEMP_COMM     ## quadruples
#TEMP_COMM     #for molD in range(molC+1,nmols):
#TEMP_COMM
#TEMP_COMM     # ntuples[3] += 1
#TEMP_COMM
#TEMP_COMM     # molecules=[molA,molB,molC,molD]
#TEMP_COMM     # if (connected_cluster(molecules, array_of_lines)):
#TEMP_COMM     #  create_new_cluster_record( molecules, array_of_lines, indivdir, snapshotdir, cp2kfname, farming_file[3] )
#TEMP_COMM     #  ntuples_kept[3] += 1
#TEMP_COMM
#TEMP_COMM     ##end molD loop
#TEMP_COMM
#TEMP_COMM    #end molC loop
#TEMP_COMM   
#TEMP_COMM   #end (Max_Cluster_Size-1)*Rcut criterion
#TEMP_COMM
#TEMP_COMM  #end molB loop 
#TEMP_COMM print ""
#TEMP_COMM
#TEMP_COMM#end molA loop
#TEMP_COMM
#TEMP_COMM############# wrap up ################
#TEMP_COMM
#TEMP_COMM# close the farming files
#TEMP_COMMfor icluster in range(0,largest_cluster):
#TEMP_COMM ngroups,ncores,ppn,wallminutes = clusters.get_schedule(ntuples_not_submitted[icluster])
#TEMP_COMM clusters.farming_file_finish_writing(farming_file[icluster],ngroups)
#TEMP_COMM print "%3d-molecule clusters: %10d --> %10d%10d%10d" % (icluster+1,ntuples_not_submitted[icluster],ngroups,ncores,wallminutes)
#TEMP_COMM clusters.pack_and_submit(snapshotdir,indivdir,icluster+1,ibatch[icluster],ncores,wallminutes)
#TEMP_COMM
#TEMP_COMMprint "=========================="
#TEMP_COMMfor icluster in range(0,largest_cluster):
#TEMP_COMM print "TOTAL %3d-molecule clusters: %10d, submitted in %03d jobs" % (icluster+1,ntuples_kept[icluster],ibatch[icluster])
#TEMP_COMMprint "=========================="
#TEMP_COMM
