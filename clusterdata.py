import clusters

########## IMPORTANT INITIAL SETTINGS ################
indivdir = "cluster"
cp2kfname = "standard"
abc = [35., 35.,35.] #[15.5356853362,15.5356853362,15.5356853362] # periodic box lengths
largest_cluster=5
atoms_per_molecule=3 # how many atoms are in a molecule (not all subrout are generalized)
Rcutoff = 4.0
doSubmit = False
only_first_N_molecules=0 # if more than zero select only the first N molecules from the central cell

########## SHARED DATA ############
# bookkeeping: farming files and counters for different cluster size
farming_file = []
ntuples = []
ntuples_kept = []
ntuples_not_submitted = []
ibatch = []
# directory/files - related data
snapshotdir = "."

array_of_lines = []
connectivity = [] 
abc_gasphase = [0.0,0.0,0.0]

def init_bookkeeping_data():
 
 for icluster in range(0,largest_cluster):
  #print "cluster=", icluster
  ntuples.append(0)
  ntuples_kept.append(0)
  ntuples_not_submitted.append(0)
  ibatch.append(1)
  farming_file.append( clusters.farming_file_start_writing( icluster+1, ibatch[icluster], snapshotdir ) )
  #print farming_file, snapshotdir

#def record_number_of_molecules(nmols_unique,nmols):
# 
# # use 0th element to keep the number of molecules
# ntuples[0] = nmols
# ntuples_kept[0] = nmols_unique

def close_submit_report():

 # close the farming files
 for icluster in range(0,largest_cluster):
  if (ntuples_not_submitted[icluster] > 0):
   ngroups,ncores,ppn,wallminutes = clusters.get_schedule(ntuples_not_submitted[icluster])
   clusters.farming_file_finish_writing(farming_file[icluster],ngroups)
   print "%3d-molecule clusters: #clusters %10d --> groups %10d cores %10d minutes %10d" % (icluster+1,ntuples_not_submitted[icluster],ngroups,ncores,wallminutes)
   clusters.pack_and_submit(snapshotdir,indivdir,icluster+1,ibatch[icluster],ncores,wallminutes)

 print "=========================="
 for icluster in range(0,largest_cluster):
  print "TOTAL %3d-molecule clusters: %10d, submitted in %03d jobs" % (icluster+1,ntuples_kept[icluster],ibatch[icluster])
 print "=========================="

