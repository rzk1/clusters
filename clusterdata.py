########## IMPORTANT INITIAL SETTINGS ################
# 0 - submit jobs, 1 - get energies
action = 1 
# Process this cluster size?
# String containing 0's and/or 1's. 0 - do not process, 1 - process
# First digit refers to 1-molecule clusters, second - to 2-molecule and so on 
docluster = [1,1]
#abc = [15.5356853362,15.5356853362,15.5356853362] # 125-molecule periodic box
#abc = [35., 35.,35.] # flat system
abc = [2*15.492205666998032,2*15.492205666998032,2*15.492205666998032] # 1000-molecule box
#Rcutoff = 7.76
Rcutoff = 11.0
#Rcutoff = 3.01
doSubmit = False # for now, relevant only for action=0
only_first_N_molecules = 1 # if more than zero select only the first N molecules from the central cell
# mode_epsilon_file tells what to do with the existing epsilon file: 1 - rewrite, 0 - do not write to it
# elements of the array refere to clusters of certain size: 1,2,3,etc molecules
mode_epsilon_file = [0,0]

atoms_per_molecule=3 # how many atoms are in a molecule (not all subrout are generalized)
indivdir = "cluster"
cp2kfname = "standard"
tempdir = "TEMPOR"
energyfile = "energies.out"
epsilonfile = "epsilon.out"
###################################

########## SHARED DATA ############
# number of molecules in the 0 cell (read from file, shifted to zero cell)
nmols_0cell = 0
# number of neigbor cells (in one dimension, one direction, not counting cell zero)
ncells = [0,0,0]
# bookkeeping: farming files and counters for different cluster size
farming_file = []
ntuples = []
ntuples_kept = []
ntuples_not_submitted = []
ibatch = []
# directory/files - related data
snapshotdir = "."
largest_cluster = 0

array_of_lines = []
array_of_charges = []
connectivity = [] 
abc_gasphase = [0.0,0.0,0.0]

action_submit = 0
action_readenergy = 1

epsilon_rewrite = 1
epsilon_donotwrite = 0

# (approximate) total energy of the system
total_energy = 0.0
##################################

