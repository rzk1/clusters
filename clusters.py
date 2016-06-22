import clusterdata
import sys, os, subprocess
from sys import argv
from math import sqrt, floor, ceil

def get_0cell_coordinates(infile,atoms_per_molecule,abc):

 nmols = 0
 linenumber = 0
 xyzStr = ["","",""]
 xyz = [0.0,0.0,0.0]
 xyz_shift = [0.0,0.0,0.0]
 array_of_lines = []
 # =============================================================
 # read and process atomic coordinates
 # loop over the lines of the coordinate file
 with open(infile) as fp:
  for line in fp:
      
   linenumber += 1
   
   # set up convenient flags 
   it_is_a_new_molecule = ( linenumber%atoms_per_molecule == 1 )
   it_is_the_last_atom = ( linenumber%atoms_per_molecule == 0 )
 
   if it_is_a_new_molecule:
    nmols += 1
 
   # split the line to extract coordinates
   symbolA,xyzStr[0],xyzStr[1],xyzStr[2],molnameA = line.split()
   
   # bring the molecule to the central cell 
   for idim in range(0, 3):
    xyz[idim]=float(xyzStr[idim])
    if (it_is_a_new_molecule):
     xyz_shift[idim] = floor( xyz[idim]/abc[idim] ) * abc[idim]
    xyz[idim] = xyz[idim] - xyz_shift[idim]
   
   # save zero-image coordinates
   array_of_lines.append(xyz[:])

 return array_of_lines

# the order of cells is
#  0  0  0
# -1 -1 -1
# -1 -1  0
# -1 -1  1
# -1  0 -1
# -1  0  0
#.....
def expand_0cell_to_multiple_cells(array_of_lines,ncells,abc):

 xyz = [0.0,0.0,0.0]
 natoms = len(array_of_lines)

 # copy the 0th cell first
 array_of_lines_final = array_of_lines[:]

 #======= multiple the cell in all directions =======
 for kc in range(-ncells[0],ncells[0]+1):
  for lc in range(-ncells[1],ncells[1]+1):
   for mc in range(-ncells[2],ncells[2]+1):
 
    # skip 0,0,0 cell - its already in the array
    if (kc==0 and lc==0 and mc==0 ): continue
 
    for iatom in range (0,natoms):
 
     # modify the coordinates
     xyz[0]=array_of_lines[iatom][0]+kc*abc[0]
     xyz[1]=array_of_lines[iatom][1]+lc*abc[1]
     xyz[2]=array_of_lines[iatom][2]+mc*abc[2]
     
     array_of_lines_final.append(xyz[:])
 
 return array_of_lines_final

def pack_and_submit(snapshotdir,indivdir,cluster_size,ibatch,ncores,wallminutes):

 csdir = "/%03d" % cluster_size
 batchdir = "/%03d" % ibatch
 submitdir = snapshotdir + csdir + batchdir 
 if (clusterdata.doSubmit):
  command = "cd %s; tar --remove-files -zcf packed.tar.gz %s*; ../../../submit_clusters.pl -1 %d 16 sandybridge -1 00 %d farm%d-%03d farming.inp farming.out; cd ../.." % (submitdir,indivdir,ncores,wallminutes,cluster_size,ibatch)
 else: 
  command = "cd %s; tar --remove-files -zcf packed.tar.gz %s*; echo fake-submission -1 %d 16 sandybridge -1 00 %d farm%d-%03d farming.inp farming.out; cd ../.." % (submitdir,indivdir,ncores,wallminutes,cluster_size,ibatch)
 process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
 proc_stdout = process.communicate()[0].strip()
 print proc_stdout
 #print command

 return

def get_abc_from_template(filename):

 command = "grep ABC %s | grep -v \"^\s*#\"" % filename
 process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
 proc_stdout = process.communicate()[0].strip()
 tmpstr,a1,b1,c1 = proc_stdout.split()
 aaa=float(a1)
 bbb=float(b1)
 ccc=float(c1)

 return (aaa,bbb,ccc) 

def get_schedule(nclusters):

 ppn = 16
 minutes_per_job = 15
 optimal_number_of_jobs_per_group = 12
 min_processes_per_group = 8 # must satisfy "ppn / min_processes_per_group is integer"
 min_nodes = 1
 max_nodes = 16

 suggested_ngroups = ceil ( nclusters / optimal_number_of_jobs_per_group )

 # suggested ngroups should request between min-max nodes
 if (suggested_ngroups*min_processes_per_group < min_nodes*ppn ):
  ngroups = min_nodes*ppn // min_processes_per_group
 elif (suggested_ngroups*min_processes_per_group > max_nodes*ppn):
  ngroups = max_nodes*ppn // min_processes_per_group
 else:
  ngroups = suggested_ngroups
 
 # make sure that ncores will be divisible by ppn
 ncores = ceil( ngroups * min_processes_per_group / ppn ) * ppn
 ngroups = ncores // min_processes_per_group
 wallminutes = minutes_per_job * ceil ( nclusters / ngroups )
  
 return (ngroups,ncores,ppn,wallminutes)

def farming_file_start_writing( cluster_size, ibatch, snapshotdir ):
 
 CSdir = "/%03d" % cluster_size
 BatchStr = "/%03d" % ibatch
 
 farmingdir = snapshotdir + CSdir
 
 if not os.path.exists(farmingdir):
  os.makedirs(farmingdir)
 else:
  if (ibatch==1):
   print "Directory %s exists!" % farmingdir
   sys.exit(0)

 # add batch number
 farmingdir = farmingdir + BatchStr

 if not os.path.exists(farmingdir):
  os.makedirs(farmingdir)
 else:
  print "Directory %s exists!" % farmingdir
  sys.exit(0)
 
 farmingfname = farmingdir + "/farming.inp" 
 fileID = open( farmingfname, 'w' )

 str = """
&GLOBAL
  PROJECT farming-job
  PROGRAM FARMING
  RUN_TYPE NONE
&END GLOBAL

&FARMING
"""
 fileID.write(str)
 
 return fileID

def farming_file_finish_writing( fileID, ngroups ):

 str = "\n  NGROUPS %d" % ngroups
 fileID.write(str)
 str = """

&END FARMING

"""
 fileID.write(str)
 fileID.close()
 return

def write_molecule(coords,cluster_center,molA,abc_gasphase,file):
  
 new_coords = [0.0,0.0,0.0]

 for iatom in range(0,clusterdata.atoms_per_molecule):
  symbol="H"
  if (iatom==0): symbol="O"
  for icoord in range(0,3):
   new_coords[icoord] = coords[clusterdata.atoms_per_molecule*molA+iatom][icoord] - cluster_center[icoord] + abc_gasphase[icoord] / 2.0
  line = "%5s%20.10f%20.10f%20.10f%6s\n" % (symbol,new_coords[0],new_coords[1],new_coords[2],"H2O")
  file.write(line)

 return 
  
def get_cluster_center(coords,mols):
  
 cluster_size = len(mols)

 center = [0.0,0.0,0.0]

 for icoord in range(0,3):
  for imol in range(0,cluster_size):
   for iatom in range(0,clusterdata.atoms_per_molecule):
    center[icoord] += coords[clusterdata.atoms_per_molecule*mols[imol]+iatom][icoord]
  center[icoord] /= cluster_size*clusterdata.atoms_per_molecule

 return center 
  
def create_new_cluster_record( mols, coords, abc_gasphase, indivdir, snapshotdir, cp2kfname, farming_file, nclusters, ibatch ):

 # progress indicatior
 #print "-",
 #sys.stdout.flush()

 cluster_size = len(mols)

 # create identification string for cluster
 strID=""
 for imol in range(0, cluster_size):
  strID = strID + "-%05d" % (mols[imol]+1)
  
 strCS = "%03d" % cluster_size
 strBatch = "%03d" % ibatch
 
 # make a new dir for this molecule
 newdirname0 = indivdir + strID
 newdirname = snapshotdir + "/" + strCS + "/" + strBatch + "/" + newdirname0
 if not os.path.exists(newdirname):
  os.makedirs(newdirname)
 else:
  print "Directory %s exists!" % newdirname
  sys.exit(0)

 # create a new coordinate file for this molecules
 newfilename = newdirname + "/" + cp2kfname + ".x"
 new_file = open( newfilename, 'w' )
 
 # get the geometric center of the cluster
 cluster_center = [0.0,0.0,0.0]
 cluster_center = get_cluster_center(coords,mols)
 
 # write cluster coordinates 
 for imol in range(0, cluster_size):
  write_molecule(coords,cluster_center,mols[imol],abc_gasphase,new_file)
 
 new_file.close()

 # add molecule to the farming job input
 newfilename0 = cp2kfname + ".inp"
 farming_lines="\n  &JOB\n    DIRECTORY ./" + newdirname0 + "\n    INPUT_FILE_NAME " + newfilename0 + "\n  &END JOB\n"
 farming_file.write(farming_lines)

 # copy the cp2k input file
 subprocess.call(["cp",cp2kfname+".inp",newdirname])

 # if there are too many clusters added close the current farming file and open a new one
 nclusters += 1 # add the current cluster
 #print "CS = %03d, Ncluster = %05d" % (cluster_size, nclusters)
 if (nclusters >= 150):
  ngroups,ncores,ppn,wallminutes = get_schedule(nclusters)
  farming_file_finish_writing(farming_file,ngroups)
  print "%3d-molecule clusters: #clusters %10d --> groups %10d cores %10d minutes %10d" % (cluster_size,nclusters,ngroups,ncores,wallminutes)
  pack_and_submit(snapshotdir,indivdir,cluster_size, ibatch, ncores,wallminutes)
  # restart the counter and open new farming files
  nclusters = 0
  ibatch += 1 
  farming_file = farming_file_start_writing( cluster_size, ibatch, snapshotdir )

 return farming_file, nclusters, ibatch

def get_OO_distance(array_of_lines,molA,molB):
 
 rAB = [0.0,0.0,0.0]

 for idim in range(0, 3):
  rAB[idim] = abs(array_of_lines[clusterdata.atoms_per_molecule*molA][idim]-array_of_lines[clusterdata.atoms_per_molecule*molB][idim])
  # modify dx, dy, dz so they are the MINIMUM distance
  #rAB[idim] -= abc[idim] * floor(rAB[idim]/abc[idim])
  #if ( rAB[idim] > (0.5 * abc[idim]) ):
  # rAB[idim] = abc[idim] - rAB[idim]
 
 ROO = sqrt ( rAB[0]**2 + rAB[1]**2 + rAB[2]**2 )

 return ROO

def read_connectivity_matrix(file):

 print "Reading the connectivity matrix"

 connectivity = []
 array1D = []
 imol = 0
 # =============================================================
 with open(file) as fp:
  for line in fp:
     
   if ( imol%1000==0 ):
    print "Reading connectivity matrix: %10d" % (imol)
   imol+=1

   line=line.rstrip(os.linesep) 
   #array1D = line.split()
   array1D = list(line)
   array1D = map(int,array1D)
   connectivity.append(array1D[:])

 return connectivity

def write_connectivity_matrix(file,connectivity,nmols):

 with open(file,'w') as fp:
  for imol in range(0, nmols):
   fp.writelines(["%1d"%connectivity[imol][jmol] for jmol in range(0,nmols)])  
   fp.write("\n")

def create_connectivity_matrix(array_of_lines,nmols,Rcutoff):

  # create the connectivity matrix
  #connectivity=[ [ 0 ]*nmols ]*nmols
  connectivity=[ [ 0 for imol in range(nmols) ] for jmol in range(nmols) ]
  
  for imol in range(0, nmols): 
   if ( imol%100==0 ):
    print "Creating connectivity matrix: %10d" % (imol) 
   for jmol in range(imol+1, nmols): 
    Rij = get_OO_distance(array_of_lines,imol,jmol) 


     #DOESN'T SEEM TO PROCESS THIS PROPERLY#
          #________________________________________________#
     #if ( Rij < Rcutoff ):
      # DO NOT WRITE NUMBERS MORE THAN 9 INTO CONNECTIVITY TABLE
      #Double Check what this means, should logically be identical to the previous code of =1#
    if Rij>9:
      connectivity[imol][jmol] = -1
      connectivity[jmol][imol] = -1
    else:
      connectivity[imol][jmol] = Rij
      connectivity[jmol][imol] = Rij

     #DOESN'T SEEM TO PROCESS THIS PROPERLY#
          #________________________________________________#

  return connectivity

def delete_epsilon_files(docluster,snapshotdir,epsilonfile):

 for isize in range(len(docluster)):
  if (docluster[isize]==1):
 
   strCS = "/%03d" % (isize+1)
   outEfile = snapshotdir + strCS + "/" + epsilonfile
   
   if (os.path.isfile(outEfile)):
 
    print "Deleting epsilon: %5s" % strCS
 
    command = "rm %s" % outEfile
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
 
 return

def extract_DFT_energies_to_one_file(docluster,snapshotdir,tempdir,indivdir,energyfile):

 # determine the largest cluster
 max_cluster_size = -1
 isize=0
 for icluster in docluster:
  isize += 1
  if (icluster==1):
   max_cluster_size = isize

 for isize in range(1,max_cluster_size+1): 
 
  if (docluster[isize-1]==1):
 
   strCS = "/%03d" % isize
   outEfile = snapshotdir + strCS + "/" + energyfile
   
   if (not (os.path.isfile(outEfile))):
 
    print "Extracting energies: %3d" % isize
 
    dirname = snapshotdir + strCS
    command = "ls -d %s/??? | wc -l" % dirname
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    ndirs = int(proc_stdout)
 
    for idir in range(1,ndirs+1):
     
     strIDir = "/%03d" % idir
     command = "cd " + snapshotdir + strCS + strIDir + "; mkdir %s; cd %s; tar -zxf ../packed.tar.gz" % (tempdir,tempdir)
     process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
     proc_stdout = process.communicate()[0].strip()
    
     dirname = snapshotdir + strCS + strIDir + "/" + tempdir + "/" + indivdir
     command = "grep 'ENERGY|' %s*/OT.out | awk '{print $1,$10}' >> %s" % (dirname,outEfile)
     process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
     proc_stdout = process.communicate()[0].strip()
  
     command = "cd " + snapshotdir + strCS + strIDir + "; rm -rf %s" % tempdir
     process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
     proc_stdout = process.communicate()[0].strip()
 
 return

def unpack_arc(max_cluster_size,snapshotdir,tempdir):

 for isize in range(1,max_cluster_size+1): 
   
  strCS = "/%03d" % isize
  
  dirname = snapshotdir + strCS
  command = "ls -d %s/??? | wc -l" % dirname
  process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
  proc_stdout = process.communicate()[0].strip()
  ndirs = int(proc_stdout)

  for idir in range(1,ndirs+1):
    
   strIDir = "%03d" % idir
   command = "cd " + snapshotdir + strCS + "; mkdir %s; cd %s; tar -zxf ../%s/packed.tar.gz" % (tempdir,tempdir,strIDir)
   process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
   proc_stdout = process.communicate()[0].strip()
 
 return

def delete_temp_dir(max_cluster_size,snapshotdir,tempdir):

 for isize in range(1,max_cluster_size+1): 
  strCS = "/%03d" % isize
  command = "cd " + snapshotdir + strCS + "; rm -rf %s" % tempdir
  process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
  proc_stdout = process.communicate()[0].strip()
 
 return

def get_atomic_charges(nmols,snapshotdir,tempdir,indivdir):

 qatom=["","",""]
 charges=[]

 for imol in range(0,nmols):
  
  command = "grep -A 5 'Mulliken Population Analysis' %s/001/%s/%s-%05d/OT.out | tail -%d | awk '{print $5}'" % (snapshotdir,tempdir,indivdir,imol+1,clusterdata.atoms_per_molecule)
  #print command 

  process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
  proc_stdout = process.communicate()[0].strip()
  #print proc_stdout
 
  if (clusterdata.atoms_per_molecule!=3):
   print "The next line of code works only for molecules with 3 atoms"
   sys.exit(5) 
  qatom[0],qatom[1],qatom[2] = proc_stdout.split()

  for iatom in range (0,clusterdata.atoms_per_molecule):
   charges.append(float(qatom[iatom]))

 # end loop over molecules
 
 #print charges
 return charges

# get cluster energy or epsilon, depending what file we look in
def get_cluster_e(mols,snapshotdir,indivdir,energyfile):

 cluster_size = len(mols)
 strCS = "%03d" % cluster_size
 energypath = snapshotdir + "/" + strCS + "/" + energyfile

 # we need to find a periodic image of the cluster energy of which is stored in the database
 # we know that database stores clusters in which at least one molecule is in the zero cell
 # search strategy: loop over molecules in the cluster, find a zero-cell index of the current molecule
 #                  shift indices of the other molecules, try to find this shifted-index cluster in the database
 entry_found=False
 for imol in range(cluster_size):

  # stop because the code was not tested for first failed search
  # make sure that shifting molecules works correctly
  if (imol>0):
   exit(555)
 
  centered_mols=center_cluster_by_molecule(imol,mols) 

  # create identification string for cluster
  strID=""
  for kmol in range(cluster_size):
   strID = strID + "-%05d" % (centered_mols[kmol]+1)
  
  #print "Centered cluster: ", centered_mols, ", original cluster: ", mols
 
  # check if the energy of this cluster has been computed
  dirname0 = indivdir + strID
  command = "grep '%s' %s | awk '{print $2}'" % (dirname0,energypath)
  process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
  proc_stdout = process.communicate()[0].strip()
  if (len(proc_stdout)>0):
   entry_found=True
   break

 if (entry_found):
  energy = float(proc_stdout)
 else:
  print "Cluster not found! ", mols
  exit(7)

 return energy

# the order of cells is
#  0  0  0 - central cell first, then it is skipped
# -1 -1 -1
# -1 -1  0
# -1 -1  1
# -1  0 -1
# -1  0  0
#.....
def center_cluster_by_molecule(imol,mols): 

 # cell coordinates of the molecule that will be shifted to the zero cell
 cell_vector = get_cell_vector(mols[imol])
 
 #print "Original cell of MTBC: ", cell_vector

 jcell_vector_new = [0,0,0]
 centered_mols = mols[:]

 for jmol in range(len(mols)):

  jcell_vector = get_cell_vector(mols[jmol])
 
  #print "           cell of mols: ", jcell_vector

  for index in range(3):
   jcell_vector_new[index] = jcell_vector[index] - cell_vector[index]
  
  #print "       New cell of mols: ", jcell_vector_new
 
  jmol_0cell_index = mols[jmol] - (mols[jmol]//clusterdata.nmols_0cell) * clusterdata.nmols_0cell

  #print "       0-cell molecule:  ", jmol_0cell_index

  centered_mols[jmol]=get_mol_index_from_cell_vector(jmol_0cell_index,jcell_vector_new)

  #print "       new molecule ind: ", centered_mols[jmol]

 return centered_mols
 
def get_cell_vector(molecule):

 K = 2*clusterdata.ncells[0]+1
 L = 2*clusterdata.ncells[1]+1
 M = 2*clusterdata.ncells[2]+1

 imol_box = molecule // clusterdata.nmols_0cell

 cell_vector = [0,0,0]
 if (imol_box != 0):

  # take into account that the central box is the first
  if ( imol_box <= K*L*M//2 ):
   box_zero_index = imol_box - 1
  else:
   box_zero_index = imol_box
  
  k_zero_index = box_zero_index // (K*L)
  l_remainder = box_zero_index - K * L * k_zero_index
  l_zero_index = l_remainder // L
  m_zero_index = l_remainder - L * l_zero_index
 
  cell_vector[0] = k_zero_index - clusterdata.ncells[0] 
  cell_vector[1] = l_zero_index - clusterdata.ncells[1] 
  cell_vector[2] = m_zero_index - clusterdata.ncells[2] 

 return cell_vector

# take cell0 index end convert it into the index in the given cell
def get_mol_index_from_cell_vector(cell0_index,cell_vector):

 K = 2*clusterdata.ncells[0]+1
 L = 2*clusterdata.ncells[1]+1
 M = 2*clusterdata.ncells[2]+1

 if (cell_vector[0]==0 and cell_vector[1]==0 and cell_vector[2]==0):

  new_index = cell0_index

 else:

  cell_vector_zero_index=cell_vector[:]
  for idim in range(3):
   cell_vector_zero_index[idim] += clusterdata.ncells[idim]

  box = K * L * cell_vector_zero_index[0] + L * cell_vector_zero_index[1] + cell_vector_zero_index[2]

  if ( box <= K*L*M//2 ):
   box += 1

  new_index = cell0_index + box * clusterdata.nmols_0cell

 return new_index


#cluster-00001-00006 -0.00423091323773548
def write_epsilon_file(mols,epsilons,snapshotdir,indivdir,epsilonfile):
 
 cluster_size = len(mols)
 strCS = "%03d" % cluster_size

 # create identification string for cluster
 strID=""
 for imol in range(0, cluster_size):
  strID = strID + "-%05d" % (mols[imol]+1)
  
 dirname0 = indivdir + strID
 energypath = snapshotdir + "/" + strCS + "/" + epsilonfile
 command = "echo %s %20.10f >> %s" % (dirname0,epsilons[cluster_size-1],energypath)
 process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
 proc_stdout = process.communicate()[0].strip()
 
