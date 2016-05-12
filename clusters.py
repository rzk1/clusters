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

 for iatom in range(0,3):
  symbol="H"
  if (iatom==0): symbol="O"
  for icoord in range(0,3):
   new_coords[icoord] = coords[3*molA+iatom][icoord] - cluster_center[icoord] + abc_gasphase[icoord] / 2.0
  line = "%5s%20.10f%20.10f%20.10f%6s\n" % (symbol,new_coords[0],new_coords[1],new_coords[2],"H2O")
  file.write(line)

 return 
  
def get_cluster_center(coords,mols):
  
 cluster_size = len(mols)

 center = [0.0,0.0,0.0]

 for icoord in range(0,3):
  for imol in range(0,cluster_size):
   for iatom in range(0,3):
    center[icoord] += coords[3*mols[imol]+iatom][icoord]
  center[icoord] /= cluster_size*3

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
  rAB[idim] = abs(array_of_lines[3*molA][idim]-array_of_lines[3*molB][idim])
  # modify dx, dy, dz so they are the MINIMUM distance
  #rAB[idim] -= abc[idim] * floor(rAB[idim]/abc[idim])
  #if ( rAB[idim] > (0.5 * abc[idim]) ):
  # rAB[idim] = abc[idim] - rAB[idim]
 
 ROO = sqrt ( rAB[0]**2 + rAB[1]**2 + rAB[2]**2 )

 return ROO

def create_connectivity_matrix(array_of_lines,nmols,Rcutoff):

  # create the connectivity matrix
  connectivity=[ [ 0 for imol in range(nmols) ] for jmol in range(nmols) ]
  
  for imol in range(0, nmols): 
   if ( imol%500==0 ):
    print "Creating connectivity matrix: %10d" % (imol) 
   for jmol in range(imol+1, nmols): 
     Rij = get_OO_distance(array_of_lines,imol,jmol) 
     if ( Rij < Rcutoff ):
      connectivity[imol][jmol] = 1
      connectivity[jmol][imol] = 1

  return connectivity

def extract_DFT_energies_to_one_file(max_cluster_size,snapshotdir,tempdir,indivdir,energyfile):

 for isize in range(1,max_cluster_size+1): 
  
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
  
  command = "grep -A 5 'Mulliken Population Analysis' %s/001/%s/%s-%05d/OT.out | tail -3 | awk '{print $5}'" % (snapshotdir,tempdir,indivdir,imol+1)
  #print command 

  process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
  proc_stdout = process.communicate()[0].strip()
  #print proc_stdout
  
  qatom[0],qatom[1],qatom[2] = proc_stdout.split()

  for iatom in range (0,3):
   charges.append(float(qatom[iatom]))

 # end loop over molecules
 
 #print charges
 return charges

def get_DFT_energy(mols,nmols_unique,monomerE,snapshotdir,tempdir,indivdir,energyfile):

 cluster_size = len(mols)
 strCS = "%03d" % cluster_size

 # create identification string for cluster
 strID=""
 for imol in range(0, cluster_size):
  strID = strID + "-%05d" % (mols[imol]+1)
  
 # check if the energy of this cluster has been computed
 dirname0 = indivdir + strID
 energypath = snapshotdir + "/" + strCS + "/" + energyfile
 command = "grep '%s' %s | awk '{print $2}'" % (dirname0,energypath)
 process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
 proc_stdout = process.communicate()[0].strip()
 #print proc_stdout
 #ENstr = proc_stdout.split()
#return zero energy if the cluster is not in the records
#danger: consider a triple cell0 - cell1 - cell2. cell1-cell2 pair will not have a well defined energy - such pairs are not calculated

 energy = float(proc_stdout)
 # subtract molecular energies one by one
 for imol in range(0, cluster_size):
  mol_0index = mols[imol] - int(floor( mols[imol] / nmols_unique )) * nmols_unique
  energy -= monomerE[mol_0index]

 # subtract the lower order terms from the interaction energy
 # the equation is given in  Molecular Physics, 84:1, 105-114
 # http://dx.doi.org/10.1080/00268979500100071
 for ituple in range(cluster_size-1,1,-1):
  ecorr = 0.0
  # loop over all possible ituples within the current cluster
  # now it is implemented only for ituple==2 (allows to handle triples)
  # for larger values use recursive nested loops
  if (ituple == 2):
   for imol in range(0, cluster_size):
    keep1 = imol+1 - int(floor(( imol+1 ) / 3)) * 3
    keep2 = imol+2 - int(floor(( imol+2 ) / 3)) * 3
    newmols=[mols[keep1],mols[keep2]]
    ecorr += get_DFT_energy(newmols,nmols_unique,monomerE,snapshotdir,tempdir,indivdir,energyfile)
  else:
   print "Cannot handle quadruples or larger clusters"
   sys.exit(1)
  energy -= ecorr

 return energy

def get_molecular_energies(nmols,snapshotdir,tempdir,indivdir,energyfile):

 energies=[]

 for imol in range(1,nmols+1):
  
  strID = "-%05d" % imol 
  dirname0 = indivdir + strID
  energypath = snapshotdir + "/001/" + energyfile
  command = "grep '%s' %s | awk '{print $2}'" % (dirname0,energypath)
  process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
  proc_stdout = process.communicate()[0].strip()
  molecule_energy = float(proc_stdout)
  energies.append(molecule_energy)

 return energies


