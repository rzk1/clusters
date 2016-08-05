from sys import exit
import random
import clusters, clusterdata
import math, numpy

def loop_over_all_clusters(cluster_index, start_index, nmols_central, nmols, cluster_size_max):

 if ( len(cluster_index) < 1 ):
  max_index = nmols_central
 else:
  max_index = nmols

 for imol in range(start_index, max_index):

  cluster_index.append(imol)
  
  cluster_size = len(cluster_index)
  if (clusterdata.action==clusterdata.action_submit):
   clusterdata.ntuples[cluster_size-1] += 1

  if (clusterdata.docluster[cluster_size-1]==1):
   if (clusterdata.action==clusterdata.action_submit): 
    process_cluster(cluster_index,clusterdata.connectivity)
   elif (clusterdata.action==clusterdata.action_readenergy):
    if (cluster_size==cluster_size_max):
     process_cluster(cluster_index,clusterdata.connectivity)
   elif (clusterdata.action==clusterdata.action_fileprep):
    if (cluster_size==cluster_size_max):
     process_cluster(cluster_index,clusterdata.connectivity)


  # can this cluster be a part of a larger structure?
  grow_it_further = connectable(cluster_index,clusterdata.connectivity)
  if ( len(cluster_index)!=cluster_size_max and grow_it_further ):
   loop_over_all_clusters(cluster_index, imol+1, nmols_central, nmols, cluster_size_max)
  
  del cluster_index[-1]

  
def process_cluster(cluster_index,connectivity):

 # consider creating a record for this cluster
 # create a record only if all molecules are connected!
 all_molecules_are_connected = is_connected(cluster_index,connectivity)

 # create a record if necessary
 if (all_molecules_are_connected):
   # how many molecules are in the cluster
   cluster_size = len(cluster_index)
   size_index=cluster_size-1
   
   if (clusterdata.action==clusterdata.action_submit):

    clusterdata.ntuples_kept[size_index] += 1
    clusterdata.farming_file[size_index], clusterdata.ntuples_not_submitted[size_index], clusterdata.ibatch[size_index] = clusters.create_new_cluster_record( cluster_index, clusterdata.array_of_lines, clusterdata.abc_gasphase, clusterdata.indivdir, clusterdata.snapshotdir, clusterdata.cp2kfname, clusterdata.farming_file[size_index], clusterdata.ntuples_not_submitted[size_index], clusterdata.ibatch[size_index] )
   
   elif (clusterdata.action==clusterdata.action_readenergy):

    energy=cluster_index[:]
    # loop over all subcluster sizes and accumulate epsilons into energy[itarget]
    # do not include the last cluster because we need to accumulate energy not epsilon
    for itarget in range(size_index):
     subcluster=[]
     energy[itarget]=loop_over_all_subclusters_in_cluster(subcluster, itarget+1, 0, cluster_index)
     #print "itarget = ", itarget, " energy = ", energy[itarget]
    # the last element of the array stores energy of the cluster
    energy[size_index] = clusters.get_cluster_e(cluster_index,clusterdata.snapshotdir,clusterdata.indivdir,clusterdata.energyfile)
    # subtract the lower order terms (epsilons) from the energy
    # the equation is given in  Molecular Physics, 84:1, 105-114
    # http://dx.doi.org/10.1080/00268979500100071
    for itarget in range(size_index):
     energy[size_index]-=energy[itarget]
    #print "Final cluster epsilon: ", energy[size_index]
    #### sum up epsilons to get the total energy
    #### DANGER: use appropriate coefficients
    coef = 1.0
    clusterdata.total_energy += coef*energy[size_index]
    clusters.write_epsilon_file(cluster_index,energy,clusterdata.snapshotdir,clusterdata.indivdir,clusterdata.epsilonfile)

    ####WRITING FILE####
   elif (clusterdata.action==clusterdata.action_fileprep):
    energy=cluster_index[:]
    for itarget in range(size_index):
      # subcluster=[]
      # energy[itarget]=loop_over_all_subclusters_in_cluster(subcluster, itarget+1, 0, cluster_index)]
      energy[size_index] = clusters.get_cluster_e(cluster_index,clusterdata.snapshotdir,clusterdata.indivdir,clusterdata.epsilonfile)
    # for itarget in range(size_index):
    #  energy[size_index]-=energy[itarget]
    # coef = 1.0
    # clusterdata.total_energy += coef*energy[size_index]
    
#Write CLUSTER-0000X-0000Y- parameter1-parameter12-energy
    print cluster_index
    write_file(cluster_index)
    #Writing energy into file
    print energy
    test = "%15f"%((energy[-1]))
    # clusterdata.filename.write("\t Energy: ")
    clusterdata.filename.write(test)
    clusterdata.filename.write("\n")
    
def is_connected(cluster_index,connectivity):

 # Use the connectivity matrix to determine if the molecules are connected
 cluster_size = len(cluster_index)
 is_connected = False
 
 if (cluster_size==1):
  
   is_connected=True
 
 elif (cluster_size==2):
  
  if (connectivity[cluster_index[0]][cluster_index[1]] == 0):
   is_connected=True
 
 # Utilizes DFS in order to get a list of all the connected molecules, 
 # if it matches cluster_size, that means that all of them are connected and returns true.
 else:
  visited=[]
  DFS(visited,cluster_index[0],connectivity,cluster_index)
  if(len(visited) == cluster_size):
    is_connected = True
  else:
    is_connected = False
 
 return is_connected

def connectable(cluster_index,connectivity):

 # Use the connectivity matrix to determine if the molecules are connected
 cluster_size = len(cluster_index)
 is_connectable = True
 
 # cluster of size 1 are connectable
 # check clusters of size 2 and 3 only
 # larger clusters are assumed to be connectable because we have not constructed a good algorithm to check them yet
 if (cluster_size==2):
  
  # check condition: r12 < (largest_cluster-1.0)*Rcutoff
  is_connectable = ( connectivity[cluster_index[0]][cluster_index[1]] < (clusterdata.largest_cluster-1) )
 
 elif (cluster_size==3 and clusterdata.largest_cluster < 6):
 
  xandy = connetableThroughNPoint(cluster_index[0],cluster_index[1],1,connectivity)
  xandz = connetableThroughNPoint(cluster_index[0],cluster_index[2],1,connectivity)
  yandz = connetableThroughNPoint(cluster_index[1],cluster_index[2],1,connectivity)
  nconnections = 0
  if xandy:
   nconnections += 1
   Connected1 = cluster_index[0]
   Connected2 = cluster_index[1]
   FurtherAway = cluster_index[2]
  if xandz:
   nconnections += 1
   Connected1 = cluster_index[0]
   Connected2 = cluster_index[2]
   FurtherAway = cluster_index[1]
  if yandz:
   nconnections += 1
   Connected1 = cluster_index[1]
   Connected2 = cluster_index[2]
   FurtherAway = cluster_index[0]

  if ( nconnections == 3 or nconnections == 2):
   is_connectable = True
  elif ( nconnections == 0 ):
   is_connectable = False
  else:
   if(connectivity[Connected1][Connected2] == 0 ) :
    if(connetableThroughNPoint(FurtherAway,Connected1,clusterdata.largest_cluster-3,connectivity) or connetableThroughNPoint(FurtherAway,Connected2,clusterdata.largest_cluster-3,connectivity)):
     is_connectable = True
    else:
     is_connectable = False
   else:
    is_connectable = True
 
 return is_connectable

# checks whether distance(a,b) < (N+1)*Rcutoff
def connetableThroughNPoint(a,b,N,connectivity):
    return ( connectivity[a][b] < (N+1) )    

def DFS(visited,node,connectivity,cluster_index):
  #Add the node being visited to the visited list
  visited.append(node)
  for i in cluster_index:
    #If i doesn't equal the present node, or is a previously visited node.
    if(node != i):
      if(i not in visited):
        #Check if the two are connected, if it is it performs depth first search
        copy = int (i)
        copy2 = int (node)
        if(connectivity[copy][copy2] == 0):
          DFS(visited,i,connectivity,cluster_index)

# initial call with subcluster=[], start_index=0
def loop_over_all_subclusters_in_cluster(subcluster, target_size, start_index, cluster):

 energy=0.0

 cluster_size = len(cluster)

 for imol in range(start_index, cluster_size):

  subcluster.append(cluster[imol])
  
  subcluster_size = len(subcluster)

  if ( subcluster_size != target_size ):
   # increase subcluster recursively
   loop_over_all_subclusters_in_cluster(subcluster, target_size, imol+1, cluster)
  else:
   # get subcluster's energy, accumulate
   energy += clusters.get_cluster_e(subcluster,clusterdata.snapshotdir,clusterdata.indivdir,clusterdata.epsilonfile)
   # print " ",
   # print subcluster

  del subcluster[-1]

 # return the accumulated epsilon
 return energy

#Uses law of cosine,
def distance_2_coordinates(Coor1,Coor2):
 return math.sqrt((Coor1[0]-Coor2[0])**2+(Coor1[1]-Coor2[1])**2+(Coor1[2]-Coor2[2])**2)

def average_between_distances(Listofdistances):
 Numpy_Array = numpy.asarray(Listofdistances) 
 return (numpy.sum(Numpy_Array)/len(Numpy_Array))

def symmetric_differences_of2_distances(Distanceij,Distanceik):
 return ((Distanceij-Distanceik)**2)

def getting_cos(bisector, molA, molB):
 distanceA=float(distance_2_coordinates(bisector,molA))
 distanceB=float(distance_2_coordinates(molA,molB))
 distanceC=float(distance_2_coordinates(bisector,molB))
 cos = ((distanceB**2)+(distanceA**2)-(distanceC**2))/(2*(distanceA*distanceB))
 return cos

def recursivesearch(Ratio,Coor2,Coor3,OriginCoor2,OriginCoo3):
 Mid = (float(Coor2[0]+Coor3[0])/2,float(Coor2[1]+Coor3[1])/2,float(Coor2[2]+Coor3[2])/2)
 D12 = distance_2_coordinates(OriginCoor2,Mid)
 D23 = distance_2_coordinates(Mid,OriginCoo3)
 CalculatedRatio = float(float(D12)/float(D23))
  
 if (CalculatedRatio-0.0001)<Ratio and (CalculatedRatio+0.00001)>Ratio:
   return Mid
 elif CalculatedRatio<Ratio:
 	return recursivesearch(Ratio,Mid,Coor3,OriginCoor2,OriginCoo3)
 elif CalculatedRatio>Ratio:
 	return recursivesearch(Ratio,Coor2,Mid,OriginCoor2,OriginCoo3)

def angle_bisector_coordinates(Coor1,Coor2,Coor3):
  #Coor1 is the apex
 D12 = distance_2_coordinates(Coor1,Coor2)
 D13 = distance_2_coordinates(Coor1,Coor3)
 Ratio = float(D12/D13)
 CoorBisect = recursivesearch(Ratio,Coor2,Coor3,Coor2,Coor3)
 return CoorBisect

def planeDquation(Coor1,Coor2,Coor3):
 v1 = (Coor3[0]-Coor1[0],Coor3[1]-Coor1[1],Coor3[2]-Coor1[2])
 v2 = (Coor2[0]-Coor1[0],Coor2[1]-Coor1[1],Coor2[2]-Coor1[2])
 crossv1v2 = (v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0])
 z = crossv1v2[0]*Coor1[0]+crossv1v2[1]*Coor1[1]+crossv1v2[2]*Coor1[2]
 # print (crossv1v2[0],crossv1v2[1],crossv1v2[2],z)
 return(crossv1v2[0],crossv1v2[1],crossv1v2[2],z)
    
def dihedralAngle(Coor1,Coor2,Coor3,Coor4):
 Equ1 = planeDquation(Coor1,Coor2,Coor3)
 Equ2 = planeDquation(Coor2,Coor3,Coor4)
 a1 , a2 = Equ1[0],Equ2[0]
 b1 , b2 = Equ1[1],Equ2[1]
 c1 , c2 = Equ1[2],Equ2[2]
 d1 , d2 = Equ1[3],Equ2[3]
 cos = (a1*a2+b1*b2+c1*c2)/(math.sqrt(a1**2+b1**2+c1**2)*math.sqrt(a2**2+b2**2+c2**2))
 angle = math.acos(cos)
 return angle

def intersperse(lst, item):
 result = [item] * (len(lst) * 2 - 1)
 for i in range(len(lst)):
 	result[2*i] = str(lst[i]).zfill(6)
 return result

def write_file(cluster_index):
 clone_index = intersperse(cluster_index,'-')
 clusterdata.filename.write("CLUSTER-")
 for i in clone_index:
 	clusterdata.filename.write(i)
 clusterdata.filename.write(" ")
 coordinates = []
 #Get all of the coordinates into a list
 for i in cluster_index:
 	 for x in range(3):
		 for y in range(3):
     #Adjust formatting here
     #0=Ox,y,z,1=H,2=H
			 coordinates.append((clusterdata.array_of_lines[(clusterdata.atoms_per_molecule*i)+x][y]))
 Coor1 = []
 Coor2 = []
 Coor3 = []
 print cluster_index
 for i in range(len(cluster_index)):
 	for x in range(3):
 		Coor1.append(coordinates[3*(i)+x])
 		Coor2.append(coordinates[3*(i+1)+x])	
 		Coor3.append(coordinates[3*(i+2)+x])

 #Get the intramolecular coordinates
 Listofdistances = []
 Listofdistances.append(distance_2_coordinates(Coor1,Coor2))
 Listofdistances.append(distance_2_coordinates(Coor1,Coor3))
 line = "%15f"%(average_between_distances(Listofdistances))
 clusterdata.filename.write(line)
 line = "%15f"%(distance_2_coordinates(Coor2,Coor3))
 clusterdata.filename.write(line)
 line = "%20f"%(symmetric_differences_of2_distances(Listofdistances[0],Listofdistances[1]))

 #Intermolecular interactions specifically for cluster groups of 2
 if len(cluster_index) == 2:
	O1 = []
 	H1 = []
 	H2 = []
 	O2 = []
 	H3 = []
 	H4 = []
	for x in range(3):
	 O1.append(coordinates[x])
	 H1.append(coordinates[3+x])
	 H2.append(coordinates[6+x])

	 O2.append(coordinates[9+x])
	 H3.append(coordinates[12+x])
	 H4.append(coordinates[15+x])

	#Distance between O molecules
	line = "%20f"%(distance_2_coordinates(O1,O2))
	clusterdata.filename.write(line)


	#Obtaining the angle bisector
	#First H2O molecule
	Xa = list(angle_bisector_coordinates(O1,H1,H2))
	#Second H2O molecule
	Xb = list(angle_bisector_coordinates(O2,H3,H4))
	#Getting cos of points
	#molA is the first O, molB is second O
	#Cos of Xa,O4,O1
	line = "%12f" %(getting_cos(Xa, O1, O2))
	clusterdata.filename.write(line)
	line = "%12f" %(getting_cos(Xb, O2, O1))
	clusterdata.filename.write(line) 
	#Getting the dihedral angle
	#Dihedral angle of Xb,O4,O1,Xa
	line ="%12f" % (dihedralAngle(Xb,O2,O1,Xa))
	clusterdata.filename.write(line)


	#Dihedral angle of Hm-Xa-O1-O4
	dista = distance_2_coordinates(H1,O1)
	distb = distance_2_coordinates(H2,O2)
	if dista>distb:
	 Hm = H1
	else:
	 Hm = H2
	d = (dihedralAngle(Hm,Xa,O1,O2))
	line = "%12f"%(d)
	clusterdata.filename.write(line)
	#Dihedral angle of Hn-Xb-O4-O1
	dista = distance_2_coordinates(H3,O1)
	distb = distance_2_coordinates(H4,O2)
	if dista>distb:
	 Hn = H3
	else:
	 Hn = H4
	d=(dihedralAngle(Hn,Xb,O2,O1))
	line = "%12f" %(d)
	clusterdata.filename.write(line)		
