from sys import exit
import random
import clusters, clusterdata

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

  # can this cluster be a part of a larger structure?
  grow_it_further = potential_subsystem(cluster_index)
  if ( len(cluster_index)!=cluster_size_max and grow_it_further ):
   loop_over_all_clusters(cluster_index, imol+1, nmols_central, nmols, cluster_size_max)
  
  del cluster_index[-1]

  
def process_cluster(cluster_index,connectivity):

 # consider creating a record for this cluster
 # create a record only if all molecules are connected!
 all_molecules_are_connected = is_connected(cluster_index,connectivity)

 # create a record if necessary
 if (all_molecules_are_connected):

  #SHOULD ONLY BE PRINTING THIS,WHICH IS CONTROLLED BY IS_CONNECTED, WHICH WAS MODIFIED?
   print cluster_index

   # how many molecules are in the cluster
   cluster_size = len(cluster_index)
   size_index=cluster_size-1
   
   if (clusterdata.action==clusterdata.action_submit):

    clusterdata.ntuples_kept[size_index] += 1
    clusterdata.farming_file[size_index], clusterdata.ntuples_not_submitted[size_index], clusterdata.ibatch[size_index] = clusters.create_new_cluster_record( cluster_index, clusterdata.array_of_lines, clusterdata.abc_gasphase, clusterdata.indivdir, clusterdata.snapshotdir, clusterdata.cp2kfname, clusterdata.farming_file[size_index], clusterdata.ntuples_not_submitted[size_index], clusterdata.ibatch[size_index] )
   
   elif (clusterdata.action==clusterdata.action_readenergy):

    print ">>>"
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
    #### use appropriate coefficients
    coef = 1.0
    clusterdata.total_energy += coef*energy[size_index]
    clusters.write_epsilon_file(cluster_index,energy,clusterdata.snapshotdir,clusterdata.indivdir,clusterdata.epsilonfile)
    print "<<<"
    
def is_connected(cluster_index,connectivity):

 # TODO: use the connectivity matrix to determine if the molecules are connected
 cluster_size = len(cluster_index)
 is_connected = False
 
 if (cluster_size==1):
  
   is_connected=True
 
 elif (cluster_size==2):
  
  if (connectivity[cluster_index[0]][cluster_index[1]]<=clusterdata.Rcutoff):
   is_connected=True
 
 #Utilizes DFS in order to get a list of all the connected molecules, if it matches cluster_size, that means that all of them are connected and returns true.
 else:
  if(cluster_size>2 and cluster_size<5):
    is_connected = cycle(cluster_index,connectivity,cluster_size)
  else:
    test = False
    test = connectable(connectivity, cluster_index,1)
    if test == True:
      is_connected = cycle(cluster_index,connectivity,cluster_size)
 return is_connected

def cycle(cluster_index,connectivity,cluster_size):
  visited=[]
  DFS(visited,cluster_index[0],connectivity,cluster_index)
  if(len(visited) == cluster_size):
    return True
  else:
    return False


def connectable(connectivity,cluster_index,num_mols):
    true = 0
    false = 0
    #x,y,z are the molecules being tested
    xandy = connectedThroughNPoint(cluster_index[0],cluster_index[1],clusterdata.Rcutoff,num_mols,connectivity)
    xandz = connectedThroughNPoint(cluster_index[0],cluster_index[2],clusterdata.Rcutoff,num_mols,connectivity)
    yandz = connectedThroughNPoint(cluster_index[1],cluster_index[2],clusterdata.Rcutoff,num_mols,connectivity)
    short_long = 0
    hashtable = {}
    if xandy:
        true = true+1
        hashtable["true"] = "xandy"
    if xandz:
        true = true+1
        hashtable["true"] = "xandz"
    if yandz:
        true = true+1
        hashtable["true"] = "yandz"
    if true == 3:
        return True
    elif true ==0:
        return False
    elif true==2:
        return True
    else:
        Connected1 = 0
        Connected2 = 0
        FurtherAway = 0
        if hashtable["true"] == "xandy":
            Connected1 = cluster_index[0]
            Connected2 = cluster_index[1]
            FurtherAway = cluster_index[2]
        elif hashtable["true"] == "xandz":
            Connected1 = cluster_index[0]
            Connected2 = cluster_index[2]
            FurtherAway = cluster_index[1]
        elif hashtable["true"] == "yandz":
            Connected1 = cluster_index[1]
            Connected2 = cluster_index[2]
            FurtherAway = cluster_index[0]
        if(connectivity[Connected1][Connected2] <= clusterdata.Rcutoff) :
            if(connectedThroughNPoint(FurtherAway,Connected1,clusterdata.Rcutoff,num_mols+1,connectivity) or connectedThroughNPoint(FurtherAway,Connected2,clusterdata.Rcutoff,num_mols+1,connectivity)):
                return True
            else:
                return False
        else:
          # if(connectivity[Connected1][Connected2] <= (2*clusterdata.Rcutoff)):
          #   if(connectedThroughNPoint(FurtherAway,Connected1,clusterdata.Rcutoff,num_mols,connectivity) or connectedThroughNPoint(FurtherAway,Connected2,clusterdata.Rcutoff,num_mols,connectivity)):
          #     return True
          #   else:
          #     return False
          return True

def connectedThroughNPoint(x,y,r,n,connectivity):
    return (connectivity[x][y]<=(n+1)*r)    

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

        #DOESN'T SEEM TO PROCESS THIS PROPERLY#
          #________________________________________________#
        if(connectivity[copy][copy2] <=clusterdata.Rcutoff):

          #________________________________________________#
          DFS(visited,i,connectivity,cluster_index)


def potential_subsystem(cluster_index):
  
  # can this cluster be a part of a larger structure?
  # even if not all molecules are connected this cluster can
  # be a part of larger clusters
  # how do we know if the current cluster is needed for larger clusters?
  # use a single distance criterion
  #TEMP_COMM   if ( rAB < ((largest_cluster-1.0+0.001)*Rcutoff) and do_this_cluster_size[2]==1):
  
  grow_it = True

  return grow_it

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
   print " ",
   print subcluster

  del subcluster[-1]

 # return the accumulated epsilon
 return energy

