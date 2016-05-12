from sys import exit
import clusters, clusterdata

def loop_over_all_clusters(cluster_index, start_index, nmols_unique, nmols, cluster_size_max):

 if ( len(cluster_index) < 1 ):
  max_index = nmols_unique
 else:
  max_index = nmols

 for imol in range(start_index, max_index):

  cluster_index.append(imol)
  
  cluster_size = len(cluster_index)
  clusterdata.ntuples[cluster_size-1] += 1

  # consider creating a record for this cluster
  # create a record only if all molecules are connected
  process_cluster(cluster_index,clusterdata.connectivity)

  # can this cluster be a part of a larger structure?
  grow_it_further = potential_subsystem(cluster_index)
  # even if not all molecules are connected in this cluster
  # larger clusters formed later can be connected
  # how do we know if the current cluster is needed for larger clusters?
  # use a single distance criterion

  if ( len(cluster_index)!=cluster_size_max and grow_it_further ):
   loop_over_all_clusters(cluster_index, imol+1, nmols_unique, nmols, cluster_size_max)
  
  del cluster_index[-1]

  
def process_cluster(cluster_index,connectivity):

 # how many molecules are in the cluster
 cluster_size = len(cluster_index)
 #print cluster_index
 if (cluster_size==1):
  print "%05d" % (cluster_index[0])

 # consider creating a record for this cluster
 # create a record only if all molecules are connected!
 all_molecules_are_connected = is_connected(cluster_index,connectivity)

 # create a record if necessary
 if (all_molecules_are_connected):
   size_index=cluster_size-1
   clusterdata.ntuples_kept[size_index] += 1
   clusterdata.farming_file[size_index], clusterdata.ntuples_not_submitted[size_index], clusterdata.ibatch[size_index] = clusters.create_new_cluster_record( cluster_index, clusterdata.array_of_lines, clusterdata.abc_gasphase, clusterdata.indivdir, clusterdata.snapshotdir, clusterdata.cp2kfname, clusterdata.farming_file[size_index], clusterdata.ntuples_not_submitted[size_index], clusterdata.ibatch[size_index] )
    
def is_connected(cluster_index,connectivity):

 # TODO: use the connectivity matrix to determine if the molecules are connected
 cluster_size = len(cluster_index)
 is_connected = False
 
 if (cluster_size==1):
  
   is_connected=True
 
 elif (cluster_size==2):
  
  if (connectivity[cluster_index[0]][cluster_index[1]]==1):
   is_connected=True
 
 elif (cluster_size==3):
  
  print "Not yet implemented"
  sys.exit(123)

 elif (cluster_size==4):

  print "Not yet implemented"
  sys.exit(123)

 else:

  print "Cluster is too large to process"
  sys.exit(123)

 return is_connected


def potential_subsystem(cluster_index):
  
  # can this cluster be a part of a larger structure?
  # even if not all molecules are connected this cluster can
  # be a part of larger clusters
  # how do we know if the current cluster is needed for larger clusters?
  # use a single distance criterion
  
  grow_it = True

  return grow_it

