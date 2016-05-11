import clusters

def loop_over_all_clusters(cluster_index, start_index, nmols_unique, nmols, cluster_size_max, connectivity):

 if ( len(cluster_index) < 1 ):
  max_index = nmols_unique
 else:
  max_index = nmols

 for imol in range(start_index, max_index):

  cluster_index.append(imol)
  
  # make a record for the cluster if necessary 
  process_cluster(cluster_index,connectivity)
  # consider creating a record for this cluster
  # create a record only if distance criteria are appropriate
  # what distance criteria? all molecules must be connected!

  # can this cluster be a part of a larger structure?
  grow_it_further = potential_subsystem(cluster_index)
  # even if not all molecules are connected in this cluster
  # larger clusters formed later can be connected
  # how do we know if the current cluster is needed for larger clusters?
  # use a single distance criterion

  if ( len(cluster_index)!=cluster_size_max and grow_it_further ):
   loop_over_all_clusters(cluster_index, imol+1, nmols_unique, nmols, cluster_size_max, connectivity)
  
  del cluster_index[-1]

  
def process_cluster(cluster_index,connectivity):

  # consider creating a record for this cluster
  # create a record only if distance criteria are appropriate
  # what distance criteria? all molecules must be connected!

  # TODO: use the connectivity matrix to determine if the molecules are connected
 
  # how many molecules are in the cluster
  nmols_cluster = len(cluster_index)
  print cluster_index
 

def potential_subsystem(cluster_index):
  
  # can this cluster be a part of a larger structure?
  # even if not all molecules are connected this cluster can
  # be a part of larger clusters
  # how do we know if the current cluster is needed for larger clusters?
  # use a single distance criterion
  
  grow_it = True

  return grow_it

#nmols = 4
#cluster_size_max = 3
#cluster_index=[]
#loop_over_all_clusters(cluster_index, 0, nmols, cluster_size_max)

