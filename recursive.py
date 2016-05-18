from sys import exit
import random
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
 
 #Utilizes DFS in order to get a list of all the connected molecules, if it matches cluster_size, that means that all of them are connected and returns true.
 else:
  visited=[]
  DFS(visited,cluster_index[0],connectivity,cluster_index)
  #Performs DFS at the first index
  # DFS(visited,cluster_index[0],connectivity,cluster_index)

  if(len(visited) == cluster_size):
    is_connected = True
    print (visited)
  else:
    is_connected = False
 return is_connected
 
def BFS(start,connectivity,cluster_index):
  visited = []
  queue = [start]
  while queue:
    check = queue.pop(0)
    if check not in visited:
      visited.append(check)
      for i in cluster_index:
        clone= int(check)
        clone2 = int(i)
        if (i != check & connectivity[clone][clone2]== 1):
          queue.append(i)
  return visited

def binarysearch(visited,search):
  if len(visited) == 0:
    return False
  else:
    midpoint = len(visited)//2
    if visited[midpoint] == search:
      return True
    else:
      if search<visited[midpoint]:
        return binarysearch(visited[:midpoint],search)
      else:
        return binarysearch(visited[midpoint+1:],search)


def rndmqsort(visited): 
    if len(visited)<2: 
      return visited
    pivot_element = random.choice(visited)
    small = [i for i in visited if i< pivot_element]
    medium = [i for i in visited if i==pivot_element]
    large = [i for i in visited if i> pivot_element]
    return rndmqsort(small) + medium + rndmqsort(large)

def DFS(visited,node,connectivity,cluster_index):
  #Add the node being visited to the visited list
  visited.append(node)
  #Quicksorts the list
  rndmqsort(visited)
  #Searches through the cluster_index
  for i in cluster_index:
    #If i doesn't equal the present node, or is a previously visited node.
    if(node != i):
      if(i not in visited):
        #Check if the two are connected, if it is it performs depth first search
        copy = int (i)
        copy2 = int (node)
        if(connectivity[copy][copy2] ==1):
          DFS(visited,i,connectivity,cluster_index)


def potential_subsystem(cluster_index):
  
  # can this cluster be a part of a larger structure?
  # even if not all molecules are connected this cluster can
  # be a part of larger clusters
  # how do we know if the current cluster is needed for larger clusters?
  # use a single distance criterion
  
  grow_it = True

  return grow_it

