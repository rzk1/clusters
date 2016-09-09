import clusterdata, clusters
import math, numpy

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

def collect_analyze_write(cluster_index):

 cluster_size=len(cluster_index) 

 if (cluster_size != 2):
  print "Analysis is implemented only for dimers!"
  return 

 clone_index = intersperse(cluster_index,'-')
 clusterdata.dbfhandle.write("CLUSTER-")
 for i in clone_index:
 	clusterdata.dbfhandle.write(i)
 clusterdata.dbfhandle.write(" ")
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
 for i in range(cluster_size):
 	for x in range(3):
 		Coor1.append(coordinates[3*(i)+x])
 		Coor2.append(coordinates[3*(i+1)+x])	
 		Coor3.append(coordinates[3*(i+2)+x])

 #Get the intramolecular coordinates
 Listofdistances = []
 Listofdistances.append(distance_2_coordinates(Coor1,Coor2))
 Listofdistances.append(distance_2_coordinates(Coor1,Coor3))
 line = "%15f"%(average_between_distances(Listofdistances))
 clusterdata.dbfhandle.write(line)
 line = "%15f"%(distance_2_coordinates(Coor2,Coor3))
 clusterdata.dbfhandle.write(line)
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
	clusterdata.dbfhandle.write(line)


	#Obtaining the angle bisector
	#First H2O molecule
	Xa = list(angle_bisector_coordinates(O1,H1,H2))
	#Second H2O molecule
	Xb = list(angle_bisector_coordinates(O2,H3,H4))
	#Getting cos of points
	#molA is the first O, molB is second O
	#Cos of Xa,O4,O1
	line = "%12f" %(getting_cos(Xa, O1, O2))
	clusterdata.dbfhandle.write(line)
	line = "%12f" %(getting_cos(Xb, O2, O1))
	clusterdata.dbfhandle.write(line) 
	#Getting the dihedral angle
	#Dihedral angle of Xb,O4,O1,Xa
	line ="%12f" % (dihedralAngle(Xb,O2,O1,Xa))
	clusterdata.dbfhandle.write(line)


	#Dihedral angle of Hm-Xa-O1-O4
	dista = distance_2_coordinates(H1,O1)
	distb = distance_2_coordinates(H2,O2)
	if dista>distb:
	 Hm = H1
	else:
	 Hm = H2
	d = (dihedralAngle(Hm,Xa,O1,O2))
	line = "%12f"%(d)
	clusterdata.dbfhandle.write(line)
	#Dihedral angle of Hn-Xb-O4-O1
	dista = distance_2_coordinates(H3,O1)
	distb = distance_2_coordinates(H4,O2)
	if dista>distb:
	 Hn = H3
	else:
	 Hn = H4
	d=(dihedralAngle(Hn,Xb,O2,O1))
	line = "%12f" %(d)
	clusterdata.dbfhandle.write(line)	

 # writing energy into file
 energy = clusters.get_cluster_e(cluster_index,clusterdata.snapshotdir,clusterdata.indivdir,clusterdata.epsilonfile)
 test = "%15f" % energy
 clusterdata.dbfhandle.write(test)
 clusterdata.dbfhandle.write("\n")

