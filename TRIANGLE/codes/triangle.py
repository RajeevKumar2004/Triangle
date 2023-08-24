import numpy as np
import matplotlib.pyplot as plt
import subprocess
import shlex
import sys                                     
sys.path.insert(0, '/home/peter/TRIANGLE/CoordGeo')
import matplotlib.image as mpimg
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen
import math
y = np.array([[-3, -1], [ 5, -1], [-2, 3]])


A = y[0]
B = y[1]
C = y[2]
print(f"A = {A}")
print(f"B = {B}")
print(f"C = {C}")
print("\n1.1 - Vector\n")


#Question1.1.1
print("Solution-1.1.1")
d = B- A
e = C - B
f = A - C
print("The direction vector of AB is ",d)
print("The direction vector of BC is ",e)
print("The direction vector of CA is ",f)


#Question1.1.2
print("Solution-1.1.2")
V1 = B - C
V2 = V1.reshape(-1,1)
print(f"The length of BC is:{math.sqrt(V1@V2)}")


#Question1.1.3
print("Solution-1.1.3")

np.set_printoptions(precision=2)
Mat = np.array([[1,1,1],[A[0],B[0],C[0]],[A[1],B[1],C[1]]])
rank = np.linalg.matrix_rank(Mat)
if (rank<=2):
	print("Hence proved that points A,B,C in a triangle are collinear")
else:
	print("The given points are not collinear")
def fig1_1_3(A,B,C):
    x_AB = line_gen(A,B)
    x_BC = line_gen(B,C)
    x_CA = line_gen(C,A)
    plt.figure(1)
    #Plotting all lines
    plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
    plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
    plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
    #Labeling the coordinates
    A = A.reshape(-1,1)
    B = B.reshape(-1,1)
    C = C.reshape(-1,1)
    tri_coords = np.block([[A, B, C]])
    plt.scatter(tri_coords[0, :], tri_coords[1, :])
    vert_labels = ['A', 'B', 'C']
    for i, txt in enumerate(vert_labels):
	    offset = 10 if txt == 'C' else -10
	    plt.annotate(txt,
		         (tri_coords[0, i], tri_coords[1, i]),
		         textcoords="offset points",
		         xytext=(0, offset),
		         ha='center')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.legend(loc='best')
    plt.grid() 
    plt.axis('equal')
    plt.savefig('/home/peter/TRIANGLE/figs/fig1.1.3.pdf')

fig1_1_3(A,B,C)
#Question1.1.4
print("Solution1.1.4")
A = y[0]
B = y[1]
C = y[2]
m1=(B-A)
m2=(C-B)
m3=(A-C)
print(f"parametric of AB form is x:{A} + k{m1}")
print(f"parametric of BC form is x:{B} + k{m2}")
print(f"parametric of CA form is x:{C} + k{m3}")


#Question1.1.5
print("Solution1.1.5")
#Orthogonal matrix
omat = np.array([[0,1],[-1,0]])
def dir_vec(C,B):
  return C-B
def norm_vec(C,B):
    return omat@dir_vec(C,B)
n=norm_vec(C,B)
pro=n@B
print(n,"x=",pro)



#Question1.1.6
print("Solution1.1.6")
def AreaCalc(A, B, C):
    AB = A - B
    AC = A - C
#cross_product calculation
    cross_product = np.cross(AB,AC)
#magnitude calculation
    magnitude = np.linalg.norm(cross_product)
    area = 0.5 * magnitude
    return area
area_ABC = AreaCalc(A, B, C)
print("Area of triangle ABC:", area_ABC)


#Question1.1.7
print("Solution1.1.7")
dotA=((B-A).T)@(C-A)
NormA=(np.linalg.norm(B-A))*(np.linalg.norm(C-A))
print('value of angle A: ', np.degrees(np.arccos((dotA)/NormA)))
dotB=(A-B).T@(C-B)
NormB=(np.linalg.norm(A-B))*(np.linalg.norm(C-B))
print('value of angle B: ', np.degrees(np.arccos((dotB)/NormB)))
dotC=(A-C).T@(B-C)
NormC=(np.linalg.norm(A-C))*(np.linalg.norm(B-C))
print('value of angle C: ', np.degrees(np.arccos((dotC)/NormC)))



print("\n1.2 - Median\n")


#Question1.2.1
print("Solution-1.2.1")
D = (B + C)/2
E = (A + C)/2
F = (A + B)/2

print("D:", list(D))
print("E:", list(E))
print("F:", list(F))


#Question1.2.2
print("Solution-1.2.2")
# Orthogonal matrix
omat = np.array([[0, 1], [-1, 0]])
D = (B + C) / 2
E = (C + A) / 2
F = (A + B) / 2
def find_median_equations(A,B,C,D,E,F):
    # Calculate slopes of medians
    slope_median_AB = (C[1] - F[1]) / (C[0] - F[0])
    slope_median_BC = (A[1] - D[1]) / (A[0] - D[0])
    slope_median_CA = (B[1] - E[1]) / (B[0] - E[0])
    
    # Calculate y-intercepts of medians
    intercept_median_AB = F[1] - slope_median_AB * F[0]
    intercept_median_BC = D[1] - slope_median_BC * D[0]
    intercept_median_CA = E[1] - slope_median_CA * E[0]
    
    return (slope_median_AB, intercept_median_AB), (slope_median_BC, intercept_median_BC), (slope_median_CA, intercept_median_CA)
(slope_median_AB, intercept_median_AB), (slope_median_BC, intercept_median_BC), (slope_median_CA, intercept_median_CA)= find_median_equations(A,B,C,D,E,F)
print(f"Equation of Median AD : y = {slope_median_BC}x + {intercept_median_BC}")
print(f"Equation of Median BE : y = {slope_median_CA}x + {intercept_median_CA}")
print(f"Equation of Median CF : y = {slope_median_AB}x + {intercept_median_AB}")   
G=np.array([0,0.33])
def dir_vec(A, B):
    return B - A
def norm_vec(A, B):
    print(omat @ dir_vec(A, B))
    return omat @ dir_vec(A, B) 
def line_gen(A, B):
    len = 10
    dim = A.shape[0]
    x_AB = np.zeros((dim, len))
    lam_1 = np.linspace(0, 1, len)
    for i in range(len):
        temp1 = A + lam_1[i] * (B - A)
        x_AB[:, i] = temp1.T
    return x_AB
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
plt.figure(2)
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')

x_AD = line_gen(A, D)
plt.plot(x_AD[0, :], x_AD[1, :], label='$AD$')

x_BE = line_gen(B, E)
plt.plot(x_BE[0, :], x_BE[1, :], label='$BE$')

x_CF = line_gen(C, F)
plt.plot(x_CF[0, :], x_CF[1, :], label='$CF$')
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
D = D.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)
G = G.reshape(-1,1)
tri_coords = np.block([[A,B,C,D,E,F,G]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D','E','F','G']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(-10,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
plt.savefig('/home/peter/TRIANGLE/figs/fig1.2.2.pdf')


#Question1.2.3
print("Solution-1.2.3")
def dir_vec(A, B):
    return B - A

#for finding the norm of the vector
def norm_vec(A,B):
    N = np.array([0.0,0.0])
    N = np.linalg.norm(B-A)
    return N

#section formula
def line_section(A,B,k):
    P = np.array([0.0,0.0])
    P[0] = ((B[0]*k)+A[0])/(k+1)
    P[1] = ((B[1]*k)+A[1])/(k+1)

    return P

#intersection of lines
def line_intersect(B,E,C,F):
    m1 = dir_vec(B,E)
    m2 = dir_vec(C,F)
    p = np.zeros(2) 
    p[0] = (m2[0]*B[0]-m1[0]*C[0])/(m2[0]-m1[0])
    p[1] = (m2[1]*B[1]-m1[1]*C[1])/(m2[1]-m1[1])
    
    return p
G = line_intersect(B,E,C,F)
print("The vector G is",G)


#Question1.2.4
print("Solution-1.2.4")
A = y[0]
D = (B + C) / 2
G = line_intersect(B,E,C,F)
#finding the direction of vector
BG = dir_vec(B,G)
GE = dir_vec(G,E)
GF = dir_vec(G,F)
CG = dir_vec(C,G)
AG = dir_vec(A,G)
GD = dir_vec(G,D)
#finding the norm of vector
n_BG = norm_vec(B,G)
n_GE = norm_vec(G,E)
n_GF = norm_vec(G,F)
n_CG = norm_vec(C,G)
n_AG = norm_vec(A,G)
n_GD = norm_vec(G,D)
print("The ratio BG/GE is",round(n_BG/n_GE,2))
print("The ratio CG/GF is",round(n_CG/n_GF,2))
print("The ratio AG/GD is",round(n_AG/n_GD,2))

#Question1.2.5
print("Solution-1.2.5")
A = y[0]
B = y[1]
C = y[2]
D = (B + C) / 2
E = (C + A) / 2
F = (A + B) / 2
G = line_intersect(B,E,C,F)
Mat = np.array([[1,1,1],[A[0],D[0],G[0]],[A[1],D[1],G[1]]])
print(Mat)
rank = np.linalg.matrix_rank(Mat)
if (rank==2):
	print("Hence proved that points A,G,D in a triangle are collinear")
else:
	print("Error")
	
	
#Question1.2.6
print("Solution-1.2.6")
G = (A + B + C) / 3
print(f"centroid of the given triangle:{G} ")


#Question1.2.7
print("Solution-1.2.7")
print(f"A - F = {A-F}")
print(f"E - D = {E-D}")

print("Hence verified that A - F = E - D and AFDE is a parallelogram")


print("\n1.3 - Altitude\n")


#Question1.3.1
print("Solution-1.3.1")
bc=C-B
t=np.array([0,1,-1,0]).reshape(2,2)
#AD_1
AD_1=t@bc
#normal vector of AD_1
AD_p=t@AD_1
print(AD_p)


#Question1.3.2
print("Solution-1.3.2")
A = y[0]
B = y[1]
C = y[2]
D = alt_foot(A,B,C)

nt =t@AD_1
result = A.T@nt
print(f"The equation of AD is {nt.T}X={result}")


#Question1.3.3
print("Solution-1.3.3")
A = y[0]
B = y[1]
C = y[2]
def alt_foot(A,B,C):
  m = B-C
  n = omat@m 
  N=np.block([[m],[n]])
  p = np.zeros(2)
  p[0] = m@A 
  p[1] = n@B
  #Intersection
  P=np.linalg.inv(N.T)@p
  return P
D =  alt_foot(A,B,C)
E =  alt_foot(B,C,A)
F =  alt_foot(C,A,B)
def dir_vec(A, B):
    return B - A

def norm_vec(A, B):
    return omat @ dir_vec(A, B)

BE_norm = norm_vec(E, B)
CF_norm = norm_vec(F, C)
print(f"The equation of BE is {BE_norm}[x-B]=0" )
print(f"The equation of CF is {CF_norm}[x-C]=0" )
#Normal vectors of AD and BE
n1 = B-C
n2 = C-A
def line_gen(A,B):
  len =10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(0,1,len)
  for i in range(len):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_AD = line_gen(D,A)
x_BE = line_gen(B,E)
x_CF = line_gen(C,F)

#Plotting all lines
plt.figure(3)
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_AD[0,:],x_AD[1,:],label='$AD$')
plt.plot(x_BE[0,:],x_BE[1,:],label='$BE$')
plt.plot(x_CF[0,:],x_CF[1,:],label='$CF$')



plt.plot(A[0], A[1], 'o')
plt.text(A[0] * (1 + 0.1), A[1] * (1 - 0.1) , 'A')
plt.plot(B[0], B[1], 'o')
plt.text(B[0] * (1 - 0.2), B[1] * (1) , 'B')
plt.plot(C[0], C[1], 'o')
plt.text(C[0] * (1 + 0.03), C[1] * (1 - 0.1) , 'C')
plt.plot(D[0], D[1], 'o')
plt.text(D[0] * (1 + 0.03), D[1] * (1 - 0.1) , 'D')
plt.plot(E[0], E[1], 'o')
plt.text(E[0] * (1 + 0.03), E[1] * (1 - 0.1) , 'E')
plt.plot(F[0], F[1], 'o')
plt.text(F[0] * (1 + 0.03), F[1] * (1 - 0.1) , 'F')

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
plt.savefig("/home/peter/TRIANGLE/figs/fig1.3.3.pdf")

    
#Question1.3.4
print("Solution-1.3.4")      
import numpy as np                       #Importing numpy
import matplotlib.pyplot as plt         
A = np.array([[1,1],[5,-7]])             #Defining the vector A
B = np.array([2,20])                     #Defining the vector B
x = np.linalg.solve(A,B)                 #applying linalg.solve to find x such that Ax=B
print(x)                                 #printing the solution of equation Ax=B 
#Orthogonal matrix
omat = np.array([[0,1],[-1,0]]) 
def dir_vec(A,B):
  return B-A
def norm_vec(A,B):
  return np.matmul(omat, dir_vec(A,B))
#Generate line points
def line_gen(A,B):
  len =10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(0,1,len)
  for i in range(len):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB
#Intersection of two lines
def line_intersect(n1,A1,n2,A2):
  N=np.vstack((n1,n2))
  p = np.zeros(2)
  p[0] = n1@A1
  p[1] = n2@A2
  #Intersection
  P=np.linalg.inv(N)@p
  return P
#Intersection of two lines
def perp_foot(n,cn,P):
  m = omat@n
  N=np.block([[n],[m]])
  p = np.zeros(2)
  p[0] = cn
  p[1] = m@P
  #Intersection
  x_0=np.linalg.inv(N)@p
  return x_0

#Triangle vertices
def tri_vert(a,b,c):
  p = (a**2 + c**2-b**2 )/(2*a)
  q = np.sqrt(c**2-p**2)
  A = np.array([p,q]) 
  B = np.array([0,0]) 
  C = np.array([a,0]) 
  return  A,B,C
#Foot of the Altitude
def alt_foot(A,B,C):
  m = B-C
  n = np.matmul(omat,m) 
  N=np.vstack((m,n))
  p = np.zeros(2)
  p[0] = m@A 
  p[1] = n@B
  #Intersection
  P=np.linalg.inv(N.T)@p
  return P
A = y[0]
B = y[1]
C = y[2]
D = alt_foot(A,B,C)
E = alt_foot(B,A,C)
F = alt_foot(C,A,B)

#Finding orthocentre
H = line_intersect(norm_vec(B,E),E,norm_vec(C,F),F)
x=H



#Question1.3.5
print("Solution-1.3.5") 
A = y[0]
B = y[1]
C = y[2]
D = alt_foot(A,B,C)
H=x
# H is the point of intersection of altitudes on side AB and AC from point C and B respectively...
  # Reference from Question 1.3.4

result = int(((A - H).T) @ (B - C))    # Checking orthogonality condition...

# printing output
if result == 0:
  print("(A - H)^T (B - C) = 0\nHence Verified...")

else:
  print("(A - H)^T (B - C)) != 0\nHence the given statement is wrong...")
    #X is point of intersection of line AH and BC



print("\n1.4 - Perpendicular Bisector\n")


#Question1.4.1
print("Solution-1.4.1") 
A = y[0]
B = y[1]
C = y[2]
def midpoint(P, Q):
    return (P + Q) / 2  
#normal vector 
def norm_vec(A,B):
  omat = np.array([[0,1],[-1,0]]) 
  return omat.T@(A-B)
#to find the coefficients and constant of the equation of perpendicular bisector of BC
def perpendicular_bisector(B, C):
    midBC=midpoint(B,C)
    dir=B-C
    constant = -dir.T @ midBC
    return dir,constant
equation_coeff1,const1 = perpendicular_bisector(A, B)
equation_coeff2,const2 = perpendicular_bisector(B, C)
equation_coeff3,const3 = perpendicular_bisector(C, A)
print(f'Equation for perpendicular bisector of AB: ({equation_coeff1[0]:.2f})x + ({equation_coeff1[1]:.2f})y + ({const1:.2f}) = 0')
print(f'Equation for perpendicular bisector of  BC: ({equation_coeff2[0]:.2f})x + ({equation_coeff2[1]:.2f})y + ({const2:.2f}) = 0')
print(f'Equation for perpendicular bisector of  CA: ({equation_coeff3[0]:.2f})x + ({equation_coeff3[1]:.2f})y + ({const3:.2f}) = 0')




#Question1.4.2
print("Solution-1.4.2") 
A = y[0]
B = y[1]
C = y[2]
def line_intersect(n1,A1,n2,A2):
  N=np.vstack((n1,n2))
  p = np.zeros(2)
  p[0] = n1@A1
  p[1] = n2@A2
  #Intersection
  P=np.linalg.inv(N)@p
  return P
  
def dir_vec(A,B):
  return B-A

def line_gen(A,B):
  len =10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(0,1,len)
  for i in range(len):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB

# direction vector along line joining A & B
AB = dir_vec(A,B)
# direction vector along line joining A & C
AC = dir_vec(A,C)
# midpoint of A & B is F
F = (A+B)/2
# midpoint of A & C is E
E = (A+C)/2
D=(B+C)/2
# O is the point of intersection of perpendicular bisectors of AB and AC
O = line_intersect(AB,F,AC,E)
print(O)
# Circle parameters
center = O
radius = np.linalg.norm(A - O)

# Create circle
theta = np.linspace(0, 2*np.pi, 100)
x_circle = center[0] + radius * np.cos(theta)
y_circle = center[1] + radius * np.sin(theta)
#Generating all lines 
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_OE = line_gen(O,E)
x_OF = line_gen(O,F)
x_OD = line_gen(O,D)
x_OC = line_gen(O,C)
x_OB = line_gen(O,B)


#plotting all lines
plt.figure(4) 
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_OE[0,:],x_OE[1,:],label='$OE$')
plt.plot(x_OF[0,:],x_OF[1,:],label='$OF$')
plt.plot(x_OD[0,:],x_OD[1,:],label='$OG$')
plt.plot(x_OC[0,:],x_OC[1,:],linestyle='dotted',label='$OC$')
plt.plot(x_OB[0,:],x_OB[1,:],linestyle='dotted',label='$OB$')
plt.plot(x_circle, y_circle, label='Circle')

A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
O = O.reshape(-1,1)
D = D.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)
tri_coords = np.block([[A,B,C,O,D,E,F]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','O','D','E','F']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center


plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
plt.savefig("/home/peter/TRIANGLE/figs/fig1.4.2.pdf")

#Question1.4.3
print("Solution-1.4.3") 
A = y[0]
B = y[1]
C = y[2]
def line_intersect(n1,A1,n2,A2):
  N=np.vstack((n1,n2))
  p = np.zeros(2)
  p[0] = n1@A1
  p[1] = n2@A2
  #Intersection
  P=np.linalg.inv(N)@p
  return P
  
def dir_vec(A,B):
  return B-A

def line_gen(A,B):
  len =10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(0,1,len)
  for i in range(len):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB

AB = dir_vec(A,B)
# direction vector along line joining A & C
AC = dir_vec(A,C)
# midpoint of A & B is F
F = (A+B)/2
# midpoint of A & C is E
E = (A+C)/2
# O is the point of intersection of perpendicular bisectors of AB and AC
O = line_intersect(AB,F,AC,E)
print("O = "+str(O)+"\n Hence verified.")

#Question1.4.4
print("Solution-1.4.4") 
import numpy as np
#PPoints of triangle
A = y[0]
B = y[1]
C = y[2]
#Circumcentre
O = line_intersect(AB,F,AC,E)
print(O)
# OA, OB, OC
O_1 = O - A
O_2 = O - B
O_3 = O - C
a = np.linalg.norm(O_1)
b = np.linalg.norm(O_2)
c = np.linalg.norm(O_3)
print(" OA, OB, OC are respectively", a,",", b,",",c, ".")
print("Here, OA = OB = OC.")
print("Hence verified.")


#Question1.4.5
print("Solution-1.4.5") 
print("Figure  Generated")
#Question1.4.6
print("Solution-1.4.6") 
# Given points
A = y[0]
B = y[1]
C = y[2]
O = line_intersect(AB,F,AC,E)

#To find angle BAC
dot_pt_A = (B - A) @ ((C - A).T)
norm_pt_A = np.linalg.norm(B - A) * np.linalg.norm(C - A)
cos_theta_A = dot_pt_A / norm_pt_A
angle_BAC = round(np.degrees(np.arccos(cos_theta_A)),3)  #Round is used to round of number till 5 decimal places
print("angle BAC = " + str(angle_BAC))
#To find angle BOC
dot_pt_O = (B - O) @ ((C - O).T)
norm_pt_O = np.linalg.norm(B - O) * np.linalg.norm(C - O)
cos_theta_O = dot_pt_O / norm_pt_O
angle_BOC = round(360-np.degrees(np.arccos(cos_theta_O)),3)  #Round is used to round of number till 5 decimal places
if angle_BOC!=2*angle_BAC:
    angle_BOC = round(np.degrees(np.arccos(cos_theta_O)),3)
print("angle BOC = " + str(angle_BOC))


#To check whether the answer is correct
if angle_BOC == 2 * angle_BAC:
  print("\nangle BOC = 2 times angle BAC\nHence the give statement is correct")
else:
  print("\nangle BOC ≠ 2 times angle BAC\nHence the given statement is wrong")



#Question1.4.7
print("Solution-1.4.7") 
A = y[0]
B = y[1]
C = y[2]
O = line_intersect(AB,F,AC,E)

# Define the matrix P
def matrix_P(theta):
    return np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])

import numpy as np
for x in np.arange(-1.57,1.57,0.01):
  left_side=C-O
  right_side=matrix_P(x)@(A-O)
  if round(left_side[0],1)==round(right_side[0],1):
    theta_in_degree=round(x*360/2/np.pi,2)
    print(f"Value of theta is {theta_in_degree}")
    break
 

print("\n1.5 - Angle Bisector\n")


#Question1.5.1
print("Solution-1.5.1")    
A = y[0]
B = y[1]
C = y[2]
O = line_intersect(AB,F,AC,E)
omat = np.array([[0,1],[-1,0]]) 

#direction vector
def dir_vec(A,B):
  return B-A
  
#normal vector 
def norm_vec(A,B):
  return omat@dir_vec(A,B)
  
t = norm_vec(B,C) 
n1 = t/np.linalg.norm(t) #unit normal vector
t = norm_vec(C,A)
n2 = t/np.linalg.norm(t)
t = norm_vec(A,B)
n3 = t/np.linalg.norm(t)

#slopes of angle bisectors
m_a=norm_vec(n2,n3)
m_b=norm_vec(n1,n3)
m_c=norm_vec(n1,n2)

#generating line using slope and point
def line_dir_pt(m,A,k1,k2):
  len = 10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(k1,k2,len)
  for i in range(len):
    temp1 = A + lam_1[i]*m
    x_AB[:,i]= temp1.T
  return x_AB

def line_gen(A,B):
  len =10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(0,1,len)
  for i in range(len):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB
  
A = y[0]
B = y[1]
C = y[2]
O = line_intersect(AB,F,AC,E)
def unit_vec(A,B):
	return ((B-A)/np.linalg.norm(B-A))
#using parallelogram theorem
E= unit_vec(A,B) + unit_vec(A,C)
#point generated to create parametric form
#generating normal form
F=np.array([E[1],(E[0]*(-1))])
#matrix multiplication
C1= F@(A.T)
print("Internal Angular bisector of angle A is:",F,"*x = ",C1)


#Question1.5.2
print("Solution-1.5.2")  
A = y[0]
B = y[1]
C = y[2]
omat = np.array([[0,1],[-1,0]]) 

#direction vector
def dir_vec(A,B):
  return B-A
  
#normal vector 
def norm_vec(A,B):
  return omat@dir_vec(A,B)
  
t = norm_vec(B,C) 
n1 = t/np.linalg.norm(t) #unit normal vector
t = norm_vec(C,A)
n2 = t/np.linalg.norm(t)
t = norm_vec(A,B)
n3 = t/np.linalg.norm(t)

#Intersection of two lines
def line_intersect(n1,A1,n2,A2):
  N=np.vstack((n1,n2))
  p = np.zeros(2)
  p[0] = n1@A1
  p[1] = n2@A2
  #Intersection
  P=np.linalg.inv(N)@p
  return P

I=line_intersect(n1-n3,B,n1-n2,C) #intersection of angle bisectors B and C
print(f"I = {I}")
#Generate line points
def line_gen(A,B):
  len =10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(0,1,len)
  for i in range(len):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB



#Question1.5.3
print("Solution-1.5.3")
A = y[0]
B = y[1]
C = y[2]  
#Generating the incircle
[I,r] = icircle(A,B,C)
x_icirc= circ_gen(I,r)

#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_IA = line_gen(I,A)

#Plotting all lines
plt.figure(18)
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_IA[0,:],x_IA[1,:],label='$IA$')

#BA, CA, and IA in vector form
BA = A - B
CA = A - C
IA = A - I

def angle_btw_vectors(v1, v2):
    dot_product = v1 @ v2
    norm = np.linalg.norm(v1) * np.linalg.norm(v2)
    angle = np.arccos(dot_product / norm)
    angle_in_deg = np.degrees(angle)
    return angle_in_deg

#Calculating the angles BAI and CAI
angle_BAI = angle_btw_vectors(BA, IA)
angle_CAI = angle_btw_vectors(CA, IA)

# Print the angles
print("Angle BAI:", angle_BAI)
print("Angle CAI:", angle_CAI)

if np.isclose(angle_BAI, angle_CAI):
    print("Angle BAI is approximately equal to angle CAI.")
else:
    print("error")




#Question1.5.4
print("Solution-1.5.4")
A = y[0]
B = y[1]
C = y[2]  
omat = np.array([[0, 1], [-1, 0]])

def dir_vec(A, B):
    return B - A

def norm_vec(A, B):
    return omat @ dir_vec(A, B)

k1=1
k2=1

p = np.zeros(2)
t = norm_vec(B, C)
n1 = t / np.linalg.norm(t)
t = norm_vec(C, A)
n2 = t / np.linalg.norm(t)
t = norm_vec(A, B)
n3 = t / np.linalg.norm(t)

p[0] = n1 @ B - k1 * n2 @ C
p[1] = n2 @ C - k2 * n3 @ A

N = np.block([[n1 - k1 * n2],[ n2 - k2 * n3]])
I = np.linalg.inv(N)@p
r = abs(n1 @ (B-I))

print("Coordinates of point I:", I)
print(f"Distance from I to BC= {r}")


#Question1.5.5
print("Solution-1.5.5")
I=line_intersect(n1-n3,B,n1-n2,C)   #Here I is incentre (point of intersection of angle bisectors)
n1 = np.array([[7], [5]])   #n1 is normal vector of AB
n2 = np.array([[4], [-4]])  #n2 is normal vector of AC
n1T = n1.T   #taking transpose of n1
n2T = n2.T   #taking transpose of n2
r1 = abs((n1T @ I) - (n1T @ A))/(np.linalg.norm(n1))   #r1 is distance between I and AB (n1T.I - n1T.A=0, n1 and I are vectors)
r2 = abs((n2T @ I) - (n2T @ C))/(np.linalg.norm(n2))   #r2 is distance between I and AC (n2T.I - n2T.C=0, n2 and I are vectors)
print("Distance between I and AB is",r1)
print("Distance between I and AC is",r2)

    
#Question1.5.7
print("Solution-1.5.7")

print("Diagram Generated")


#Question1.5.8
print("Solution-1.5.8")
A = y[0]
B = y[1]
C = y[2]
D=C-B
a=np.linalg.norm(C-B)
b=np.linalg.norm(A-C)
c=np.linalg.norm(A-B)
print("the incentre coordiantes are",I)
p=pow(np.linalg.norm(C-B),2)
q=2*(D@(I-B))
r=pow(np.linalg.norm(I-B),2)-r1*r1

Discre=float(q*q-4*p*r)

print("the Value of discriminant is ",abs(round(Discre,6)))
k=((I-B)@(C-B))/((C-B)@(C-B))
print("the value of parameter k is ",k)
D3=B+k*(C-B)
print("the point of tangency of incircle by side BC is ",D3)
print("Hence we prove that side BC is tangent To incircle and also found the value of k!")


#Question1.5.9
print("Solution-1.5.9")
A = y[0]
B = y[1]
C = y[2]
b=np.linalg.norm(C-A)
c=np.linalg.norm(A-B)
#Orthogonal matrix
omat = np.array([[0,1],[-1,0]]) 
#Triangle vertices
def tri_vert(a,b,c):
  p = (a**2 + c**2-b**2 )/(2*a)
  q = np.sqrt(c**2-p**2)
  A = np.array([p,q]) 
  B = np.array([0,0]) 
  C = np.array([a,0]) 
  return  A,B,C

def dir_vec(A,B):
  return B-A

def norm_vec(A,B):
  return omat@dir_vec(A,B)
  #return np.matmul(omat, dir_vec(A,B))
#finding Incentre 
t = norm_vec(B,C) 
n1 = t/np.linalg.norm(t) #unit normal vector
t = norm_vec(C,A)
n2 = t/np.linalg.norm(t)
t = norm_vec(A,B)
n3 = t/np.linalg.norm(t)

#Intersection of two lines
def line_intersect(n1,A1,n2,A2):
  N=np.block([[n1],[n2]])
  p = np.zeros(2)
  p[0] = n1@A1
  p[1] = n2@A2
  #Intersection
  P=np.linalg.inv(N)@p
  return P
D_3=D3
I=line_intersect(n1-n3,B,n1-n2,C)  #Incentre
#Incentre
print("I = ",I)
#finding k for E_3 and F_3
k1=((I-A)@(A-B))/((A-B)@(A-B))
k2=((I-A)@(A-C))/((A-C)@(A-C))
#finding E_3 and F_3
E_3=A+(k1*(A-B))
F_3=A+(k2*(A-C))
print("k1 = ",k1)
print("k2 = ",k2)
print("E3 = ",E_3)
print("F3 = ",F_3)
#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_BI = line_gen(B,I)
x_CI = line_gen(C,I)
x_IA = line_gen(I,A)
#Generating the incircle
[I,r] = icircle(A,B,C)
x_icirc= circ_gen(I,r)
#Plotting all lines
plt.figure(5)
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
# plt.plot(x_BD[0,:],x_BD[1,:],label='$BD$')
#Plotting the incircle
plt.plot(x_icirc[0,:],x_icirc[1,:],label='$incircle$')
plt.plot(D_3[0],D_3[1],label='$D_3$')
plt.plot(E_3[0],E_3[1],label='$E_3$')
plt.plot(F_3[0],F_3[1],label='$F_3$')
plt.plot(x_CI[0,:],x_CI[1,:],label='$IC$')
plt.plot(x_BI[0,:],x_BI[1,:],label='$IB$')
plt.plot(x_IA[0,:],x_IA[1,:],label='$IA$')
A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
I = I.reshape(-1,1)
D_3=D_3.reshape(-1,1)
E_3=E_3.reshape(-1,1)
F_3=F_3.reshape(-1,1)
#Labeling the coordinates
tri_coords = np.block([[A,B,C,D_3,E_3,F_3,I]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D_3','E_3','F_3','I']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
plt.savefig("/home/peter/TRIANGLE/figs/fig1.5.9.pdf")

#Question1.5.10
print("Solution-1.5.10")

def norm(X,Y):
    magnitude=round(float(np.linalg.norm([X-Y])),3)
    return magnitude 
print("AE_3=", norm(A,E_3) ,"\nAF_3=", norm(A,F_3) ,"\nBD_3=", norm(B,D3) ,"\nBF_3=", norm(B,F_3) ,"\nCD_3=", norm(C,D3) ,"\nCE_3=",norm(C,E_3))


#Question1.5.11
print("Solution-1.5.11")

#finding sidelengths a, b & c
a = np.linalg.norm(B-C)
b = np.linalg.norm(C-A)
c = np.linalg.norm(A-B)

#creating array containing coefficients
Y = np.array([[1,1,0],[0,1,1],[1,0,1]])

#solving the equations
X = np.linalg.solve(Y,[c,a,b])

#printing output 
print(X)
