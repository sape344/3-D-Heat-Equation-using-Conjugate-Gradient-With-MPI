
from Draw3dPlot import Draw3d
import numpy as np
from scipy.sparse import linalg
import time



def CG(A,b,tol):
    x=np.random.rand(A.shape[0])
    rk = np.dot(A,x)-b
    pk=-rk
    rk_norm= np.linalg.norm(rk)

    num_iter=0
    #curve_x = [x]
    while rk_norm>tol:
        apk = np.dot(A,pk)
        rkrk=np.dot(rk,rk)
        alpha = rkrk /np.dot(pk,apk)
        x=x+alpha*pk
        rk= rk+ alpha *apk
        beta = np.dot(rk,rk)/rkrk
        pk = -rk +beta*pk
        num_iter += 1
       # curve_x.append(x)
        rk_norm = np.linalg.norm(rk)
        # print('Iteration: {} \t x = {} \t residual = {:.4f}'.
        #       format(num_iter, x, rk_norm))
    print('\nSolution: \t x = {} \t iter_num ={}'.format(x,num_iter))

    return x


    err=1/2
print("3D heat equation solver Initialize")
plate_length = 20
draw=Draw3d(plate_length)

u = np.empty((plate_length, plate_length, plate_length),dtype=float)


# Initial condition everywhere inside the grid
u_initial = 0

# Boundary conditions
u_top =  100
u_left = 0
u_bottom = 0
u_right = 100
u_front = 0
u_back = 100

## Set the initial condition
u.fill(u_initial)

# Set the boundary conditions
u[ :, :,(plate_length-1):] = u_left
u[ :, :,:1] = u_right

u[ :, :1,:] = u_front
u[ :, (plate_length-1):,:] = u_back

u[ (plate_length-1):, :,:] = u_top
u[ :1, :,:] = u_bottom



draw.SetData(u)
draw.Plot()

N=plate_length**3

A=np.zeros((N,N))
b=np.zeros((N))
t=np.zeros((N))
alpha=1 # akışkanlık katsayısı
step=1 #discreate step size
coeff2= alpha/(step**2) #
coeff=-1

for i in range(N):

    #A[i, i ]=1-6*coeff
    A[i, i ]=1-(1-6*coeff2)

    if (i%plate_length!=0): #left
        A[i,i-1]=coeff
    else:
        b[i]+=u_left


    if (i%plate_length!=plate_length-1): #right
        A[i,i+1]=coeff
    else:
        b[i]+=u_right

    if (np.floor((i%plate_length**2)/plate_length)!=0): #back
        A[i,i-plate_length]=coeff
    else:
        b[i]+=u_back

    if (np.floor((i%plate_length**2)/plate_length)!=plate_length-1): #front
        A[i,i+plate_length]=coeff
    else:
        b[i]+=u_front

    if (np.floor(i/plate_length**2)!=0):  # Top
        A[i, i - plate_length**2 ] = coeff
    else:
        b[i] += u_bottom

    if (np.floor(i/plate_length**2)!=plate_length-1):  # Bottom
        A[i, i + plate_length**2 ] = coeff
    else:
        b[i] += u_top

print("A: ",A)
print("I-A: ",np.identity(N)-A)


print("A matrix is created")
start=time.time()
t0=np.matmul(np.linalg.inv(A),b)
stop=time.time()

print("Classic method : ",stop-start, " Second")
#t0=np.matmul(np.linalg.inv(A),b)


draw.SetData(t0)
draw.Plot()
start=time.time()
t1=linalg.cg(A,b,tol=1e-12)
stop=time.time()
print("CG method : ",stop-start, " Second")
#t1=linalg.cg(A,b,tol=1e-12)

print("Error rate: ",np.sum(np.power(t1[0]-t0,2)))
draw.SetData(t1[0])
draw.Plot()

start=time.time()
t2=CG(A,b,1e-12)
stop=time.time()
print("MY CG method : ",stop-start, " Second")

print("Error rate: ",np.sum(np.power(t2-t0,2)))
draw.SetData(t2)
draw.Plot()
# print("t= {")
# for i in t0:
#     print(i ,", ")
#
# print("};")
# for i in range(10):
#     t1 = linalg.cg(A, t1[0]-b,tol=1e-9)
#
#
#
#     print("Error rate: ",np.sum(np.power(t1[0]-t0,2)))
#
#
#     draw.SetData(t1[0])
#     draw.Plot()

print("T is calculated")







