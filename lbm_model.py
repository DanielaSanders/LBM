import matplotlib.pyplot as plt
import numpy as np

#particle speed direction

ex=np.zeros(9)
ey=np.zeros(9)

ex[1] = 1
ex[3] = -1
ex[5] = 1
ex[6] = -1
ex[7] = -1
ex[8] = 1

ey[2] = 1
ey[4] = -1
ey[5] = 1
ey[6] = 1
ey[7] = -1
ey[8] = -1

LY = 10
LX = 20
tau = 1.
g = 0.00001

#set solid nodes

is_solid_node = np.zeros((LY, LX), dtype=bool)
for i in range(LX):
    is_solid_node[1,i] = 1 
    is_solid_node[LY-1,i] = 1
#print(is_solid_node)


#define initial density and fs

rho = np.ones((LY,LX))

f = np.ones((LY,LX,9))
f[:,:,0] = 4/9.*rho
f[:,:,1] = 1/9.*rho
f[:,:,2] = 1/9.*rho
f[:,:,3] = 1/9.*rho
f[:,:,4] = 1/9.*rho
f[:,:,5] = 1/36.*rho
f[:,:,6] = 1/36.*rho
f[:,:,7] = 1/36.*rho
f[:,:,8] = 1/36.*rho

#Computing macroscopic density, rho, and velocity, u=(ux,uy)

u_x = np.zeros((LY,LX))
u_y = np.zeros((LY,LX))
rho = np.zeros((LY,LX))
def update_macro(u_x, u_y, rho, f):
    for j in range(LY):
        for i in range (LX):
            u_x[j][i] = 0.0
            u_y[j][i] = 0.0
            rho[j][i] = 0.0
            if(not is_solid_node[j][i]):
                for a in range(9):
                    u_x[j][i] = u_x[j][i] + ex[a]*f[j][i][a]
                    u_y[j][i] = u_y[j][i] + ey[a]*f[j][i][a]
                    rho[j][i] = rho[j][i] + f[j][i][a]
                u_x[j][i] = u_x[j][i] / rho[j][i]
                u_y[j][i] = u_y[j][i] / rho[j][i]
    return {'u_x':u_x[j][i], 'u_y':u_y[j][i], 'rho':rho[j][i]}

u_x, u_y, rho = update_macro(u_x, u_y, rho, f)

#periodic boundaries

i = 0
j = 0

def apply_pbc(i, j):
    ip = (i < LX - 1)
    if ip:
        (i + 1)
    else:
        0

    i_n = (i > 0)      
    if i_n:
        (i - 1)
    else:
        (LX - 1)  
        
    jp = (j < LY - 1)
    if jp:
        (j + 1)
    else:
        0

    jn = (j > 0)
    if jn:
        (j - 1)
    else:
        (LY - 1) 
    return i, j

pbc = apply_pbc(i, j) 

#Compute the equilibrium distribution function, feq.

rhoij = np.zeros((LY,LX))
uxij = np.zeros((LY,LX))
uyij = np.zeros((LY,LX))
feqij = np.zeros((LY,LX,9))
                 

f1 = 3.
f2 = 9./2.
f3 = 3./2.

def equilibrium(uxij, uyij, rhoij):
    for j in range (LY):
        for i in range (LX):
            if(not is_solid_node[j][i]):
                rt0 = (4./9.)*rhoij
                rt1 = (1./9.)*rhoij
                rt2 = (1./36.)*rhoij
                ueqxij = uxij
                ueqyij = uyij
                uxsq = ueqxij * ueqxij
                uysq = ueqyij * ueqyij
                uxuy5 = ueqxij + ueqyij
                uxuy6 = -ueqxij + ueqyij
                uxuy7 = -ueqxij + -ueqyij
                uxuy8 = ueqxij + -ueqyij
                usq = uxsq + uysq
                feqij[:,:,0] = rt0 * (1. - f3 * usq)
                feqij[:,:,1] = rt1 * (1. + f1 * ueqxij + f2 * uxsq - f3 * usq)
                feqij[:,:,2] = rt1 * (1. + f1 * ueqyij + f2 * uysq - f3 * usq)
                feqij[:,:,3] = rt1 * (1. + f1 * ueqxij + f2 * uxsq - f3 * usq)
                feqij[:,:,4] = rt1 * (1. + f1 * ueqyij + f2 * uysq - f3 * usq)
                feqij[:,:,5] = rt2 * (1. + f1 * uxuy5 + f2 * uxuy5 * uxuy5 - f3 * usq)
                feqij[:,:,6] = rt2 * (1. + f1 * uxuy6 + f2 * uxuy6 * uxuy6 - f3 * usq)
                feqij[:,:,7] = rt2 * (1. + f1 * uxuy7 + f2 * uxuy7 * uxuy5 - f3 * usq)
                feqij[:,:,8] = rt2 * (1. + f1 * uxuy8 + f2 * uxuy8 * uxuy5 - f3 * usq)
    return feqij

feqij = equilibrium(uxij, uyij, rhoij)
 


#Collision step

def collision(feqij):
    for j in range (LY):
        for i in range (LX):
            if (not is_solid_node[j][i]):
                for a in range (8):
                    f[j,i,a] = f[j,i,a] + (f[j,i,a] - feqij[j,i,a] / tau)
    return f

f = collision(feqij)

#Standard bounceback

temp = f[:,:,1]; f[:,:,1] = f[:,:,3]; f[:,:,3] = temp
temp = f[:,:,2]; f[:,:,2] = f[:,:,4]; f[:,:,4] = temp
temp = f[:,:,5]; f[:,:,5] = f[:,:,7]; f[:,:,7] = temp
temp = f[:,:,6]; f[:,:,6] = f[:,:,8]; f[:,:,8] = temp


#Streaming step

is_interior_solid_node = np.zeros((LY, LX))
ftemp = np.zeros((LY, LX,9))

#print(temp.shape)

# for j in range (LY):
    
#     jn = (j - 1) if (j > 0) else (LY - 1)
#     jp = (j + 1) if (j < LY - 1) else 0

#     for i in range (LX):
#         if (not is_interior_solid_node[j][i]):
#             i_n =  (i - 1) if (i > 0) else (LX - 1)
#             i_n = int(i_n)
#             ip = (i < (LX - 1)) if (i + 1) else 0
#             #Typecasting of the boolean
#             ip = int(ip)
#             ftemp[j][i][0] = f[j][i][0]
#             ftemp[j][ip][1] = f[j][i][1]
#             ftemp[jp][i][2] = f[j][i][2]
#             ftemp[j][i_n][3] = f[j][i][3]
#             ftemp[jn][i][4] = f[j][i][4]
#             ftemp[jp][ip][5] = f[j][i][5]
#             ftemp[jp][i_n][6] = f[j][i][6]
#             ftemp[jn][i_n][7] = f[j][i][7]
#             ftemp[jn][ip][8] = f[j][i][8] 

def streaming(i, j): 
    for j in range (LY):
    
        jn = (j - 1) if (j > 0) else (LY - 1)
        jp = (j + 1) if (j < LY - 1) else 0

        for i in range (LX):
            if (not is_interior_solid_node[j][i]):
                i_n =  (i - 1) if (i > 0) else (LX - 1)
                i_n = int(i_n)
                ip = (i < (LX - 1)) if (i + 1) else 0
                #Typecasting of the boolean
                ip = int(ip)
                ftemp[j][i][0] = f[j][i][0]
                ftemp[j][ip][1] = f[j][i][1]
                ftemp[jp][i][2] = f[j][i][2]
                ftemp[j][i_n][3] = f[j][i][3]
                ftemp[jn][i][4] = f[j][i][4]
                ftemp[jp][ip][5] = f[j][i][5]
                ftemp[jp][i_n][6] = f[j][i][6]
                ftemp[jn][i_n][7] = f[j][i][7]
                ftemp[jn][ip][8] = f[j][i][8]   
    return ftemp

ftemp = streaming(i, j)

# #Function Definitions

# # def equilibrium(uxij, uyij, rhoij):
# #     feqij = np.zeros((LY,LX,9))
# #     return feq

# def equilibrium(LY, LX, is_solid_node, uxij, uyij, rhoij):
#     feqij = np.zeros((LY,LX,9))
#     return feq

# #feqij = equilibrium(uxij, uyij, rhoij)
# feqij = equilibrium(LY, LX, is_solid_node, uxij, uyij, rhoij)


# def collision(feqij):
#     #f = collision(feqij)
#     f[j,i,a] = f[j,i,a] + (f[j,i,a] - feqij[j,i,a] / tau)
#     return f

# f = collision(feqij)

# def streaming(ip, i_n, jp, jn):
#     ftemp[j][i][0] = f[j][i][0]
#     ftemp[j][ip][1] = f[j][i][1]
#     ftemp[jp][i][2] = f[j][i][2]
#     ftemp[j][i_n][3] = f[j][i][3]
#     ftemp[jn][i][4] = f[j][i][4]
#     ftemp[jp][ip][5] = f[j][i][5]
#     ftemp[jp][i_n][6] = f[j][i][6]
#     ftemp[jn][i_    u_x[j][i] = u_x[j][i] + ex[a]*f[j][i][a]
                
#     ftemp[jn][ip    u_x[j][i] = u_x[j][i] + ex[a]*f[j][i][a]
                
#     return ftemp    u_x[j][i] = u_x[j][i] + ex[a]*f[j][i][a]
                

# ftemp = streaming(ip, i_n, jp, jn)


#Main time loop

maxIter = 4000 #total number of iterations

for time in range(maxIter):
    macro = update_macro(u_x, u_y, rho, f)
    pbc = apply_pbc(i, j)
    feqij = equilibrium(uxij, uyij, rhoij)
    f = collision(feqij)
    ftemp = streaming(i, j)

