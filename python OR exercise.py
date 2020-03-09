#Maximise z = 2x1 + 3x2 + 4x3 + x4 + 8x5 + x6
#subject to
#               x1 - x2 + 2x3 + x5 + x6 = 18
#               x2 - x3 + x4 + 3x6 <= 8
#               x1 + x2 - 3x3 + x4 + x5 <= 36
#               x1 - x2 + x5 + x6 <= 23
#               x1, x2, x3, x4, x5, x6 >= 0

#imports
import numpy as np

# Z = vx sbject to A_x_ = _b_

#1. Create function with variables of matrix A and vector b to get serplus matrix A_
def cannon(A, B):
    (n, m) = A.shape
    A_ = np.asmatrix(np.zeros((n, m+n), np.float))
    A_[:, :m] = A
    for i in range(0, n):
        A_[i, m+i] = 1
    return A_


def simplex(c, A, b):
    z = 0
    cprime = c.T
    (na, ma) = A.shape
    A_ = cannon(A, b)
    (na_, ma_) = A_.shape
    nc = c.size
    Tableau = np.zeros((na_ + 1, ma_ + 2), np.float)
    M = np.zeros((na_+ 1, ma_ + 1), np.float)
    V = np.zeros((ma_ + 1, 1), np.float)
    Tableau[0, 0] = 1
    M[0, 0] = 1
    for i in range(0, nc):
        Tableau[0, i+1] = -1 * c[i]
        M[0, i+1] = -1 * c[i]
    Tableau[1:na_+1, 1:ma_+1] = A_
    M[1:na+1, 1:ma_+1] = A_
    Tableau[1:na_+1, ma_+1:ma+2] = b
    V[1:na_+1, :] = b
    #print(Tableau)
    #Iteration 1
    D = np.asmatrix(np.identity(na)) #Let D be the identity matrix
    Dinv = D
    c_D = np.asmatrix(np.zeros((na, 1), np.float))
    #print(c_D)
    cond = np.array([[-1]], np.float) #Conditional value
    i = 0
    #ent = 0 #entering value
    #ent_prev = 0
    leave = 0
    N = A
    #prev_z = -1
    #print(all_less_than(cond, 0).all())
    while all_less_than(cond, 0).all() == True: #and i < 4:
        ent = 0
        ent_prev = 0
        prev_z = -1
        #print(D)
        #print(D.I)
        Dinv = D.I #inverse of D
        for j in range(ma):
            z_xj = c_D.T * Dinv * N[:, j] - c[j]
            #print(z_xj)
            if all_less_than(z_xj, prev_z).all() == True and all_less_than(z_xj, 0).all() == True:
                 prev_z = z_xj[0, 0]
                 ent_prev = ent
                 ent = j
            elif all_less_than(z_xj, 0).any() == False:
                 cond = z_xj
                 break;
        nxt = Dinv * b
        if all_less_than(nxt, 0).all() == False:
            prev_b = 100000000000
            for k in range(na):
                #print(nxt[k, k])
                if nxt[k, 0] < prev_b:
                    prev_b = nxt[k, 0]
                    leave = k
        else:
            break
        #Swap entering and leaving values
        if N[leave, ent] == 0:
            ent = ent_prev
        N, D = swap(N, D, ent, leave)
        c, c_D = swap(c, c_D, ent, leave)
        #temp1 = D[:, leave]
        #D[:, leave] = N[:, ent]
        #N[:, ent] = temp1
        i = i + 1
    x_D = Dinv * b
    z = c_D.T * x_D
    return z[0, 0]

def swap(N, B, x, y): #Swaps the values within a matrix
    [nn, mn] = N.shape
    [nb, mb] = B.shape
    
    N0 = np.asmatrix(np.zeros((nn, mn)))
    B0 = np.asmatrix(np.zeros((nb, mb)))
    
    if mn > 1 and mb > 1:
        N0[:, x] = B[:, y] - N[:, x]
        B0[:, y] = N[:, x] - B[:, y]
    else:
        N0[x, :] = B[y, :] - N[x, :]
        B0[y, :] = N[x, :] - B[y, :]  
    
    N = N + N0
    B = B + B0
    return N, B

def all_less_than(A, num): #A is multidimentional array, num is scalar value
    A_linear = A.flatten()
    y = lambda a: a < num
    return y(A_linear)
                 
##############################=== TEST ===######################################
A = np.matrix([[1, -1, 2, 0, 1, 1],
         [0, 1, -1, 1, 0, 3],
         [1, 1, -3, 1, 1, 0],
         [1, -1, 0, 0, 1, 1]])

b = np.matrix([[18], [8], [36], [23]])

c = np.matrix([[2], [3], [4], [1], [8], [1]])
#print(c.size)

z = simplex(c, A, b)
print("The optimal value of z is: ")
print(z)



