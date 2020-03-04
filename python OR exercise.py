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
    ent = 0 #entering value
    leave = 0
    N = A
    prev_z = -1
    #print(all_less_than(cond, 0).all())
    while all_less_than(cond, 0).all() == True and i < 4:
        print(D)
        Dinv = inv(D)
        for j in range(ma):
            z_xj = c_D.T * Dinv * N[:, j] - c[j]
            if all_less_than(z_xj, prev_z).all() == True and all_less_than(z_xj, 0).all() == True:
                 prev_z = z_xj[0, 0]
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
        #Swap entering and leaving value
        tmp1 = D[:, leave]
        tmp2 = N[:, ent]
        D[:, leave] = tmp2
        N[:, ent] = tmp1
        i = i + 1
    x_D = Dinv * b
    z = c_D * x_D
    return z


def inv(A): # Inverse of matrix A
    [na, ma] = A.shape #size of A
    Id = np.identity(na) #identity matrix of dimentions na x na
    #print(np.equal(A, Id).all())
    if np.equal(A, Id).all():
        return A
    else:
        Echilon_form = A
        for i in range(na-1):
            a_ii = Echilon_form[i, i]
            Echilon_form[i, :] = (1/a_ii) * Echilon_form[i, :]
            for j in range(i+1, na):
                a_ij = Echilon_form[i, j]
                Echilon_form[j, :] = (-1 * a_ij/a_ii) * Echilon_form[i, :] + Echilon_form[j, :]
                Id[j, :] = (-1 * a_ij/a_ii) * Id[i, :] + Id[j, :]
        Ainv = Id
        for i2 in range(na-1, 1, -1):
            a_ii = Echilon_form[i2, i2]
            Echilon_form[i2, :] = (1/a_ii) * Echilon_form[i2, :]
            for j2 in range(na-2, 0, -1):
                a_ij = Echilon_form[i2, j2]
                Echilon_form[j2, :] = (-1 * a_ij/a_ii) * Echilon_form[i2, :] + Echilon_form[j2, :]
                Ainv[j2, :] = (-1 * a_ij/a_ii) * Ainv[i2, :] + Ainv[j2, :]
        return Ainv

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
print(z)



