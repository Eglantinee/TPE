import numpy as np
from random import randint
from scipy.stats import f
from scipy.stats import t as par_t
from copy import deepcopy
import sys

x1_min = 10
x1_max = 60
x2_min = -25
x2_max = 10
x3_min = 10
x3_max=  15

y_min = 200 + sum([x1_min, x2_min, x3_min]) / len([x1_min, x2_min, x3_min])
y_max = 200 + sum([x1_max, x2_max, x3_max]) / len([x1_max, x2_max, x3_max])

m=3
p = 0.95
N = 4

# just norm matrix 

norm_matrix = [(1, -1, -1, -1),
               (1, -1, 1, 1),
               (1, 1, -1, 1),
               (1, 1, 1, -1),
               (1, -1, -1, 1),
               (1, -1, 1, -1),
               (1, 1, -1, -1),
               (1, 1, 1, 1)]

norm_matrix = norm_matrix[:N]

init_matrix  = np.delete(norm_matrix, 0, axis=1)

# generate matrix with min max values  
# TODO do smth with 

def mat(mx):
    matrix = mx
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] == 1:
                matrix[i][j] = globals()['x' + str(j + 1) + '_max']
            else:
                matrix[i][j] = globals()['x' + str(j + 1) + '_min']

mat(init_matrix)



def main_check(numerate):
    
    m = numerate
    
    
    def gen_matrix(num):
        
        matrix = init_matrix[:]
        
        for i in range(num):
            matrix = np.insert(matrix, 2, matrix[:,2], axis = 1)    # 2 is number of x columns (3-1) 
        
        for i in range(len(matrix)):
            
            tmp = len(matrix[i])
            
            for j in range(tmp - 1, tmp - num - 1, -1):
                matrix[i][j] = np.random.randint(y_min, y_max)
                
        return matrix
            
    matrix = gen_matrix(m)
    
    
    matrix  = [[ 10, -25,  10, 199, 224, 200],
                    [ 10,  10,  15, 210, 204, 217],
                    [ 60, -25,  15, 211, 218, 207],
                    [ 60,  10,  10, 201, 209, 227]]

       
    row = len(matrix[0])
        
    y1_abs = sum([matrix[0][j] for j in range(row - 1, row - m - 1, -1)]) / m
    y2_abs = sum([matrix[1][j] for j in range(row - 1, row - m - 1, -1)]) / m
    y3_abs = sum([matrix[2][j] for j in range(row - 1, row - m - 1, -1)]) / m
    y4_abs = sum([matrix[3][j] for j in range(row - 1, row - m - 1, -1)]) / m
    y_abs = (y1_abs, y2_abs, y3_abs, y4_abs)

    mx1 = sum([matrix[i][0] for i in range(N)]) / N
    mx2 = sum([matrix[i][1] for i in range(N)]) / N
    mx3 = sum([matrix[i][2] for i in range(N)]) / N
    mx = (mx1, mx2, mx3)   
    
    my = sum([y1_abs, y2_abs, y3_abs, y4_abs]) / N
    
    

    a1 = sum([matrix[i][0] * y_abs[i] for i in range(N)]) / N
    a2 = sum([matrix[i][1] * y_abs[i] for i in range(N)]) / N
    a3 = sum([matrix[i][2] * y_abs[i] for i in range(N)]) / N
    
    

    a11 = sum([matrix[i][0] ** 2 for i in range(N)]) / N
    a22 = sum([matrix[i][1] ** 2 for i in range(N)]) / N
    a33 = sum([matrix[i][2] ** 2 for i in range(N)]) / N
    
    
    def prod(lst):
        return lst[0] * lst[1]


    a12 = a21 = sum([prod(matrix[i][:2]) for i in range(N)]) / N
    a13 = a31 = sum([prod(matrix[i][:3:2]) for i in range(N)]) / N
    a23 = a32 = sum([prod(matrix[i][1:3]) for i in range(N)]) / N
    
    b = np.linalg.solve([[1, mx1, mx2, mx3],
                         [mx1, a11, a21, a31],
                         [mx2, a12, a22, a32],
                         [mx3, a13, a23, a33]],
                        [my, a1, a2, a3])
    
    # Regression equation =>  y1 = b[0] + b[1] * x1 + b[2] * x2 + b[3] * x3

    def check1():
        # Confirm that y[j]_abs = y[j]
        res = dict()
        for i in range(N):
            res.update(
                {y_abs[i]: b[0] + b[1] * matrix[i][0] + b[2] * matrix[i][1] +
                                                        b[3]
                                                        * matrix[i][2]})
        return res


    check1_res = check1()
    print(check1_res)
    print("Matrix:\n", matrix)
    print("Coefficients: ", b)
    print("Check results for truth: \n", check1_res)
    
     # Korhen check:

    def korhen(f1, f2, q=0.05):
        # Get critical value
        
        q1 = q / f1
        critical = f.ppf(q=1 - q1, dfn=f2, dfd=(f1 - 1) * f2)

        return critical / (critical + f2 - 1)


    s_y1 = sum([(matrix[0][j] - y1_abs) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y2 = sum([(matrix[1][j] - y2_abs) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y3 = sum([(matrix[2][j] - y3_abs) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y4 = sum([(matrix[3][j] - y4_abs) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y = (s_y1, s_y2, s_y3, s_y4)

    gp = max(s_y1, s_y2, s_y3, s_y4) / sum([s_y1, s_y2, s_y3, s_y4])
    f1 = m - 1
    f2 = N
    g_kr = korhen(f1, f2)

    print("Gp={}, G_kr={}".format(gp, g_kr))

    if gp < g_kr:
        print("dispersion is uniform")
        return True, s_y, y_abs, matrix, f1, f2, b
    else:
        print("dispersion isn't uniform")
        if m < 23:
            return main_check(m+1)
        else:
            return False
    
    

a, s_y, y_abs, matrix, f1, f2, b = main_check(m)

def check_for_adequate():
    if a:
        s2_b = sum(s_y) / N
        s2_beta = s2_b / (N * m)
        s_beta = np.sqrt(s2_beta)

        # Student criteria: we want to find which coefficients are valuable:
        #   to do it we compare t[i] with t_kr
        #   but t[i] = b_[i]/ s_beta
        #   so we find s2_b then s2_beta and then s_beta;  find b_[i] and then t[i]    <= using formulas

        b_0 = sum([y_abs[i] * norm_matrix[i][0] for i in range(N)]) / N
        b_1 = sum([y_abs[i] * norm_matrix[i][1] for i in range(N)]) / N
        b_2 = sum([y_abs[i] * norm_matrix[i][2] for i in range(N)]) / N
        b_3 = sum([y_abs[i] * norm_matrix[i][3] for i in range(N)]) / N
        b_sum = (b_0, b_1, b_2, b_3)
    
    
        f3 = f1 * f2
        t_kr = par_t.ppf(df=f3, q=(1 + p) / 2)

        d = 0
        t = [abs(b_sum[i]) / s_beta for i in range(N)]

        for i in range(len(t)):
            if t[i] < t_kr:
                t[i] = 0
            else:
                t[i] = 1
                d += 1

        print("b{} is unnecessary".format([i for i in range(len(t)) if t[i] == 0]))
    
        y1 = b[0] * t[0] + b[1] * t[1] * matrix[0][0] + b[2] * t[2] * matrix[0][1] + b[3] * t[3] * matrix[0][2]
        y2 = b[0] * t[0] + b[1] * t[1] * matrix[1][0] + b[2] * t[2] * matrix[1][1] + b[3] * t[3] * matrix[1][2]
        y3 = b[0] * t[0] + b[1] * t[1] * matrix[2][0] + b[2] * t[2] * matrix[2][1] + b[3] * t[3] * matrix[2][2]
        y4 = b[0] * t[0] + b[1] * t[1] * matrix[3][0] + b[2] * t[2] * matrix[3][1] + b[3] * t[3] * matrix[3][2]

    # Fisher criteria

        y = (y1, y2, y3, y4)

    s2_ad = (m / (N - d)) * sum(
        [(y[i] - y_abs[i]) ** 2 for i in range(N)])

    Fp = s2_ad / s2_beta
    f4 = N - d
    F_kr = f.ppf(dfn=f4, dfd=f3, q=0.95)

    print("Fp={}, F_kr={}".format(Fp, F_kr))

    if Fp < F_kr:
        print("Regression equation is adequate to original")
    else:
        print("Regression equation is inadequate to original")

check_for_adequate()

def second_check(number):
    
    N = 8
    global norm_matrix
    
    # norm_matrix = [[1, -1, -1, -1],
                   # [1, -1, -1, 1],
                   # [1, -1, 1, -1],
                   # [1, -1, 1, 1],
                   # [1, 1, -1, -1],
                   # [1, 1, -1, 1],
                   # [1, 1, 1, -1],
                   # [1, 1, 1, 1]]
    
    m = number
    
    matrix2 = deepcopy(norm_matrix)
        
    def mat(mx):
        matrix = mx
        for i in range(len(matrix)):
            for j in range(1,4):
                if matrix[i][j] == 1:
                    matrix[i][j] = globals()['x' + str(j) + '_max']
                else:
                    matrix[i][j] = globals()['x' + str(j) + '_min']
    
    mat(matrix2)
    
    def maker(mx):
        norm_matrix = mx
        for i in range(N):
            norm_matrix[i].append(norm_matrix[i][1]*norm_matrix[i][2])
            norm_matrix[i].append(norm_matrix[i][1]*norm_matrix[i][3])
            norm_matrix[i].append(norm_matrix[i][2]*norm_matrix[i][3])
            norm_matrix[i].append(norm_matrix[i][1]*norm_matrix[i][2]*norm_matrix[i][3])
        
    maker(norm_matrix)
    maker(matrix2)
    norm_matrix = np.array(norm_matrix)    
    
    def gen2_matrix(mx):
        for i in range(N):
            mx[i].remove(1)
            for j in range(m):
                mx[i].append(np.random.randint(y_min, y_max))
    
    gen2_matrix(matrix2)
    matrix2 = np.array(matrix2)
    # print(matrix2)
    
    row = len(matrix2[0])
    
    y1_abs = sum([matrix2[0][j] for j in range(row - 1, row - m - 1, -1)]) / m
    y2_abs = sum([matrix2[1][j] for j in range(row - 1, row - m - 1, -1)]) / m
    y3_abs = sum([matrix2[2][j] for j in range(row - 1, row - m - 1, -1)]) / m
    y4_abs = sum([matrix2[3][j] for j in range(row - 1, row - m - 1, -1)]) / m
    y5_abs = sum([matrix2[4][j] for j in range(row - 1, row - m - 1, -1)]) / m
    y6_abs = sum([matrix2[5][j] for j in range(row - 1, row - m - 1, -1)]) / m
    y7_abs = sum([matrix2[6][j] for j in range(row - 1, row - m - 1, -1)]) / m
    y8_abs = sum([matrix2[7][j] for j in range(row - 1, row - m - 1, -1)]) / m
    y_abs = (y1_abs,y2_abs,y3_abs,y4_abs,y5_abs,y6_abs,y7_abs,y8_abs)
    
    mx1 = sum([matrix2[i][0] for i in range(N)]) / N
    mx2 = sum([matrix2[i][1] for i in range(N)]) / N
    mx3 = sum([matrix2[i][2] for i in range(N)]) / N
    mx12 = sum([matrix2[i][3] for i in range(N)]) / N
    mx13 = sum([matrix2[i][4] for i in range(N)]) / N
    mx23 = sum([matrix2[i][5] for i in range(N)]) / N
    mx123 = sum([matrix2[i][6] for i in range(N)]) / N
    
    
    
    
    m00 = N / N
    m10 = sum([matrix2[i][0] for i in range(N)]) / N
    m20 = sum([matrix2[i][1] for i in range(N)]) / N
    m30 = sum([matrix2[i][2] for i in range(N)]) / N
    m40 = sum([matrix2[i][3] for i in range(N)]) / N
    m50 = sum([matrix2[i][4] for i in range(N)]) / N
    m60 = sum([matrix2[i][5] for i in range(N)]) / N
    m70 = sum([matrix2[i][6] for i in range(N)]) / N
    
    m01 = sum([matrix2[i][0] for i in range(N)]) / N
    m11 = sum([matrix2[i][0] * matrix2[i][0] for i in range(N)]) / N
    m21 = sum([matrix2[i][1] * matrix2[i][0] for i in range(N)]) / N
    m31 = sum([matrix2[i][2] * matrix2[i][0] for i in range(N)]) / N
    m41 = sum([matrix2[i][3] * matrix2[i][0] for i in range(N)]) / N
    m51 = sum([matrix2[i][4] * matrix2[i][0] for i in range(N)]) / N
    m61 = sum([matrix2[i][5] * matrix2[i][0] for i in range(N)]) / N
    m71 = sum([matrix2[i][6] * matrix2[i][0] for i in range(N)]) / N
    
    
    
    m02 = sum([matrix2[i][1] for i in range(N)]) / N
    m12 = sum([matrix2[i][0] * matrix2[i][1] for i in range(N)]) / N
    m22 = sum([matrix2[i][1] * matrix2[i][1] for i in range(N)]) / N
    m32 = sum([matrix2[i][2] * matrix2[i][1] for i in range(N)]) / N
    m42 = sum([matrix2[i][3] * matrix2[i][1] for i in range(N)]) / N
    m52 = sum([matrix2[i][4] * matrix2[i][1] for i in range(N)]) / N
    m62 = sum([matrix2[i][5] * matrix2[i][1] for i in range(N)]) / N
    m72 = sum([matrix2[i][6] * matrix2[i][1] for i in range(N)]) / N
    
    
    m03 = sum([matrix2[i][2] for i in range(N)]) / N
    m13 = sum([matrix2[i][0] * matrix2[i][2] for i in range(N)])  / N
    m23 = sum([matrix2[i][1] * matrix2[i][2] for i in range(N)]) / N
    m33 = sum([matrix2[i][2] * matrix2[i][2] for i in range(N)]) / N
    m43 = sum([matrix2[i][3] * matrix2[i][2] for i in range(N)]) / N
    m53 = sum([matrix2[i][4] * matrix2[i][2] for i in range(N)]) / N
    m63 = sum([matrix2[i][5] * matrix2[i][2] for i in range(N)]) / N
    m73 = sum([matrix2[i][6] * matrix2[i][2] for i in range(N)]) / N
    
    
    m04 = sum([matrix2[i][3] for i in range(N)]) / N
    m14 = sum([matrix2[i][0] * matrix2[i][3] for i in range(N)]) / N
    m24 = sum([matrix2[i][1] * matrix2[i][3] for i in range(N)]) / N
    m34 = sum([matrix2[i][2] * matrix2[i][3] for i in range(N)]) / N
    m44 = sum([matrix2[i][3] * matrix2[i][3] for i in range(N)]) / N
    m54 = sum([matrix2[i][4] * matrix2[i][3] for i in range(N)]) / N
    m64 = sum([matrix2[i][5] * matrix2[i][3] for i in range(N)]) / N
    m74 = sum([matrix2[i][6] * matrix2[i][3] for i in range(N)]) / N
    
    
    
    m05 = sum([matrix2[i][4] for i in range(N)]) / N
    m15 = sum([matrix2[i][0] * matrix2[i][4] for i in range(N)]) / N
    m25 = sum([matrix2[i][1] * matrix2[i][4] for i in range(N)]) / N
    m35 = sum([matrix2[i][2] * matrix2[i][4] for i in range(N)]) / N
    m45 = sum([matrix2[i][3] * matrix2[i][4] for i in range(N)]) / N
    m55 = sum([matrix2[i][4] * matrix2[i][4] for i in range(N)]) / N
    m65 = sum([matrix2[i][5] * matrix2[i][4] for i in range(N)]) / N
    m75 = sum([matrix2[i][6] * matrix2[i][4] for i in range(N)]) / N
    
    m06 = sum([matrix2[i][5] for i in range(N)]) / N
    m16 = sum([matrix2[i][0] * matrix2[i][5] for i in range(N)]) / N
    m26 = sum([matrix2[i][1] * matrix2[i][5] for i in range(N)]) / N
    m36 = sum([matrix2[i][2] * matrix2[i][5] for i in range(N)]) / N
    m46 = sum([matrix2[i][3] * matrix2[i][5] for i in range(N)]) / N
    m56 = sum([matrix2[i][4] * matrix2[i][5] for i in range(N)]) / N
    m66 = sum([matrix2[i][5] * matrix2[i][5] for i in range(N)]) / N
    m76 = sum([matrix2[i][6] * matrix2[i][5] for i in range(N)]) / N
    
    m07 = sum([matrix2[i][6] for i in range(N)]) / N
    m17 = sum([matrix2[i][0] * matrix2[i][6] for i in range(N)]) / N
    m27 = sum([matrix2[i][1] * matrix2[i][6] for i in range(N)]) / N
    m37 = sum([matrix2[i][2] * matrix2[i][6] for i in range(N)]) / N
    m47 = sum([matrix2[i][3] * matrix2[i][6] for i in range(N)]) / N
    m57 = sum([matrix2[i][4] * matrix2[i][6] for i in range(N)]) / N
    m67 = sum([matrix2[i][5] * matrix2[i][6] for i in range(N)]) / N
    m77 = sum([matrix2[i][6] * matrix2[i][6] for i in range(N)]) / N
    
    
    k0 = sum([y_abs[i] for i in range(N)]) / N
    k1 = sum([y_abs[i] * matrix2[i][0] for i in range(N)]) / N
    k2 = sum([y_abs[i] * matrix2[i][1] for i in range(N)]) / N
    k3 = sum([y_abs[i] * matrix2[i][2] for i in range(N)]) / N
    k4 = sum([y_abs[i] * matrix2[i][3] for i in range(N)]) / N
    k5 = sum([y_abs[i] * matrix2[i][4] for i in range(N)]) / N
    k6 = sum([y_abs[i] * matrix2[i][5] for i in range(N)]) / N
    k7 = sum([y_abs[i] * matrix2[i][6] for i in range(N)]) / N
    
      
    
    kkk  = np.linalg.solve([[m00, m10, m20, m30, m40, m50, m60, m70],
                            [m01, m11, m21, m31, m41, m51, m61, m71],
                            [m02, m12, m22, m32, m42, m52, m62, m72],
                            [m03, m13, m23, m33, m43, m53, m63, m73],
                            [m04, m14, m24, m34, m44, m54, m64, m74],
                            [m05, m15, m25, m35, m45, m55, m65, m75],
                            [m06, m16, m26, m36, m46, m56, m66, m76],
                            [m07, m17, m27, m37, m47, m57, m67, m77]],
                            [k0, k1, k2, k3, k4, k5, k6, k7])
    
        
    b = kkk
    
    
    def check1():
        # Confirm that y[j]_abs = y[j]
        res = dict()
        for i in range(N):
            res.update(
                {y_abs[i]: b[0] + b[1] * matrix2[i][0] + b[2] * matrix2[i][1] +
                                                        b[3]
                                                        * matrix2[i][2] + b[4] * matrix2[i][3] + b[5] * matrix2[i][4]+ b[6] *matrix2[i][5] + b[7] * matrix2[i][6]})
        return res
        
    check1 = check1()
    
    def korhen(f1, f2, q=0.05):
        
        # Get critical value
        
        q1 = q / f1
        critical = f.ppf(q=1 - q1, dfn=f2, dfd=(f1 - 1) * f2)

        return critical / (critical + f2 - 1)


    s_y1 = sum([(matrix2[0][j] - y1_abs) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y2 = sum([(matrix2[1][j] - y2_abs) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y3 = sum([(matrix2[2][j] - y3_abs) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y4 = sum([(matrix2[3][j] - y4_abs) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y5 = sum([(matrix2[4][j] - y5_abs) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y6 = sum([(matrix2[5][j] - y6_abs) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y7 = sum([(matrix2[6][j] - y7_abs) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y8 = sum([(matrix2[7][j] - y8_abs) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y = (s_y1, s_y2, s_y3, s_y4, s_y5, s_y6, s_y7, s_y8)

    gp = max(s_y1, s_y2, s_y3, s_y4, s_y5, s_y6, s_y7, s_y8) / sum(s_y)
    f1 = m - 1
    f2 = N
    g_kr = korhen(f1, f2)

    print("Gp={}, G_kr={}".format(gp, g_kr))

    if gp < g_kr:
        print("dispersion is uniform")
        return True, s_y, y_abs, matrix2, f1, f2, b, check1, norm_matrix
    else:
        print("dispersion isn't uniform")
        if m < 23:
            return second_check(m+1)
        else:
            return False
    
    
z, s_y, y_abs, matrix2, f1, f2, b, check1, norm_matrix  = second_check(m)
N = 8
if z:
    s2_b = sum(s_y) / N
    s2_beta = s2_b / (N * m)
    s_beta = np.sqrt(s2_beta)
        
    # Student criteria: we want to find which coefficients are valuable:
    #   to do it we compare t[i] with t_kr
    #   but t[i] = b_[i]/ s_beta
    #   so we find s2_b then s2_beta and then s_beta;  find b_[i] and then t[i]    <= using formulas
    
    
    b_0 = sum([y_abs[i] * norm_matrix[i][0] for i in range(N)]) / N
    b_1 = sum([y_abs[i] * norm_matrix[i][1] for i in range(N)]) / N
    b_2 = sum([y_abs[i] * norm_matrix[i][2] for i in range(N)]) / N
    b_3 = sum([y_abs[i] * norm_matrix[i][3] for i in range(N)]) / N
    b_4 = sum([y_abs[i] * norm_matrix[i][4] for i in range(N)]) / N
    b_5 = sum([y_abs[i] * norm_matrix[i][5] for i in range(N)]) / N
    b_6 = sum([y_abs[i] * norm_matrix[i][6] for i in range(N)]) / N
    b_7 = sum([y_abs[i] * norm_matrix[i][7] for i in range(N)]) / N
    
    
    
    
    b_sum = (b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7)
    
        # # print("********",b_sum)
    
    f3 = f1 * f2
    t_kr = par_t.ppf(df=f3, q=(1 + p) / 2)
    d = 0
    t = [abs(b_sum[i]) / s_beta for i in range(N)]
    for i in range(len(t)):
        if t[i] < t_kr:
            t[i] = 0
        else:
            t[i] = 1
            d += 1
    
    b = [b[i] * t[i] for i in range(len(t))]
    print(b)
    print("b{} is unnecessary".format([i for i in range(len(t)) if t[i] == 0]))
    
    
# Fisher criteria

    y = [b[0] + b[1] * matrix2[i][0] + b[2] * matrix2[i][1] + b[3] * matrix2[i][2] + b[4] * matrix2[i][3] + b[5] * matrix2[i][4]+ b[6] *matrix2[i][5] + b[7] * matrix2[i][6] for i in range(N)]
    
    
    s2_ad = (m / (N - d)) * sum(
        [(y[i] - y_abs[i]) ** 2 for i in range(N)])

    Fp = s2_ad / s2_beta
    f4 = N - d
    F_kr = f.ppf(dfn=f4, dfd=f3, q=0.95)

    print("Fp={}, F_kr={}".format(Fp, F_kr))

    if Fp < F_kr:
        print("Regression equation is adequate to original")
    else:
        print("Regression equation is inadequate to original")
        sys.exit()




