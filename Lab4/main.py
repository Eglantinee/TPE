import numpy as np
from random import randint
from scipy.stats import f
from scipy.stats import t as par_t
from copy import deepcopy
import sys
from random import uniform

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

initial_matrix = [[1, -1, -1, -1],
                   [1, -1, 1, 1],
                   [1, 1, -1, 1],
                   [1, 1, 1, -1],
                   [1, -1, -1, 1],
                   [1, -1, 1, -1],
                   [1, 1, -1, -1],
                   [1, 1, 1, 1]]

def abs_maker(mx, column, start, stop, step, val):
        return sum([mx[column][j] for j in range(start, stop, step)]) / val
        
        
def beta(lst, mx, n):
        ac = 0
        res = []
        for j in range(n):
            for i in range(n):
                ac += lst[i] * mx[i][j]
                if i == n - 1:
                    res.append(ac / n)
                    ac = 0
        return res


def add_y_vals(func):
    def inner(*args):
        returned_values = func(*args)
        mx = returned_values[0]
        m = returned_values[1]
        for i in range(len(mx)):
            for j in range(m):
                mx[i].append(round(uniform(y_min, y_max), ndigits=3))
        return mx
    return inner
    
    
@add_y_vals
def gen_matrix(mx, num, mins, maxes):
    for i in range(len(mx)):
        mx[i] = list(mx[i])
        mx[i].pop(0)
        for j in range(len(mx[i])):
            if mx[i][j] == -1:
                mx[i][j] = mins[j]
            elif mx[i][j] == 1:
                mx[i][j] = maxes[j]
    return mx, num


def second_check(s_y, y_abs, matrix2, f1, f2, b, check1, norm_matrix):
    N = 8
    m=3
    s2_b = sum(s_y) / N
    s2_beta = s2_b / (N * m)
    s_beta = np.sqrt(s2_beta)
        
    b_sum = beta(y_abs, norm_matrix, N)
        
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
    
    
    s2_ad = (m / (N - d)) * sum([(y[i] - y_abs[i]) ** 2 for i in range(N)])

    Fp = s2_ad / s2_beta
    f4 = N - d
    F_kr = f.ppf(dfn=f4, dfd=f3, q=0.95)

    print("Fp={}, F_kr={}".format(Fp, F_kr))

    if Fp < F_kr:
        print("Regression equation is adequate to original")
    else:
        print("Regression equation is inadequate to original")
        start(4)

def finish(number):
    
   # define static variables
    N = 8
    m = 3
    x_mins = [x1_min, x2_min, x3_min]
    x_maxes = [x1_max, x2_max, x3_max]
    norm_matrix = deepcopy(initial_matrix[:N])
    matrix2 = deepcopy(norm_matrix)
    solve_mx = [[] for i in range(N)]
            
    matrix2 = gen_matrix(matrix2, m, x_mins, x_maxes)
               
    def cooperation(mx, pos, zero):
        for i in range(N):
            mx[i].insert(pos, mx[i][zero] * mx[i][zero+1])
            mx[i].insert(pos + 1, mx[i][zero] * mx[i][zero+2])
            mx[i].insert(pos + 2, mx[i][zero+1] * mx[i][zero+2])
            mx[i].insert(pos + 3, mx[i][zero] * mx[i][zero+1] * mx[i][zero+2])
        
    cooperation(norm_matrix, 4, 1)
    cooperation(matrix2, 3, 0)
    norm_matrix = np.array(norm_matrix)    
    
    row = len(matrix2[0])
    y_abs = [abs_maker(matrix2, i, row - 1, row - m - 1, -1, m) for i in range(N)]
        
    
    def mk_solve(mx1, res_mx):
        vector0 = [1]
        tmp = [sum(mx1[i][j] for i in range(N)) / N for j in range(N-1)]
        for i in tmp:
            vector0.append(i)
        k = 0
        
        def mk_vector():
            nonlocal k
            nonlocal mx1
            mx = mx1
            vector = []
            
            for z in range(N):
                vector.clear()
                flag = True
                
                for j in range(N):
                    suma = 0
                    if flag:
                        suma = sum(mx[i][k] for i in range(N)) 
                        flag = False
                    else:
                        j -= 1
                        for i in range(N):
                            suma += mx[i][j] * mx[i][k]
                    vector.append(suma / N)
                k += 1
                return  vector
        
        res_mx[0] = vector0
        for i in range(1, len(res_mx)):
            res_mx[i] = mk_vector()
                
                
    mk_solve(matrix2, solve_mx)
        
    
    
    
    k0 = sum([y_abs[i] for i in range(N)]) / N
    k1 = sum([y_abs[i] * matrix2[i][0] for i in range(N)]) / N
    k2 = sum([y_abs[i] * matrix2[i][1] for i in range(N)]) / N
    k3 = sum([y_abs[i] * matrix2[i][2] for i in range(N)]) / N
    k4 = sum([y_abs[i] * matrix2[i][3] for i in range(N)]) / N
    k5 = sum([y_abs[i] * matrix2[i][4] for i in range(N)]) / N
    k6 = sum([y_abs[i] * matrix2[i][5] for i in range(N)]) / N
    k7 = sum([y_abs[i] * matrix2[i][6] for i in range(N)]) / N
    
    b = np.linalg.solve(solve_mx,
                            [k0, k1, k2, k3, k4, k5, k6, k7])
    
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
    print(check1)
    
    def korhen(f1, f2, q=0.05):
        
        # Get critical value
        
        q1 = q / f1
        critical = f.ppf(q=1 - q1, dfn=f2, dfd=(f1 - 1) * f2)

        return critical / (critical + f2 - 1)


    s_y1 = sum([(matrix2[0][j] - y_abs[0]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y2 = sum([(matrix2[1][j] - y_abs[1]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y3 = sum([(matrix2[2][j] - y_abs[2]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y4 = sum([(matrix2[3][j] - y_abs[3]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y5 = sum([(matrix2[4][j] - y_abs[4]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y6 = sum([(matrix2[5][j] - y_abs[5]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y7 = sum([(matrix2[6][j] - y_abs[6]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y8 = sum([(matrix2[7][j] - y_abs[7]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y = (s_y1, s_y2, s_y3, s_y4, s_y5, s_y6, s_y7, s_y8)

    gp = max(s_y1, s_y2, s_y3, s_y4, s_y5, s_y6, s_y7, s_y8) / sum(s_y)
    f1 = m - 1
    f2 = N
    g_kr = korhen(f1, f2)

    print("Gp={}, G_kr={}".format(gp, g_kr))

    if gp < g_kr:
        print("dispersion is uniform")
        second_check(s_y, y_abs, matrix2, f1, f2, b, check1, norm_matrix)
    else:
        print("dispersion isn't uniform")
        if m < 23:
            return finish(m+1)
        else:
            return False





def check_for_adequate(s_y, y_abs, matrix, f1, f2, b, norm_matrix, N):
    
    #Define variables
    norm_matrix = norm_matrix
    matrix = matrix
    N = N
    s_y = s_y
    y_abs = y_abs
    f1 = f1
    f2 = f2
    b = b
    
    # Main action
     
    s2_b = sum(s_y) / N
    s2_beta = s2_b / (N * m)
    s_beta = np.sqrt(s2_beta)

    # Student criteria: we want to find which coefficients are valuable:
    #   to do it we compare t[i] with t_kr
    #   but t[i] = b_[i]/ s_beta
    #   so we find s2_b then s2_beta and then s_beta;  find b_[i] and then t[i]    <= using formulas
    
    b_sum = beta(y_abs, norm_matrix, N)
    
    # b_0 = sum([y_abs[i] * norm_matrix[i][0] for i in range(N)]) / N    
    
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
    print("b{} is unnecessary".format([i for i in range(len(t)) if t[i] == 0]))
    
    y = [b[0] + b[1] * matrix[i][0] + b[2] * matrix[i][1] + b[3] * matrix[i][2] for i in range(4)]
    
    # Fisher criteria    
    
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
        print('*'*50)
        print(" ")
        finish(8)


def start(number):
    
    # define static variables
    N = 4
    m = number
    x_mins = [x1_min, x2_min, x3_min]
    x_maxes = [x1_max, x2_max, x3_max]
    norm_matrix = deepcopy(initial_matrix[:N])
    matrix = deepcopy(norm_matrix)
    
    matrix = gen_matrix(matrix, m, x_mins, x_maxes)
    row = len(matrix[0])
    
    # Equals to y1_abs = sum([matrix[0][j] for j in range(row - 1, row - m - 1, -1)]) / m
    def abs_maker(mx, column, start, stop, step, val):
        return sum([mx[column][j] for j in range(start, stop, step)]) / val
        
            
    y_abs = [abs_maker(matrix, i, row - 1, row - m - 1, -1, m) for i in range(N)]
    
    # mmx = manual_maker(matrix, 3, N)
    
    mx1 = sum([matrix[i][0] for i in range(N)]) / N
    mx2 = sum([matrix[i][1] for i in range(N)]) / N
    mx3 = sum([matrix[i][2] for i in range(N)]) / N
    mx = (mx1, mx2, mx3)  
    
    # print("***"*30)
    # print(mmx)
    # print(mx) 
    
    my = sum(y_abs) / N
    
    

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


    s_y1 = sum([(matrix[0][j] - y_abs[0]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y2 = sum([(matrix[1][j] - y_abs[1]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y3 = sum([(matrix[2][j] - y_abs[2]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y4 = sum([(matrix[3][j] - y_abs[3]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    s_y = (s_y1, s_y2, s_y3, s_y4)

    gp = max(s_y1, s_y2, s_y3, s_y4) / sum([s_y1, s_y2, s_y3, s_y4])
    f1 = m - 1
    f2 = N
    g_kr = korhen(f1, f2)

    print("Gp={}, G_kr={}".format(gp, g_kr))

    if gp < g_kr:
        print("dispersion is uniform")
        check_for_adequate(s_y, y_abs, matrix, f1, f2, b, norm_matrix, N)
    else:
        print("dispersion isn't uniform")
        if m < 23:
            return start(m+1)
        else:
            return False
    


start(m)
