from math import sqrt
from scipy import linalg
from scipy.stats import f
from scipy.stats import t as par_t
from copy import deepcopy
from random import uniform
from MyPrettyTable import table as tb

x1_min = 10
x1_max = 60
x2_min = -25
x2_max = 10
x3_min = 10
x3_max = 15

y_min = 200 + sum([x1_min, x2_min, x3_min]) / len([x1_min, x2_min, x3_min])
y_max = 200 + sum([x1_max, x2_max, x3_max]) / len([x1_max, x2_max, x3_max])

p = 0.95
N = 15

initial_matrix = [[-1, -1, -1],
                  [-1, -1, 1],
                  [-1, 1, -1],
                  [-1, 1, 1],
                  [1, -1, -1],
                  [1, -1, 1],
                  [1, 1, -1],
                  [1, 1, 1],
                  [-1.215, 0, 0],
                  [1.215, 0, 0],
                  [0, -1.215, 0],
                  [0, 1.215, 0],
                  [0, 0, -1.215],
                  [0, 0, 1.215],
                  [0, 0, 0]]


def abs_maker(mx, row, start, stop, step, val):
    return sum([mx[row][j] for j in range(start, stop, step)]) / val


def check(mx, coef, y):
    y_abs = y
    b = coef
    matrix = mx
    # Confirm that y[j]_abs = y[j]
    res = dict()
    for i in range(15):
        res.update(
            {y_abs[i]: b[0] + b[1] * matrix[i][0] + b[2] * matrix[i][1] + b[3] * matrix[i][2] + b[4] * matrix[i][3] + b[
                5] * matrix[i][4] + b[6] * matrix[i][5] + b[7] * matrix[i][6] + b[8] * matrix[i][7] + b[9] * matrix[i][
                           8] + b[10] * matrix[i][9]})
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
def gen_matrix(mx, num):
    for i in range(len(mx)):
        st = 0
        ac = st + 1
        fixed = len(mx[i]) - 1
        while st < fixed:
            mx[i].append(mx[i][st] * mx[i][ac])
            ac += 1
            if ac > fixed:
                st += 1
                ac = st + 1
        mx[i].append(mx[i][0] * mx[i][1] * mx[i][2])
    for i in range(len(mx)):
        mx[i].append(mx[i][0] ** 2)
        mx[i].append(mx[i][1] ** 2)
        mx[i].append(mx[i][2] ** 2)
    for i in range(len(mx)):  # add zero column for student criteria
        mx[i].insert(0, 1)
    return mx, num


@add_y_vals
def naturalize_matrix(mx, num, mins, maxes):
    for i in range(len(mx)):
        for j in range(3):  # 3 means x1, x2, x3 only
            x_zero = (maxes[j] + mins[j]) / 2
            delta_x = maxes[j] - x_zero
            if mx[i][j] == -1:
                mx[i][j] = mins[j]
            elif mx[i][j] == 1:
                mx[i][j] = maxes[j]
            elif abs(mx[i][j]) == 1.215:
                mx[i][j] = mx[i][j] * delta_x + x_zero
            elif mx[i][j] == 0:
                mx[i][j] = x_zero

    for i in range(len(mx)):
        st = 0
        ac = st + 1
        fixed = len(mx[i]) - 1
        while st < fixed:
            mx[i].append(mx[i][st] * mx[i][ac])
            ac += 1
            if ac > fixed:
                st += 1
                ac = st + 1
        mx[i].append(mx[i][0] * mx[i][1] * mx[i][2])
    for i in range(len(mx)):
        mx[i].append(mx[i][0] ** 2)
        mx[i].append(mx[i][1] ** 2)
        mx[i].append(mx[i][2] ** 2)
    return mx, num


def korhen_check(f1, f2, q=0.05):
    # Get critical value
    q1 = q / f1
    critical = f.ppf(q=1 - q1, dfn=f2, dfd=(f1 - 1) * f2)
    return critical / (critical + f2 - 1)


def square_equation(val):
    m = val
    x_min = [x1_min, x2_min, x3_min]
    x_maxes = [x1_max, x2_max, x3_max]
    # matrix -> matrix with naturalized values
    matrix = naturalize_matrix(deepcopy(initial_matrix), m, x_min, x_maxes)
    norm_matrix = gen_matrix(deepcopy(initial_matrix), m)

    num_x = 11  # number of columns with x values x1 -> x3^2 + 1 for zero column
    solve_mx = [[] for i in range(num_x)]

    def mk_solve(mx1, res_mx):
        vector0 = [1]
        tmp = [sum(mx1[i][j] for i in range(N)) / N for j in range(num_x - 1)]
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

                for j in range(num_x):
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
                return vector

        res_mx[0] = vector0
        for i in range(1, len(res_mx)):
            res_mx[i] = mk_vector()

    mk_solve(matrix, solve_mx)
    row = len(matrix[0])
    y_abs = [abs_maker(matrix, i, row - 1, row - m - 1, -1, m) for i in range(N)]

    def res_vector():
        res = [sum([y_abs[i] for i in range(N)]) / N]
        ac = 0
        for i in range(num_x - 1):
            for j in range(N):
                ac += y_abs[j] * matrix[j][i]
                if j == N - 1:
                    res.append(ac)
                    ac = 0
        res = [res[0]] + [i / N for i in res[1:]]
        return res

    koef = res_vector()

    # useless now
    # k0 = sum([y_abs[i] for i in range(N)]) / N
    # k1 = sum([y_abs[i] * matrix[i][0] for i in range(N)]) / N
    # k2 = sum([y_abs[i] * matrix[i][1] for i in range(N)]) / N
    # k3 = sum([y_abs[i] * matrix[i][2] for i in range(N)]) / N
    # k4 = sum([y_abs[i] * matrix[i][3] for i in range(N)]) / N
    # k5 = sum([y_abs[i] * matrix[i][4] for i in range(N)]) / N
    # k6 = sum([y_abs[i] * matrix[i][5] for i in range(N)]) / N
    # k7 = sum([y_abs[i] * matrix[i][6] for i in range(N)]) / N
    # k8 = sum([y_abs[i] * matrix[i][7] for i in range(N)]) / N
    # k9 = sum([y_abs[i] * matrix[i][8] for i in range(N)]) / N
    # k10 = sum([y_abs[i] * matrix[i][9] for i in range(N)]) / N

    b = linalg.solve(solve_mx, koef)

    check1 = check(matrix, b, y_abs)

    f1 = m - 1
    f2 = N
    Gt = korhen_check(f1, f2)

    def s_y_maker():
        return [sum((matrix[i][j] - y_abs[i])**2 for j in range(row - 1, row -m -1, -1)) for i in range(N)]


    # s_y1 = sum([(matrix[0][j] - y_abs[0]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y2 = sum([(matrix[1][j] - y_abs[1]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y3 = sum([(matrix[2][j] - y_abs[2]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y4 = sum([(matrix[3][j] - y_abs[3]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y5 = sum([(matrix[4][j] - y_abs[4]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y6 = sum([(matrix[5][j] - y_abs[5]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y7 = sum([(matrix[6][j] - y_abs[6]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y8 = sum([(matrix[7][j] - y_abs[7]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y9 = sum([(matrix[8][j] - y_abs[8]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y10 = sum([(matrix[9][j] - y_abs[9]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y11 = sum([(matrix[10][j] - y_abs[10]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y12 = sum([(matrix[11][j] - y_abs[11]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y13 = sum([(matrix[12][j] - y_abs[12]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y14 = sum([(matrix[13][j] - y_abs[13]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    # s_y15 = sum([(matrix[14][j] - y_abs[14]) ** 2 for j in range(row - 1, row - m - 1, -1)]) / m
    #
    # s_y = (s_y1, s_y2, s_y3, s_y4, s_y5, s_y6, s_y7, s_y8, s_y9, s_y10, s_y11, s_y12, s_y13, s_y14, s_y15)
    s_y = s_y_maker()
    s_y = [i / m for i in s_y]

    gp = max(s_y) / sum(s_y)
    if gp < Gt:
        tb(matrix)
        print("\nCoefficients of regression are: \n", b)
        print("\nCheck is all is norm:")
        for i in check1.items():
            print(i)
        print()
        print('G_kr', Gt)
        print('G_p', gp)
        # print('ok')
    else:
        m = m + 1
        square_equation(m)

    s2_b = sum(s_y) / N
    s2_beta = s2_b / (N * m)
    s_beta = sqrt(s2_beta)

    def beta(lst, mx, n, nn):
        ac = 0
        res = []
        try_n = nn
        for i in range(try_n):
            for j in range(n):
                ac += lst[j] * mx[j][i]
                if j == n - 1:
                    res.append(ac / n)
                    ac = 0
        return res

    b_sum = beta(y_abs, matrix, N, num_x)
    f3 = f1 * f2
    t_kr = par_t.ppf(df=f3, q=(1 + p) / 2)
    d = 0
    t = [abs(b_sum[i]) / s_beta for i in range(num_x)]

    for i in range(len(t)):
        if t[i] < t_kr:
            t[i] = 0
        else:
            t[i] = 1
            d += 1

    b = [b[i] * t[i] for i in range(len(t))]
    print("\nStudent criteria:")
    print("Coefficients:")
    print(b)
    print("b{} is unnecessary".format([i for i in range(len(t)) if t[i] == 0]))

    check2 = check(mx=matrix, coef=b, y=y_abs)
    print("\nCheck if all is norm: ")
    for i in check2.items():
        print(i)
    print()

    y = [
        b[0] + b[1] * matrix[i][0] + b[2] * matrix[i][1] + b[3] * matrix[i][2] + b[4] * matrix[i][3] + b[5] * matrix[i][
            4] + b[6] * matrix[i][5] + b[7] *
        matrix[i][6] + b[8] * matrix[i][7] + b[9] * matrix[i][8] + b[10] * matrix[i][9] for i in range(N)]

    s2_ad = (m / (N - d)) * sum([(y[i] - y_abs[i]) ** 2 for i in range(N)])

    Fp = s2_ad / s2_beta
    f4 = N - d
    F_kr = f.ppf(dfn=f4, dfd=f3, q=0.95)

    print("Fp={}, F_kr={}".format(Fp, F_kr))

    if Fp < F_kr:
        print("Regression equation is adequate to original")
    else:
        print("Regression equation is inadequate to original")


square_equation(val=3)
