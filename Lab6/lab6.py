import sys
from math import sqrt
from copy import deepcopy
from random import uniform
from random import randint
from myprettytable import table as tb
from scipy import linalg
from scipy.stats import f
from scipy.stats import t as par_t

x1_min = -40
x1_max = 20
x2_min = 30
x2_max = 80
x3_min = -40
x3_max = -25


def make_y(x1, x2, x3):
    return (1.7 + 2.8 * x1 + 9.7 * x2 +
            2.9 * x3 + 4.9 * x1 * x1 +
            0.2 * x2 * x2 + 0.9 * x3 * x3 + 3.2 * x1 * x2 +
            1.0 * x1 * x3 + 6.8 * x2 * x3 +
            2.8 * x1 * x2 * x3) / 10 + randint(0, 10) - 5


y_min = make_y(x1_min, x2_min, x3_min)
y_max = make_y(x1_max, x2_max, x3_max)

p = 0.95
N = 14

initial_matrix = [[-1, -1, -1],
                  [-1, -1, 1],
                  [-1, 1, -1],
                  [-1, 1, 1],
                  [1, -1, -1],
                  [1, -1, 1],
                  [1, 1, -1],
                  [1, 1, 1],
                  [-1.73, 0, 0],
                  [1.73, 0, 0],
                  [0, -1.73, 0],
                  [0, 1.73, 0],
                  [0, 0, -1.73],
                  [0, 0, 1.73]]


def fisher_criteria(matrix, m, y_abs, s2_beta, d, b, f3, num_x, N, step):
    y = [b[0] + sum([b[j] * matrix[i][j - 1] for j in range(1, num_x)]) for i in range(N)]
    if N - d == 0:
        N += 1
        matrix, norm_matrix = pre_run(m)
        equations(matrix, norm_matrix, m, num_x, step)
    s2_ad = (m / (N - d)) * sum([(y[i] - y_abs[i]) ** 2 for i in range(N)])
    Fp = s2_ad / s2_beta
    f4 = N - d
    F_kr = f.ppf(dfn=f4, dfd=f3, q=0.95)
    return Fp, F_kr


def student_criteria(matrix, num_x, y_abs, s_y, b, f3, m):
    s2_b = sum(s_y) / N
    s2_beta = s2_b / (N * m)
    s_beta = sqrt(s2_beta)

    b_sum = beta(y_abs, matrix, N, num_x)
    t_kr = par_t.ppf(df=f3, q=(1 + p) / 2)
    d = 0
    t = [abs(b_sum[i]) / s_beta for i in range(num_x)]

    for i, _value in enumerate(t):
        if t[i] < t_kr:
            t[i] = 0
        else:
            t[i] = 1
            d += 1

    b = [b[i] * t[i] for i in range(len(t))]
    return t, s2_beta, d


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


def mk_solve(mx1, res_mx, num_x):
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

        for _ in range(N):
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


def res_vector(matrix, y_abs, num_x):
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


def abs_maker(mx, row, start, stop, step, val):
    return sum([mx[row][j] for j in range(start, stop, step)]) / val


def my_check(mx, b, y_abs, num_x):
    res = dict()
    for i in range(N):
        res.update({y_abs[i]: b[0] + sum([b[j] * mx[i][j - 1] for j in range(1, num_x)])})
    return res


def s_y_maker(matrix, y_abs, row, m):
    return [sum((matrix[i][j] - y_abs[i]) ** 2 for j in range(row - 1, row - m - 1, -1)) for i in range(N)]


def add_y_vals(func):
    def inner(*args):
        returned_values = func(*args)
        mx = returned_values[0]
        m = returned_values[1]
        for i, _value in enumerate(mx):
            for _ in range(m):
                mx[i].append(round(uniform(y_min, y_max), ndigits=3))
        return mx

    return inner


@add_y_vals
def gen_matrix(mx, num):
    for i, _value in enumerate(mx):
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
    for i, _value in enumerate(mx):
        mx[i].append(mx[i][0] ** 2)
        mx[i].append(mx[i][1] ** 2)
        mx[i].append(mx[i][2] ** 2)
    for i, _value in enumerate(mx):  # add zero column for student criteria
        mx[i].insert(0, 1)
    return mx, num


@add_y_vals
def naturalize_matrix(mx, num, mins, maxes):
    for i, _value in enumerate(mx):
        for j in range(3):  # 3 means x1, x2, x3 only
            x_zero = (maxes[j] + mins[j]) / 2
            delta_x = maxes[j] - x_zero
            if mx[i][j] == -1:
                mx[i][j] = mins[j]
            elif mx[i][j] == 1:
                mx[i][j] = maxes[j]
            elif abs(mx[i][j]) == 1.73:
                mx[i][j] = mx[i][j] * delta_x + x_zero
            elif mx[i][j] == 0:
                mx[i][j] = x_zero

    for i, _value in enumerate(mx):
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
    for i, _value in enumerate(mx):
        mx[i].append(mx[i][0] ** 2)
        mx[i].append(mx[i][1] ** 2)
        mx[i].append(mx[i][2] ** 2)
    return mx, num


def korhen_check(f1, f2, q=0.05):
    # Get critical value
    q1 = q / f1
    critical = f.ppf(q=1 - q1, dfn=f2, dfd=(f1 - 1) * f2)
    return critical / (critical + f2 - 1)


def pre_run(m):
    x_min = [x1_min, x2_min, x3_min]
    x_maxes = [x1_max, x2_max, x3_max]

    matrix = naturalize_matrix(deepcopy(initial_matrix), m, x_min, x_maxes)
    norm_matrix = gen_matrix(deepcopy(initial_matrix), m)
    return matrix, norm_matrix


def equations(matrix, norm_matrix, val, num_x, step):
    m = val

    # num_x = 4  # number of columns with x values x1 -> x3^2 + 1 for zero column
    solve_mx = [[] for i in range(num_x)]

    row = len(matrix[0])
    y_abs = [abs_maker(matrix, i, row - 1, row - m - 1, -1, m) for i in range(N)]
    mk_solve(matrix, solve_mx, num_x)
    koef = res_vector(matrix, y_abs, num_x)

    b = linalg.solve(solve_mx, koef)

    check1 = my_check(matrix, b, y_abs, num_x)

    f1 = m - 1
    f2 = N
    Gt = korhen_check(f1, f2)

    s_y = s_y_maker(matrix, y_abs, row, m)
    s_y = [i / m for i in s_y]

    gp = max(s_y) / sum(s_y)

    if gp > Gt:
        m = m + 1
        matrix, norm_matrix = pre_run(m)
        equations(matrix, norm_matrix, m, num_x, step)

    check2 = my_check(matrix, b, y_abs, num_x)

    f3 = f1 * f2
    t, s2_beta, d = student_criteria(matrix, num_x, y_abs, s_y, b, f3, m)
    Fp, F_kr = fisher_criteria(matrix, m, y_abs, s2_beta, d, b, f3, num_x, N, step)

    if Fp < F_kr:
        tb(matrix)
        print("\nCoefficients of regression are: \n", b)
        print('G_kr', Gt)
        print('G_p', gp)
        print("\nStudent criteria:")
        print("Coefficients:")
        print(b)
        out = [i for i in range(len(t)) if t[i] == 0]
        if out:
            print("b{} is unnecessary".format([i for i in range(len(t)) if t[i] == 0]))
        else:
            print("all coefficients are necessary")
        print("Fp={}, F_kr={}".format(Fp, F_kr))
        print("Regression equation is adequate to original")
        sys.exit()
    else:
        step += 1
        if step == 1:
            num_x = 8
            equations(matrix, norm_matrix, 3, num_x, step)
        elif step == 2:
            num_x = 11
            equations(matrix, norm_matrix, 3, num_x, step)
        else:
            step = 0
            num_x = 4
            matrix, norm_matrix = pre_run(3)
            equations(matrix, norm_matrix, 3, num_x, step)


matrix, norm_matrix = pre_run(3)
equations(matrix, norm_matrix, 3, 4, 0)
