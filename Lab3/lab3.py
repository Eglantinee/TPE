import sys
import numpy as np
from scipy.stats import f, t

# 227 Variant

x1_min = -40
x1_max = 20
x2_min = 30
x2_max = 80
x3_min = -40
x3_max = 25
y_max = 200 + sum([x1_max, x2_max, x3_max]) / len([x1_max, x2_max, x3_max])
y_min = 200 + sum([x1_min, x2_min, x3_min]) / len([x1_min, x2_min, x3_min])

# Global variables

k_x = 4
m = 3
N = 4
s_y1 = 0
s_y2 = 0
s_y3 = 0
s_y4 = 0
f1 = 0
f2 = 0
b = []

# Create initial matrix

matrix = np.zeros((N, k_x + m))
norm_matrix = matrix[:]
val_x = [(1, -1, -1, -1),
         (1, -1, 1, 1),
         (1, 1, -1, 1),
         (1, 1, 1, -1),
         (1, -1, -1, 1),
         (1, -1, 1, -1),
         (1, 1, -1, -1),
         (1, 1, 1, 1)]

for i in range(N):
    for j in range(len(val_x[i])):
        matrix[i][j] = val_x[i][j]

while m < 24:
    if m == 23:
        sys.exit()
    for i in range(len(matrix)):
        for j in range(len(matrix[i]) - 1, len(matrix[i]) - m - 1, -1):
            matrix[i][j] = np.random.randint(y_min, y_max)

    # Matrix without x0 and min, max values => norm_matrix

    norm_matrix = np.delete(matrix, 0, axis=1)

    for i in range(len(norm_matrix)):
        for j in range(k_x - 1):
            norm_matrix[i][j] = globals()['x' + str(j + 1) + '_max'] if norm_matrix[i][j] == 1 \
                else globals()['x' + str(j + 1) + '_min']

    y1_abs = sum([norm_matrix[0][j] for j in range(len(norm_matrix[0]) - 1, len(norm_matrix[0]) - m - 1, -1)]) / m
    y2_abs = sum([norm_matrix[1][j] for j in range(len(norm_matrix[1]) - 1, len(norm_matrix[1]) - m - 1, -1)]) / m
    y3_abs = sum([norm_matrix[2][j] for j in range(len(norm_matrix[2]) - 1, len(norm_matrix[2]) - m - 1, -1)]) / m
    y4_abs = sum([norm_matrix[3][j] for j in range(len(norm_matrix[3]) - 1, len(norm_matrix[3]) - m - 1, -1)]) / m

    mx1 = sum([norm_matrix[i][0] for i in range(len(norm_matrix))]) / N
    mx2 = sum([norm_matrix[i][1] for i in range(len(norm_matrix))]) / N
    mx3 = sum([norm_matrix[i][2] for i in range(len(norm_matrix))]) / N

    my = sum([y1_abs, y2_abs, y3_abs, y4_abs]) / N

    a1 = sum([norm_matrix[i][0] * globals()['y' + str(i + 1) + '_abs'] for i in range(len(norm_matrix))]) / N
    a2 = sum([norm_matrix[i][1] * globals()['y' + str(i + 1) + '_abs'] for i in range(len(norm_matrix))]) / N
    a3 = sum([norm_matrix[i][2] * globals()['y' + str(i + 1) + '_abs'] for i in range(len(norm_matrix))]) / N

    a11 = sum([norm_matrix[i][0] ** 2 for i in range(len(norm_matrix))]) / N
    a22 = sum([norm_matrix[i][1] ** 2 for i in range(len(norm_matrix))]) / N
    a33 = sum([norm_matrix[i][2] ** 2 for i in range(len(norm_matrix))]) / N


    def prod(lst):
        return lst[0] * lst[1]


    a12 = a21 = sum([prod(norm_matrix[i][:2]) for i in range(len(norm_matrix))]) / N
    a13 = a31 = sum([prod(norm_matrix[i][:3:2]) for i in range(len(norm_matrix))]) / N
    a23 = a32 = sum([prod(norm_matrix[i][1:3]) for i in range(len(norm_matrix))]) / N

    b = np.linalg.solve([[1, mx1, mx2, mx3],
                         [mx1, a11, a21, a31],
                         [mx2, a12, a22, a32],
                         [mx3, a13, a23, a33]],
                        [my, a1, a2, a3])

    # Regression equation =>  y1 = b[0] + b[1] * x1 + b[2] * x2 + b[3] * x3

    def check1():
        # Confirm that y[j]_abs = y[j]
        res = dict()
        for it in range(len(norm_matrix)):
            res.update(
                {globals()['y' + str(it + 1) + '_abs']: b[0] + b[1] * norm_matrix[it][0] + b[2] * norm_matrix[it][1] +
                                                        b[3]
                                                        * norm_matrix[it][2]})
        return res


    check1_res = check1()

    print("Matrix1\n", matrix)
    print("Matrix2\n", norm_matrix)
    print("Coefficients:\t", b)
    print("Check results for truth\n", check1_res)

    # Korhen check:

    def korhen(ff1, ff2, q=0.05):
        # Get critical value
        q1 = q / ff1
        critical = f.ppf(q=1 - q1, dfn=ff2, dfd=(ff1 - 1) * ff2)

        return critical / (critical + ff2 - 1)


    s_y1 = sum([(matrix[0][i] - y1_abs) ** 2 for i in range(4, 7)]) / m
    s_y2 = sum([(matrix[1][i] - y2_abs) ** 2 for i in range(4, 7)]) / m
    s_y3 = sum([(matrix[2][i] - y3_abs) ** 2 for i in range(4, 7)]) / m
    s_y4 = sum([(matrix[3][i] - y4_abs) ** 2 for i in range(4, 7)]) / m

    gp = max(s_y1, s_y2, s_y3, s_y4) / sum([s_y1, s_y2, s_y3, s_y4])
    f1 = m - 1
    f2 = N
    g_kr = korhen(f1, f2)

    print("Gp={}, G_kr={}".format(gp, g_kr))

    if gp < g_kr:
        print("dispersion is uniform")
        break
    else:
        print("dispersion isn't uniform")
        m += 1

# Student criteria

s2_b = sum([s_y1, s_y2, s_y3, s_y4]) / N

s2_beta = s2_b / (N * m)
s_beta = np.sqrt(s2_beta)

# Student criteria: we want to find which coefficients are valuable:
#   to do it we compare t[i] with t_kr
#   but t[i] = b_[i]/ s_beta
#   so we find s2_b then s2_beta and then s_beta;  find b_[i] and then t[i]    <= using formulas

b_0 = sum([globals()['y' + str(i + 1) + '_abs'] * matrix[i][0] for i in range(len(matrix))]) / N
b_1 = sum([globals()['y' + str(i + 1) + '_abs'] * matrix[i][1] for i in range(len(matrix))]) / N
b_2 = sum([globals()['y' + str(i + 1) + '_abs'] * matrix[i][2] for i in range(len(matrix))]) / N
b_3 = sum([globals()['y' + str(i + 1) + '_abs'] * matrix[i][3] for i in range(len(matrix))]) / N

f3 = f1 * f2
t_kr = t.ppf(df=f3, q=(1 + 0.95) / 2)

d = 0
t = [abs(globals()['b_' + str(i)]) / s_beta for i in range(N)]

for i in range(len(t)):
    if t[i] < t_kr:
        t[i] = 0
    else:
        t[i] = 1
        d += 1

print("b{} is unnecessary".format([i for i in range(len(t)) if t[i] == 0]))

# Equations with valuable coefficients

y1 = b[0] * t[0] + b[1] * t[1] * norm_matrix[0][0] + b[2] * t[2] * norm_matrix[0][1] + b[3] * t[3] * norm_matrix[0][2]
y2 = b[0] * t[0] + b[1] * t[1] * norm_matrix[1][0] + b[2] * t[2] * norm_matrix[1][1] + b[3] * t[3] * norm_matrix[1][2]
y3 = b[0] * t[0] + b[1] * t[1] * norm_matrix[2][0] + b[2] * t[2] * norm_matrix[2][1] + b[3] * t[3] * norm_matrix[2][2]
y4 = b[0] * t[0] + b[1] * t[1] * norm_matrix[3][0] + b[2] * t[2] * norm_matrix[3][1] + b[3] * t[3] * norm_matrix[3][2]


# Fisher criteria

s2_ad = (m / (N - d)) * sum(
    [(globals()['y' + str(i + 1)] - globals()['y' + str(i + 1) + '_abs']) ** 2 for i in range(N)])

Fp = s2_ad / s2_beta
f4 = N - d
F_kr = f.ppf(dfn=f4, dfd=f3, q=0.95)

print("Fp={}, F_kr={}".format(Fp, F_kr))

if Fp < F_kr:
    print("Regression equation is adequate to original")
else:
    print("Regression equation is inadequate to original")

#sms