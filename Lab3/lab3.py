from statistics import mean
import numpy as np

x1_min = -25
x1_max = 75
x2_min = 5
x2_max = 40
x3_min = 15
x3_max = 25
# y_max = 200 + mean([x1_max, x2_max, x3_max])
# y_min = 200 + mean([x1_min, x2_min, x3_min])
y_max = 20
y_min = 10
m = 3
N = 4

matrix = np.zeros((4, 7))
val_x = [(1, -1, -1, -1), (1, -1, 1, 1), (1, 1, -1, 1), (1, 1, 1, -1)]

for i in range(len(val_x)):
    for j in range(len(val_x[i])):
        matrix[i][j] = val_x[i][j]

for i in range(len(matrix)):
    for j in range(len(matrix[i]) - 1, -1, -1):
        if not matrix[i][j]:
            matrix[i][j] = np.random.randint(y_min, y_max)


matrix = [[1, -1,  -1,  -1,  15, 18,  16],
        [1, -1,  1,  1,  10,  19,  13],
        [1,  1,   -1,  1,  11,  14,  12],
        [ 1, 1,  1,  -1,  16,  19,  16]]

norm_matrix = np.delete(matrix, 0, axis=1)

for i in range(len(norm_matrix)):
    for j in range(3):
        norm_matrix[i][j] = globals()['x' + str(j + 1) + '_max'] if norm_matrix[i][j] == 1 \
            else globals()['x' + str(j + 1) + '_min']

y1_abs = sum([norm_matrix[0][j] for j in range(len(norm_matrix[0]) - 1, 2, -1) ]) / 3
y2_abs = sum([norm_matrix[1][j] for j in range(len(norm_matrix[1]) - 1, 2, -1) ]) / 3
y3_abs = sum([norm_matrix[2][j] for j in range(len(norm_matrix[2]) - 1, 2, -1) ]) / 3
y4_abs = sum([norm_matrix[3][j] for j in range(len(norm_matrix[3]) - 1, 2, -1) ]) / 3

mx1 = sum([norm_matrix[i][0] for i in range(len(norm_matrix))]) / 4
mx2 = sum([norm_matrix[i][1] for i in range(len(norm_matrix))]) / 4
mx3 = sum([norm_matrix[i][2] for i in range(len(norm_matrix))]) / 4

my = sum([y1_abs, y2_abs, y3_abs, y4_abs]) / 4

a1 = sum([norm_matrix[i][0] * globals()['y' + str(i + 1) + '_abs'] for i in range(len(norm_matrix))]) / 4
a2 = sum([norm_matrix[i][1] * globals()['y' + str(i + 1) + '_abs'] for i in range(len(norm_matrix))]) / 4
a3 = sum([norm_matrix[i][2] * globals()['y' + str(i + 1) + '_abs'] for i in range(len(norm_matrix))]) / 4

a11 = sum([norm_matrix[i][0] ** 2 for i in range(len(norm_matrix))]) / 4
a22 = sum([norm_matrix[i][1] ** 2 for i in range(len(norm_matrix))]) / 4
a33 = sum([norm_matrix[i][2] ** 2 for i in range(len(norm_matrix))]) / 4


def prod(lst):
    return lst[0] * lst[1]


a12 = a21 = sum([prod(norm_matrix[i][:2]) for i in range(len(norm_matrix))]) / 4
a13 = a31 = sum([prod(norm_matrix[i][:3:2]) for i in range(len(norm_matrix))]) / 4
a23 = a32 = sum([prod(norm_matrix[i][1:3]) for i in range(len(norm_matrix))]) / 4

b = np.linalg.solve([[1, mx1, mx2, mx3],
                     [mx1, a11, a21, a31],
                     [mx2, a12, a22, a32],
                     [mx3, a13, a23, a33]],
                    [my, a1, a2, a3])


# y1 = b[0] + b[1] * x1 + b[2] * x2 + b[3] * x3


def check1():
    res = dict()
    for it in range(len(norm_matrix)):
        res.update(
            {globals()['y' + str(it + 1) + '_abs']: b[0] + b[1] * norm_matrix[it][0] + b[2] * norm_matrix[it][1] + b[3]
                                                    * norm_matrix[it][2]})
    # print(res)


check1()

# Korhen

s_y1 = sum([(matrix[0][i] - y1_abs) ** 2 for i in range(4, 7)]) / 3
s_y2 = sum([(matrix[1][i] - y2_abs) ** 2 for i in range(4, 7)]) / 3
s_y3 = sum([(matrix[2][i] - y3_abs) ** 2 for i in range(4, 7)]) / 3
s_y4 = sum([(matrix[3][i] - y4_abs) ** 2 for i in range(4, 7)]) / 3

gp = max(s_y1, s_y2, s_y3, s_y4) / sum([s_y1, s_y2, s_y3, s_y4])
# print(gp)

# Student

s2_b = sum([s_y1, s_y2, s_y3, s_y4]) / N

s2_beta = s2_b / (N * m)
s_beta = np.sqrt(s2_beta)

b_0 = sum([globals()['y' + str(i + 1) + '_abs'] * matrix[i][0] for i in range(len(matrix))]) / N
b_1 = sum([globals()['y' + str(i + 1) + '_abs'] * matrix[i][1] for i in range(len(matrix))]) / N
b_2 = sum([globals()['y' + str(i + 1) + '_abs'] * matrix[i][2] for i in range(len(matrix))]) / N
b_3 = sum([globals()['y' + str(i + 1) + '_abs'] * matrix[i][3] for i in range(len(matrix))]) / N

t_kr = 2.306
d = 0
t = [abs(globals()['b_' + str(i)]) / s_beta for i in range(4)]
print(t)
for i in range(len(t)):
    if t[i] < t_kr:
        t[i] = 0
    else:
        t[i] = 1
        d += 1

print(t)

print('enter', b[0])
y1 = b[0] * t[0] + b[1] * t[1] * norm_matrix[0][0] + b[2] * t[2] * norm_matrix[0][1] + b[3] * t[3] * norm_matrix[0][2]
y2 = b[0] * t[0] + b[1] * t[1] * norm_matrix[1][0] + b[2] * t[2] * norm_matrix[1][1] + b[3] * t[3] * norm_matrix[1][2]
y3 = b[0] * t[0] + b[1] * t[1] * norm_matrix[2][0] + b[2] * t[2] * norm_matrix[2][1] + b[3] * t[3] * norm_matrix[2][2]
y4 = b[0] * t[0] + b[1] * t[1] * norm_matrix[3][0] + b[2] * t[2] * norm_matrix[3][1] + b[3] * t[3] * norm_matrix[3][2]

f_4 = N - d

s2_ad = (m / (N - d)) * sum(
    [(globals()['y' + str(i + 1)] - globals()['y' + str(i + 1) + '_abs']) ** 2 for i in range(4)])

Fp = s2_ad / s2_beta
# print(Fp)

# print(s2_ad)
# print(t)
# print(b_0)
# print(t)
# print(norm_matrix)
# print(matrix)
# print(norm_matrix)
# print(y1_abs, y2_abs, y3_abs, y4_abs)
# print(my)
# print(a1, a2, a3)
# print(mx1, mx2, mx3)
# print(a11, a22, a33)
# print(a12, a13, a23)
# print(b)
# print(s_y1, s_y2, s_y3, s_y4)
# print(gp)
# print(s2_b)
# print(s2_beta)
# print(b_0, b_1, b_2, b_3)
print(y1, y2, y3, y4)
print(s2_ad)
print(Fp)