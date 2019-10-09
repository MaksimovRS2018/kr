import control.matlab as c
import numpy as np
import matplotlib.pyplot as plt
import control as con
import sympy as sp

# W0 = c.tf(1, 1)
# W1 = c.tf([50, 2], [1, 0])
# W2 = c.tf(3, [0.2, 1])
# W3 = c.tf(0.3, [4, 1])

s = sp.symbols("s", real=True)
# для упрощенной
W0 = 1
W1 = (50*s+2)/s
W2 = 3/(0.2*s+1)
W3 = 0.3/(4*s+1)

# W1 = c.tf([50, 2], [1, 14]) #пример неустойчивой системы
# W2 = c.tf([1, 25], [2, 1])
# W3 = c.tf(15, [0.1, 1])

W4 = W1 * W2
# W5 = W0.feedback(W1,sign=-1) #sign = -1 ->  indicates negative feedback
W5 = 1/(1+W1)
# W6 = W4.feedback(1 / W1,sign=-1)
W6 = W4/(1+W4/W1)
W7 = 1 - (W3 / W1)
W8 = W7 * W6 * W5
W9 = W3 * W5
W10 = W9 + W8
# W = W10.feedback(W0,sign=1) #sign = -1 ->  indicates positive feedback
W=W10/(1-W10)
W1 = sp.sympify(sp.factor(W))
W_upr =2.63484675507157*c.tf([1,0.29047095290471,0.00999900009999],[0.179151664178449,1,0.272679371212786,0.00878194432247299])
print(W_upr)

t = np.linspace(0, stop=10, num=1000)
plt.figure(1)
y1, t1 = c.step(W_upr, t)
lines = [y1]
lines[0] = plt.plot(t, y1, "r")
plt.title('Переходная характеристика', fontsize=10)
plt.ylabel('h')
plt.xlabel('t, c')
plt.grid()
plt.show()

# f = con.root_locus(W_upr)
# g = f[0]  # возвращает список всех корней, 1 список - корни
# D = True
# u = 0
# G = False
# while D and u < len(g):  # перебор списка
#     h = g[u]
#     j = 0
#     while D and j < len(h):
#         r = h[j].real
#         if float(r) <= 0.000000:
#             print("Левый корень = " + str(r))
#             D = True
#             if float(r)==0.0:
#                 G=True
#         else:
#             D = False
#             print("Первый правый  корень = "+ str(r))
#         j=j+1
#     u+=1
# if D and G:
#     print("Система на границе устойчивости, есть нулевые корни")
# elif D:
#     print("Система устойчива, все корни левые")
# else:
#     print("Cистема не устойчива, есть хотя бы один правый корень")
