import control.matlab as c
import numpy as np
import matplotlib.pyplot as plt
import control as con
import sympy as sp

W0 = c.tf(1, 1)
W1 = c.tf([50, 2], [1, 0])
W2 = c.tf(3, [0.2, 1])
W3 = c.tf(0.3, [4, 1])

# s = sp.symbols("s", real=True)
# для упрощенной
# W0 = 1
# W1 = (50*s+2)/s
# W2 = 3/(0.2*s+1)
# W3 = 0.3/(4*s+1)

# W1 = c.tf([50, 2], [1, 14]) #пример неустойчивой системы
# W2 = c.tf([1, 25], [2, 1])
# W3 = c.tf(15, [0.1, 1])

W4 = W1 * W2
W5 = W0.feedback(W1,sign=-1) #sign = -1 ->  indicates negative feedback
# W5 = 1/(1+W1)
W6 = W4.feedback(1 / W1,sign=-1)
# W6 = W4/(1+W4/W1)
W7 = 1 - (W3 / W1)
W8 = W7 * W6 * W5
W9 = W3 * W5
W10 = W9 + W8
W = W10.feedback(W0,sign=1) #sign = -1 ->  indicates positive feedback
# W=W10/(1-W10)
# W1_num = sp.sympify(sp.factor(W))
# W_upr =2.63484675507157*c.tf([1,0.29047095290471,0.00999900009999],[0.179151664178449,1,0.272679371212786,0.00878194432247299])
print(W)
# W_upr = c.tf([204832.,21664095.,1201854.,210754.,11619.,257,2.], [256.,215180.,2247575.,1235015.,214670.,11758.,258.,2.])
t = np.linspace(0, stop=10, num=1000)
plt.figure(1)
y1, t1 = c.step(W, t)
lines = [y1]
lines[0] = plt.plot(t, y1, "r")
plt.title('Переходная характеристика', fontsize=10)
plt.ylabel('h')
plt.xlabel('t, c')
plt.grid()
plt.show()

def nyquis(W):
    plt.figure(1)
    con.nyquist_plot(W)
    plt.grid()
    plt.show()

def pzmap(W):
    f = con.pzmap(W)
    plt.title("Расположение нулей на комплексной плоскости")
    plt.grid(True)
    plt.show()
    g = f[0]  # возвращает список всех полюсов - корней знаменателя
    D = True
    j = 0
    G = False
    while D and j < len(g):
        r = g[j].real
        if float(r) <= 0.000000:
            print("Левый корень = " + str(r))
            D = True
            if float(r) == 0.0:
                G = True
        else:
            D = False
            print("Первый правый  корень = " + str(r))
        j = j + 1
    if D and G:
        print("Система на границе устойчивости, есть нулевые корни")
    elif D:
        print("Система устойчива, все корни левые")
    else:
        print("Cистема не устойчива, есть хотя бы один правый корень")

def mihalych(W1):
    P2 = (W1.den[0])
    a = P2[0]  # коэфициенты
    plt.figure(1)
    w = sp.symbols("w", real=True)  # Годограф Михайлова
    z = 0
    for i in range(0,len(a),1):
        z = z + a[i]*(1j*w)**((len(a)-1)-i)
    z = sp.factor(z)
    zR = sp.re(z)
    zIm = sp.im(z)
    x = [zR.subs({w: q}) for q in np.arange(0, 100, 0.01)]
    y = [zIm.subs({w: q}) for q in np.arange(0, 100, 0.01)]
    plt.title("Годограф Михайлова")
    # plt.axis([-500000, 5000000, -5000000, 50000000])
    plt.plot(x, y)
    plt.grid()
    plt.show()

def matrica_from_spisok(z,k): #k - kolichestvo strok=stolbcov
    matr = ''
    z=z
    k=k
    h=1
    for i in range(0, len(z), 1):
        if ((h*k)-i) == 0 and i!=0 and i!=len(z)-1: # последний?
            matr = matr + ';' + str(z[i])+' '
            h=h+1
        else:
            if i == len(z) - 1: # проверка на самый посл элемент
                matr = matr + str(z[i])
            else:
                matr = matr + str(z[i]) + ' '
    Mat = np.matrix(matr)
    return Mat

def matrica_spisok(z1,z2,k):

    """ к - количество столбцов и строк в матрице, ранг матрицы
    z1 - список четных коэф
    z2 - список нечетных коэф
    """
    z =[]
    h1 = 0
    h2 = 0
    for g in range(0, k, 1):
        index_nach_nech = h1
        index_nach_ch = h2
        index_konec_nech = k - len(z1) - h1
        index_konec_ch = k - len(z2) - h2
        if g % 2 == 0:  # no четный
            for i in range(0, index_nach_nech, 1):
                z.append("0")
            for i in range(0, len(z1), 1):
                z.append(z1[i])
            for i in range(0, index_konec_nech, 1):
                z.append("0")
            h1 = h1 + 1
        else:
            for i in range(0, index_nach_ch, 1):
                z.append("0")
            for i in range(0, len(z2), 1):
                z.append(z2[i])
            for i in range(0, index_konec_ch, 1):
                z.append("0")
            h2 = h2 + 1
    return z

def gurych(W1):
    T = False
    P2 = (W1.den[0])
    a = P2[0]  # коэфициенты
    z1 = []
    z2 = []
    for i in range(0, len(a), 1):
        if (i % 2 != 0):
            z1.append(str(a[i]))
        else:
            z2.append(str(a[i]))

    # print(z1) #нечетные коэф
    # print(z2) #четные коэф
    kol_str_stol = len(a) - 1  # количество строк и столбцов
    # print(kol_str_stol)
    z = matrica_spisok(z1,z2,kol_str_stol) #создаем матрицу n порядка
    Mat_gl = matrica_from_spisok(z,kol_str_stol)
    # print(Mat_gl)
    Del_mat = np.linalg.det(Mat_gl)
    # print("определитель главной Matrica del = " +str(Del_mat))

    """ необходимо рассмотреть определители главных миноров матрицы"""
    Del_minors = []
    for i in range(1,kol_str_stol,1):
        Mat__minor_gl = Mat_gl[:i, :i]
        Del_mat_minor = np.linalg.det(Mat__minor_gl)
        Del_minors.append(Del_mat_minor)

    """Проверка определителей главных миноров"""
    g = 1
    for i in Del_minors:
        g = g*i

    if (g*Del_mat > 0):
        T = True
        print("Ситема устойчива")
    elif (g*Del_mat == 0):
        print("Система на границе устойчивости устойчива")
    else:
        print("Ситема не устойчива")
    return T

W_raz = W10.feedback(W0,sign=-1)
nyquis(W_raz)
pzmap(W)
mihalych(W)
gurych(W)