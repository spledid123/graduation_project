import iapws
from math import pi, log10, exp
import matplotlib.pyplot as plt
import numpy as np

m = 18e-3 / 6.02e23
def a(t):  
    temp_p = iapws._iapws._Sublimation_Pressure(t)
    return temp_p, temp_p * m ** 0.5 / iapws._iapws._Ice(t, temp_p)["rho"]/(2*1.38e-23*t*pi)**0.5 * (3600*24*365)


def b(r, t):
    return iapws._iapws._Sublimation_Pressure(t) * exp(2 * 0.138 * pi * 4/3 * (4e-10) ** 3 / (r*1.38e-23*t))
def c(t):
    E =  iapws._iapws._Sublimation_Pressure(t) * 1e6 / 2 / pi / 1.38e-23 / t / m
    th = ((iapws._Ice(t, 1e-9)['rho']) / 0.018) ** (2/3)
    return th / E


x = list(np.arange(50, 200, 1))     # 50-200, 间隔1K
#r = list(np.arange(1e-6, 1e-3, (1e-3 - 1e-6) / 100))
log_r = []
s_p = []    # 饱和蒸气压
sub_speed = []      # 升华速率
c = []  # 随曲率变化的饱和蒸汽压
t_num = 5
delta_t = 10
for i in range(150):
    temp_a = c(x[i])
    s_p.append(temp_a)

for i in range(t_num):
    temp_c = []
    for j in range(100):
        temp_c.append(b(r[j], 100 + i * delta_t))
    c.append(temp_c)
for i in range(t_num):
    temp_t = 100 + i * delta_t
    plt.plot(r, c[i], label=temp_t)
# plt.xscale("log")
# plt.yscale("log")
plt.title("sublimation_pressure-r我")
# plt.xlabel("r")
# plt.ylabel("pressure_pressure")
# plt.legend()
# plt.plot(x, s_p)
plt.xlim(50, 200)
plt.yscale("log")
# plt.title("sublimation_pressure-t")
# plt.xlabel("t")
# plt.ylabel("pressure_pressure")
#plt.plot(x, sub_speed)
plt.plot(x, s_p)
plt.show()
