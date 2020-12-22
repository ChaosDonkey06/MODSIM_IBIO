import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import sympy as sym

c = "c"
a = "a"
k = 2


alphas = [0.01, 0.02, 0.03, 0.04, 0.05] #1, 2.1, 2.2, 3.1, 3.2
sigma_c = 1/3
gamma_c = 1/4

#params adults
sigma_a = 1/3
gamma_a = 1/4



def beta_val(coh):
    if coh == 1:
        B11 = 0.01
        B12 = 0.025
        B21 = 0.04
        B22 = 0.05
    elif coh == 2.1:
        B11 = 0.01
        B12 = 0.025
        B21 = 0.04
        B22 = 0.05
    elif coh == 2.21 or coh ==2.22:
        B11 = 0.01
        B12 = 0.025
        B21 = 0.04
        B22 = 0.05
    elif coh == 3.11 or coh ==3.12:
        B11 = 0.01
        B12 = 0.025
        B21 = 0.04
        B22 = 0.05
    elif coh == 3.21 or coh == 3.22 or coh == 3.23:
        B11 = 0.01
        B12 = 0.025
        B21 = 0.04
        B22 = 0.05

    beta = np.array([[B11, B12],[B21, B22]])
    return beta

def beta_t(coh, beta, t):
    if coh ==1 or coh ==2.1:
        if t%7==0 or t%7==6:
            return beta
        else:
            return beta*k
    elif coh == 2.21 or coh==3.11:
        if t%14<=5 and t%14>=1:
            return beta*k
        else:
            return beta
    elif coh == 2.22 or coh == 3.12:
        if t%14<=12 and t%14>=8:
            return beta*k
        else:
            return beta
    elif coh == 3.21:
        if t % 21 <= 5 and t % 21 >= 1:
            return beta * k
        else:
            return beta
    elif coh == 3.22:
        if t % 21 <= 12 and t % 21 >= 8:
            return beta * k
        else:
            return beta
    elif coh == 3.23:
        if t % 21 <= 19 and t % 21 >= 15:
            return beta * k
        else:
            return beta


def inf_force(popu, coh, alpha, t):
    betita = beta_t(coh, beta_val(coh), t)


    if popu == c:
        cc = alpha*(betita[0,0])
        ca = alpha*(betita[0,1])
        both = [cc, ca]
        return both

    if popu == a:
        ac = betita[1, 0]
        aa = betita[1, 1]
        both = [aa, ac]
        return both


def model(vars, t, param): #param = coh, alpha #vars =x, y
    x = vars[0:4]
    y = vars[4:8]
    coh = param[0]
    alpha = param [1]
    SC=x[0]
    EC = x[1]
    IC = x[2]
    RC = x[3]

    SA = y[0]
    EA = y[1]
    IA = y[2]
    RA = y[3]

    dSC= -inf_force(c, coh, alpha, t)[0]*SC*IC -inf_force(c, coh, alpha,  t)[1]*SC*IA
    dEC = inf_force(c, coh, alpha, t)[0]*SC*IC +inf_force(c, coh, alpha,  t)[1]*SC*IA - sigma_c*EC
    dIC = sigma_c * EC - gamma_c*IC
    dRC = gamma_c*IC
    systC = [dSC, dEC, dIC, dRC]

    dSA= -inf_force(a, coh, alpha, t)[0]*SA*IA -inf_force(a, coh, alpha,  t)[1]*SA*IC
    dEA = inf_force(a, coh, alpha, t)[0]*SA*IA +inf_force(a, coh, alpha, t)[1]*SA*IC - sigma_a*EA
    dIA = sigma_a * EA - gamma_a*IA
    dRA = gamma_a*IA
    systA = [dSA, dEA, dIA, dRA]
    systfin = systC + systA
    return systfin


#ESCENARIO 1 => 1 COHORTE

ts1 = np.linspace(0, 200, 1000)
x01 = [200, 70, 30, 0]
y01= [8, 4, 1, 0]
z01 = x01 + y01
sol1 = odeint(model, z01, ts1, args=([1, alphas[0]],))

plt.plot(ts1, sol1[:, 0], label="SC")
plt.plot(ts1, sol1[:, 1], label="EC")
plt.plot(ts1, sol1[:, 2], label="IC")
plt.plot(ts1, sol1[:, 3], label="RC")

plt.plot(ts1, sol1[:, 4], label="SA")
plt.plot(ts1, sol1[:, 5], label="EA")
plt.plot(ts1, sol1[:, 6], label="IA")
plt.plot(ts1, sol1[:, 7], label="RA")
plt.legend(loc="best")
plt.ylabel("Individuos")
plt.xlabel("Tiempo(días)")
plt.suptitle("Simulación para 1 solo cohorte")
plt.title("$CI: SC_0: 200, EC_0: 70, IC_0:30, RC_0:0$" " & " "$SA_0: 8, EA_0: 4, IA_0:1, RA_0:0$")
plt.show()

#ESCENARIO 2.1 => 2 COHORTES:ONLINE/PRESENCIAL

ts21 = np.linspace(0, 200, 1000)
x021 = [200, 40, 10, 0]
y021= [8, 4, 1, 0]
z021 = x021 + y021
sol21 = odeint(model, z021, ts21, args=([2.1, alphas[1]],))

plt.plot(ts1, sol21[:, 0], label="SC")
plt.plot(ts1, sol21[:, 1], label="EC")
plt.plot(ts1, sol21[:, 2], label="IC")
plt.plot(ts1, sol21[:, 3], label="RC")

plt.plot(ts1, sol21[:, 4], label="SA")
plt.plot(ts1, sol21[:, 5], label="EA")
plt.plot(ts1, sol21[:, 6], label="IA")
plt.plot(ts1, sol21[:, 7], label="RA")
plt.legend(loc="best")
plt.ylabel("Individuos")
plt.xlabel("Tiempo(días)")
plt.suptitle("Simulacion para 2 cohortes: 1 online/ 1 presencial")
plt.title("$CI: SC_0: 200, EC_0: 40, IC_0:10, RC_0:0$" " & " "$SA_0: 8, EA_0: 4, IA_0:1, RA_0:0$")
plt.show()

#ESCENARIO 2.2 => 2 COHORTES:ROTACION

ts22 = np.linspace(0, 200, 1000)
x022 = [100, 40, 10, 0]
y022= [8, 4, 1, 0]
z022 = x022 + y022
sol22_1 = odeint(model, z022, ts22,  args=([2.21, alphas[2]],))

plt.plot(ts1, sol22_1[:, 0], label="SC_G1")
plt.plot(ts1, sol22_1[:, 1], label="EC_G1")
plt.plot(ts1, sol22_1[:, 2], label="IC_G1")
plt.plot(ts1, sol22_1[:, 3], label="RC_G1")

plt.plot(ts1, sol22_1[:, 4], label="SA_G1")
plt.plot(ts1, sol22_1[:, 5], label="EA_G1")
plt.plot(ts1, sol22_1[:, 6], label="IA_G1")
plt.plot(ts1, sol22_1[:, 7], label="RA_G1")


sol22_2 = odeint(model, z022, ts22,  args=([2.22, alphas[2]],))

plt.plot(ts1, sol22_2[:, 0], label="SC_G2")
plt.plot(ts1, sol22_2[:, 1], label="EC_G2")
plt.plot(ts1, sol22_2[:, 2], label="IC_G2")
plt.plot(ts1, sol22_2[:, 3], label="RC_G2")

plt.plot(ts1, sol22_2[:, 4], label="SA_G2")
plt.plot(ts1, sol22_2[:, 5], label="EA_G2")
plt.plot(ts1, sol22_2[:, 6], label="IA_G2")
plt.plot(ts1, sol22_2[:, 7], label="RA_G2")

plt.legend(loc="best")
plt.ylabel("Individuos")
plt.xlabel("Tiempo(días)")
plt.suptitle("Simulacion para 2 cohortes: rotando")
plt.title("$CI: SC_0: 100, EC_0: 40, IC_0:10, RC_0:0$" " & " "$SA_0: 8, EA_0: 4, IA_0:1, RA_0:0$")
plt.show()


#ESCENARIO 3.1 => 3 COHORTES: 1 ONLINE/2 ROTANDO

ts31 = np.linspace(0, 200, 1000)
x031 = [80, 40, 5, 0]
y031 = [8, 4, 1, 0]
z031 = x031 + y031
sol31_1 = odeint(model, z031, ts31,  args=([3.11, alphas[3]],))

plt.plot(ts1, sol31_1[:, 0], label="SC_G1")
plt.plot(ts1, sol31_1[:, 1], label="EC_G1")
plt.plot(ts1, sol31_1[:, 2], label="IC_G1")
plt.plot(ts1, sol31_1[:, 3], label="RC_G1")

plt.plot(ts1, sol31_1[:, 4], label="SA_G1")
plt.plot(ts1, sol31_1[:, 5], label="EA_G1")
plt.plot(ts1, sol31_1[:, 6], label="IA_G1")
plt.plot(ts1, sol31_1[:, 7], label="RA_G1")

sol31_2 = odeint(model, z031, ts31,  args=([3.12, alphas[3]],))

plt.plot(ts1, sol31_2[:, 0], label="SC_G2")
plt.plot(ts1, sol31_2[:, 1], label="EC_G2")
plt.plot(ts1, sol31_2[:, 2], label="IC_G2")
plt.plot(ts1, sol31_2[:, 3], label="RC_G2")

plt.plot(ts1, sol31_2[:, 4], label="SA_G2")
plt.plot(ts1, sol31_2[:, 5], label="EA_G2")
plt.plot(ts1, sol31_2[:, 6], label="IA_G2")
plt.plot(ts1, sol31_2[:, 7], label="RA_G2")

plt.legend(loc="best")
plt.ylabel("Individuos")
plt.xlabel("Tiempo(días)")
plt.suptitle("Simulacion para 3 cohortes: 1 online/ 2 rotando")
plt.title("$CI: SC_0: 80, EC_0: 40, IC_0:5, RC_0:0$" " & " "$SA_0: 8, EA_0: 4, IA_0:1, RA_0:0$")
plt.show()

#ESCENARIO 3.2 => 3 COHORTES: 3 ROTANDO

ts = np.linspace(0, 200, 1000)
x032 = [70, 25, 5, 0]
y032= [8, 4, 1, 0]
z032 = x032 + y032
sol32_1 = odeint(model, z032, ts,  args=([3.21, alphas[4]],))

plt.plot(ts1, sol32_1[:, 0], label="SC_G1")
plt.plot(ts1, sol32_1[:, 1], label="EC_G1")
plt.plot(ts1, sol32_1[:, 2], label="IC_G1")
plt.plot(ts1, sol32_1[:, 3], label="RC_G1")

plt.plot(ts1, sol32_1[:, 4], label="SA_G1")
plt.plot(ts1, sol32_1[:, 5], label="EA_G1")
plt.plot(ts1, sol32_1[:, 6], label="IA_G1")
plt.plot(ts1, sol32_1[:, 7], label="RA_G1")

sol32_2 = odeint(model, z032, ts, args=([3.22, alphas[4]],))

plt.plot(ts1, sol32_2[:, 0], label="SC_G2")
plt.plot(ts1, sol32_2[:, 1], label="EC_G2")
plt.plot(ts1, sol32_2[:, 2], label="IC_G2")
plt.plot(ts1, sol32_2[:, 3], label="RC_G2")

plt.plot(ts1, sol32_2[:, 4], label="SA_G2")
plt.plot(ts1, sol32_2[:, 5], label="EA_G2")
plt.plot(ts1, sol32_2[:, 6], label="IA_G2")
plt.plot(ts1, sol32_2[:, 7], label="RA_G2")

sol32_3 = odeint(model, z032, ts, args=([3.23, alphas[4]],))

plt.plot(ts1, sol32_3[:, 0], label="SC_G3")
plt.plot(ts1, sol32_3[:, 1], label="EC_G3")
plt.plot(ts1, sol32_3[:, 2], label="IC_G3")
plt.plot(ts1, sol32_3[:, 3], label="RC_G3")

plt.plot(ts1, sol32_3[:, 4], label="SA_G3")
plt.plot(ts1, sol32_3[:, 5], label="EA_G3")
plt.plot(ts1, sol32_3[:, 6], label="IA_G3")
plt.plot(ts1, sol32_3[:, 7], label="RA_G3")

plt.legend(loc="best")
plt.ylabel("Individuos")
plt.xlabel("Tiempo(días)")
plt.suptitle("Simulacion para 3 cohortes: 3 rotando")
plt.title("$CI: SC_0: 70, EC_0: 25, IC_0:5, RC_0:0$" " & " "$SA_0: 8, EA_0: 4, IA_0:1, RA_0:0$")
plt.show()
