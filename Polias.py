''' Gabriel de Oliveira Frais
	Henrique dos Santos Flor
	Wender James Ferreira Gomes'''

import math
import pandas as pd


def Cabos(CMT):
    tabela = pd.read_excel("Tabelas/Cabos_v2.xlsx")
    for i in range(len(tabela)):
        if CMT < tabela["CMT"][i]:
            d = tabela["mm"][i]
            m = tabela["Massa"][i]
            modelo = tabela["Modelo"][i]
            break
    return d, m, modelo


def MenorValor(lst):
    i = float("inf")
    x = 0
    for nr in lst:
        if nr < i:
            i = nr
            Indice = x
        x = x + 1
    return Indice


def MenorValorLista(lst):
    y = float("inf")
    for i in range(len(lst)):
        if lst[i][-1] < y:
            y = lst[i][-1]
            indice = i
    return indice


def Motor(Wt, Er, Ttm_min):  # Motor(Wt, Er, Ttm)
    tabela = pd.read_excel("Tabelas/Motores.xlsx")
    for i in range(len(tabela)):
        P = tabela["Potência KW"][i]
        RPM = tabela["RPM"][i]
        Peso = tabela["Peso"][i]
        Polos = tabela["Polos"][i]
        Modelo = tabela["Modelo"][i]
        P_W = P * 1000
        Ttm = (P_W / Wt) * Er
        # print("Ttm_min: ", Ttm_min)
        # print("Ttm: ", Ttm)

        # print("Ttw: ", Ttw)

        if Ttm > Ttm_min:
            break

    return Modelo, Polos, RPM, P, Peso, Ttm


# Valores Iniciais
W = 50  # tonf
W_lbf = W / 0.00044642857142857  # lbf
Vl = 20  # ft/min
Vl_ms = Vl * 0.00508  # m/s
H = 5  # m
Ne = 2
CS_cabo = 5
pi = math.pi
den = 7800  # densidade do metal
Sy_MPa = 1300  # MPa
Sy = 1300 * (10 ** 6)

DAF = 1.4 + (0.1 * math.sqrt(50 / W))

Valores = []
Massa_total = []
x = 0

for Np in range(5, 11):
    print(x)

    N = Np + 1
    Ep = 0.99 ** Np

    Pc = (W * DAF * CS_cabo) / (N * Ep)  # CARGA DE RUPTURA
    Pc_N = Pc * 9964.016384

    dc, mc, modelo_c = Cabos(Pc)  # diametro e massa da corda mm e kg/m
    dc_m = dc / 1000
    dp = 20 * dc  # mm
    dt = dp
    dt_m = dp / 1000  # m
    dp_m = dt_m

    Lc = H * N  # comprimento do cabo

    Vc = (Vl * N) / Ne  # velocidade do cabo
    Vc_m = (Vl_ms * N) / Ne

    Wt = (2 * Vc_m) / dt_m  # rad/s
    Nt = (60 * Vc_m) / (pi * dt_m)  # rpm

    massa_polia = ((pi * (dp_m ** 2)) / 4) * dc_m * den  # kg
    MT_polia = Np * massa_polia

    MT_cabo = Lc * mc

    # TAMBOR

    Nran = int(math.ceil((Lc / Ne) - (1 / (pi * dt_m))))  # numero de ranhuras

    p = 1.2 * dc  # passe de ranhura
    p_m = 1.2 * dc_m

    Ttw = (Pc_N * dt_m * Ne) / 2  # torque no tambor que vem da carga

    Lt = (Nran + 2) * Ne * p_m  # comprimento minimo do tambor

    dti_max = math.ceil(dt - ((7 / 8) * dc))

    I = (pi * (dt_m ** 4)) / 64  # momento de inercia
    J = I * 2  # momento polar de inercia

    Sy_raiz3 = Sy / math.sqrt(3)

    for i in range(dti_max, 0, -1):
        dti = i  # diametro interno do tambor
        dti_m = dti / 1000

        h = (dt_m / 2) - (dti_m / 2) - ((7 * dc_m) / 16)
        Volume_tambor = (((dt_m ** 2) - (dti_m ** 2)) * (pi / 4)) * Lt
        Pt = Volume_tambor * den  # Pt é o peso do tambor

        Ac = h * p_m
        Sc = Pc_N / Ac  # tensão de compressão  Sc <= Sy

        F = (Pc_N * Ne) + Pt
        Mf = (F * Lt) / 4  # momento fletor
        Smf = (Mf * (dt_m / 2)) / I  # Smf <= Sy
        Smt = (Ttw * (dt_m / 2)) / J  # Smt <= Sy / RAIZ(3)

        As = (dt_m ** 2 - dti_m ** 2) * pi / 4
        Ts = Pc / As  # Tensão na solda Ts <= Sy / RAIZ(3)    Pc_N e As_metro

        if h > 0 and Sc <= Sy and Smf <= Sy and Smt <= Sy_raiz3 and Ts <= Sy_raiz3:
            break

    valor_alto_motor = float('inf')

    for Nenge in range(4, 50):  # numero de engrenamentos
        Npe = int(math.floor(Np / 2))  # numero de polias entre o tambor e a polia equalizadora

        # norma AISI 6 / 91
        Kv = 1.38  # fator de tensão
        Ks = 1.40  # fator de serviço

        Ec = (0.97 ** Nenge) * (0.99 ** Npe)  # eficiencia de acoplamento

        Pm = (Vl * W_lbf * Ks * Kv) / (33000 * Ec)  # potencia do motor [hp] dps conveter pra cv e para W
        Pm_W = Pm * 745.6999

        Er = 0.97 ** Nenge

        Ttm_min = (Pm_W / Wt) * Er  # torque no tambor que vem do motor

        Valores_Motor = Motor(Wt, Er, Ttm_min)

        if Valores_Motor[4] < valor_alto_motor:
            valor_alto_motor = Valores_Motor[4]
            Nenge_final = Nenge
            Npe_final = Npe
            Ec_final = Ec
            Pm_W_final = Pm_W
            Er_final = Er
            Ttm_min_final = Ttm_min
            Modelo = Valores_Motor[0]
            Polos = Valores_Motor[1]
            Nm = Valores_Motor[2]
            Potencia_motor = Valores_Motor[3]
            Peso_motor = Valores_Motor[4]
            Ttm = Valores_Motor[5]

    i = Nm / Nt  # taxa de redução

    MT = MT_cabo + MT_polia + Pt + Peso_motor

    Resultado = [Np, Ep, Pc_N, dc, mc, modelo_c, dp, Lc, MT_cabo, Vc_m, Wt, Nt, massa_polia, MT_polia, Nran, p_m, Ttw, Lt, dti_max, I, J, dti, h, Volume_tambor, Pt, Ac, Sc, F, Mf, Smf, Smt, As, Ts, Nenge_final, Npe_final, Ec_final, Pm_W_final, Er_final, Ttm_min_final, Modelo, Polos, Nm, Potencia_motor, Peso_motor, Ttm, i, MT]

    Massa_total.append(MT)
    Valores.append(Resultado)
    x = x + 1

Menor = MenorValorLista(Valores)

Np, Ep, Pc_N, dc, mc, modelo_c, dp, Lc, MT_cabo, Vc_m, Wt, Nt, massa_polia, MT_polia, Nran, p_m, Ttw, Lt, dti, I, J, dti_m, h, Volume_tambor, Pt, Ac, Sc, F, Mf, Smf, Smt, As, Ts, Nenge, Npe, Ec, Pm_W, Er, Ttm_min, Modelo, Polos, Nm, Potencia_motor, Peso_motor, Ttm, i, MT = Valores[Menor]

df = pd.DataFrame(Valores, columns=['Numero de Polias', 'Ep', 'Carga de Ruptura [N]', 'Diametro do Cabo [mm]', 'Densidade do Cabo [kg/m]', 'Modelo do Cabo', 'Diametro da Polia/Tambor [mm]', 'Comprimento do Cabo [m]', 'Massado Cabo [kg]', 'Velocidade do Cabo [m/s]', 'Wt [rot/s]', 'Nt [RPM]', 'Massa da Polia [kg]', 'Soma de todas as Polias [kg]', 'Nran', 'Passo de Ranhura [m]', 'Ttw [Nm]', 'Comprimento do Tambor [m]', 'Diametro Interno Maximo do Tambor [mm]', 'Momento de Inercia [kg.m2]', 'Momento Polar de Inercia [kg.m2]', 'Diametro Interno do Tambor [mm]', 'Altura do Diametro Interno até o rasgo do enrolamento [m]', 'Volume do Tambor [m3]', 'Massa do Tambor [kg]', 'Ac', 'Sc', 'F', 'Mf', 'Smf', 'Smt', 'As', 'Ts', 'Número de Engrenamentos', 'Npe', 'Ec', 'Potência do Motor Necessário [W]', 'Er', 'Ttm Necessário [Nm]', 'Modelo do Motor', 'Polos', 'Rotação do Motor [RPM]', 'Potência do Motor [kW]', 'Massa do Motor [kg]', 'Ttm fornecido pelo motor [Nm]', 'Taxa de Redução', 'Massa Total do Projeto [kg]'])

with pd.option_context('display.max_columns', None,
    'display.precision', 3,
                       ):
    print(df.iloc[[Menor]])




