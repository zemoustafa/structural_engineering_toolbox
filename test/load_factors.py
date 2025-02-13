import numpy as np

# -------------------------------------------------------------------------------------------------------------
# STRESS IN WALLS


def EQCompStress(N_G, N_Q, N_RS, M_G, M_Q, M_RS, tw, Lw):
    Awall = tw * Lw
    Zwall = tw * np.power(Lw, 2) / 6
    eq_comp = (N_G + 0.3 * N_Q + N_RS) * 1000 / Awall + (
        M_G + 0.3 * M_Q + M_RS
    ) * 1000000 / Zwall
    return eq_comp


def EQTensStress(N_G, N_Q, N_RS, M_G, M_Q, M_RS, tw, Lw):
    Awall = tw * Lw
    Zwall = tw * np.power(Lw, 2) / 6
    eq_tens = (N_G + 0.3 * N_Q - N_RS) * 1000 / Awall + (
        M_G + 0.3 * M_Q - M_RS
    ) * 1000000 / Zwall
    return eq_tens


# wind X compression
def WindCompStress1(N_G, N_Q, N_WX, M_G, M_Q, M_WX, tw, Lw):
    Awall = tw * Lw
    Zwall = tw * np.power(Lw, 2) / 6
    wind = (1.2 * N_G + 0.4 * N_Q + N_WX) * 1000 / Awall + (
        1.2 * M_G + 0.4 * M_Q + M_WX
    ) * 1000000 / Zwall
    return wind


# wind X tension
def WindTensStress1(N_G, N_Q, N_WX, M_G, M_Q, M_WX, tw, Lw):
    Awall = tw * Lw
    Zwall = tw * np.power(Lw, 2) / 6
    wind = (1.2 * N_G + 0.4 * N_Q - N_WX) * 1000 / Awall + (
        1.2 * M_G + 0.4 * M_Q - M_WX
    ) * 1000000 / Zwall
    return wind


# wind Y compression
def WindCompStress2(N_G, N_Q, N_WY, M_G, M_Q, M_WY, tw, Lw):
    Awall = tw * Lw
    Zwall = tw * np.power(Lw, 2) / 6
    wind = (1.2 * N_G + 0.4 * N_Q + N_WY) * 1000 / Awall + (
        1.2 * M_G + 0.4 * M_Q + M_WY
    ) * 1000000 / Zwall
    return wind


# wind Y tension
def WindTensStress2(N_G, N_Q, N_WY, M_G, M_Q, M_WY, tw, Lw):
    Awall = tw * Lw
    Zwall = tw * np.power(Lw, 2) / 6
    wind = (1.2 * N_G + 0.4 * N_Q - N_WY) * 1000 / Awall + (
        1.2 * M_G + 0.4 * M_Q - M_WY
    ) * 1000000 / Zwall
    return wind


# -------------------------------------------------------------------------------------------------------------
# SHEAR IN WALLS


def EQShear(V_G, V_Q, V_RS):
    shear = V_G + 0.3 * V_Q + np.abs(V_RS)
    return shear


def WXShear(V_G, V_Q, V_WX):
    shear = V_G + 0.4 * V_Q + np.abs(V_WX)
    return shear


def WYShear(V_G, V_Q, V_WY):
    shear = V_G + 0.4 * V_Q + np.abs(V_WY)
    return shear
