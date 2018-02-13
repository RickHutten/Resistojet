import math
from properties import *


def Vandenkerckhove(y):
    return math.sqrt(y) * (2 / (y + 1)) ** ((y + 1) / (2 * y - 2))


def ExitPressure(y, Ae, At, Pc):
    def f(Pe):
        return Vandenkerckhove(y) / math.sqrt(
            ((2 * y) / (y - 1)) * ((Pe / Pc) ** (2 / y)) * (1 - (Pe / Pc) ** ((y - 1) / y)))

    step = 0.25 * Pc
    direction = 0
    Pe = 0.5 * Pc
    while step / Pe > 0.00001:  # When the change in stepsize is less then 0.001%, Pe is very accurate
        difference = f(Pe) - Ae / At
        if difference > 0:
            if direction == -1:
                step /= 2
            if direction == 1:
                step *= 1.5
            direction = 1
            Pe += step
        elif difference < 0:
            if direction == -1:
                step *= 1.5
            if direction == 1:
                step /= 2
            direction = -1
            Pe -= step
        if Pe <= 0:  # Pe can't be negative
            Pe = 1e-20
    return Pe


def AreaThroat(mass_flow, gamma, Pc, Tc):
    return mass_flow * math.sqrt(R / M * Tc) / (Vandenkerckhove(gamma) * Pc)


def ChamberPressure(At, gamma, Tc):
    return mass_flow * math.sqrt(R * Tc) / (Vandenkerckhove(gamma) * At)


def ExhaustVelocity(gamma, Pe, Pc, Tc):
    v_max = (2 * gamma * R * Tc / ((gamma - 1) * M)) ** 0.5
    return v_max * (1 - (Pe / Pc) ** ((gamma - 1) / gamma)) ** 0.5


def RocketThrustEquation(mass_flow, Ue, Pe, Pa, Ae):
    return mass_flow * Ue + (Pe - Pa) * Ae
