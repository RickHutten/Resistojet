from ChemicalProperties import *
from functions import *
from annotations import *
import convert
import numpy as np


def performance_by_chamber_pressure_area_ratio(Pc, AreaRatio, print_results=True):
    data = reader.FileReader(reader.Files.isothermal_water_500_K)
    gamma = get_specific_heat_ratio(Pc, data)
    At = AreaThroat(gamma, Pc, Tc)
    Ae = AreaRatio * At
    Pe = ExitPressure(gamma, Ae, At, Pc)
    Dt = 2 * (At / math.pi) ** 0.5
    De = 2 * (Ae / math.pi) ** 0.5
    Ue = ExhaustVelocity(gamma, Pe, Pc, Tc)
    Thrust = RocketThrustEquation(Ue, Pe, Pa, Ae)
    c_star = Pc * At / mass_flow
    CF = Thrust / (Pc * At)
    Isp = CF * c_star / g0

    if print_results:
        print "Thrust:", Thrust, "N"
        print "Isp:", Isp, "s"
        print "Specific heat ratio:", gamma
        print "Exhaust Velocity:", Ue, "m/s"
        print "c*:", c_star, "m/s"
        print "Thrust coefficient:", CF
        print "Chamber Temperature:", Tc, "K"
        print "Chamber pressure:", Pc / 100000., "bar"
        print "Exit Pressure:", Pe / 100000., "bar"
        print "Throat diameter:", Dt * 1000, "mm"
        print "Exit diameter:", De * 1000, "mm"
        print "Mass flow: ", mass_flow * 1E6, "mg/s"
    return Dt, Isp

def f(a, b):
    return round(100 * (b / float(a)) - 100, 2)

def performance_by_throat_diameter_area_ratio(power, mass_flow, Dt, De=None, AreaRatio=None, print_results=True):
    if not De and not AreaRatio:
        raise StandardError("Either 'De' or 'AreaRatio' has to be defined!")
    if De is not None and AreaRatio is not None:
        raise StandardError("Only 'De' OR 'AreaRatio' may be defined!")

    enthalpy = (power / mass_flow) / 1000. + H_sto  # Enthalpy of the propellant in the chamber

    def calculate_pc(T):
        # Guess a chamber pressure
        Pc = convert.bar_to_pa(1)  # 1 bar
        step = convert.bar_to_pa(0.5)  # 0.5 bar
        direction = 0
        while step / Pc > 0.0001:  # When the change in stepsize is less then 0.01%, Pc is very accurate
            gamma = get_specific_heat_ratio_pressure_temperature(Pc, T)
            At = AreaThroat(mass_flow, gamma, Pc, T)
            temp_Dt = convert.area_to_diameter(At)
            difference = temp_Dt - Dt

            if difference > 0:
                if direction == -1:
                    step /= 2
                if direction == 1:
                    step *= 1.5
                direction = 1
                Pc += step
            elif difference < 0:
                if direction == -1:
                    step *= 1.5
                if direction == 1:
                    step /= 2
                direction = -1
                Pc -= step
            if Pc <= 0:  # Pe can't be negative
                Pc = 1e-20
        return Pc

    Tc = Pc = 0
    H_min = 1E100
    for t in range(temperature_min, temperature_max):
        Pc_temp = calculate_pc(t)
        if Pc_temp > convert.bar_to_pa(pressure_max):  # Can't calculate pressures above pressure_max
            continue
        enthalpy_temp = get_enthalpy_pressure_temperature(Pc_temp, t)
        enthalpy_diff = abs(enthalpy - enthalpy_temp)
        if enthalpy_diff < H_min:
            Pc = Pc_temp
            Tc = t
            H_min = enthalpy_diff
    if Pc == 0: raise ValueError("Pressure too high!")
    if Tc == temperature_min: raise ValueError("Temperature too low!")
    if Tc == temperature_max - 1: raise ValueError("Temperature too high!")

    gamma = get_specific_heat_ratio_pressure_temperature(Pc, Tc)
    At = convert.diameter_to_area(Dt)
    if De is not None:
        Ae = convert.diameter_to_area(De)
    elif AreaRatio is not None:
        Ae = AreaRatio * At
        De = convert.area_to_diameter(Ae)

    Pe = ExitPressure(gamma, Ae, At, Pc)
    Ue = ExhaustVelocity(gamma, Pe, Pc, Tc)
    Thrust = RocketThrustEquation(mass_flow, Ue, Pe, Pa, Ae)
    c_star = Pc * At / mass_flow
    CF = Thrust / (Pc * At)
    Isp = CF * c_star / g0

    if print_results:
        print "Thrust:", Thrust * 1000, "mN", '\t\t', f(2.31386612509, Thrust * 1000)
        print "Isp:", Isp, "s", '\t\t', f(131.082600134, Isp)
        print "Specific Heat Ratio:", gamma
        print "Exhaust Velocity:", Ue, "m/s"
        print "Effective Exhaust Velocity", Thrust / mass_flow, "m/s"
        print "c*:", c_star, "m/s", '\t\t', f(691.043939954, c_star)
        print "Thrust Coefficient:", CF, '\t\t', f(1.86020179945, CF)
        print "Chamber Temperature:", Tc, "K", '\t\t', f(465, Tc)
        print "Chamber Pressure:", convert.pa_to_bar(Pc), "bar", '\t\t', f(0.703891577199, convert.pa_to_bar(Pc))
        print "Exit Pressure:", convert.pa_to_bar(Pe), "bar"
        print "Pressure Ratio: 1/" + str(round(Pc / Pe, 2))
        print "Throat Diameter:", Dt * 1000, "mm"
        print "Exit Diameter:", De * 1000, "mm"
        print "Area Ratio: 1/" + str(Ae / At)
        print "Mass flow: ", mass_flow * 1E6, "mg/s"

    return Thrust * 1000, Isp, gamma, convert.pa_to_bar(Pc)


performance_by_throat_diameter_area_ratio(power=5, mass_flow=1.8E-6, Dt=0.15E-3, De=3.0E-3, print_results=True)


tlist = []
ilist = []
tlim = None
ilim = None
ylabel1 = None
ylabel2 = None

""" Power change """
# xlist = np.arange(4.78, 5.48, 0.01)
# for p in xlist:
#     t, i, _, _ = performance_by_throat_diameter_area_ratio(power=p, mass_flow=1.8E-6, Dt=0.15E-3, AreaRatio=400,
#                                                            print_results=False)
#     tlist.append(t)
#     ilist.append(i)
#     xlabel = "$P_{in}$ [W]"
#     title = 'Performance as a function of power input to the propellant'
#     tlim = [2.0, 2.8]
#     ilim = [120, 155]

""" Mass flow change """
# xlist = np.arange(1.645, 1.89, 0.005)
# for p in xlist:
#     t, i, _, _ = performance_by_throat_diameter_area_ratio(power=5, mass_flow=p * 1E-6, Dt=0.15E-3, AreaRatio=400,
#                                                            print_results=False)
#     tlist.append(t)
#     ilist.append(i)
#     xlabel = "Mass flow [mg/s]"
#     title = 'Performance as a function of mass flow'
#     tlim = [2.0, 2.8]
#     ilim = [120, 155]

""" Dt change, constant AreaRatio """
# xlist = np.arange(0.04, 0.3, 0.01)
# for p in xlist:
#     t, i, _, _ = performance_by_throat_diameter_area_ratio(power=5, mass_flow=1.8E-6, Dt=p * 1E-3, AreaRatio=400,
#                                                            print_results=False)
#     tlist.append(t)
#     ilist.append(i)
#     xlabel = "$D_t$ [mm]"
#     title = 'Performance as a function of throat diameter, with $\\frac{A_e}{A_t} = 400$'
#     tlim = [2.24, 2.34]
#     ilim = [126, 132]

""" Area Ratio change, constant Dt """
# xlist = range(2, 20, 2) + range(20, 100, 10) + range(100, 501, 25)
# for p in xlist:
#     t, i, _, _ = performance_by_throat_diameter_area_ratio(power=5, mass_flow=1.8E-6, Dt=0.15E-3, AreaRatio=p,
#                                                            print_results=False)
#     tlist.append(t)
#     ilist.append(i)
#     xlabel = "Area Ratio, $\\frac{A_e}{A_t}$"
#     title = 'Performance as a function of area ratio'
#     tlim = [1.7, 2.4]
#     ilim = [100, 135]

""" Dt change, constant De """
# xlist = np.arange(0.2, 10, 0.1)
# for p in xlist:
#     t, i, _, _ = performance_by_throat_diameter_area_ratio(power=5, mass_flow=1.8E-6, Dt=0.15 * 1E-3, De=p * 1E-3,
#                                                            print_results=False)
#     tlist.append(t)
#     ilist.append(i)
#     xlabel = "$D_t$ [mm]"
#     title = 'Performance as a function of throat diameter'
#     tlim = [1.7, 2.4]
#     ilim = [100, 135]
#
""" Dt change, plot gamma and Pc """
# xlist = np.arange(0.03, 0.3, 0.005)
# for p in xlist:
#     _, _, y, pc = performance_by_throat_diameter_area_ratio(power=5, mass_flow=1.8E-6, Dt=p * 1E-3, De=3.0E-3,
#                                                             print_results=False)
#     tlist.append(y)
#     ilist.append(pc)
#     xlabel = "$D_t$ [mm]"
#     title = '$\gamma$ and $P_c$ as a function of throat diameter'
#     tlim = [1.31, 1.43]
#     ilim = [0, 19]
#     ylabel1 = "Specific Heat Ratio"
#     ylabel2 = "Chamber Pressure [bar]"
#
# import matplotlib.pyplot as plt
#
# fig, ax1 = plt.subplots()
# ax1.plot(xlist, tlist, 'b', label='Thrust [mN]')
# ax1.set_xlabel(xlabel)
# if ylabel1:
#     ax1.set_ylabel(ylabel1)
# else:
#     ax1.set_ylabel("Thrust [mN]")
# if tlim:
#     ax1.set_ylim(tlim)
# ax2 = ax1.twinx()
# ax2.plot(xlist, ilist, 'r', label='$I_{sp}$')
# if ylabel2:
#     ax2.set_ylabel(ylabel2)
# else:
#     ax2.set_ylabel("Specific Impulse [s]")
# if ilim:
#     ax2.set_ylim(ilim)
# fig.legend(loc="upper center", bbox_to_anchor=(0.5, 0.8))
# plt.title(title)
# fig.tight_layout()
# plt.show()
