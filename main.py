from ChemicalProperties import *
from functions import *
from annotations import *
import convert
import numpy as np
import sys

sys.path.append('C:\\Users\\Rick\\Dropbox\\Thesis\\Code\\HeatingChamberModel\\')  # Add module to path
import HeatingChamberModel.Functions.Power_budget as pb


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


def performance_difference(a, b):
    """
    Returns how much b differs from a in percentage
    """
    return "Change: " + str(round(100 * (b / float(a)) - 100, 2)) + "%"


def performance_by_throat_diameter_area_ratio(power, mass_flow, Dt, De=None, Tc=None, AreaRatio=None, print_results=False):
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

    if Tc:
        Pc = calculate_pc(Tc)
    else:
        Tc = Pc = 0
        H_min = 1E100
        for t in range(temperature_min, temperature_max):
            Pc_temp = calculate_pc(t)
            if Pc_temp > convert.bar_to_pa(pressure_max):  # Can't calculate pressures above pressure_max
                pass
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
        print "Thrust:", Thrust * 1000, "mN", '\t\t', performance_difference(2.31386612509, Thrust * 1000)
        print "Isp:", Isp, "s", '\t\t', performance_difference(131.082600134, Isp)
        print "Specific Heat Ratio:", gamma
        print "Exhaust Velocity:", Ue, "m/s"
        print "Effective Exhaust Velocity", Thrust / mass_flow, "m/s"
        print "c*:", c_star, "m/s", '\t\t', performance_difference(691.043939954, c_star)
        print "Thrust Coefficient:", CF, '\t\t', performance_difference(1.86020179945, CF)
        print "Chamber Temperature:", Tc, "K", '\t\t', performance_difference(465, Tc)
        print "Chamber Pressure:", convert.pa_to_bar(Pc), "bar", '\t\t', performance_difference(0.703891577199,
                                                                                                convert.pa_to_bar(Pc))
        print "Exit Pressure:", convert.pa_to_bar(Pe), "bar"
        print "Pressure Ratio: 1/" + str(round(Pc / Pe, 2))
        print "Throat Diameter:", Dt * 1000, "mm"
        print "Exit Diameter:", De * 1000, "mm"
        print "Area Ratio: 1/" + str(Ae / At)
        print "Mass flow:", mass_flow * 1E6, "mg/s"
        print "Enthalpy:", enthalpy

    return Thrust * 1000, Isp, gamma, convert.pa_to_bar(Pc)


def calculate_efficiency(P, massflow, print_debug=False):
    [eff, P_in, T_struc, h, Tc, p] = pb.HCM_mflow(m_flow=massflow, D_char=14.6691E-3, h_in=H_sto*1000, p_in=0.8E5,
                                                 P_in=P, k_m=0.000, k_c=0.0826, A_ex=0.00354, e_ex=0.07,
                                                 L_cham=35E-3, prop='H2O', eps=0.952, PPI=50)
    if print_debug:
        print "Efficiency:", eff, "\tPin =>", eff * P_in
        print "Enthalpy:", h[-1]
        print "Chamber Temperature:", Tc[-1]
        print "Structure Temperature:", T_struc
        print "Pressure:", p[-1], performance_difference(p[0], p[-1])
    return eff, Tc[-1]

P_in = 10
massflow = 3.6E-6
eff, Tc = calculate_efficiency(P_in, massflow, print_debug=False)
thrust, isp, _, _ = performance_by_throat_diameter_area_ratio(power=P_in*eff, mass_flow=massflow, Dt=0.15E-3, De=3.0E-3, Tc=Tc, print_results=False)
print "Data:", thrust, isp, Tc, eff
1/0

y1list = []
y2list = []
xlim = None
ylim = None
ylabel1 = None
ylabel2 = None


""" Efficiency """
# xlist = [6.5] + range(7, 18) + [17.4]
# for P_in in xlist:
#     massflow = 2.2E-6
#     eff, Tc = calculate_efficiency(P_in, massflow)
#     print Tc, P_in, eff
#     Thrust, Isp, _, _ = performance_by_throat_diameter_area_ratio(power=P_in*eff, mass_flow=massflow, Dt=0.15E-3, De=3.0E-3, Tc=Tc)
#     y1list.append(Thrust)
#     y2list.append(Isp)
#     xlabel = "$P_{in}$ [W]"
#     title = 'Heating efficiency and chamber temperature'
#     ylabel1 = "Chamber Temperature [K]"
#     ylabel2 = "Heating efficiency"
#
# print y1list
# print y2list

f_list_14 = [1.60779513209, 1.64390536092, 1.6725452514172388, 1.76964098997, 1.9411085225443268, 2.1144180971044984, 2.2458716471308366, 2.353370681356975, 2.4457112235001137, 2.527605714641428, 2.6016976896252832, 2.669577254504918, 2.732524919359868, 2.7914142196209597]
i_list_14 = [117.10676881, 119.736924938, 121.822964986386, 128.895114902, 141.38427368486305, 154.00760104509973, 163.5822664890542, 171.4121599233301, 178.13795620465368, 184.10289480253766, 189.49952252700854, 194.44365772969786, 199.0285687598188, 203.3178812060446]
t_list_14 = [378.58762894, 394.253951659, 406.82791886571215, 450.459585679, 531.2229628080277, 617.3451237703897, 685.4119868385337, 742.6861834224286, 792.9679510313081, 838.3729537224158, 880.0976574949385, 918.8639053516492, 955.2846614804573, 989.7694060023457]
e_list_14 = [0.93444817822, 0.933936520827, 0.9310723932059984, 0.894309055187, 0.8138857148766329, 0.7186981789147054, 0.64402158652483, 0.5845295482287318, 0.5362470423655332, 0.49637056547608127, 0.46289928955312826, 0.43436954270898737, 0.40976899235693093, 0.38834252456853474]
x_list_14 = [3.9, 3.95, 4, 4.3] + range(5, 15)

f_list_16 = [1.84964774187, 1.88527842288, 1.9716641616760382, 2.067775495767558, 2.305962991164661, 2.4760546319078136, 2.611314010497081, 2.7252381545661706, 2.8245668960391463, 2.913437500151662, 2.994540990792221, 3.0693186529750704, 3.1389549701437662, 3.2043276859170198]
i_list_16 = [117.882236918, 120.153060862, 125.65861951303698, 131.78401236454076, 146.96424053860525, 157.80456577346837, 166.4249521050181, 173.68559565232334, 180.01604115824128, 185.67996590015846, 190.84887492111355, 195.61462457714094, 200.05270467895295, 204.21905581397695]
t_list_16 = [383.736165247, 397.259045226, 430.726424018, 469.086316025465, 568.9950127395884, 644.2513239313516, 706.1976910293872, 759.7037482019339, 807.3011895332288, 850.6097751173015, 890.723367406005, 928.2122635401702, 963.5631153664002, 997.1441448982532]
e_list_16 = [0.939742845139, 0.939046787716, 0.921800619423, 0.8907523336224532, 0.7952855801957249, 0.7165979917849762, 0.6526908348502838, 0.6002257748746878, 0.5565154511026839, 0.5195980229100932, 0.4880742500855107, 0.46080967298035475, 0.43700448821505467, 0.41602784979381335]
x_list_16 = [4.45, 4.5, 4.7] + range(5, 16)

f_list_18 = [2.08180901954, 2.12229995989, 2.15775503097, 2.3069770265409133, 2.4531624385867237, 2.6747692971670785, 2.844073539256185, 2.983005860044771, 3.1027254212238087, 3.208373640543391, 3.303853205595211, 3.3913225003135214, 3.472593949524054, 3.5485564640186023, 3.6061092829180064]
i_list_18 = [117.936356086, 120.230204328, 122.23876604, 130.69232649618786, 138.97386176085297, 151.52808990538043, 161.1193277156139, 168.98996882754568, 175.77218980226053, 181.75725658628423, 187.16625994691873, 192.12147427814128, 196.72557507881868, 201.02891994596652, 204.28933897565014]
t_list_18 = [384.606250909, 398.228219842, 410.309752935, 462.534203312, 515.8621219029533, 600.5145224372201, 668.0320542813419, 725.1118322857219, 775.4275621160477, 820.6580591779788, 862.1898138586012, 900.7853967196442, 937.1101752671298, 971.4828128337085, 997.79226775]
e_list_18 = [0.941559532333, 0.942062057709, 0.94135429619, 0.906764309201, 0.8628500143698253, 0.7832431235693716, 0.7164286542815349, 0.6606174139703492, 0.6137579251210808, 0.5738845758258164, 0.5396395038525504, 0.5099392496378228, 0.4839555721308482, 0.46101203821590214, 0.444502628333]
x_list_18 = [5.0, 5.05, 5.1, 5.5] + range(6, 16) + [15.8]

f_list_20 = [2.28726230828, 2.42715419484, 2.543259424055738, 2.835640888590843, 3.047899300707304, 3.2174857393708867, 3.3601452387119806, 3.4849954320186893, 3.5968148212939695, 3.6983324392826793, 3.791746130156672, 3.8793657078500443, 3.9611268474929333, 4.0003596129375305]
i_list_20 = [116.617922954, 123.750424194, 129.6701434259272, 144.5774494139611, 155.39961662276642, 164.0461186730885, 171.3197288937599, 177.68531720917386, 183.38651941763854, 188.56247746593792, 193.3252502208538, 197.79260541826437, 201.96126340253466, 203.96157775272547]
t_list_20 = [377.399625108, 419.93093156, 456.4740679352886, 553.3485934605402, 627.6796302280151, 689.2471731018736, 742.4171283466867, 789.9201644407983, 833.2125556197797, 873.1131697174276, 910.3311106541812, 945.6843577918542, 979.063821824, 995.216414894]
e_list_20 = [0.94569826599, 0.942583703089, 0.9195566647462061, 0.8431280534060229, 0.7753043614162876, 0.7173522489004216, 0.6679202281623018, 0.625596391255147, 0.5890453742824869, 0.5571817084731561, 0.5291627681084033, 0.5044521073375562, 0.48237706422, 0.472224238947]
x_list_20 = [5.5, 5.7] + range(6, 17) + [16.5]

f_list_22 = [2.52128045797, 2.60211506347, 2.67021838689, 2.776735140362681, 2.9502188737154076, 3.21861749052261, 3.4253466024822545, 3.5959612642568, 3.742211577380597, 3.872027698204446, 3.9892079967318645, 4.096686758027306, 4.196057074014907, 4.288754966952326, 4.376537761247238, 4.410182627601315]
i_list_22 = [116.863207293, 120.60995083, 123.766590059, 128.7037200806077, 136.74481897095836, 149.18529265728358, 158.76733934636442, 166.67545465469397, 173.45426444766824, 179.47133731630373, 184.90273050948565, 189.88445030238321, 194.49033767971486, 198.7869533314048, 202.85575053881516, 204.4152148890301]
t_list_22 = [379.406488116, 401.50592808, 420.482686604, 450.80852468336315, 501.91480629751015, 584.8322984272843, 651.6178941351704, 708.477016415042, 758.3788364098128, 803.5274321355589, 844.94878918107, 883.4943158214851, 919.6012992575677, 953.6962575104405, 986.359003758057, 998.9714433001623]
e_list_22 = [0.947194291141, 0.947811103017, 0.945986464313, 0.9299209763042608, 0.8952561566045649, 0.8289096945235561, 0.7700693884689548, 0.7190079521476364, 0.6746895843105803, 0.6361795464590146, 0.60243511432186, 0.572713192041142, 0.5463057449438419, 0.5226949171672925, 0.5015546448735537, 0.4936655362842654]
x_list_22 = [6.05, 6.15, 6.25, 6.5] + range(7, 18) + [17.4]

import matplotlib.pyplot as plt

plt.figure(0)
plt.plot(x_list_14, t_list_14, 'k', label='$\dot{m}$ = 1.4mg/s')
plt.plot(x_list_16, t_list_16, 'k-.', label='$\dot{m}$ = 1.6mg/s')
plt.plot(x_list_18, t_list_18, 'k--', label='$\dot{m}$ = 1.8mg/s')
plt.plot(x_list_20, t_list_20, 'k:', label='$\dot{m}$ = 2.0mg/s')
plt.plot(x_list_22, t_list_22, 'k', alpha=0.2, label='$\dot{m}$ = 2.2mg/s')
plt.xlabel("$P_{tot}$ [W]")
plt.ylabel("Chamber temperature [K]")
plt.legend(loc="best")
plt.title('Chamber temperature')
plt.tight_layout()

plt.figure(1)
plt.plot(x_list_14, e_list_14, 'k', label='$\dot{m}$ = 1.4mg/s')
plt.plot(x_list_16, e_list_16, 'k-.', label='$\dot{m}$ = 1.6mg/s')
plt.plot(x_list_18, e_list_18, 'k--', label='$\dot{m}$ = 1.8mg/s')
plt.plot(x_list_20, e_list_20, 'k:', label='$\dot{m}$ = 2.0mg/s')
plt.plot(x_list_22, e_list_22, 'k', alpha=0.2, label='$\dot{m}$ = 2.2mg/s')
plt.xlabel("$P_{tot}$ [W]")
plt.ylabel("Heating efficiency")
plt.legend(loc="best")
plt.title('Heating efficiency')
plt.tight_layout()

plt.show()

import matplotlib.pyplot as plt

plt.figure(0)
plt.plot(x_list_14, f_list_14, 'k', label='$\dot{m}$ = 1.4mg/s')
plt.plot(x_list_16, f_list_16, 'k-.', label='$\dot{m}$ = 1.6mg/s')
plt.plot(x_list_18, f_list_18, 'k--', label='$\dot{m}$ = 1.8mg/s')
plt.plot(x_list_20, f_list_20, 'k:', label='$\dot{m}$ = 2.0mg/s')
plt.plot(x_list_22, f_list_22, 'k', alpha=0.2, label='$\dot{m}$ = 2.2mg/s')
plt.xlabel("$P_{tot}$ [W]")
plt.ylabel("Thrust [mN]")
plt.legend(loc="best")
plt.title('Thrust')
plt.tight_layout()

plt.figure(1)
plt.plot(x_list_14, i_list_14, 'k', label='$\dot{m}$ = 1.4mg/s')
plt.plot(x_list_16, i_list_16, 'k-.', label='$\dot{m}$ = 1.6mg/s')
plt.plot(x_list_18, i_list_18, 'k--', label='$\dot{m}$ = 1.8mg/s')
plt.plot(x_list_20, i_list_20, 'k:', label='$\dot{m}$ = 2.0mg/s')
plt.plot(x_list_22, i_list_22, 'k', alpha=0.2, label='$\dot{m}$ = 2.2mg/s')
plt.xlabel("$P_{tot}$ [W]")
plt.ylabel("Specific impulse [s]")
plt.legend(loc="best")
plt.title('Specific impulse')
plt.tight_layout()

plt.show()

""" Power change """
# xlist = np.arange(4.78, 5.48, 0.01)
# for p in xlist:
#     t, i, _, _ = performance_by_throat_diameter_area_ratio(power=p, mass_flow=1.8E-6, Dt=0.15E-3, AreaRatio=400,
#                                                            print_results=False)
#     y1list.append(t)
#     y2list.append(i)
#     xlabel = "$P_{in}$ [W]"
#     title = 'Performance as a function of power input to the propellant'
#     xlim = [2.0, 2.8]
#     ylim = [120, 155]

""" Mass flow change """
# xlist = np.arange(1.645, 1.89, 0.005)
# for p in xlist:
#     t, i, _, _ = performance_by_throat_diameter_area_ratio(power=5, mass_flow=p * 1E-6, Dt=0.15E-3, AreaRatio=400,
#                                                            print_results=False)
#     y1list.append(t)
#     y2list.append(i)
#     xlabel = "Mass flow [mg/s]"
#     title = 'Performance as a function of mass flow'
#     xlim = [2.0, 2.8]
#     ylim = [120, 155]

""" Dt change, constant AreaRatio """
# xlist = np.arange(0.04, 0.3, 0.01)
# for p in xlist:
#     t, i, _, _ = performance_by_throat_diameter_area_ratio(power=5, mass_flow=1.8E-6, Dt=p * 1E-3, AreaRatio=400,
#                                                            print_results=False)
#     y1list.append(t)
#     y2list.append(i)
#     xlabel = "$D_t$ [mm]"
#     title = 'Performance as a function of throat diameter, with $\\frac{A_e}{A_t} = 400$'
#     xlim = [2.24, 2.34]
#     ylim = [126, 132]

""" Area Ratio change, constant Dt """
# xlist = range(2, 20, 2) + range(20, 100, 10) + range(100, 501, 25)
# for p in xlist:
#     t, i, _, _ = performance_by_throat_diameter_area_ratio(power=5, mass_flow=1.8E-6, Dt=0.15E-3, AreaRatio=p,
#                                                            print_results=False)
#     y1list.append(t)
#     y2list.append(i)
#     xlabel = "Area Ratio, $\\frac{A_e}{A_t}$"
#     title = 'Performance as a function of area ratio'
#     xlim = [1.7, 2.4]
#     ylim = [100, 135]

""" Dt change, constant De """
# xlist = np.arange(0.2, 10, 0.1)
# for p in xlist:
#     t, i, _, _ = performance_by_throat_diameter_area_ratio(power=5, mass_flow=1.8E-6, Dt=0.15 * 1E-3, De=p * 1E-3,
#                                                            print_results=False)
#     y1list.append(t)
#     y2list.append(i)
#     xlabel = "$D_t$ [mm]"
#     title = 'Performance as a function of throat diameter'
#     xlim = [1.7, 2.4]
#     ylim = [100, 135]
#
""" Dt change, plot gamma and Pc """
# xlist = np.arange(0.03, 0.3, 0.005)
# for p in xlist:
#     _, _, y, pc = performance_by_throat_diameter_area_ratio(power=5, mass_flow=1.8E-6, Dt=p * 1E-3, De=3.0E-3,
#                                                             print_results=False)
#     y1list.append(y)
#     y2list.append(pc)
#     xlabel = "$D_t$ [mm]"
#     title = '$\gamma$ and $P_c$ as a function of throat diameter'
#     xlim = [1.31, 1.43]
#     ylim = [0, 19]
#     xlabel = 'Thrust [mN]'
#     ylabel1 = "Specific Heat Ratio"
#     ylabel2 = "Chamber Pressure [bar]"

# import matplotlib.pyplot as plt
#
# fig, ax1 = plt.subplots()
# ax1.plot(xlist, y1list, 'b', label='$T_c$')
# ax1.set_xlabel(xlabel)
# if ylabel1:
#     ax1.set_ylabel(ylabel1)
# else:
#     ax1.set_ylabel("No Label")
# if xlim:
#     ax1.set_ylim(xlim)
# ax2 = ax1.twinx()
# ax2.plot(xlist, y2list, 'r', label='$\eta$')
# if ylabel2:
#     ax2.set_ylabel(ylabel2)
# else:
#     ax2.set_ylabel("No Label")
# if ylim:
#     ax2.set_ylim(ylim)
# fig.legend(loc="upper center", bbox_to_anchor=(0.5, 0.9))
# plt.title(title)
# fig.tight_layout()
# plt.show()
