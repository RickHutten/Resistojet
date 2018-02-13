import os
import convert
from annotations import *
import properties


class Files:
    file_dir = os.path.dirname(os.path.realpath(__file__)) + '\\'
    isobaric_water_6_bar = file_dir + 'isobaric_water_6_bar.cgi'
    isobaric_ammonia_6_bar = file_dir + 'isobaric_ammonia_6_bar.cgi'
    isobaric_butane_6_bar = file_dir + 'isobaric_butane_6_bar.cgi'
    isobaric_nitrogen_6_bar = file_dir + 'isobaric_nitrogen_6_bar.cgi'
    isothermal_water_375_K = file_dir + 'isothermal_water_375_K.cgi'
    isothermal_water_400_K = file_dir + 'isothermal_water_400_K.cgi'
    isothermal_water_425_K = file_dir + 'isothermal_water_425_K.cgi'
    isothermal_water_450_K = file_dir + 'isothermal_water_450_K.cgi'
    isothermal_water_475_K = file_dir + 'isothermal_water_475_K.cgi'
    isothermal_water_500_K = file_dir + 'isothermal_water_500_K.cgi'
    isothermal_water_525_K = file_dir + 'isothermal_water_525_K.cgi'
    isothermal_water_550_K = file_dir + 'isothermal_water_550_K.cgi'
    isothermal_water_575_K = file_dir + 'isothermal_water_575_K.cgi'
    isothermal_water_600_K = file_dir + 'isothermal_water_600_K.cgi'


class Index:
    temperature = 0
    pressure = 1
    density = 2
    volume = 3
    enthalpy = 5
    cv = 7
    cp = 8
    sound_speed = 9
    viscosity = 11
    thermal_conductivity = 12


class FileReader:
    class titles:
        temperature = ''
        pressure = ''
        density = ''
        volume = ''
        enthalpy = ''
        cv = ''
        cp = ''
        sound_speed = ''
        thermal_conductivity = ''
        gamma = ''

    def __init__(self, filename):
        # print "Reading:", filename  # Check to see if the file is being read too much
        with open(filename) as f:
            title_list = f.readline().split('\t')[:-1]
            self.titles.temperature = title_list[Index.temperature]
            self.titles.pressure = title_list[Index.pressure]
            self.titles.density = title_list[Index.density]
            self.titles.volume = title_list[Index.volume]
            self.titles.enthalpy = title_list[Index.enthalpy]
            self.titles.cv = title_list[Index.cv]
            self.titles.cp = title_list[Index.cp]
            self.titles.sound_speed = title_list[Index.sound_speed]
            self.titles.viscosity = title_list[Index.viscosity]
            self.titles.thermal_conductivity = title_list[Index.thermal_conductivity]
            self.titles.gamma = 'Specific Heat Ratio'

            data = f.readlines()
            for i in range(len(data)):
                data[i] = data[i].split('\t')
                data[i] = [float(j) for j in data[i] if j != data[i][-1]]
        self.temperature = [i[Index.temperature] for i in data]
        self.pressure = [i[Index.pressure] for i in data]
        self.density = [i[Index.density] for i in data]
        self.volume = [i[Index.volume] for i in data]
        self.enthalpy = [i[Index.enthalpy] for i in data]
        self.cv = [i[Index.cv] for i in data]
        self.cp = [i[Index.cp] for i in data]
        self.sound_speed = [i[Index.sound_speed] for i in data]
        self.viscosity = [i[Index.viscosity] for i in data]
        self.thermal_conductivity = [i[Index.thermal_conductivity] for i in data]
        self.gamma = [i[Index.cp] / i[Index.cv] for i in data]


class DataContainer:
    """
    Used to store data in memory so we only need to read the files once
    """

    def __init__(self):
        self.data = {'375': FileReader(Files.isothermal_water_375_K),
                     '400': FileReader(Files.isothermal_water_400_K),
                     '425': FileReader(Files.isothermal_water_425_K),
                     '450': FileReader(Files.isothermal_water_450_K),
                     '475': FileReader(Files.isothermal_water_475_K),
                     '500': FileReader(Files.isothermal_water_500_K),
                     '525': FileReader(Files.isothermal_water_525_K),
                     '550': FileReader(Files.isothermal_water_550_K),
                     '575': FileReader(Files.isothermal_water_575_K),
                     '600': FileReader(Files.isothermal_water_600_K)
                     }

    def get_data(self, temperature):
        if temperature < properties.temperature_min or temperature >= properties.temperature_max:
            raise ValueError('Temperature out of range:' + str(temperature))
        temp_rounded = int(temperature / properties.temperature_step) * properties.temperature_step
        data1 = self.data[str(temp_rounded)]
        data2 = self.data[str(temp_rounded + properties.temperature_step)]
        return data1, data2


dataContainer = DataContainer()


@deprecated
def get_specific_heat_ratio(pressure, data):
    pressure = convert.pa_to_bar(pressure)  # Convert to bar
    for i, p in enumerate(data.pressure):
        if p == pressure:
            return data.gamma[i]
        if p > pressure:
            # We're in between two pressures that are listed
            fraction = (pressure - data.pressure[i - 1]) / (data.pressure[i] - data.pressure[i - 1])
            return data.gamma[i - 1] + fraction * (data.gamma[i] - data.gamma[i - 1])


def get_specific_heat_ratio_pressure_temperature(pressure, temperature):
    pressure = convert.pa_to_bar(pressure)  # Convert to bar
    # Usually gets out of range because of the guessing algorithm
    if pressure < properties.pressure_min: pressure = properties.pressure_min
    if pressure > properties.pressure_max: pressure = properties.pressure_max

    if temperature < properties.temperature_min or temperature >= properties.temperature_max:
        raise ValueError('Temperature out of range:' + str(temperature))

    # Get isothermal water Cv and Cp for given temperature (interpolate)
    data_low, data_high = dataContainer.get_data(temperature)
    fraction = (temperature % properties.temperature_step) / float(properties.temperature_step)
    cv = [data_low.cv[i] + fraction * (data_high.cv[i] - data_low.cv[i]) for i in range(len(data_low.cv))]
    cp = [data_low.cp[i] + fraction * (data_high.cp[i] - data_low.cp[i]) for i in range(len(data_low.cp))]

    # Cv and Cp are known for the exact temperature, calculate at exact pressure (interpolate)
    for i, p in enumerate(data_low.pressure):
        if p == pressure:
            return cp[i] / cv[i]
        if p > pressure:
            # We're in between two pressures that are listed, i and i-1
            fraction = (pressure - data_low.pressure[i - 1]) / (data_low.pressure[i] - data_low.pressure[i - 1])
            gamma_low = cp[i - 1] / cv[i - 1]
            gamma_high = cp[i] / cv[i]
            return gamma_low + fraction * (gamma_high - gamma_low)


def get_enthalpy_pressure_temperature(pressure, temperature):
    pressure = convert.pa_to_bar(pressure)  # Convert to bar
    if pressure < properties.pressure_min or pressure > properties.pressure_max:
        raise ValueError('Pressure out of range:' + str(pressure))
    if temperature < properties.temperature_min or temperature >= properties.temperature_max:
        raise ValueError('Temperature out of range:' + str(temperature))

    # Get enthalpy for given temperature (interpolate)
    data_low, data_high = dataContainer.get_data(temperature)
    fraction = (temperature % properties.temperature_step) / float(properties.temperature_step)
    enthalpy = [data_low.enthalpy[i] + fraction * (data_high.enthalpy[i] - data_low.enthalpy[i]) for i in
                range(len(data_low.enthalpy))]

    # Enthalpy known for the exact temperature, calculate at exact pressure (interpolate)
    for i, p in enumerate(data_low.pressure):
        if p == pressure:
            return enthalpy[i]
        if p > pressure:
            # We're in between two pressures that are listed, i and i-1
            fraction = (pressure - data_low.pressure[i - 1]) / (data_low.pressure[i] - data_low.pressure[i - 1])
            return enthalpy[i - 1] + fraction * (enthalpy[i] - enthalpy[i - 1])
