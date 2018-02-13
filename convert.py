import math


def bar_to_pa(bar):
    return bar * 100000


def pa_to_bar(pa):
    return pa / 100000.


def celcius_to_kelvin(celcius):
    return celcius + 273.15


def kelvin_to_celcius(kelvin):
    return kelvin - 273.15


def diameter_to_area(diameter):
    return math.pi * (diameter / 2.) ** 2


def area_to_diameter(area):
    return (area / math.pi) ** 0.5 * 2.
