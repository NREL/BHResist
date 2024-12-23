from math import exp


def inch_to_m(x_inch: float) -> float:
    """
    Convert inches to meters

    :param x_inch: value to convert, inches
    :return: float value, meters
    """
    return x_inch * 0.0254


def smoothing_function(x: float, a: float, b: float) -> float:
    """
    Sigmoid smoothing function

    https://en.wikipedia.org/wiki/Sigmoid_function

    :param x: independent variable
    :param a: fitting parameter 1
    :param b: fitting parameter 2
    :return: float between 0-1
    """

    return 1 / (1 + exp(-(x - a) / b))

def coth(x):
    return cosh(x) / sinh(x)