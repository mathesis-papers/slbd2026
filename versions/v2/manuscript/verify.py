from fractions import Fraction
import json
from math import gcd


class Exact:
    def __init__(self, rational=0, radical=0, imaginary=0, imaginary_radical=0):
        self.coefficients = tuple(map(Fraction, (rational, radical, imaginary, imaginary_radical)))

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Exact) else Exact(value)

    def __add__(self, other):
        other = self.coerce(other)
        return Exact(*(first + second for first, second in zip(self.coefficients, other.coefficients)))

    __radd__ = __add__

    def __neg__(self):
        return Exact(*(-coefficient for coefficient in self.coefficients))

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        result = [Fraction(0) for index in range(4)]
        for first, first_coefficient in enumerate(self.coefficients):
            if first_coefficient == 0:
                continue
            for second, second_coefficient in enumerate(other.coefficients):
                if second_coefficient == 0:
                    continue
                factor = (3 if first & second & 1 else 1) * (-1 if first & second & 2 else 1)
                result[first ^ second] += factor * first_coefficient * second_coefficient
        return Exact(*result)

    __rmul__ = __mul__

    def __truediv__(self, divisor):
        return Exact(*(coefficient / Fraction(divisor) for coefficient in self.coefficients))

    def __pow__(self, exponent):
        result = Exact(1)
        for index in range(exponent):
            result = result * self
        return result

    def conjugate(self):
        rational, radical, imaginary, imaginary_radical = self.coefficients
        return Exact(rational, radical, -imaginary, -imaginary_radical)

    def __eq__(self, other):
        return self.coefficients == self.coerce(other).coefficients

    def real(self):
        return self.coefficients[2:] == (0, 0)


LABELS = [(first, second) for first in range(3) for second in range(3)]
IDENTITY = [[Exact(int(row == column)) for column in range(9)] for row in range(9)]
IMAGINARY = Exact(imaginary=1)
OMEGA = Exact(Fraction(-1, 2), imaginary_radical=Fraction(1, 2))


def phase(exponent):
    return OMEGA ** (exponent % 3)


def multiply(left, right):
    dimension = len(left)
    return [[sum((left[row][middle] * right[middle][column]
                  for middle in range(dimension)
                  if left[row][middle] != 0 and right[middle][column] != 0), Exact())
             for column in range(dimension)] for row in range(dimension)]


def combine(left, right, factor=1):
    return [[left[row][column] + factor * right[row][column]
             for column in range(len(left))] for row in range(len(left))]


def scale(matrix, scalar):
    return [[scalar * entry for entry in row] for row in matrix]


def adjoint(matrix):
    return [[matrix[column][row].conjugate() for column in range(len(matrix))]
            for row in range(len(matrix))]


def trace(matrix):
    return sum(matrix[index][index] for index in range(len(matrix)))


def symbols(matrix, fourier):
    transformed = multiply(matrix, fourier)
    dimension = len(matrix)
    return [[dimension * fourier[row][column].conjugate() * transformed[row][column]
             for column in range(dimension)] for row in range(dimension)]


def construct():
    fourier = [[phase(first * frequency + second * momentum) / 3
                for frequency, momentum in LABELS] for first, second in LABELS]
    projection = [[int(second == other_second) * phase(second * (first - other_first)) / 3
                   for other_first, other_second in LABELS] for first, second in LABELS]
    other_projection = [[int(first == other_first) * phase(first * (other_second - second)) / 3
                         for other_first, other_second in LABELS] for first, second in LABELS]
    observable = scale(combine(scale(other_projection, 3), IDENTITY, -1), Fraction(1, 2))
    commutator = combine(multiply(projection, other_projection), multiply(other_projection, projection), -1)
    return fourier, projection, other_projection, observable, commutator


def require(condition, label):
    if not condition:
        raise ValueError(label)


def square_check(order, root):
    dimension = order ** 2
    labels = [(first, second) for first in range(order) for second in range(order)]
    phases = [root ** exponent for exponent in range(order)]
    identity = [[Exact(int(row == column)) for column in range(dimension)] for row in range(dimension)]
    fourier = [[phases[(first * frequency + second * momentum) % order] / order
                for frequency, momentum in labels] for first, second in labels]
    projection = [[int(second == other_second) * phases[(second * (first - other_first)) % order] / order
                   for other_first, other_second in labels] for first, second in labels]
    other_projection = [[int(first == other_first) * phases[(first * (other_second - second)) % order] / order
                         for other_first, other_second in labels] for first, second in labels]
    require(multiply(adjoint(fourier), fourier) == identity, 'square Fourier unitarity')
    for matrix in (projection, other_projection):
        require(adjoint(matrix) == matrix and multiply(matrix, matrix) == matrix, 'square projection')
        require(trace(matrix) == order, 'square projection rank')
    product = multiply(projection, other_projection)
    reverse = multiply(other_projection, projection)
    common = gcd(2, order)
    require(trace(product) == 1, 'square trace PQ')
    require(trace(multiply(product, product)) == Fraction(common, order), 'square gcd trace')
    hermitian_commutator = scale(combine(product, reverse, -1), IMAGINARY)
    square = multiply(hermitian_commutator, hermitian_commutator)
    require(trace(square) == Fraction(2 * (order - common), order), 'square H trace')
    require(trace(multiply(projection, square)) == Fraction(order - common, order), 'square PH trace')
    for matrix, sign in ((product, 1), (reverse, -1)):
        actual = symbols(matrix, fourier)
        for row, (first, second) in enumerate(labels):
            for column, (frequency, momentum) in enumerate(labels):
                require(actual[row][column] == phases[(sign * (first + momentum) * (second - frequency)) % order] / order,
                        'square product symbol')
    density = scale(combine(identity, projection), Fraction(1, order * (order + 1)))
    observable = scale(combine(scale(other_projection, order), identity, -1), Fraction(1, order - 1))
    derivative = scale(hermitian_commutator, Fraction(2 * order, 3 * (order - 1)))
    tangent = scale(combine(multiply(density, observable), multiply(observable, density), -1), IMAGINARY)
    require(tangent == scale(combine(multiply(density, derivative), multiply(derivative, density)), Fraction(1, 2)),
            'square SLD equation')
    fisher = Fraction(4 * (order - common), 3 * (order - 1) ** 2 * (order + 1))
    require(trace(multiply(density, multiply(derivative, derivative))) == fisher, 'square QFI')
    for matrix in (density, observable, derivative):
        require(all(entry.real() for row in symbols(matrix, fourier) for entry in row), 'square KD reality')
    transformed_density = multiply(multiply(adjoint(fourier), density), fourier)
    transformed_tangent = multiply(multiply(adjoint(fourier), tangent), fourier)
    require(all(density[index][index] == Fraction(1, dimension)
                and transformed_density[index][index] == Fraction(1, dimension)
                and tangent[index][index] == 0 and transformed_tangent[index][index] == 0
                for index in range(dimension)), 'square marginal identities')
    return {'order': order, 'dimension': dimension, 'gcd_2_order': common, 'lambda': '1', 'qfi': str(fisher)}


def endpoint_checks(projection, observable, commutator, fourier):
    outcomes = []
    for parameter in (Fraction(0), Fraction(1, 3), Fraction(1)):
        density = combine(scale(projection, parameter / 3),
                          scale(combine(IDENTITY, projection, -1), (1 - parameter) / 6))
        derivative = scale(commutator, 3 * (3 * parameter - 1) * IMAGINARY / (1 + parameter))
        tangent = scale(combine(multiply(density, observable), multiply(observable, density), -1), IMAGINARY)
        require(tangent == scale(combine(multiply(density, derivative), multiply(derivative, density)), Fraction(1, 2)), 'endpoint SLD')
        fisher = (3 * parameter - 1) ** 2 / (1 + parameter)
        require(trace(multiply(density, multiply(derivative, derivative))) == fisher, 'endpoint QFI')
        density_symbols = symbols(density, fourier)
        for row, (first, second) in enumerate(LABELS):
            for column, (frequency, momentum) in enumerate(LABELS):
                expected = parameter / 3 if frequency == second else (1 - parameter) / 6
                require(density_symbols[row][column] == expected, 'endpoint positive symbol')
        outcomes.append({'t': str(parameter), 'qfi': str(fisher)})
    return outcomes


def main():
    fourier, projection, other_projection, observable, commutator = construct()
    require(OMEGA ** 3 == 1 and 1 + OMEGA + OMEGA ** 2 == 0, 'cube root identities')
    require(multiply(adjoint(fourier), fourier) == IDENTITY, 'Fourier unitarity')
    require(all(entry != 0 for row in fourier for entry in row), 'nonzero overlaps')
    for matrix in (projection, other_projection):
        require(adjoint(matrix) == matrix, 'Hermitian projection')
        require(multiply(matrix, matrix) == matrix and trace(matrix) == 3, 'rank-three projection')
    require(multiply(multiply(projection, other_projection), projection) == scale(projection, Fraction(1, 3)), 'PQP=P/3')
    require(multiply(multiply(other_projection, projection), other_projection) == scale(other_projection, Fraction(1, 3)), 'QPQ=Q/3')
    require(trace(multiply(projection, other_projection)) == 1, 'trace PQ')
    require(trace(multiply(adjoint(commutator), commutator)) == Fraction(4, 3), 'nonzero commutator')
    for matrix, sign in ((multiply(projection, other_projection), 1), (multiply(other_projection, projection), -1)):
        actual = symbols(matrix, fourier)
        for row, (first, second) in enumerate(LABELS):
            for column, (frequency, momentum) in enumerate(LABELS):
                require(actual[row][column] == phase(sign * (first + momentum) * (second - frequency)) / 3,
                        'product-symbol identity')
    square = scale(multiply(commutator, commutator), -1)
    require(trace(square) == Fraction(4, 3), 'trace H squared')
    require(trace(multiply(projection, square)) == Fraction(2, 3), 'trace P H squared')
    outcomes = []
    for parameter in (Fraction(-1, 2), Fraction(0), Fraction(1, 2), Fraction(1), Fraction(2), Fraction(10)):
        normalization = 9 + 3 * parameter
        density = scale(combine(IDENTITY, projection, parameter), 1 / normalization)
        derivative = scale(commutator, 3 * parameter * IMAGINARY / (2 + parameter))
        require(adjoint(density) == density and trace(density) == 1, 'density trace and symmetry')
        require((1 + parameter) / normalization > 0 and 1 / normalization > 0, 'positive spectral coefficients')
        require(adjoint(observable) == observable and adjoint(derivative) == derivative, 'Hermitian observables')
        for matrix in (density, observable, derivative):
            require(all(entry.real() for row in symbols(matrix, fourier) for entry in row), 'KD reality')
        density_symbols = symbols(density, fourier)
        observable_symbols = symbols(observable, fourier)
        for row, (first, second) in enumerate(LABELS):
            for column, (frequency, momentum) in enumerate(LABELS):
                require(density_symbols[row][column] == (1 + parameter * int(frequency == second)) / normalization, 'density symbol')
                require(observable_symbols[row][column] == Fraction(3 * int((momentum + first) % 3 == 0) - 1, 2), 'generator symbol')
        tangent = scale(combine(multiply(density, observable), multiply(observable, density), -1), IMAGINARY)
        jordan = scale(combine(multiply(derivative, density), multiply(density, derivative)), Fraction(1, 2))
        require(tangent == jordan, 'SLD equation')
        fisher = 2 * parameter ** 2 / ((2 + parameter) * (3 + parameter))
        require(trace(multiply(multiply(density, derivative), derivative)) == fisher, 'QFI value')
        fourier_density = multiply(multiply(adjoint(fourier), density), fourier)
        fourier_tangent = multiply(multiply(adjoint(fourier), tangent), fourier)
        require(all(density[index][index] == Fraction(1, 9) and fourier_density[index][index] == Fraction(1, 9)
                    and tangent[index][index] == 0 and fourier_tangent[index][index] == 0 for index in range(9)), 'marginals and tangent diagonals')
        outcomes.append({'lambda': str(parameter), 'qfi': str(fisher)})
    obstruction = combine(multiply(observable, scale(commutator, IMAGINARY)),
                          multiply(scale(commutator, IMAGINARY), observable), -1)
    transformed_obstruction = multiply(multiply(adjoint(fourier), obstruction), fourier)
    require(all(obstruction[index][index] == 0 and transformed_obstruction[index][index] == 0
                for index in range(9)), 'basis-diagonal trace obstruction')
    endpoints = endpoint_checks(projection, observable, commutator, fourier)
    squares = [square_check(2, Exact(-1)), square_check(4, IMAGINARY),
               square_check(6, Exact(Fraction(1, 2), imaginary_radical=Fraction(1, 2)))]
    print(json.dumps({'arithmetic': 'exact Q(sqrt(3), i)', 'dimension': 9,
                      'identities': 'Fourier, projection, product symbols, trace formulas',
                      'finite_parameter_checks': outcomes, 'closed_interval_checks': endpoints,
                      'square_checks': squares, 'trace_obstruction': 'zero diagonals in both bases',
                      'result': 'ALL_CHECKS_OK'}, indent=2))


if __name__ == '__main__':
    main()
