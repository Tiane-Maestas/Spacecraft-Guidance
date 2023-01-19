from orbits import TwoBodyKeplerOrbit
import numpy

test_position = numpy.array([12426, 0, 0])
test_velocity = numpy.array([0, 4.78, 0])

test_orbit = TwoBodyKeplerOrbit(test_position, test_velocity)

print(test_orbit)
