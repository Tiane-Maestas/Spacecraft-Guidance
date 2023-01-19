from orbits import TwoBodyKeplerOrbit
import numpy

test_position = numpy.array([8250, 0, 0])
test_velocity = numpy.array([1.2054, 7.0263, 0])

test_orbit = TwoBodyKeplerOrbit(test_position, test_velocity)

print(test_orbit)
