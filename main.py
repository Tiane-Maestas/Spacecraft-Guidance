from orbits import TwoBodyKeplerOrbit

# test_position = [12426, 0, 0]
# test_velocity = [0, 4.78, 0]
# test_orbit = TwoBodyKeplerOrbit(test_position, test_velocity)
# print(test_orbit)

test_orbit = TwoBodyKeplerOrbit.build_from_known_params(9649.6, 0.2877, 180)
print(test_orbit)
