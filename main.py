from orbits import TwoBodyKeplerOrbit

# test_position = [12426, 0, 0]
# test_velocity = [0, 4.78, 0]
# test_orbit = TwoBodyKeplerOrbit(test_position, test_velocity)
# print(test_orbit)

# test_orbit = TwoBodyKeplerOrbit.build_from_known_params(9649.6, 0.2877, 180)
# print(test_orbit)

# test_orbit = TwoBodyKeplerOrbit.build_from_inertial_frame([-6796, 4025, 3490], [-3.7817, -6.0146, 1.1418])
# print(test_orbit)

test_state_vectors = TwoBodyKeplerOrbit.find_state_vectors_from_orbit(6922.3, 0.001143, 50.75, 193.89, 85.41, 10.23)
print(test_state_vectors)

# This re-building from position and velocity vectors isn't exactly right yet!
test_orbit = TwoBodyKeplerOrbit.build_from_inertial_frame(test_state_vectors[0], test_state_vectors[1])
print(test_orbit)
