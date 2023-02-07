from orbits import TwoBodyKeplerOrbit
from orbits import OrbitUtilities


# test_position = [-7208.2, -3822.1, 142.3]
# test_velocity = [4.0643, -6.0305, 0.8915]
# test_orbit = TwoBodyKeplerOrbit(test_position, test_velocity, angle_type='deg')
# print(test_orbit)

# test_perifocal = TwoBodyKeplerOrbit.convert_position_and_velocity_to_perifocal_frame(test_orbit)
# print(test_perifocal)

# test_orbit = TwoBodyKeplerOrbit.build_from_known_orbital_params([9056, 0.142, 7.2, 200, 60, 320], angle_type='deg', time_of_flight=0)
# print(test_orbit)

# velocity = OrbitUtilities.calculate_velocity_gibbs([[-11052.902, -12938.738, 8505.244], 
#                                                     [-10378.257, -15955.205, 14212.351], 
#                                                     [-9366.222, -17747.079, 18337.068]])
# print(velocity)

test_position = [63942.2, 37491.1, -23787.6]
test_measured_positions = [[208668.2, 94880.1, -83019.1], test_position, [35211.9, 25745.8, -12098.5]]
test_velocity = OrbitUtilities.calculate_velocity_gibbs(test_measured_positions)
test_orbit = TwoBodyKeplerOrbit(test_position, test_velocity, angle_type='deg')
print(test_orbit)
