from orbits import TwoBodyKeplerOrbit
from orbits import OrbitUtilities


# test_position = [-7208.2, -3822.1, 142.3]
# test_velocity = [4.0643, -6.0305, 0.8915]
# test_orbit = TwoBodyKeplerOrbit(test_position, test_velocity, angle_type='deg')
# print(test_orbit)

# old_test_params = [9649.6, 0.2877, 180]
# test_orbit = TwoBodyKeplerOrbit.build_from_known_orbital_params([9056, 0.142, 7.2, 200, 60, 320], angle_type='deg', time_of_flight=0)
# print(test_orbit)

velocity = OrbitUtilities.calculate_velocity_gibbs([[-11052.902, -12938.738, 8505.244], 
                                         [-10378.257, -15955.205, 14212.351], 
                                         [-9366.222, -17747.079, 18337.068]])
print(velocity)
