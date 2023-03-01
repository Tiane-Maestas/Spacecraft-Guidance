from orbits import TwoBodyKeplerOrbit
from orbits import OrbitUtilities

# -----Testing Building Basic Orbit-----
# test_position = [-7208.2, -3822.1, 142.3]
# test_velocity = [4.0643, -6.0305, 0.8915]
# test_orbit = TwoBodyKeplerOrbit(test_position, test_velocity, angle_type='deg')
# print(test_orbit)

# -----Testing Perifocal Transformation-----
# test_perifocal = TwoBodyKeplerOrbit.convert_position_and_velocity_to_perifocal_frame(test_orbit)
# print(test_perifocal)

# -----Testing Build from Known Kepler Orbital Elements-----
# test_orbit = TwoBodyKeplerOrbit.build_from_known_orbital_params([9056, 0.142, 7.2, 200, 60, 320], angle_type='deg', time_of_flight=0)
# print(test_orbit)

# -----Testing Gibbs Method of Orbital Determination-----
# velocity = OrbitUtilities.calculate_velocity_gibbs([[-11052.902, -12938.738, 8505.244],
#                                                     [-10378.257, -15955.205, 14212.351],
#                                                     [-9366.222, -17747.079, 18337.068]])
# print(velocity)

# -----Example Gibbs Method to Full Orbit and Verifcation of Hyperbolic Reconstruction from Known Kepler Orbital Elements-----
# test_position = [63942.2, 37491.1, -23787.6]
# test_measured_positions = [[208668.2, 94880.1, -83019.1], test_position, [35211.9, 25745.8, -12098.5]]
# test_velocity = OrbitUtilities.calculate_velocity_gibbs(test_measured_positions)
# test_orbit = TwoBodyKeplerOrbit(test_position, test_velocity, angle_type='deg')
# print(test_orbit)

# test_hyperbolic_reconstruction = TwoBodyKeplerOrbit.build_from_known_orbital_params([-9664.2, 1.663, 27.7059, 68.0541, 76.9129, 346.3288], angle_type='deg', time_of_flight=0)
# print(test_hyperbolic_reconstruction)

# -----Testing SEZ to ECI Transform-----
# test_position = OrbitUtilities.transform_position_SEZ_to_ECI([19.8, 283.5], [1298.4, 62.7, 158.2])
# print(test_position)

# -----Testing True Anomaly Propagation-----
# test_orbit = TwoBodyKeplerOrbit([-6796, 4025, 3490], [-3.7817, -6.0146, 1.1418], angle_type='deg')
# print(str(test_orbit.propagate_true_anomaly(200)) + ' [s]')
# print(test_orbit)

# -----Testing True Anomaly Propagation Lagrange-----
# test_orbit = TwoBodyKeplerOrbit([-6796, 4025, 3490], [-3.7817, -6.0146, 1.1418], angle_type='deg')
# print(str(test_orbit.propagate_true_anomaly_lagrange(200)) + ' [s]')
# print(test_orbit)

# -----Testing True Anomaly Propagation Lagrange in Perifocal-----(Unfinished) I assume transformtion to perifocal is wrong...
# test_orbit = TwoBodyKeplerOrbit([-6796, 4025, 3490], [-3.7817, -6.0146, 1.1418], angle_type='deg')
# test_pos_vel = OrbitUtilities.propagate_true_anomaly_lagrange_perifocal(test_orbit, 200)
# print(test_pos_vel)
# # Compare to my propagation (Doesn't give same answer yet...)
# test_orbit.propagate_true_anomaly_lagrange(200)
# print(TwoBodyKeplerOrbit.convert_position_and_velocity_to_perifocal_frame(test_orbit))

#-----Testing Lambert's Problem by P-Iteration-----
test_position = [3108.9, 6491.1, 428.2]
test_params = OrbitUtilities.find_velocities_from_lambert_problem_p_iteration(test_position, [-2865, 6610.2, 2234.2], 840)
print(test_params)
test_orbit = TwoBodyKeplerOrbit(test_position, test_params[0], angle_type='deg')
print(test_orbit)