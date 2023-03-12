from orbits import TwoBodyKeplerOrbit
from orbits import OrbitUtilities


site_lat = 32.881191 # N (deg)
site_lon = -117.233614 # W (deg)
# lst = [297.113238841384, 297.301252188220, 297.489265534848]
alt = 111 # (m)
RAs  = [11.199735719497, 333.281391660113, 321.343138567912]
DECs = [81.554921692864, 71.437660976025, 54.089413168164]
JD = [2454871.264894439, 2454871.265415273, 2454871.265936106]
JD_prop = 2454873.205555555






line_of_sights = OrbitUtilities.line_of_sights_from_ra_and_dec(RAs, DECs)

r_sites = OrbitUtilities.site_positions(site_lat, site_lon, alt, JD)

test_pos = OrbitUtilities.positions_from_line_of_sight_gauss(JD, r_sites, line_of_sights)
test_vel = OrbitUtilities.calculate_velocity_gibbs(test_pos)
test_orbit = TwoBodyKeplerOrbit(test_pos[1], test_vel, angle_type='deg')
delta_t = JD_prop - JD[1] # days
delta_t = delta_t * 24 * 3600 # s
test_orbit.propagate_mean_anomaly(test_orbit.period)
print(test_orbit)



# import pandas
# from IPython.display import display

# headers = ['Orbital Element', 'Value']
# # Build Data List from orbit here.
# data_lists = [["Test", "V"], [1234, 1.0012]]
# orbit_dataframe = pandas.DataFrame(data_lists, columns=headers)
# orbit_dataframe
# display(orbit_dataframe)