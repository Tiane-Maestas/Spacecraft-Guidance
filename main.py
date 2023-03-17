from orbits import TwoBodyKeplerOrbit
from orbits import OrbitUtilities

# --Manually set observations in these variables--
# Observation location.
site_lat = 32.881191 # N (deg)
site_lon = -117.233614 # W (deg)
alt = 0.111 # (km)
# Observation information. (3 observations)
lsts = [297.113238841384, 297.301252188220, 297.489265534848] # (deg)
RAs  = [11.199735719497, 333.281391660113, 321.343138567912] # (deg)
DECs = [81.554921692864, 71.437660976025, 54.089413168164] # (deg)
JD = [2454871.264894439, 2454871.265415273, 2454871.265936106] # (days)
# What date to propogate the orbit to.
JD_prop = 2454873.205555555

# --Run Orbit Calculations--
# Get 3 line of sight unit vectors from measurments and the 3 corresponding site vectors.
line_of_sights = OrbitUtilities.line_of_sights_from_ra_and_dec(RAs, DECs)
r_sites = OrbitUtilities.site_positions(site_lat, alt, lsts)

# Estimate the three measured positions using gauss method.
measured_positions = OrbitUtilities.positions_from_line_of_sight_gauss(JD, r_sites, line_of_sights, printRoots=False)

# Find a velocity of the middle measurment from gibbs method.
measured_velocity = OrbitUtilities.calculate_velocity_gibbs(measured_positions)

# Build the orbit and propagate it the final julian date.
orbit = TwoBodyKeplerOrbit(measured_positions[1], measured_velocity, angle_type='deg')

delta_t = JD_prop - JD[1] # days
delta_t = delta_t * 24 * 3600 # s
time_of_flight = orbit.propagate_mean_anomaly(delta_t) # 'time_of_flight' should equal 'delta_t'

# Print a clean output of only the important orbit attributes.
TwoBodyKeplerOrbit.print_clean_output(orbit)