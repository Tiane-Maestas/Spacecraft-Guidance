"""This module holds tools for use with orbit calculations and representations."""
import numpy as np
from prettytable import PrettyTable

class KeplerOrbitTypes:
    """Different characterizations of Kepler Orbits."""
    CIRCULAR = "Circular Orbit"
    ELIPTICAL = "Eliptical Orbit"
    PARABOLIC = "Parabolic Orbit"
    HYPERBOLIC = "Hyperbolic Orbit"

    @staticmethod
    def get_orbit_type(eccentricity) -> str:
        """Characterizes orbit from eccentricity."""
        if eccentricity == 0:
            return KeplerOrbitTypes.CIRCULAR
        if eccentricity == 1:
            return KeplerOrbitTypes.PARABOLIC
        if eccentricity > 1:
            return KeplerOrbitTypes.HYPERBOLIC

        return KeplerOrbitTypes.ELIPTICAL


class TwoBodyKeplerOrbit:
    """A classic representation of two-body keplerian motion."""

    def __init__(self, position_vector, velocity_vector, **options):
        """Initialize from a measured position and velocity in the cartesian ECI (Earth Centered Inertial) frame."""
        options = TwoBodyKeplerOrbit.process_optional_params(options)
        self.angle_type = options[0] # This specifies what unit to print in.
        self.time_of_flight = options[1] # This variable doesn't do anything yet here.

        self.position_vector = np.array(position_vector)
        self.velocity_vector = np.array(velocity_vector)

        self.calculate_orbit_info()

    @classmethod
    def build_from_known_orbital_params(cls, orbital_params_list, **options):
        """Initialize from a list of known orbital parameters. [semimajor axis (a), eccentricity (ec), inclination (i), right ascension of the ascending node (omega), argument of periapsis (w), mean or hyperbolic anomaly at epoch (theta)]"""
        options = TwoBodyKeplerOrbit.process_optional_params(options)
        angle_type = options[0]
        time_of_flight = options[1]

        a = orbital_params_list[0]
        ec = orbital_params_list[1]
        i = orbital_params_list[2]
        omega = orbital_params_list[3]
        w = orbital_params_list[4]
        theta = orbital_params_list[5]

        if (angle_type == 'deg'):
            i = np.radians(i); omega = np.radians(omega); w = np.radians(w); theta = np.radians(theta)
        i = np.mod(i, 2 * np.pi); omega = np.mod(omega, 2 * np.pi); w = np.mod(w, 2 * np.pi); theta = np.mod(theta, 2 * np.pi)

        orbit_type = KeplerOrbitTypes.get_orbit_type(ec)

        # Orbital mean motion [ rad/s ]
        nm = np.sqrt(OrbitUtilities.EARTH_MU / np.power(float(np.abs(a)), 3))
        # Mean anomaly [ rad ]
        M = theta + nm * time_of_flight

        if (orbit_type == KeplerOrbitTypes.ELIPTICAL) or (orbit_type == KeplerOrbitTypes.CIRCULAR):
            # Converged eccentric anomaly [ rad ]
            E = OrbitUtilities.eccentric_anomaly_from_mean(M, ec)
            # True anomaly [ rad ]
            f = np.mod(2 * np.arctan(np.sqrt((1 + ec)/(1 - ec)) * np.tan(E / 2)), 2 * np.pi)
        elif orbit_type == KeplerOrbitTypes.HYPERBOLIC:
            # Converged hyperbolic anomaly [ rad ]
            H = OrbitUtilities.hyperbolic_anomaly_from_mean(M, ec)
            # True anomaly [ rad ]
            f = np.mod(2 * np.arctan(np.sqrt((ec + 1)/(ec - 1)) * np.tanh(H / 2)), 2 * np.pi)

        # Argument of Latitude
        u = w + f

        r = a * (1 - ec**2) / (1 + ec * np.cos(f))
        v = np.sqrt(OrbitUtilities.EARTH_MU * (2 / r - 1 / a))

        # Ascending node vector
        nhat = np.array([ np.cos(omega) , np.sin(omega), 0 ])
        rT = np.array([ -1 * np.cos(i) * np.sin(omega), np.cos(i) * np.cos(omega), np.sin(i) ])

        gamma = np.arctan2(ec * np.sin(f), 1 + ec * np.cos(f))

        rhat = np.cos(u) * nhat + np.sin(u) * rT
        vhat = np.sin(gamma - u) * nhat + np.cos(gamma - u) * rT

        position_vector = r * rhat
        velocity_vector = v * vhat

        return cls(position_vector, velocity_vector, angle_type=angle_type, time_of_flight=time_of_flight)

    def propagate_mean_anomaly(self, delta_t):
        """This will propagate the orbit forward by a given time amount. (delta_t > 0)"""
        if self.orbit_type != KeplerOrbitTypes.ELIPTICAL:
            return None

        final_mean_anomaly = np.sqrt(OrbitUtilities.EARTH_MU / self.semi_major_axis**3) * delta_t + self.mean_anomaly
        num_of_full_revs = int((final_mean_anomaly) / (2 * np.pi))

        final_eccentric_anomaly = OrbitUtilities.eccentric_anomaly_from_mean(final_mean_anomaly, self.eccentricity)

        final_true_anomaly = 2 * np.arctan((np.sqrt(1 + self.eccentricity) / np.sqrt(1 - self.eccentricity)) * np.tan(final_eccentric_anomaly / 2))
        delta_true_anomaly = final_true_anomaly - self.true_anomaly

        if self.angle_type == 'deg':
            delta_true_anomaly = np.degrees(delta_true_anomaly)

        return self.propagate_true_anomaly(delta_true_anomaly, num_of_revolutions=num_of_full_revs) # self.propagate_true_anomaly_lagrange(delta_true_anomaly)

    def propagate_true_anomaly(self, delta_true_anomaly, num_of_revolutions=0) -> float:
        """This will increase the true anomaly by the given amount and re-calculate the changed orbital elements. This function returns the time of flight this propagation took."""
        if self.angle_type == 'deg':
            delta_true_anomaly = np.radians(delta_true_anomaly)

        self.true_anomaly = np.mod(self.true_anomaly + delta_true_anomaly, 2 * np.pi)

        # Fix position given new true anomaly.
        f = self.true_anomaly
        ec = self.eccentricity
        a = self.semi_major_axis
        i = self.inclination
        w = self.argument_of_periapsis
        omega = self.right_ascension_of_ascending_node

        u = w + f
        r = a * (1 - ec**2) / (1 + ec * np.cos(f))
        v = np.sqrt(OrbitUtilities.EARTH_MU * (2 / r - 1 / a))

        nhat = np.array([ np.cos(omega) , np.sin(omega), 0 ])
        rT   = np.array([ -1 * np.cos(i) * np.sin(omega), np.cos(i) * np.cos(omega), np.sin(i) ])
        gamma = np.arctan2(ec * np.sin(f), 1 + ec * np.cos(f))

        rhat = np.cos(u) * nhat + np.sin(u) * rT
        vhat = np.sin(gamma - u) * nhat + np.cos(gamma - u) * rT

        self.position_vector = r * rhat
        self.position_hat = rhat
        self.velocity_vector = v * vhat

        # Fix Mean Anomaly and Find TOF.
        eccentric_anomaly = 2 * np.arctan((np.sqrt(1 - self.eccentricity) / np.sqrt(1 + self.eccentricity)) * np.tan(self.true_anomaly / 2))
        old_mean_anomaly = self.mean_anomaly
        self.mean_anomaly = eccentric_anomaly - self.eccentricity * np.sin(eccentric_anomaly)

        if self.mean_anomaly < 0:
            self.mean_anomaly += 2 * np.pi

        time_of_flight = num_of_revolutions * self.period + (self.mean_anomaly - old_mean_anomaly) / (self.mean_motion)

        return time_of_flight

    def propagate_true_anomaly_lagrange(self, delta_true_anomaly, num_of_revolutions=0):
        """This will use lagrange method to propagate the changed orbital elements given change in true anomaly.  This function returns the time of flight this propagation took."""
        if self.angle_type == 'deg':
            delta_true_anomaly = np.radians(delta_true_anomaly)

        position = self.position_vector
        velocity = self.velocity_vector

        # Constants needed for Lagrange Coefficients
        r_0 = np.linalg.norm(position)
        v_0 = np.linalg.norm(velocity)
        v_r_0 = np.dot(position, velocity) / r_0
        h = r_0 * np.sqrt(v_0**2 - v_r_0**2)
        b = h**2 / OrbitUtilities.EARTH_MU
        r = b * (1 / (1 + (b / r_0 - 1) * np.cos(delta_true_anomaly) - (b / h) * v_r_0 * np.sin(delta_true_anomaly)))

        # Lagrange Coefficients -> Transformation Matrix -> Propagated Position and Velocity
        f = 1 - (1 / b) * r * (1 - np.cos(delta_true_anomaly))
        g = r * r_0 / h * np.sin(delta_true_anomaly)
        f_dot = (h / b) * ((1 - np.cos(delta_true_anomaly)) / np.sin(delta_true_anomaly)) * (1 / b * (1 - np.cos(delta_true_anomaly)) - 1 / r_0 - 1 / r)
        g_dot = 1 - 1 / b * r_0 * (1 - np.cos(delta_true_anomaly))

        # Update position and velocity
        self.position_vector = f * position + g * velocity
        self.velocity_vector = f_dot * position + g_dot * velocity

        # Update True Anomaly and Mean Anomaly
        self.true_anomaly = np.mod(self.true_anomaly + delta_true_anomaly, 2 * np.pi)
        eccentric_anomaly = 2 * np.arctan((np.sqrt(1 - self.eccentricity) / np.sqrt(1 + self.eccentricity)) * np.tan(self.true_anomaly / 2))
        old_mean_anomaly = self.mean_anomaly
        self.mean_anomaly = eccentric_anomaly - self.eccentricity * np.sin(eccentric_anomaly)
        if self.mean_anomaly < 0:
            self.mean_anomaly += 2 * np.pi

        time_of_flight = num_of_revolutions * self.period + (self.mean_anomaly - old_mean_anomaly) / (self.mean_motion)

        return time_of_flight

    def calculate_orbit_info(self) -> None:
        """This calculates all of the info that can be displayed via this 'to string' method."""
        # using this weird extra definition because numpy is bugged and it says code isn't reachable otherwise.
        def cross(x, y): return np.cross(x, y)

        #  ECI Coordinate frame unit vectors
        I = np.array([1, 0, 0]); J = np.array([0, 1, 0]); K = np.array([0, 0, 1])

        # Angular momentum is orthogonal to orbital plane.
        self.angular_momentum = cross(self.position_vector, self.velocity_vector)

        self.position_hat = self.position_vector / np.linalg.norm(self.position_vector)
        self.angular_momentum_hat = self.angular_momentum / np.linalg.norm(self.angular_momentum)

        # Normalized ascending node vector.
        self.ascending_node_hat = cross(K, self.angular_momentum) / np.linalg.norm(cross(K, self.angular_momentum))

        self.total_energy = 0.5 * np.dot(self.velocity_vector, self.velocity_vector) - OrbitUtilities.EARTH_MU / np.linalg.norm(self.position_vector)

        # Eccentricity
        self.eccentricity_vector = 1 / OrbitUtilities.EARTH_MU * cross(self.velocity_vector, self.angular_momentum) - self.position_hat
        self.eccentricity =  np.linalg.norm(self.eccentricity_vector)
        self.orbit_type = KeplerOrbitTypes.get_orbit_type(self.eccentricity)

        # If orbit_type is 'parabolic' energy is 0.
        if self.orbit_type == KeplerOrbitTypes.PARABOLIC:
            self.semi_major_axis = np.Inf
            self.parameter = np.dot(self.angular_momentum, self.angular_momentum) / OrbitUtilities.EARTH_MU
        else:
            self.semi_major_axis = -1 * OrbitUtilities.EARTH_MU / (2 * self.total_energy)
            self.parameter = self.semi_major_axis * (1 - np.power(self.eccentricity, 2))
            self.perigee = self.semi_major_axis * (1 - self.eccentricity) # Closest point

        if self.orbit_type == KeplerOrbitTypes.HYPERBOLIC:
            self.semi_minor_axis = None
            self.period = None
            self.apogee = None # Furthest point
        else:
            self.semi_minor_axis = self.semi_major_axis * np.sqrt(1 - self.eccentricity**2)
            self.period = 2 * np.pi * np.sqrt(np.power(self.semi_major_axis, 3)/OrbitUtilities.EARTH_MU)
            self.apogee = self.semi_major_axis * (1 + self.eccentricity) # Furthest point

        # If the inclination is less than 90 degrees, hte elliptical orbit is a direct (prograde) orbit.
        self.inclination = np.arccos(np.dot(K, self.angular_momentum_hat))

        self.right_ascension_of_ascending_node = np.mod(np.arctan2(np.dot(J, self.ascending_node_hat), np.dot(I, self.ascending_node_hat)), 2 * np.pi)

        self.argument_of_periapsis = np.mod(np.arctan2(np.dot(self.angular_momentum_hat, cross(self.ascending_node_hat, self.eccentricity_vector)), np.dot(self.ascending_node_hat, self.eccentricity_vector)), 2 * np.pi)

        self.true_anomaly = np.mod(np.arctan2(np.dot(self.angular_momentum_hat, cross(self.eccentricity_vector, self.position_vector)), np.dot(self.eccentricity_vector, self.position_vector)), 2 * np.pi)

        self.flight_path_angle = np.degrees(np.arctan(self.eccentricity * np.sin(self.true_anomaly) / (1 + self.eccentricity * np.cos(self.true_anomaly))))

        if (self.orbit_type == KeplerOrbitTypes.ELIPTICAL) or (self.orbit_type == KeplerOrbitTypes.CIRCULAR):
            # Eccentric and Mean anomaly at epoch.
            self.eccentric_anomaly = 2 * np.arctan2(np.sqrt(1 - self.eccentricity) * np.tan(self.true_anomaly / 2), np.sqrt(1 + self.eccentricity))
            self.mean_anomaly = np.mod((self.eccentric_anomaly - self.eccentricity * np.sin(self.eccentric_anomaly)), 2 * np.pi)
        elif self.orbit_type == KeplerOrbitTypes.HYPERBOLIC:
            self.hyperbolic_anomaly = 2 * np.arctanh(np.sqrt(self.eccentricity - 1) * np.tan(self.true_anomaly / 2) / np.sqrt(self.eccentricity + 1))
            self.mean_anomaly = np.mod((self.eccentricity * np.sinh(self.hyperbolic_anomaly) - self.hyperbolic_anomaly), 2 * np.pi)

        # Orbital mean motion [ rad/s ] (The same as {2pi / Period})
        self.mean_motion = np.sqrt(OrbitUtilities.EARTH_MU / np.power(float(np.abs(self.semi_major_axis)), 3))

    @staticmethod
    def convert_position_and_velocity_to_perifocal_frame(orbit):
        """This will take in a pre-constructed orbit and return the position and velocity of the orbit from the perifocal frame."""
        H = np.linalg.norm(orbit.angular_momentum)
        e = orbit.eccentricity
        theta = orbit.mean_anomaly

        P = orbit.semi_major_axis * (1 - e**2)
        gamma = P / (1 + e * np.cos(theta))

        perifocal_position = np.array([gamma * np.cos(theta), gamma * np.sin(theta)])
        perifocal_velocity = np.array([-1 * (OrbitUtilities.EARTH_MU / H) * (np.sin(theta)), (OrbitUtilities.EARTH_MU / H) * (e + np.cos(theta))])

        return (perifocal_position, perifocal_velocity)


    ORBIT_INFO_FORMAT = "\n-----Orbit INFO-----\nOrbit Type: {orbit_type}\nPosition: {position} [km]\nVelocity: {velocity} [km/s]\nAngular Momentum(H): {angular_momentum} [km^3/s]\nTotal Energy(E): {total_energy} [km^2/s^2]\nSemi-Major Axis(a): {semi_major_axis} [km]\nSemi-Minor Axis(b): {semi_minor_axis} [km]\nParameter(p): {parameter}\nEccentricity(e): {eccentricity}\nPeriod(T): {period} [s]\nPerigee: {perigee} [km]\nApogee: {apogee} [km]\nTrue Anomaly(f): {true_anomaly} [{a_unit}]\nFlight Path Angle(gamma): {flight_path_angle} [{a_unit}]\nMean Anomaly(M): {mean_anomaly} [{a_unit}]\n \
                         \n-----Orientation INFO-----\nRight Ascension of Ascending Node(Omega): {right_ascension_of_ascending_node} [{a_unit}]\nInclination(i): {inclination} [{a_unit}]\nArgument of Periapsis(w): {argument_of_periapsis} [{a_unit}]\n"
    def __str__(self) -> str:
        true_anomaly = self.true_anomaly; flight_path_angle = self.flight_path_angle; mean_anomaly = self.mean_anomaly; right_ascension_of_ascending_node = self.right_ascension_of_ascending_node; inclination = self.inclination; argument_of_periapsis = self.argument_of_periapsis
        a_unit = 'rad'
        if(self.angle_type == 'deg'):
            true_anomaly = np.degrees(true_anomaly); flight_path_angle = np.degrees(flight_path_angle); mean_anomaly = np.degrees(mean_anomaly); right_ascension_of_ascending_node = np.degrees(right_ascension_of_ascending_node); inclination = np.degrees(inclination); argument_of_periapsis = np.degrees(argument_of_periapsis)
            a_unit = 'deg'
        return TwoBodyKeplerOrbit.ORBIT_INFO_FORMAT.format(orbit_type=self.orbit_type, position=self.position_vector, velocity=self.velocity_vector, angular_momentum=self.angular_momentum, total_energy=self.total_energy, semi_major_axis=self.semi_major_axis, semi_minor_axis=self.semi_minor_axis, parameter=self.parameter, eccentricity=self.eccentricity, period=self.period, perigee=self.perigee, apogee=self.apogee, true_anomaly=true_anomaly, flight_path_angle=flight_path_angle, mean_anomaly=mean_anomaly, right_ascension_of_ascending_node=right_ascension_of_ascending_node, inclination=inclination, argument_of_periapsis=argument_of_periapsis, a_unit=a_unit)

    @staticmethod
    def process_optional_params(options):
        """This will return a list of all the optional paramater values. [angle_type, time_of_flight]"""
        return_list = ['rad', 0]
        for key, value in options.items(): # Process optional paramaters.
            if key == "angle_type":
                return_list[0] = value  # Can be 'deg'.
            elif key == "time_of_flight":
                return_list[1] = value
        return return_list
    
    @staticmethod
    def print_clean_output(orbit):
        """This prints the orbit infomrmation in a clean grid (unlike __str__) with only key information."""
        orbit_table = PrettyTable()
        orbit_table.field_names = ['Orbital Element', 'Value']
        orbit_table.add_row(["Orbit Type", orbit.orbit_type])
        orbit_table.add_row(["Position [km]", orbit.position_vector])
        orbit_table.add_row(["Velocity [km/s]", orbit.velocity_vector])
        orbit_table.add_row(["Semi-Major Axis [km]", orbit.semi_major_axis])
        orbit_table.add_row(["Eccentricity [None]", orbit.eccentricity])
        orbit_table.add_row(["True Anomaly [deg]", np.degrees(orbit.true_anomaly)])
        orbit_table.add_row(["Inclination [deg]", np.degrees(orbit.inclination)])
        orbit_table.add_row(["R.A. of Ascending Node [deg]", np.degrees(orbit.right_ascension_of_ascending_node)])
        orbit_table.add_row(["Argument of Periapsis [deg]", np.degrees(orbit.argument_of_periapsis)])
        print(orbit_table)
        

class OrbitUtilities:
    """This is a collection of utitity functions and constants for orbital calculations."""
    # Known Constants
    EARTH_RADIUS = 6378.1366  # km
    EARTH_MU = 3.986004354e5  # km^3/s^2
    EARTH_ROTATION_RATE = ωE = 0.000072921159 # rad/s

    @staticmethod
    def eccentric_anomaly_from_mean(M, ec, tolerance=1e-14):
        """This will return an estimated Eccentric Anomaly given a Mean Anomaly (M) and eccentricity (ec)"""
        # First make an initial guess at the Eccentric anomaly
        E = 1
        if ((-np.pi < M) and (M < 0)) or (M > np.pi):
            E = M - ec
        else:
            E = M + ec

        # Newton iteration to find eccentric anomaly [goal: find E so that f = 0 (i.e. find the roots for f)]
        def f(Ea):
            return Ea - ec * np.sin(Ea) - M # Kepler's Equ.
        def df(Ea):
            return 1 - ec * np.cos(Ea) # Derivative of Kepler's Equ.

        while np.abs(f(E)) > tolerance:
            E = E - f(E) / df(E)

        return E

    @staticmethod
    def hyperbolic_anomaly_from_mean(M, ec, tolerance=1e-14):
        """This will return an estimated Hyperbolic Anomaly given a Mean Anomaly (M) and eccentricity (ec)"""
        # First make an initial guess at the Eccentric anomaly
        H = M

        # Newton iteration to find hyperbolic anomaly [goal: find H so that f = 0 (i.e. find the roots for f)]
        def f(Ha):
            return ec * np.sinh(Ha) - M - Ha # Kepler's Equ.
        def df(Ha):
            return ec * np.cosh(Ha) - 1 # Derivative of Kepler's Equ.

        while np.abs(f(H)) > tolerance:
            H = H - f(H) / df(H)

        return H

    @staticmethod
    def calculate_velocity_gibbs(measured_positions):
        """This function uses Gibbs method of orbital determination to calculate orbital velocity in the ECI frame given three known positions in the ECI frame. This returns the velocity of the middle position measurement."""
        position1 = np.array(measured_positions[0])
        position2 = np.array(measured_positions[1])
        position3 = np.array(measured_positions[2])

        # using this weird extra definition because numpy is bugged and it says code isn't reachable otherwise.
        def cross(x, y): return np.cross(x, y)

        c1 = np.linalg.norm(cross(position2, position3)) / np.linalg.norm(cross(position1, position3))
        c3 = np.linalg.norm(cross(position2, position1)) / np.linalg.norm(cross(position3, position1))
        r1 = np.linalg.norm(position1)
        r2 = np.linalg.norm(position2)
        r3 = np.linalg.norm(position3)

        H = np.sqrt(OrbitUtilities.EARTH_MU * (r2 - c1 * r1 - c3 * r3) / (1 - c1 - c3))
        h_hat = cross(position1, position3) / np.linalg.norm(cross(position1, position3))
        angular_momentum_vector = H * h_hat

        eccentricity_vector = (((H**2 / OrbitUtilities.EARTH_MU - r1) * cross(position3, h_hat)) - ((H**2 / OrbitUtilities.EARTH_MU - r3) * cross(position1, h_hat))) / np.linalg.norm(cross(position1, position3))

        calculated_velocity = cross((OrbitUtilities.EARTH_MU / H**2) * angular_momentum_vector, eccentricity_vector + position2 / r2)

        return calculated_velocity

    @staticmethod
    def find_velocities_from_lambert_problem_p_iteration(position_1, position_2, TOF, short_direction=True, tolerance=1e-6):
        """This uses p-iteration to solve lambert's problem. Given two position vectors and the time of flight between them this will return the two coresponding velocities and the converged 'p' parameter."""
        position_1 = np.array(position_1)
        position_2 = np.array(position_2)
        r1 = np.linalg.norm(position_1)
        r2 = np.linalg.norm(position_2)
        delta_theta = np.arccos(np.dot(position_1, position_2) / (r1 * r2))

        if delta_theta < np.pi:
            # Short Way
            if not short_direction:
                delta_theta = np.pi + (np.pi - delta_theta)
        else:
            # Long Way
            if short_direction:
                delta_theta = np.pi - (delta_theta - np.pi)
            
        
        k = r1 * r2 * (1 - np.cos(delta_theta))
        m = r1 * r2 * (1 + np.cos(delta_theta))
        l = r1 + r2
        
        # Lagrange Functions from delta_theta
        def lagrange_f(current_p):
            return 1 - r2 / current_p * (1 - np.cos(delta_theta))
        def lagrange_f_dot(current_p):
            return np.sqrt(OrbitUtilities.EARTH_MU / current_p) * np.tan(delta_theta / 2) * ((1 - np.cos(delta_theta)) / current_p - 1 / r1 - 1 / r2)
        def lagrange_g(current_p):
            return r1 * r2 * np.sin(delta_theta) / np.sqrt(current_p * OrbitUtilities.EARTH_MU)
        def lagrange_g_dot(current_p):
            return 1 - r1 / current_p * (1 - np.cos(delta_theta))

        # Functions Used Directly in P-Iteration.
        def semi_major_axis(current_p):
            return (m * k * current_p) / ((2 * m - l**2) * current_p**2 + (2 * k * l * current_p) - k**2)
        def delta_E(current_p):
            d_E = np.arctan2(-1 * r1 * r2 * lagrange_f_dot(current_p), (1 - r1 / semi_major_axis(current_p) * (1 - lagrange_f(current_p))) * np.sqrt(OrbitUtilities.EARTH_MU * semi_major_axis(current_p)))
            if d_E < 0:
                return d_E + 2 * np.pi
            return d_E
        def TOF_error(TOF_i):
            return TOF_i - TOF
        def next_TOF(current_p):
            return lagrange_g(current_p) + np.sqrt(semi_major_axis(current_p)**3 / OrbitUtilities.EARTH_MU) * (delta_E(current_p) - np.sin(delta_E(current_p)))

        # Set Initial Conditions
        p_min = k / (l + np.sqrt(2 * m))
        p_max = k / (l - np.sqrt(2 * m))
        previous_p = 0.7 * p_min + 0.3 * p_max
        p = 0.3 * p_min + 0.7 * p_max
        previous_TOF = next_TOF(previous_p)
        current_TOF = next_TOF(p)

        # P-Iteration
        while np.abs(TOF_error(current_TOF)) > tolerance:
            tmp_previous_p = previous_p
            previous_p = p
            p = p - TOF_error(current_TOF) * ((p - tmp_previous_p) / (TOF_error(current_TOF)  - TOF_error(previous_TOF)))
            previous_TOF = current_TOF
            current_TOF = next_TOF(p)
        
        v1 = (position_2 - lagrange_f(p) * position_1) / lagrange_g(p)
        v2 = lagrange_f_dot(p) * position_1 + lagrange_g_dot(p) * v1
        
        return (v1, v2, p)

    @staticmethod
    def transform_position_SEZ_to_ECI(observation_lat_long, line_of_sight_elements):
        """This will take a location on earth and a position measurment of an orbit in South-East-Zenith(SEZ) coordinates and convert the position to ECI coordinates.
        'observation_lat_long' and 'line_of_sight_elements' are lists [latitude, longitude] (not earth lat, long) and [distance, elevation, azimuth] (Note: all angles are given in degrees and longitude is from I axis)."""
        lat = np.radians(observation_lat_long[0])
        long = np.radians(observation_lat_long[1])

        dist = line_of_sight_elements[0]
        ele = np.radians(line_of_sight_elements[1])
        azi = np.radians(line_of_sight_elements[2])

        r_site_hat = np.array([np.cos(lat) * np.cos(long), np.cos(lat) * np.sin(long), np.sin(lat)])
        r_site_vector = OrbitUtilities.EARTH_RADIUS * r_site_hat

        p_vector_hat = np.array([-1 * np.cos(ele) * np.cos(azi), np.cos(ele) * np.sin(azi), np.sin(ele)])
        p_vector = dist * p_vector_hat

        sez_to_eci_transformation_matrix = np.matrix([[np.sin(lat) * np.cos(long), -1 * np.sin(long), np.cos(lat) * np.cos(long)],
                                                      [np.sin(lat) * np.sin(long), np.cos(long), np.cos(lat) * np.sin(long)],
                                                      [-1 * np.cos(lat), 0, np.sin(lat)]])

        return r_site_vector + np.array(np.matmul(sez_to_eci_transformation_matrix, p_vector))[0]

    @staticmethod
    def propagate_true_anomaly_lagrange_perifocal(orbit, delta_true_anomaly) -> None:
        """This will use lagrange method to propagate position and velocity in perifocal frame given change in true anomaly."""
        if orbit.angle_type == 'deg':
            delta_true_anomaly = np.radians(delta_true_anomaly)

        r_v = TwoBodyKeplerOrbit.convert_position_and_velocity_to_perifocal_frame(orbit) # This is probably wrong.
        
        position = r_v[0]
        velocity = r_v[1]

        # Constants needed for Lagrange Coefficients
        r_0 = np.linalg.norm(position)
        v_0 = np.linalg.norm(velocity)
        v_r_0 = np.dot(position, velocity) / r_0
        h = r_0 * np.sqrt(v_0**2 - v_r_0**2)
        b = h**2 / OrbitUtilities.EARTH_MU
        r = b * (1 / (1 + (b / r_0 - 1) * np.cos(delta_true_anomaly) - (b / h) * v_r_0 * np.sin(delta_true_anomaly)))

        # Lagrange Coefficients -> Transformation Matrix -> Propagated Position and Velocity
        f = 1 - (1 / b) * r * (1 - np.cos(delta_true_anomaly))
        g = r * r_0 / h * np.sin(delta_true_anomaly)
        f_dot = (h / b) * ((1 - np.cos(delta_true_anomaly)) / np.sin(delta_true_anomaly)) * (1 / b * (1 - np.cos(delta_true_anomaly)) - 1 / r_0 - 1 / r)
        g_dot = 1 - 1 / b * r_0 * (1 - np.cos(delta_true_anomaly))

        new_position = f * position + g * velocity
        new_velocity = f_dot * position + g_dot * velocity

        return (new_position, new_velocity)
    
    @staticmethod
    def positions_from_line_of_sight_gauss(julian_dates, r_sites, line_of_sights, rootIndex=0, printRoots=True):
        """This returns a set of 3 positions given 3 line of site measurments."""
        delta_t1 = julian_dates[0] - julian_dates[1] # Days
        delta_t3 = julian_dates[2] - julian_dates[1] # Days
        delta_t1 = delta_t1 * 24 * 3600 # s
        delta_t3 = delta_t3 * 24 * 3600 # s

        a1 = delta_t3 / (delta_t3 - delta_t1)
        a3 = -1 * delta_t1 / (delta_t3 - delta_t1)
        a1_u = delta_t3 * (np.power(delta_t3 - delta_t1, 2) - np.power(delta_t3, 2)) / (6 * (delta_t3 - delta_t1))
        a3_u = -1 * delta_t1 * (np.power(delta_t3 - delta_t1, 2) - np.power(delta_t1, 2)) / (6 * (delta_t3 - delta_t1))

        line_of_sights = np.matrix(line_of_sights)        
        line_of_sights_inv = np.linalg.inv(line_of_sights)

        r_sites = np.array(r_sites)
        M = np.array(np.matmul(line_of_sights_inv, r_sites))
    
        d1 = M[1][0] * a1 - M[1][1] + M[1][2] * a3
        d2 = M[1][0] * a1_u + M[1][2] * a3_u
        
        C = np.array(np.dot(np.transpose(line_of_sights)[1], np.transpose(r_sites)[1]))[0][0]
        
        alpha = d1**2 + 2 * C * d1 + np.dot(np.transpose(r_sites)[1], np.transpose(r_sites)[1])
        beta = 2 * OrbitUtilities.EARTH_MU * (C * d2 + d1 * d2) 
        gamma = OrbitUtilities.EARTH_MU**2 * d2**2

        # Find the real roots of polynomial.
        coeff = [1, 0, -1 * alpha, 0, 0, -1 * beta, 0, 0, -1 * gamma]
        roots = np.roots(coeff)
        real_roots = []
        for root in roots:
            if np.isreal(root):
                real_roots.append(root.real)
        real_roots = np.array(real_roots)
        
        if printRoots:
            print("Gauss Roots: " + str(real_roots))

        u = OrbitUtilities.EARTH_MU / real_roots[rootIndex]**3

        c = np.array([a1 + a1_u * u, -1, a3 + a3_u * u])

        roe = np.array(np.matmul(M, -1 * c)) / c
        
        p1 = np.array(roe[0] * np.transpose(line_of_sights)[0] + np.transpose(r_sites)[0])[0]
        p2 = np.array(roe[1] * np.transpose(line_of_sights)[1] + np.transpose(r_sites)[1])[0]
        p3 = np.array(roe[2] * np.transpose(line_of_sights)[2] + np.transpose(r_sites)[2])[0]
                
        return p1, p2, p3
    
    @staticmethod # Unfinished.
    def positions_from_line_of_sight_laplace(julian_dates, r_sites, line_of_sights):
        """This returns a set of 3 positions given 3 line of site measurments."""
        delta_t1 = julian_dates[0] - julian_dates[1] # Days
        delta_t3 = julian_dates[2] - julian_dates[1] # Days
        delta_t1 = delta_t1 * 24 * 3600 # s
        delta_t3 = delta_t3 * 24 * 3600 # s

        roe2_dot = delta_t3**2 * line_of_sights[0] + (delta_t1**2 - delta_t3**2) * line_of_sights[1] - delta_t1**2 * line_of_sights[2] / (delta_t1 * delta_t3 * (delta_t3 - delta_t1))
        roe2_dot_dot = 2 * (-1 * delta_t3 * line_of_sights[0] + (delta_t3 - delta_t1) * line_of_sights[1] + delta_t1 * line_of_sights[2]) / (delta_t1 * delta_t3 * (delta_t3 - delta_t1)) 


    @staticmethod
    def line_of_sights_from_ra_and_dec(RAs, DECs):
        """This returns the corresponding line of site vectors in ECI frame from RAs and DECs."""

        line_of_sites = []
        for i in range(len(RAs)):
            # Assume both inputs are in deg.
            alpha = np.radians(RAs[i])
            delta = np.radians(DECs[i])
            line_of_sites.append([np.cos(delta) * np.cos(alpha), np.cos(delta) * np.sin(alpha), np.sin(delta)])
        line_of_sites = np.array(line_of_sites)

        return np.transpose(line_of_sites)
    
    @staticmethod
    def site_positions(lat, alt, LSTs):
        """This will return the corresponding site positions in ECI frame given site lat, long, and altitude with the time of each measurement."""
        # Assume inputs are in deg.
        phi = np.radians(lat)
        
        site_positions = []
        for lst in LSTs:
            lamb = np.radians(lst)
            r_site = (OrbitUtilities.EARTH_RADIUS + alt) * np.array([np.cos(phi) * np.cos(lamb), np.cos(phi) * np.sin(lamb), np.sin(phi)])
            site_positions.append(r_site)
        #site_positions = np.array(site_positions)
        #print(site_positions)
        return np.transpose(site_positions)

        
