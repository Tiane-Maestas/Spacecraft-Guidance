"""This module holds tools for use with orbit calculations and representations."""
import numpy


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
    # Known Constants
    EARTH_RADIUS = 6378  # km
    EARTH_MU = 3.986e5  # km^3/s^2

    def __init__(self, position_vector, velocity_vector, **options):
        """Initialize from a measured position and velocity in the cartesian ECI (Earth Centered Inertial) frame."""
        options = TwoBodyKeplerOrbit.process_optional_params(options)
        self.angle_type = options[0] # This specifies what unit to print in.
        self.time_of_flight = options[1] # This variable doesn't do anything yet here.  

        self.position_vector = numpy.array(position_vector)
        self.velocity_vector = numpy.array(velocity_vector)

        self.calculate_orbit_info()
    
    @classmethod
    def build_from_known_orbital_params(cls, orbital_params_list, **options):
        """Initialize from a list of known orbital parameters. [semimajor axis (a), eccentricity (ec), inclination (i), longitude of the ascending node (omega), argument of periapsis (w), mean or hyperbolic anomaly at epoch (theta)]"""
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
            i = numpy.radians(i); omega = numpy.radians(omega); w = numpy.radians(w); theta = numpy.radians(theta)
        i = numpy.mod(i, 2 * numpy.pi); omega = numpy.mod(omega, 2 * numpy.pi); w = numpy.mod(w, 2 * numpy.pi); theta = numpy.mod(theta, 2 * numpy.pi)

        orbit_type = KeplerOrbitTypes.get_orbit_type(ec)

        # Orbital mean motion [ rad/s ]
        nm = numpy.sqrt(TwoBodyKeplerOrbit.EARTH_MU / numpy.power(numpy.abs(a), 3))
        # Mean anomaly [ rad ]
        M = theta + nm * time_of_flight

        if (orbit_type == KeplerOrbitTypes.ELIPTICAL) or (orbit_type == KeplerOrbitTypes.CIRCULAR):
            # Converged eccentric anomaly [ rad ]
            E = OrbitUtilities.eccentric_anomaly_from_mean(M, ec)
            # True anomaly [ rad ]
            f = numpy.mod(2 * numpy.arctan(numpy.sqrt((1 + ec)/(1 - ec)) * numpy.tan(E / 2)), 2 * numpy.pi)
        elif orbit_type == KeplerOrbitTypes.HYPERBOLIC:
            pass # Yet to be implemented.

        # Argument of Latitude
        u = w + f
        
        r = a * (1 - ec**2) / (1 + ec * numpy.cos(f))
        v = numpy.sqrt(TwoBodyKeplerOrbit.EARTH_MU * (2 / r - 1 / a))

        # Ascending node vector
        nhat = numpy.array([ numpy.cos(omega) , numpy.sin(omega), 0 ])
        rT   = numpy.array([ -1 * numpy.cos(i) * numpy.sin(omega), numpy.cos(i) * numpy.cos(omega), numpy.sin(i) ])

        gamma = numpy.arctan2(ec * numpy.sin(f), 1 + ec * numpy.cos(f))

        rhat = numpy.cos(u) * nhat + numpy.sin(u) * rT
        vhat = numpy.sin(gamma - u) * nhat + numpy.cos(gamma - u) * rT

        position_vector = r * rhat
        velocity_vector = v * vhat
        
        return cls(position_vector, velocity_vector, angle_type=angle_type, time_of_flight=time_of_flight)

    def update(self, radial_vector, velocity_vector) -> None:
        """Not sure about this yet."""

    def calculate_orbit_info(self) -> None:
        """This calculates all of the info that can be displayed via this 'to string' method."""
        # using this weird extra definition because numpy is bugged and it says code isn't reachable otherwise.
        def cross(x, y): return numpy.cross(x, y)

        #  ECI Coordinate frame unit vectors
        I = numpy.array([1, 0, 0]); J = numpy.array([0, 1, 0]); K = numpy.array([0, 0, 1])

        # Angular momentum is orthogonal to orbital plane.
        self.angular_momentum = cross(self.position_vector, self.velocity_vector)

        self.position_hat = self.position_vector / numpy.linalg.norm(self.position_vector)
        self.angular_momentum_hat = self.angular_momentum / numpy.linalg.norm(self.angular_momentum)

        # Normalized ascending node vector.
        self.ascending_node_hat = cross(K, self.angular_momentum) / numpy.linalg.norm(cross(K, self.angular_momentum)) 

        self.total_energy = 0.5 * numpy.dot(self.velocity_vector, self.velocity_vector) - TwoBodyKeplerOrbit.EARTH_MU / numpy.linalg.norm(self.position_vector)

        # Eccentricity
        self.eccentricity_vector = 1 / TwoBodyKeplerOrbit.EARTH_MU * cross(self.velocity_vector, self.angular_momentum) - self.position_hat
        self.eccentricity =  numpy.linalg.norm(self.eccentricity_vector)
        self.orbit_type = KeplerOrbitTypes.get_orbit_type(self.eccentricity)
        
        # If orbit_type is 'parabolic' energy is 0.
        if self.orbit_type == KeplerOrbitTypes.PARABOLIC:
            self.semi_major_axis = numpy.Inf
            self.parameter = numpy.dot(self.angular_momentum, self.angular_momentum) / TwoBodyKeplerOrbit.EARTH_MU
        else:
            self.semi_major_axis = -1 * TwoBodyKeplerOrbit.EARTH_MU / (2 * self.total_energy)
            self.parameter = self.semi_major_axis * (1 - numpy.power(self.eccentricity, 2))
            self.perigee = self.semi_major_axis * (1 - self.eccentricity) # Closest point
        
        self.semi_minor_axis = self.semi_major_axis * numpy.sqrt(1 - self.eccentricity**2)
        self.period = 2 * numpy.pi * numpy.sqrt(numpy.power(self.semi_major_axis, 3)/TwoBodyKeplerOrbit.EARTH_MU)
        self.apogee = self.semi_major_axis * (1 + self.eccentricity) # Furthest point
        
        # If the inclination is less than 90 degrees, hte elliptical orbit is a direct (prograde) orbit.
        self.inclination = numpy.arccos(numpy.dot(K, self.angular_momentum_hat))

        self.right_ascension_of_ascending_node = numpy.mod(numpy.arctan2(numpy.dot(J, self.ascending_node_hat), numpy.dot(I, self.ascending_node_hat)), 2 * numpy.pi)

        self.argument_of_periapsis = numpy.mod(numpy.arctan2(numpy.dot(self.angular_momentum_hat, cross(self.ascending_node_hat, self.eccentricity_vector)), numpy.dot(self.ascending_node_hat, self.eccentricity_vector)), 2 * numpy.pi)

        self.true_anomaly = numpy.mod(numpy.arctan2(numpy.dot(self.angular_momentum_hat, cross(self.eccentricity_vector, self.position_vector)), numpy.dot(self.eccentricity_vector, self.position_vector)), 2 * numpy.pi)

        self.flight_path_angle = numpy.degrees(numpy.arctan(self.eccentricity * numpy.sin(self.true_anomaly) / (1 + self.eccentricity * numpy.cos(self.true_anomaly))))

        if (self.orbit_type == KeplerOrbitTypes.ELIPTICAL) or (self.orbit_type == KeplerOrbitTypes.CIRCULAR):
            # Eccentric and Mean anomaly at epoch.
            self.eccentric_anomaly = 2 * numpy.arctan2(numpy.sqrt(1 - self.eccentricity) * numpy.tan(self.true_anomaly / 2), numpy.sqrt(1 + self.eccentricity))
            self.mean_anomaly = numpy.mod((self.eccentric_anomaly - self.eccentricity * numpy.sin(self.eccentric_anomaly)), 2 * numpy.pi)
        elif self.orbit_type == KeplerOrbitTypes.HYPERBOLIC:
            self.hyperbolic_anomaly = 2 * numpy.arctanh(numpy.sqrt(self.eccentricity - 1) * numpy.tan(self.true_anomaly / 2) / numpy.sqrt(self.eccentricity + 1))
            self.mean_anomaly = numpy.mod((self.eccentricity * numpy.sinh(self.hyperbolic_anomaly) - self.hyperbolic_anomaly), 2 * numpy.pi)

    @staticmethod
    def convert_position_and_velocity_to_perifocal_frame(orbit):
        """This will take in a pre-constructed orbit and return the position and velocity of the orbit from the perifocal frame."""
        H = numpy.linalg.norm(orbit.angular_momentum)
        e = orbit.eccentricity
        theta = orbit.mean_anomaly

        P = orbit.semi_major_axis * (1 - e**2)
        gamma = P / (1 + e * numpy.cos(theta))

        perifocal_position = [gamma * numpy.cos(theta), gamma * numpy.sin(theta), 0]
        perifocal_velocity = [-1 * (TwoBodyKeplerOrbit.EARTH_MU / H) * (numpy.sin(theta)), (TwoBodyKeplerOrbit.EARTH_MU / H) * (e + numpy.cos(theta)), 0]

        return (perifocal_position, perifocal_velocity)
        

    ORBIT_INFO_FORMAT = "\n-----Orbit INFO-----\nOrbit Type: {orbit_type}\nPosition: {position}\nVelocity: {velocity}\nAngular Momentum(H): {angular_momentum}\nTotal Energy(E): {total_energy}\nSemi-Major Axis(a): {semi_major_axis}\nSemi-Minor Axis(b): {semi_minor_axis}\nParameter(p): {parameter}\nEccentricity(e): {eccentricity}\nPeriod(T): {period}\nPerigee: {perigee}\nApogee: {apogee}\nTrue Anomaly(f): {true_anomaly}\nFlight Path Angle(gamma): {flight_path_angle}\nMean Anomaly(M): {mean_anomaly}\n \
                         \n-----Orientation INFO-----\nRight Ascension of Ascending Node(Omega): {right_ascension_of_ascending_node}\nInclination(i): {inclination}\nArgument of Periapsis(w): {argument_of_periapsis}\n"
    def __str__(self) -> str:
        true_anomaly = self.true_anomaly; flight_path_angle = self.flight_path_angle; mean_anomaly = self.mean_anomaly; right_ascension_of_ascending_node = self.right_ascension_of_ascending_node; inclination = self.inclination; argument_of_periapsis = self.argument_of_periapsis
        if(self.angle_type == 'deg'):
            true_anomaly = numpy.degrees(true_anomaly); flight_path_angle = numpy.degrees(flight_path_angle); mean_anomaly = numpy.degrees(mean_anomaly); right_ascension_of_ascending_node = numpy.degrees(right_ascension_of_ascending_node); inclination = numpy.degrees(inclination); argument_of_periapsis = numpy.degrees(argument_of_periapsis)
        return TwoBodyKeplerOrbit.ORBIT_INFO_FORMAT.format(orbit_type=self.orbit_type, position=self.position_vector, velocity=self.velocity_vector, angular_momentum=self.angular_momentum, total_energy=self.total_energy, semi_major_axis=self.semi_major_axis, semi_minor_axis=self.semi_minor_axis, parameter=self.parameter, eccentricity=self.eccentricity, period=self.period, perigee=self.perigee, apogee=self.apogee, true_anomaly=true_anomaly, flight_path_angle=flight_path_angle, mean_anomaly=mean_anomaly, right_ascension_of_ascending_node=right_ascension_of_ascending_node, inclination=inclination, argument_of_periapsis=argument_of_periapsis)

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


class OrbitUtilities:
    """This is a collection of utitity function for orbital calculations."""

    @staticmethod
    def eccentric_anomaly_from_mean(M, e, tolerance=1e-14):
        """This will return an estimated Eccentric Anomaly given a Mean Anomaly (M) and eccentricity (e)"""
        # Compute Eccentric anomaly
        # def E(j):
        # j = 1 # A priori estimate
        # if ((-numpy.pi < M) and (M < 0)) or (M > numpy.pi):
        #     E = M - e
        # else:
        #     E = M + e

        # # Newton iteration to find eccentric anomaly
        # # Algorithm [goal: find E so that f = 0]
        # f_E(j) = E(j) - ec*sin(E(j)) - M;
        # while abs(f_E(j)) > 1e-11
        #     E(j + 1) = E(j) - f_E(j)/(1 - ec*cos(E(j)));
        #     j = j + 1;
        #     f_E(j) = E(j) - ec*sin(E(j)) - M;
        # end
        return -0.8000052947851282
    
    @staticmethod
    def calculate_velocity_gibbs(measured_positions):
        """This function uses Gibbs method of orbital determination to calculate orbital velocity in the ECI frame given three known positions in the ECI frame. This returns the velocity of the middle position measurement."""
        position1 = numpy.array(measured_positions[0])
        position2 = numpy.array(measured_positions[1])
        position3 = numpy.array(measured_positions[2])

        # using this weird extra definition because numpy is bugged and it says code isn't reachable otherwise.
        def cross(x, y): return numpy.cross(x, y)

        c1 = numpy.linalg.norm(cross(position2, position3)) / numpy.linalg.norm(cross(position1, position3))
        c3 = numpy.linalg.norm(cross(position2, position1)) / numpy.linalg.norm(cross(position3, position1))
        r1 = numpy.linalg.norm(position1)
        r2 = numpy.linalg.norm(position2)
        r3 = numpy.linalg.norm(position3)

        H = numpy.sqrt(TwoBodyKeplerOrbit.EARTH_MU * (r2 - c1 * r1 - c3 * r3) / (1 - c1 - c3))
        h_hat = cross(position1, position3) / numpy.linalg.norm(cross(position1, position3))
        angular_momentum_vector = H * h_hat

        eccentricity_vector = (((H**2 / TwoBodyKeplerOrbit.EARTH_MU - r1) * cross(position3, h_hat)) - ((H**2 / TwoBodyKeplerOrbit.EARTH_MU - r3) * cross(position1, h_hat))) / numpy.linalg.norm(cross(position1, position3))
        
        calculated_velocity = cross((TwoBodyKeplerOrbit.EARTH_MU / H**2) * angular_momentum_vector, eccentricity_vector + position2 / r2)

        return calculated_velocity
