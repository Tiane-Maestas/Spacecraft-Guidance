"""This library holds tools for use with orbit calculations and representations."""
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

    def __init__(self, position_vector, velocity_vector):
        """Initialize from a measured position and velocity in an orbital frame. {r, theta, h}"""
        # only in radial direction (r, 0, 0)
        self.position_vector = numpy.array(position_vector)
        # only in radial and theta (r, theta, 0)
        self.velocity_vector = numpy.array(velocity_vector)

        self.calculate_orbit_info()

        # Orientation variables can manually be set if desired. Since we are in the 'in-plane' orbital frame these tell orientation in 3D.
        self.longitude_of_ascending_node = None
        self.inclination = None
        self.argument_of_periapsis = None
    
    @classmethod
    def build_from_known_params(cls, semi_major_axis, eccentricity, true_anomaly):
        """Initialize from known orbital parameters."""
        # Note: maybe add angle type so they can pass in deg or rad.
        true_anomaly = numpy.radians(true_anomaly)
        H = numpy.sqrt(TwoBodyKeplerOrbit.EARTH_MU * semi_major_axis * (1 - numpy.power(eccentricity, 2)))
        
        # Calculate position vector.
        r = semi_major_axis * (1 - numpy.power(eccentricity, 2)) / (1 + eccentricity * numpy.cos(true_anomaly))
        position_vector = [r, 0, 0]

        # Calculate velocity vector.
        v_r = (eccentricity * TwoBodyKeplerOrbit.EARTH_MU / H) * numpy.sin(true_anomaly)
        v_theta = (TwoBodyKeplerOrbit.EARTH_MU / H) * (1 + eccentricity * numpy.cos(true_anomaly))
        velocity_vector = [v_r, v_theta, 0]

        return cls(position_vector, velocity_vector)

    @classmethod
    def build_from_inertial_frame(cls, position_vector, velocity_vector):
        """Initialize from measured position and velocity in cartesian coordinates."""

        # Calculate orbital elements in from cartesian frame.
        r = numpy.sqrt(numpy.dot(position_vector, position_vector))
        v_squared = numpy.dot(velocity_vector, velocity_vector)
        h_vector = numpy.cross(position_vector, velocity_vector)
        h = numpy.sqrt(numpy.dot(h_vector, h_vector))

        semi_major_axis = 1 / (2 / r - v_squared / TwoBodyKeplerOrbit.EARTH_MU)

        eccentricity_vector = numpy.cross(velocity_vector, h_vector) / TwoBodyKeplerOrbit.EARTH_MU - position_vector / r
        eccentricity = numpy.linalg.norm(eccentricity_vector)

        # Calcualte orientation of perifocal frame.
        e_hat = eccentricity_vector / eccentricity
        h_hat = h_vector / h
        e_perpendicular_hat = numpy.cross(h_hat, e_hat)

        true_anomaly = numpy.arctan2(numpy.dot(numpy.cross(e_hat, e_perpendicular_hat), h_hat), numpy.dot(e_hat, e_perpendicular_hat))

        orbit = TwoBodyKeplerOrbit.build_from_known_params(semi_major_axis, eccentricity, true_anomaly)

        orbit.longitude_of_ascending_node = numpy.degrees(numpy.arctan2(h_hat[0], -h_hat[1]))
        orbit.inclination = numpy.degrees(numpy.arccos(h_hat[2]))
        orbit.argument_of_periapsis = numpy.degrees(numpy.arctan2(e_hat[2], e_perpendicular_hat[2]))

        return orbit

    @staticmethod
    def find_state_vectors_from_orbit(semi_major_axis, eccentricity, inclination, longitude_of_ascending_node, argument_of_periapsis, true_anomaly):
        """Returns a tuple containing the position and velocity vectors coresponding to these orbital parameters in the inertial frame.
           Can be used to re-build orbit in this 'build_from_inertial_frame' class method."""
        # Note: maybe add angle type so they can pass in deg or rad.
        inclination = numpy.radians(inclination)
        longitude_of_ascending_node = numpy.radians(longitude_of_ascending_node)
        argument_of_periapsis = numpy.radians(argument_of_periapsis)
        true_anomaly = numpy.radians(true_anomaly)

        # Calculate position vector.
        r = semi_major_axis * (1 - numpy.power(eccentricity, 2)) / (1 + eccentricity * numpy.cos(true_anomaly))
        r_x = numpy.cos(longitude_of_ascending_node) * numpy.cos(true_anomaly) - numpy.sin(longitude_of_ascending_node) * numpy.sin(true_anomaly) * numpy.cos(inclination)
        r_y = numpy.sin(longitude_of_ascending_node) * numpy.cos(true_anomaly) + numpy.cos(longitude_of_ascending_node) * numpy.sin(true_anomaly) * numpy.cos(inclination)
        r_z = numpy.sin(true_anomaly) * numpy.sin(inclination)
        r_hat = numpy.array([r_x, r_y, r_z])
        position_vector = r * r_hat

        # Calculate velocity vector.
        H = numpy.sqrt(TwoBodyKeplerOrbit.EARTH_MU * semi_major_axis * (1 - numpy.power(eccentricity, 2)))
        v = (-TwoBodyKeplerOrbit.EARTH_MU / H)
        v_x = numpy.cos(longitude_of_ascending_node) * (numpy.sin(true_anomaly) + eccentricity * numpy.sin(argument_of_periapsis)) + numpy.sin(longitude_of_ascending_node) * (numpy.cos(true_anomaly) + eccentricity * numpy.cos(argument_of_periapsis)) * numpy.cos(inclination)
        v_y = numpy.sin(longitude_of_ascending_node) * (numpy.sin(true_anomaly) + eccentricity * numpy.sin(argument_of_periapsis)) - numpy.cos(longitude_of_ascending_node) * (numpy.cos(true_anomaly) + eccentricity * numpy.cos(argument_of_periapsis)) * numpy.cos(inclination)
        v_z = -(numpy.cos(true_anomaly) + eccentricity * numpy.cos(argument_of_periapsis)) * numpy.sin(inclination)
        v_hat = numpy.array([v_x, v_y, v_z])
        velocity_vector = v * v_hat

        return (position_vector, velocity_vector)

    def update(self, radial_vector, velocity_vector) -> None:
        """Not sure about this yet."""

    def calculate_orbit_info(self) -> None:
        """This calculates all of the info that can be displayed via this 'to string' method."""
        # using this weird extra definition because numpy is bugged and it says code isn't reachable otherwise.
        def cross(x, y): return numpy.cross(x, y)
        # the final dimmension (0, 0, h). Angular momentum is orthogonal to orbital plane.
        self.angular_momentum = cross(self.position_vector, self.velocity_vector)

        self.total_energy = 0.5 * numpy.dot(self.velocity_vector, self.velocity_vector) - TwoBodyKeplerOrbit.EARTH_MU / numpy.linalg.norm(self.position_vector)

        self.semi_major_axis = -TwoBodyKeplerOrbit.EARTH_MU / (2 * self.total_energy)

        self.parameter = numpy.dot(self.angular_momentum, self.angular_momentum) / TwoBodyKeplerOrbit.EARTH_MU

        self.eccentricity_vector = 1 / TwoBodyKeplerOrbit.EARTH_MU * cross(self.velocity_vector, self.angular_momentum) - self.position_vector / numpy.linalg.norm(self.position_vector)
        self.eccentricity = numpy.sqrt(1 - self.parameter / self.semi_major_axis) # this is the same as the norm of the eccentricity vector.
        self.orbit_type = KeplerOrbitTypes.get_orbit_type(self.eccentricity)

        self.semi_minor_axis = self.semi_major_axis * numpy.sqrt(1 - self.eccentricity**2)

        self.period = 2 * numpy.pi * numpy.sqrt(numpy.power(self.semi_major_axis, 3)/TwoBodyKeplerOrbit.EARTH_MU)

        # Closest point
        self.perigee = self.semi_major_axis * (1 - self.eccentricity)

        # Furthest point
        self.apogee = self.semi_major_axis * (1 + self.eccentricity)

        argument = numpy.clip(1 / self.eccentricity * (self.parameter / numpy.linalg.norm(self.position_vector) - 1), -1, 1) # in case of rounding errors.
        self.true_anomaly = numpy.degrees(numpy.arccos(argument))

        self.flight_path_angle = numpy.degrees(numpy.arctan(self.eccentricity * numpy.sin(self.true_anomaly) / (1 + self.eccentricity * numpy.cos(self.true_anomaly))))

    ORBIT_INFO_FORMAT = "\n-----Orbit INFO-----\nOrbit Type: {orbit_type}\nPosition: {position}\nVelocity: {velocity}\nAngular Momentum(H): {angular_momentum}\nTotal Energy(E): {total_energy}\nSemi-Major Axis(a): {semi_major_axis}\nSemi-Minor Axis(b): {semi_minor_axis}\nParameter(p): {parameter}\nEccentricity(e): {eccentricity}\nPeriod(T): {period}\nPerigee: {perigee}\nApogee: {apogee}\nTrue Anomaly(f): {true_anomaly}\nFlight Path Angle(gamma): {flight_path_angle}\n \
                         \n-----Orientation INFO-----\nLongitude of Ascending Node(Omega): {longitude_of_ascending_node}\nInclination(i): {inclination}\nArgument of Periapsis(w): {argument_of_periapsis}\n"
    def __str__(self) -> str:
        return TwoBodyKeplerOrbit.ORBIT_INFO_FORMAT.format(orbit_type=self.orbit_type,
                                                           position=self.position_vector,
                                                           velocity=self.velocity_vector,
                                                           angular_momentum=self.angular_momentum,
                                                           total_energy=self.total_energy,
                                                           semi_major_axis=self.semi_major_axis,
                                                           semi_minor_axis=self.semi_minor_axis,
                                                           parameter=self.parameter,
                                                           eccentricity=self.eccentricity,
                                                           period=self.period,
                                                           perigee=self.perigee,
                                                           apogee=self.apogee,
                                                           true_anomaly=self.true_anomaly,
                                                           flight_path_angle=self.flight_path_angle,
                                                           longitude_of_ascending_node=self.longitude_of_ascending_node,
                                                           inclination=self.inclination,
                                                           argument_of_periapsis=self.argument_of_periapsis)
