"""This library holds tools for use with orbit calculations and representations."""
import numpy


class KeplerOrbitTypes:
    """Different characterizations of Kepler Orbits."""
    CIRCULAR = 'Circular Orbit'
    ELIPTICAL = 'Eliptical Orbit'
    PARABOLIC = 'Parabolic Orbit'
    HYPERBOLIC = 'Hyperbolic Orbit'
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
    def __init__(self, radial_vector, velocity_vector) -> None:
        self.radial_vector = radial_vector
        self.velocity_vector = velocity_vector
        self.angular_momentum = 0
        self.total_energy = 0
        self.semi_major_axis = 0
        self.semi_minor_axis = 0
        self.eccentricity = 0

    def update(self, radial_vector, velocity_vector) -> None:
        """Not sure about this yet."""

    def __str__(self) -> str:
        pass
