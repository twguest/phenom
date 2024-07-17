# -*- coding: utf-8 -*-
import numpy as np


class Mesh:
    """
    Generalized Mesh class.

    This class can be constructed via two methods: from single-value field definitions or from
    numpy arrays defining exact coordinates.

    Parameters
    ----------
    nx : int, optional
        Number of points in the x direction, by default 512.
    ny : int, optional
        Number of points in the y direction, by default 512.
    nt : int, optional
        Number of points in the t direction, by default 10.
    dx : float, optional
        Sampling interval in the x direction, by default 1.
    dy : float, optional
        Sampling interval in the y direction, by default 1.
    dt : float, optional
        Sampling interval in the t direction, by default 1.
    x : np.ndarray, optional
        Array of x coordinates, by default None.
    y : np.ndarray, optional
        Array of y coordinates, by default None.
    t : np.ndarray, optional
        Array of t coordinates, by default None.
    """

    def __init__(self, nx=512, ny=512, nt=10, dx=1, dy=1, dt=1, x=None, y=None, t=None):
        self.x = x
        self.y = y
        self.t = t

        self.build_from_array()

    def sampling(self):
        """
        Returns and records the field sampling in each of the spatial
        dimensions.
        """
        if self.ndims >= 2:
            self.dx = (self.xMax - self.xMin) / self.nx
            self.dy = (self.yMax - self.yMin) / self.ny

            if self.ndims == 3:
                self.dt = (self.tMax - self.tMin) / self.nt

    def get_sampling(self):
        """
        Returns the two or three dimensional sampling of the field.

        Returns
        -------
        tuple
            Sampling intervals in each dimension.
        """
        if self.ndims >= 2:
            return self.dx, self.dy
        elif self.ndims == 3:
            return self.dx, self.dy, self.dt

    def build_from_array(self):
        """
        Constructs a mesh object from linear-space variables.
        """
        self.nx = len(self.x)
        self.xMin = np.min(self.x)
        self.xMax = np.max(self.x)
        self.dx = (self.xMax - self.xMin) / 2

        self.ny = len(self.y)
        self.yMin = np.min(self.y)
        self.yMax = np.max(self.y)
        self.dy = (self.yMax - self.yMin) / 2

        self.nt = len(self.t)
        self.tMin = np.min(self.t)
        self.tMax = np.max(self.t)
        self.dt = (self.tMax - self.tMin) / 2

    def get_extent(self):
        """
        Returns wpg style extent with respect to the transverse axis.

        Returns
        -------
        list
            Extent of the mesh in the form [xMin, xMax, yMin, yMax].
        """
        return [self.xMin, self.xMax, self.yMin, self.yMax]

    def get_array(self, axes=0):
        """
        Returns the mesh in array form.

        Parameters
        ----------
        axes : int or str, optional
            Axis to retrieve the array for, by default 0.

        Returns
        -------
        np.ndarray
            Array of coordinates for the specified axis.
        """
        self.check_axes_options(axes)

        if axes in ["x", 0]:
            return np.linspace(self.xMin, self.xMax, self.nx)
        if axes in ["y", 1]:
            return np.linspace(self.yMin, self.yMax, self.ny)
        if axes in ["t", 2]:
            return np.linspace(self.tMin, self.tMax, self.nt)

    def check_axes_options(self, axes):
        """
        Checks the axis options for validity.

        Parameters
        ----------
        axes : int or str
            Axis to check.

        Raises
        ------
        AssertionError
            If the axis is not valid.
        """
        options = ["x", "y", "t", 0, 1, 2]

        assert axes in options, f"the supplied axes should be in {options}"

    @property
    def meshgrid(self):
        """
        Returns a 2D numpy style meshgrid of spatial coordinates.

        Returns
        -------
        tuple of np.ndarray
            Meshgrid arrays for the x and y coordinates.
        """
        return np.meshgrid(self.get_array(0), self.get_array(1))

    def __update__(self, **kwargs):
        """
        Updates the attributes of the mesh via various methods.

        Parameters
        ----------
        **kwargs
            Arbitrary keyword arguments.
        """
        if "wfr" in kwargs:
            for o in [x for x in list(self.__dict__.keys()) if x in dir(kwargs["wfr"].params.Mesh)]:
                setattr(self, o, getattr(kwargs["wfr"].params.Mesh, o))


if __name__ == "__main__":
    m = Mesh(x=np.linspace(10, 10, 10), y=np.linspace(10, 10, 10))
    print(m)
