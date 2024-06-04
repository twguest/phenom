import numpy as np
import pandas as pd
from pandas import DataFrame
from felpy.utils.vis_utils import Grids
from matplotlib.pyplot import cm
from skimage.measure import EllipseModel
from matplotlib.patches import Ellipse
from numpy import sin, cos

def solution_for_y(x, h, w, cx, cy, theta):
    """
    Solve for the y-coordinates of an ellipse given x.

    :param x: x-coordinate
    :param h: height of the ellipse
    :param w: width of the ellipse
    :param cx: horizontal center
    :param cy: vertical center
    :param theta: rotation angle
    :return: tuple of y-coordinates (lower, upper)
    """
    x -= cx
    a = ((sin(theta)**2) / w**2 + (cos(theta)**2) / h**2)
    b = 2 * sin(theta) * cos(theta) * ((1 / w**2) - (1 / h**2))
    c = sin(theta)**2 / w**2 + cos(theta)**2 / h**2
    discriminant = b**2 * x**2 - 4 * a * (c * x**2 - 1)
    if discriminant < 0:
        return np.nan, np.nan
    return (cy + ((-b * x) - np.sqrt(discriminant)) / (2 * a),
            cy + ((-b * x) + np.sqrt(discriminant)) / (2 * a))

def solution_from_df(df, pulse, val, method=None, dim='x'):
    """ 
    Return the (n+1)th pulse value from the value and pulse number of the nth pulse.

    :param df: DataFrame mapping the 95% confidence ellipse parameter for each pulse
    :param pulse: the nth pulse (int)
    :param val: the value of the nth pulse
    :param method: method for distribution (default is None, assumes normal distribution)
    :param dim: dimension ('x' or 'y')
    :return: value of the (n+1)th pulse
    """
    lower, upper = solution_for_y(val,
                                  h=df['Height'][pulse] / 2,
                                  w=df['Width'][pulse] / 2,
                                  cx=df['Centre'][pulse][0],
                                  cy=df['Centre'][pulse][1],
                                  theta=df['Angle'][pulse])
    
    while np.isnan(lower) and np.isnan(upper):
        val = df['Centre'][pulse][0]
        lower, upper = solution_for_y(val,
                                      h=df['Height'][pulse] / 2,
                                      w=df['Width'][pulse] / 2,
                                      cx=df['Centre'][pulse][0],
                                      cy=df['Centre'][pulse][1],
                                      theta=df['Angle'][pulse])

    std = abs(upper - lower) / (2 * df['Sigma'][pulse])
    
    if method is None:
        sol = np.random.normal(loc=(upper + lower) / 2, scale=std)

    return sol

def statengine_covariance_ellipse(arr, M=1, nstd=2, constraint=None, initial_condition=0):
    """
    Generate a statistical model based on covariance ellipses.

    :param arr: input data array
    :param M: number of pulse trains to generate
    :param nstd: number of standard deviations for ellipse
    :param constraint: constraint for generation (default None)
    :param initial_condition: initial condition for the first pulse
    :return: tuple of input and output DataFrames
    """
    df_input = extract_elliptical_stats(arr, nstd=nstd)
    outputs = {}

    if constraint is None:
        for m in range(M):
            output = []
            if isinstance(initial_condition, (float, int)):
                output.append(initial_condition)
            elif isinstance(initial_condition, list):
                output.append(initial_condition[m])

            for n in range(arr.shape[0] - 1):
                output.append(solution_from_df(df_input, pulse=n, val=output[n], method=None, dim='x'))

            outputs[f'{m}'] = output

    df_output = pd.DataFrame.from_dict(outputs)
    return df_input, df_output

def get_cov_ellipse_params(cov, centre, nstd, **kwargs):
    """
    Return a matplotlib Ellipse patch representing the covariance matrix.

    :param cov: covariance matrix
    :param centre: center of the ellipse
    :param nstd: number of standard deviations for scaling
    :return: list of ellipse parameters
    """
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    vx, vy = eigvecs[:, 0][0], eigvecs[:, 0][1]
    theta = np.arctan2(vy, vx)
    width, height = 2 * nstd * np.sqrt(eigvals)
    params = [centre, width, height, theta]
    return params

def sequential_pairs(arr, n=0):
    """
    Generate sequential pairs from an array.

    :param n: starting index
    :return: array of sequential pairs
    """
    M = arr.shape[1]
    return np.array([(arr[n, m], arr[n + 1, m]) for m in range(M - 1)])

def extract_elliptical_stats(arr, nstd=2):
    """
    Extract elliptical statistics from an array.

    :param arr: input data array
    :param nstd: number of standard deviations for ellipse
    :return: DataFrame of elliptical statistics
    """
    pre_ = []
    cov = confidence_ellipse(arr, nstd=nstd)
    for n in range(len(cov)):
        pre_.append([n, *cov[n], nstd])
    df = DataFrame(pre_, columns=['Pulse', 'Centre', 'Width', 'Height', 'Angle', 'Sigma'])
    df.set_index('Pulse')
    return df

def confidence_ellipse(arr, nstd=2):
    """
    Compute covariance ellipses from an array.

    :param arr: input data array
    :param nstd: number of standard deviations for ellipse
    :return: list of covariance ellipse parameters
    """
    N = arr.shape[0] - 1
    M = arr.shape[1]
    cov = []

    for n in range(N):
        coordinates = sequential_pairs(arr, n=n)
        N0, N1 = coordinates[:, 0], coordinates[:, 1]
        cov_matrix = np.cov(N0, N1)
        cov.append(get_cov_ellipse_params(cov_matrix, (N0.mean(), N1.mean()), nstd=nstd))

    return cov 

def plot_covariance_ellipse(arr, cov=None, grid=None, N=0, axis=0, link_dots=False, add_colorbar=True):
    """
    Plot covariance ellipse of a given array.

    :param arr: input data array
    :param cov: covariance ellipse parameters
    :param grid: grid object for plotting
    :param N: index of the ellipse to plot
    :param axis: axis index for plotting
    :param link_dots: flag to link dots with lines
    :param add_colorbar: flag to add colorbar
    :return: grid object
    """
    color = cm.rainbow(np.linspace(0, 1, arr.shape[1] - 1))

    if cov is None:
        cov = confidence_ellipse(arr)
    
    if grid is None:
        grid = Grids(scale=2, global_aspect=1)
        grid.create_grid(n=1, m=1, share_x=False, share_y=False)
        grid.pad(4)
    
    if add_colorbar:
        grid.add_global_colorbar(vmin=0, vmax=arr.shape[-1], cmap='rainbow', clabel="Train No.", fontsize=22)
    
    ax = grid.axes if isinstance(grid.axes, list) else grid.axes[list(grid.axes.keys())[axis]] if isinstance(grid.axes, dict) else grid.axes
    grid.set_fontsize(22)
    
    for m in range(arr.shape[-1] - 1):
        c = color[m]
        coordinates = sequential_pairs(arr, n=N)
        N0, N1 = coordinates[:, 0], coordinates[:, 1]
        ax.scatter(N0[m], N1[m], color=c)
    
    if link_dots:
        ax.plot(coordinates[:, 0], coordinates[:, 1], linewidth=0.25)

    ex = Ellipse(xy=cov[N][0], width=cov[N][1], height=cov[N][2], angle=np.degrees(cov[N][3]), fill=None)
    ax.add_artist(ex)

    return grid
