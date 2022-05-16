import pandas as pd


def show_versions():
    """Method to print the version of dependencies.

    """
    from pyet import __version__ as ps_version
    from pandas import __version__ as pd_version
    from numpy import __version__ as np_version
    from scipy import __version__ as sc_version
    from matplotlib import __version__ as mpl_version
    from sys import version as os_version

    msg = (
        f"Python version: {os_version}\n"
        f"Numpy version: {np_version}\n"
        f"Scipy version: {sc_version}\n"
        f"Pandas version: {pd_version}\n"
        f"Matplotlib version: {mpl_version}\n"
        f"Pyet version: {ps_version}"
    )
    return print(msg)


def get_index_shape(df):
    """Method to return the index and shape of the input data.

    """
    try:
        index = pd.DatetimeIndex(df.index)
<<<<<<< HEAD
    except AttributeError:
=======
    except:
>>>>>>> 9dd04566d670bf640b4a10154e8430ae75f9c615
        index = pd.DatetimeIndex(df.time)
    return index, df.shape
