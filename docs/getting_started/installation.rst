Installing and Updating Phydrus
===============================

Installing Python
-----------------
To install Phydrus, a working version of Python 3.7 or higher has to be
installed. We recommend using the `Anaconda Distribution
<https://www.continuum.io/downloads>`_ of Python. This Python distribution
includes most of the python package dependencies and the Jupyter Lab
software to run the notebooks. Moreover, it includes the Graphical User
Interface (GUI) Spyder to start scripting in Python. However, you are free
to install any Python distribution you want.

Installing the Phydrus package
------------------------------
The latest stable version of the Phydrus package is available from the Pypi
package index.

>>> pip install phydrus

To install in developer mode, clone the GitHub repository and use the
following syntax:

>>> pip install -e .

Updating the Phydrus package
----------------------------
If you have already installed Phydrus, it is possible to update Phydrus
easily. To update, open a Windows command screen or a Mac terminal and type:

>>> pip install phydrus --upgrade

Dependencies
------------
Phydrus depends on a number of Python packages, which are all automatically
installed when using the pip install manager. The following packages are
necessary for the installation of Phydrus:

.. include:: ../../requirements.txt
    :literal: