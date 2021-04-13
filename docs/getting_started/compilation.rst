Compiling the Fortran Source Code
=================================

After installing the Phydrus package, you need to obtain an executable for
the Hydrus-1D model. There after two options to obtain this executable: **1)**
download a pre-compiled version, or **2)** compiled the fortan code locally.
Below both approaches are described. Compiling the Fortran source code can
be a bit more challenging depending on your environment. Please carefully read
the instructions below and read through the dedicated `GitHub Discussions category
<https://github.com/phydrus/phydrus/discussions/categories/compiling-fortran
-code>`_ on this topic before opening a new Discussion.

Compiling the source code
-------------------------
The recommended option is to compile the adapted Hydrus-1D Fortran77 files for
your own environment and computer. The following steps should be taken:

1. Download the Phydrus-optimized fortran code from `this
dedicated repository <https://github.com/phydrus/source_code>`_.

2. Open a Windown command line or Linux / MacOS Terminal and move into the
`source` folder (e.g., `cd path/to/directory`).

3. and use the following syntax in your terminal or windows command line to
compile the source code:

>>> make

4. This should create a Windows or Unix Executable that can be used to run the
HYDRUS-1D simulation. In the Python code, you have to reference to the
location of the executable, so you can store it anywhere you want.

Troubleshooting
~~~~~~~~~~~~~~~
Depending on your operation system (e.g., Windows/MacOS/Linux) you may need
to install additional tools to compile:

- The compile the source code for MacOS/Linux, Gfortran needs to be installed. Instructions can be found here: https://gcc.gnu.org/wiki/GFortranBinaries.
- Please let us know when you find other requirements to add to this list!

Using pre-compiled executables
------------------------------
It is in principle possible to use pre-compiled versions of Hydrus-1D. Note
that the Phydrus software has been developed based on the the Hydrus-1D 4.08
Fortran Code and other versions may not be suppported.

If you have the Graphical User Interface for Hydrus-1D for Windows installed
(which can be obtained from https://www.pc-progress.com/ ), you may directly
point to that exectuable. You can also download the pre-compiled versions
from the `source code repository <https://github.com/phydrus/source_code>`_.
However, it can not be guaranteed that these executables work for your
personal computing environment.
