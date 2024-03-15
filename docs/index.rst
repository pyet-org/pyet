*pyet* - Estimation of Potential Evapotranspiration
===================================================

*pyet* is an open source Python package for the estimation of reference and potential evapotranspiration (PET) from
time series data (`Pandas <https://pandas.pydata.org>`_) and gridded data (`Xarray <https://xarray.dev>`_). This
allows users to estimate potential evapotranspiration and calibrate the models with just a few lines of python code.

.. grid::

    .. grid-item-card:: Getting started
        :link: userguide/index
        :link-type: doc

        User guide on the basic concepts of Pastas.

    .. grid-item-card:: Examples
        :link: examples/index
        :link-type: doc

        Examples of *pyet* usage.

    .. grid-item-card:: Code Reference
        :link: api/index
        :link-type: doc

        *pyet* code reference.

.. grid::

    .. grid-item-card:: Publications
        :link: publications
        :link-type: doc

        Overview of publications that use *pyet*

    .. grid-item-card:: References
        :link: references
        :link-type: doc

        References used in the package.


Currently, 18 methods are implemented for calculating daily PET
-----------------------------

.. list-table::
   :widths: 20 20 5 5 5 5 5 5 15
   :header-rows: 1

   * - Method name
     - pyet function
     - T
     - RH
     - R
     - u2
     - Lat.
     - El.
     - Benchmarked?
   * - Penman
     - penman
     - ✓ `a`_
     - ✓ `b,c`_
     - ✓ `d`_
     - ✓
     - ✓ `d`_
     - ✓ `e`_
     - ✓
   * - Penman-Monteith
     - pm
     - ✓ `a`_
     - ✓ `b,c`_
     - ✓ `d`_
     - ✓
     - ✓ `d`_
     - ✓ `e`_
     - ✓
   * - ASCE-PM
     - pm_asce
     - ✓ `a`_
     - ✓ `b,c`_
     - ✓ `d`_
     - ✓
     - ✓ `d`_
     - ✓ `e`_
     - ✓
   * - FAO-56
     - pm_fao56
     - ✓ `a`_
     - ✓ `b,c`_
     - ✓ `d`_
     - ✓
     - ✓ `d`_
     - ✓ `e`_
     - ✓
   * - Priestley-Taylor
     - priestley_taylor
     - ✓
     - ✓ `h`_
     - ✓ `h`_
     - -
     - ✓ `h`_
     - ✓ `e`_
     - ✓
   * - Kimberly-Penman
     - kimberly_penman
     - ✓ `a`_
     - ✓ `b,c`_
     - ✓ `d`_
     - ✓
     - ✓ `d`_
     - ✓ `e`_
     - -
   * - Thom-Oliver
     - thom_oliver
     - ✓ `a`_
     - ✓ `b,c`_
     - ✓ `d`_
     - ✓
     - ✓ `d`_
     - ✓ `e`_
     - -
   * - Blaney-Criddle
     - blaney_criddle
     - ✓
     - - `i`_
     - - `i`_
     - - `i`_
     - ✓
     - -
     - ✓
   * - Hamon
     - hamon
     - ✓
     - -
     - -
     - -
     - ✓
     - -
     - ✓
   * - Romanenko
     - romanenko
     - ✓
     - ✓
     - -
     - -
     - -
     - -
     - ✓
   * - Linacre
     - linacre
     - ✓ `j`_
     - -
     - -
     - -
     - -
     - ✓
     - ✓
   * - Haude
     - haude
     - ✓
     - ✓ `k`_
     - -
     - -
     - -
     - -
     - ✓
   * - Turc
     - turc
     - ✓
     - ✓
     - ✓
     - -
     - -
     - -
     - ✓
   * - Jensen-Haise
     - jensen_haise
     - ✓
     - -
     - ✓ `l`_
     - -
     - ✓ `l`_
     - -
     - ✓
   * - McGuinness-Bordne
     - mcguinness_bordne
     - ✓
     - -
     - -
     - -
     - ✓
     - -
     - ✓
   * - Hargreaves
     - hargreaves
     - ✓ `m`_
     - -
     - -
     - -
     - ✓
     - -
     - ✓
   * - FAO-24 radiation
     - fao_24
     - ✓
     - ✓
     - ✓
     - ✓
     - -
     - ✓ `e`_
     - -
   * - Abtew
     - abtew
     - ✓
     - -
     - ✓
     - -
     - -
     - -
     - ✓
   * - Makkink
     - makkink
     - ✓
     - -
     - ✓
        - -
     - -
     - ✓ `e`_
     - ✓
   * - Oudin
     - oudin
     - ✓
     - -
     - -
     - -
     - ✓
     - -
     - -

.. rubric:: Footnotes

.. _`a`: T_max and T_min can also be provided.
.. _`b`: RH_max and RH_min can also be provided.
.. _`c`: If actual vapor pressure is provided, RH is not needed.
.. _`d`: Input for radiation can be (1) Net radiation, (2) solar radiation, or (3) sunshine hours. If (1), then latitude is not needed. If (1, 3) then latitude and elevation are needed.
.. _`e`: One must provide either the atmospheric pressure or elevation.
.. _`h`: If net radiation is provided, RH and Lat are not needed.
.. _`i`: If method==2, u2, RH_min, and sunshine hours are required.
.. _`j`: Additional input of Tmax and Tmin, or Tdew.
.. _`k`: Input can be RH or actual vapor pressure.
.. _`l`: If method==1, latitude is needed instead of Rs.
.. _`m`: Tmax and Tmin also needed.
 

Using *pyet*? Show your support by citing us!
-----------------------------

If you find *pyet* useful and use it in your research or project, we kindly ask you to cite 
the *pyet* preprint published in Hydrology and Earth System Sciences (HESS) as follows:

- Vremec, M., Collenteur, R. A., and Birk, S.: Technical note: Improved handling of potential 
   evapotranspiration in hydrological studies with PyEt, Hydrol. Earth Syst. Sci. Discuss. 
   [preprint], https://doi.org/10.5194/hess-2022-417, 2023.