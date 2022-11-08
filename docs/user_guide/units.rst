Notation and units
------------------

Many issues and errors in the evapotranspiration estimation come from the wrong units of the input data. Throughout
PyET we have tried to be consistent in the notation of the variables and their units. Table 1 provides an
overview of the different variables, their units, and python function argument name. When providing arguments to any
the evapotranspiration methods, it is important to make sure the units of each variable are as listed in Table 1.

.. list-table::
   :widths: 25 50 25
   :header-rows: 1

   * - .. math:: Variable
     - .. math:: Description
     - .. math:: Units
   * - .. math:: PE
     - Potential evaporation
     - .. math:: mm d^{-1}
   * - .. math:: \Delta
     - Slope of vapor pressure curve
     - .. math:: kPa °C^{-1}
   * - .. math:: \gamma
     - Latent heat of vaporization
     - .. math:: MJ kg^{-1}
   * - .. math::  \rho_w
     - Water density (=1000)
     - .. math:: kg L^{-1}
   * - .. math:: \rho_a
     - Air density
     - .. math::  kg m^{-3}
   * - .. math:: \gamma
     - Psychrometric constant
     - .. math:: kPa °C^{-1}
   * - .. math:: e_s
     - Saturation vapour pressure
     - .. math::  kPa
   * - .. math::  e_a
     - Actual vapour pressure
     - .. math::  kPa
   * - .. math:: r_a
     - Aerodynamic resistance
     - .. math:: s m^{-1}
   * - .. math::  r_s
     - Surface resistance of reference crop (=70)
     - .. math:: s m^{-1}
   * - .. math:: u_2
     - Wind speed 2 m above soil surface
     - .. math:: m s^{-1}
   * - .. math:: T_a
     - Air temperature
     - .. math:: °C
   * - .. math:: T_d
     - Dew point temperature
     - .. math:: °C
   * - .. math::  T_{max}
     - Maximum air temperature
     - .. math:: °C
   * - .. math::  T_{min}
     - Minimum air temperature
     - .. math:: °C
   * - .. math::  R_a
     - Extraterrestrial radiation
     - .. math:: MJ m^{-2} d^{-1}
   * - .. math::  R_s
     - Global short-wave radiation
     - .. math:: MJ m^{-2} d^{-1}
   * - .. math::  R_n
     - Net incoming solar radiation
     - .. math:: MJ m^{-2} d^{-1}
   * - .. math::  RH
     - Relative humidity
     - .. math::  \%
   * - .. math::  DL
     - Day length
     - .. math:: h d^{-1}
   * - .. math::  \alpha
     - Surface albedo
     - .. math::  \%
   * - .. math::  J_D
     - Julian day
     - ordinal date
   * - .. math::  c_p
     - Air specific heat capacity
     - .. math:: MJ kg^{-1} °C^{-1}