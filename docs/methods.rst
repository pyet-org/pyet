Penman (1948)
-------------

.. math::

      PET = \frac{\Delta (R_n - G) + \gamma 2.6 (1 + 0.536 u_2)(e_s-e_a)}{\lambda (\Delta +\gamma)}

FAO-56 Penman-Monteith (Allen, 1998)
------------------------------------

.. math::

    PET = \frac{0.408 \Delta (R_n - G)+\gamma \frac{900}{T+273}(e_s-e_a) u_2}{\Delta+\gamma(1+0.34 u_2)}

Penman-Monteith (Monteith, 1965)
--------------------------------

.. math::

    PET = \frac{\Delta (R_n-G)+ \rho_a c_p 86400 \frac{e_s-e_a}{r_a}}{\lambda(\Delta +\gamma(1+\frac{r_s}{r_a}))}

Kimberly-Penman (Wright, 1982)
------------------------------

.. math::

    PET = \frac{\Delta (R_n-G)+ \gamma (e_s-e_a) w}{\lambda(\Delta+\gamma)}

.. math::

    w =  u_2 * (0.4 + 0.14 * exp(-(\frac{J_D-173}{58})^2))+(0.605 + 0.345 * exp(-(\frac{J_D-243}{80})^2))

Doorenbos–Pruitt (FAO-24) (Jensen et al., 1990)
-----------------------------------------------

.. math::

    PET = \frac{-0.3 \Delta + R_s (1-\alpha) w}{\lambda(\Delta +\gamma)}

.. math::

    w = 1.066-0.13*\frac{rh}{100}+0.045*u_2-0.02*\frac{rh}{100}*u_2-3.15*(\frac{rh}{100})^2-0.0011*u_2

    
Thom and Oliver (1977)
----------------------

.. math::

    PET = \frac{\Delta (R_n-G)+ 2.5 \gamma (e_s-e_a) w}{\lambda(\Delta+\gamma(1+\frac{r_s}{r_a}))}

.. math::

    w=2.6(1+0.53u_2)

Priestley and Taylor (1972)
---------------------------

.. math::

    PET = \frac{\alpha_{PT} \Delta (R_n-G)}{\lambda(\Delta + \gamma)}

Makkink (1957)
--------------

.. math::

    PET = \frac{0.65 \Delta (R_s)}{\lambda(\Delta+\gamma)}

Hamon (1961)
------------

.. math::

    PET = (\frac{DL}{12})^2 exp(\frac{T_a}{16})

Hargreaves (Hargreaves and Samani, 1982; 1985)
----------------------------------------------

.. math::

    PET = 0.0023 \frac{R_a (T_a+17.8)\sqrt{(T_{max}-T_{min})}}{\lambda}$

Jensen and Haise (2005)
-----------------------

.. math::

    PET = \frac{C_r (T - T_x) R_s}{\lambda}

.. math::

    PET = \frac{R_a (T_a+5)}{\lambda 68}

Oudin (1961)
-------------

.. math::

    PET = \frac{R_a (T_a+K_2)}{\lambda K_1}; if T_a+K_2>0

.. math::

    else: P = 0

Abtew (1996)
------------

.. math::

    PET = \frac{k R_s}{\lambda}

Turc (1961)
-----------

.. math::

    PET=k(\frac{T_a}{T_a+15})(R_s/4.184 + 50)*4.184; for RH>50

.. math::

    PET=k(\frac{T_a}{T_a+15})(R_s/4.184 + 50)(1+\frac{50-RH}{70})*4.184;for RH<50
   
McGuinness and Bordne (1972)
----------------------------

.. math::

    PET = \frac{0.0147 R_a (T_a + 5)}{\lambda}

Linacre (1977)
--------------

.. math::

    PET = \frac{\frac{500 T_m}{(100-A)}+15 (T_a-T_d)}{80-T_a}

.. math::

    T_m = T_a + 0.006 * elevation

Blaney-Criddle (1950)
---------------------

.. math::

    PET=kp(0.46 * T_a + 8.13)

Romanenko (1961)
----------------

.. math::

    PET=4.5(1 + (\frac{T_a}{25})^2 (1  \frac{e_a}{e_s})

Notation and units
------------------

.. math::

    Variable ---- Description --- Units

    PET --- Potential evaporation --- mm d^{-1}

    \Delta --- Slope of vapor pressure curve --- kPa °C^{-1}

    \gamma --- Latent heat of vaporization --- MJ kg^{-1}

    \rho_w --- Water density (=1000) --- kg L^{-1}

    \rho_a --- Air density --- kg m^{-3}