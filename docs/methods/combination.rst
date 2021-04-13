Combination
===============================

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


Priestley and Taylor (1972)
---------------------------

.. math::

    PET = \frac{\alpha_{PT} \Delta (R_n-G)}{\lambda(\Delta + \gamma)}

Thom and Oliver (1977)
----------------------

.. math::

    PET = \frac{\Delta (R_n-G)+ 2.5 \gamma (e_s-e_a) w}{\lambda(\Delta+\gamma(1+\frac{r_s}{r_a}))}

.. math::

    w=2.6(1+0.53u_2)


    else: P = 0