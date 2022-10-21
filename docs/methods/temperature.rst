 Temperature
===============================


Blaney-Criddle (1950)
---------------------

.. math::

    PE=kp(0.46 * T_a + 8.13)

Hamon (1961)
------------

Method = 0: After Oudin et al. 2005

.. math::

    PE = k (\frac{DL}{12})^2 exp(\frac{T_a}{16})

Method = 1: After Ansorge et al. 2019, eq. 7

.. math::

    pt = 4.95 exp(0.062 T_{mean}) / 100

.. math::

    PE = 13.97 (\frac{DL}{12})^2 pt

Method = 2: After Ansorge et al. 2019, eq. 12

.. math::

    PE = 218.527 (\frac{DL}{12})^2 / (T_{mean}+273.3) exp(17.26939 T_{mean} / (T_{mean+273.3}))

Romanenko (1961)
----------------

.. math::

    PE=4.5(1 + (\frac{T_a}{25})^2 (1  \frac{e_a}{e_s})

Linacre (1977)
--------------

.. math::

    PE = \frac{\frac{500 T_m}{(100-A)}+15 (T_a-T_d)}{80-T_a}

.. math::

    T_m = T_a + 0.006 * elevation




