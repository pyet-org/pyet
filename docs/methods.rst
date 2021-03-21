Penman (1948)
-------------

.. math::

      ET = \frac{\Delta (R_n - G) + \gamma 2.6 (1 + 0.536 u_2)(e_s-e_a)}{\lambda (\Delta +\gamma)}

FAO-56 Penman-Monteith (Allen, 1998)
------------------------------------

.. math::

   \begin{eqnarray}
      ET_0 = \frac{0.408 \Delta (R_n - G)+\gamma \frac{900}{T+273}(e_s-e_a) u_2}{\Delta+\gamma(1+0.34 u_2)}
   \end{eqnarray}
   
Penman-Monteith (Monteith, 1965)
--------------------------------

.. math::

   \begin{eqnarray}
      ET = \frac{\Delta (R_n-G)+ \rho_a c_p K_{min} \frac{e_s-e_a}{r_a}}{\lambda(\Delta +\gamma(1+\frac{r_s}{r_a}))}
   \end{eqnarray}
   
Kimberly-Penman (Wright, 1982)
------------------------------

.. math::

   \begin{eqnarray}
       ET = \frac{\Delta (R_n-G)+ \gamma (e_s-e_a) w}{\lambda(\Delta+\gamma)};
   \end{eqnarray}
     where:

    \begin{eqnarray}
        w =  u_2 * (0.4 + 0.14 * exp(-(\frac{J_D-173}{58})^2))+(0.605 + 0.345 * exp(-(\frac{J_D-243}{80})^2))
   \end{eqnarray}
   
Doorenbos–Pruitt (FAO-24) (Jensen et al., 1990)
-----------------------------------------------

.. math::

   \begin{eqnarray}
       ET = \frac{-0.3 \Delta + R_s (1-\alpha) w}{\lambda(\Delta +\gamma)}
   \end{eqnarray}
     where:

    \begin{eqnarray}
        w = 1.066-0.13*\frac{rh}{100}+0.045*u_2-0.02*\frac{rh}{100}*u_2-3.15*(\frac{rh}{100})^2-0.0011*u_2
    \end{eqnarray}   
    
Thom and Oliver (1977)
----------------------

.. math::

   \begin{eqnarray}
       ET = \frac{\Delta (R_n-G)+ 2.5 \gamma (e_s-e_a) w}{\lambda(\Delta+\gamma(1+\frac{r_s}{r_a}))}
   \end{eqnarray}
    where:

    \begin{eqnarray}
        w=2.6(1+0.53u_2)
    \end{eqnarray} 

Priestley and Taylor (1972)
---------------------------

.. math::

   \begin{eqnarray}
     ET = \frac{\alpha_{PT} \Delta (R_n-G)}{\lambda(\Delta + \gamma)}
   \end{eqnarray}

Makkink (1957)
--------------

.. math::

   \begin{eqnarray}
     ET = \frac{0.65 \Delta (R_s)}{\lambda(\Delta+\gamma)}
   \end{eqnarray}

Hamon (1961)
------------

.. math::

   \begin{eqnarray}
     ET = (\frac{DL}{12})^2 exp(\frac{T_a}{16})
   \end{eqnarray}
   
Hargreaves (Hargreaves and Samani, 1982; 1985)
----------------------------------------------

.. math::

   \begin{eqnarray}
     ET = 0.0023 \frac{R_a (T_a+17.8)\sqrt{(T_{max}-T_{min})}}{\lambda}$
   \end{eqnarray}

Jensen and Haise (2005)
-----------------------

.. math::

   \begin{eqnarray}
     ET = \frac{C_r (T - T_x) R_s}{\lambda}
   \end{eqnarray}

    or

   \begin{eqnarray}
     ET = \frac{R_a (T_a+5)}{\lambda 68}
   \end{eqnarray}   

Oudin (1961)
-------------

.. math::

   \begin{eqnarray}
        ET = \frac{R_a (T_a+K_2)}{\lambda K_1};
   \end{eqnarray}
        if $T_a+K_2>0$

    else: P = 0

Abtew (1996)
------------

.. math::

   \begin{eqnarray}
       ET = \frac{k R_s}{\lambda}
   \end{eqnarray}

Turc (1961)
-----------

.. math::

   \begin{eqnarray}
       ET=k(\frac{T_a}{T_a+15})(R_s/4.184 + 50)*4.184; for RH>50
   \end{eqnarray}
   \begin{eqnarray}
       ET=k(\frac{T_a}{T_a+15})(R_s/4.184 + 50)(1+\frac{50-RH}{70})*4.184;for RH<50
   \end{eqnarray}
   
McGuinness and Bordne (1972)
----------------------------

.. math::

   \begin{eqnarray}
       ET = \frac{0.0147 R_a (T_a + 5)}{\lambda}
   \end{eqnarray}

Linacre (1977)
--------------

.. math::

   \begin{eqnarray}
       ET = \frac{\frac{500 T_m}{(100-A)}+15 (T_a-T_d)}{80-T_a}
   \end{eqnarray}

    \begin{eqnarray}
       T_m = T_a + 0.006 * elevation
   \end{eqnarray}

Blaney-Criddle (1950)
---------------------

.. math::

   \begin{eqnarray}
       ET=kp(0.46 * T_a + 8.13)
   \end{eqnarray}
   
Romanenko (1961)
----------------

.. math::

   \begin{eqnarray}
       ET=4.5(1 + (\frac{T_a}{25})^2 (1  \frac{e_a}{e_s})
   \end{eqnarray}