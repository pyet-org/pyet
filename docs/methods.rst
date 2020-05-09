FAO-56 Penman-Monteith (Allen, 1998)
-----------------

.. math::

   \begin{eqnarray}
      ET_0 = \frac{0.408 \Delta (R_{n}-G)+ \gamma \frac{900}{T+273} u_2 (e_{s}-e_{a})}{\Delta +\gamma (1+0.34u_2)}
   \end{eqnarray}
   
Penman-Monteith (ASCE, 2005)
-----------------

.. math::

   \begin{eqnarray}
      \lambda ET = \frac{\Delta (R_{n}-G)+ \rho_a c_p \frac{(e_{s}-e_{a})}{r_a}}{\Delta +\gamma (1+\frac{r_s}{r_a})}
   \end{eqnarray}

Penman-Monteith (FAO, 1990)
-----------------

Penman (1948)
-----------------

.. math::

   \begin{eqnarray}
      \lambda \ET = \frac{\Delta (R_{n}-G)+ K_w \frac{\gamma}{\Delta + \gamma}(a_w+b_w u_2)(e_{s}-e_{a})}{\Delta +\gamma)}
   \end{eqnarray}
where:
	k:sub:`w`				= unit constant
	a:sub:`w` and b:sub:`w` = wind function coefficients
	u:sub:`2`				= wind speed at 2m [ms :sup:`-1`]
	
Hargreaves (Hargreaves and Samani, 1982; 1985)
-----------------

.. math::

   \begin{eqnarray}
     \lambda ET_0 = 0.0023(T_{max}-T_{min})^{0.5} (T_{mean}+17.8)R_a
   \end{eqnarray}
   
Hamon (1961)
-----------------

Priestley and Taylor (1972)
-----------------

.. math::

   \begin{eqnarray}
     ET = 1.26 \frac{\Delta}{\Delta + gamma} (R_n - G)
   \end{eqnarray}
   
Jensen and Haise (1972)
-----------------

.. math::

   \begin{eqnarray}
     \lambda ET = C_r (T - T_x) R_s
   \end{eqnarray}

Makkink (1957)
-----------------

.. math::

   \begin{eqnarray}
     ET = 0.61 \frac{\Delta}{\Delta + gamma} \frac{R_s}{2.45} -0.12)
   \end{eqnarray}
   
Corrected Penman-Monteith by Schymanski (2007)
-----------------

.. math::

   \begin{eqnarray}
      \lambda ET = \frac{\Delta (R_{n}-G)+ \rho_a c_p \frac{(e_{s}-e_{a}) a_{sh}}{r_a}}{\Delta +\gamma \frac{a_{sh}}{a_s}(1+\frac{r_s}{r_a})}
   \end{eqnarray}
