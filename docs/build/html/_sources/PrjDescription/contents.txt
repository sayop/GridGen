Project description
===================


In this exercise you will generate an inviscid, 2-D computational grid around a modified NACA 00xx series airfoil in a channel. The thickness distribution of a modified NACA 00xx series airfoil is given by:

.. math::
   y(x) = \pm 5t [0.2969 \sqrt{x_{int}x} - 0.126 x_{int} x - 0.3516 (x_{int}x)^{2} + 0.2843(x_{int}x)^{3} - 0.1015 (x_{int}x)^{4}]

where the ":math:`+`" sign is used for the upper half of the airfoil, the ":math:`-`" sign is used for the lower half and :math:`x_{int} = 1.008930411365`. Note that in the expression above :math:`x`, :math:`y`, and :math:`t` represent values which have been normalized by the airfoil by the airf
oil chord.

A sketch of the computational domain is shown below:

.. image:: ./images/CP1_compute_domain.png
   :scale: 60%

Each grid point can be described by :math:`(x,y)` location or :math:`(i,j)` location where :math:`i` is the index in the :math:`\xi` direction and the :math:`j` index is in the :math:`\eta` direction. The grid should have imax=41 points in the :math:`\xi` direction and jmax=19 points in the :math:`\eta` direction. The coordinates of points A-F shown in the figure are given in the following table:

