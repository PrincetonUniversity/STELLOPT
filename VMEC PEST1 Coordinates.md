VMEC PEST1 Coordinates
======================

------------------------------------------------------------------------

The coordinates output by the VMEC code are not straight field line
coordinates. In order to convert to straight field line coordinates
(PEST1) the VMEC theta coordinate must be transformed to the straight
field line coordinate. In general the VMEC coordinates are written
[math](math)
R\_{m,n}\\left(s\\right)=\\frac{2\*N}{\\left(2\\pi\\right)\^2}\\int\^{2\\pi/L}\_0
dv \\int\^{2\\pi}\_0 d\\theta\_\* R\\left(s,\\theta\_\*,v\\right)
\\cos\\left(m\\theta\_\*-nv\\right) [math](math) where [math](math)
\\theta\_\*= u + \\lambda\\left(s,u,v\\right) [math](math) In VMEC, the
lambda quantity is in terms of normalized flux. Using these definitions,
the following identity can be written [math](math)
R\_{m,n}\\left(s\\right)=\\frac{2\*N}{\\left(2\\pi\\right)\^2}\\int\^{2\\pi/L}\_0
dv \\int\^{2\\pi}\_0 du \\left(1+\\frac{\\partial \\lambda}{\\partial
u}\\right) R\\left(s,u,v\\right)
\\cos\\left(m\\left(u+\\lambda\\right)-nv\\right) -
g\_{00}\\frac{2\*N}{\\left(2\\pi\\right)\^2}\\int\^{2\\pi/L}\_0 dv
\\int\^{2\\pi}\_0 du \\left(1+\\frac{\\partial \\lambda}{\\partial
u}\\right) R\\left(s,u,v\\right) [math](math) As a result the quantities
in VMEC coordinates must first be transformed to suv space, then Fourier
decomposed to the PEST1 coordinates.
