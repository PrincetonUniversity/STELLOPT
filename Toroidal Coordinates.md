{% include head.html %}

Toroidal Coordinate Systems
===========================

------------------------------------------------------------------------

Introduction
------------

In order to represent the fields in a compact form toroidal coordinates
are often used. In these coordinates the various coordinates and
quantities are represented by periodic functions in the poloidal (theta)
and toroidal (phi) directions. This allows the representation of
quantities on a given radial surface in terms of trigonometric function
(sine and cosine). The choice of radial coordinate can vary, but is
almost always normalized to some edge value (say toroidal flux in VMEC).
This allows quantities in 3D to be represented by continuous functions
in the poloidal and toroidal direction and a discrete set of points in
the radial direction. Thus on any radial surface the value of a
parameter is known to machine accuracy at any point. Of course the
accuracy of this representation is a function of the ability to
represent a spatially varying quantity with a truncated trigonometric
series. In general any quantity can be represented by

\$$ f\left(r,\theta,\zeta\right)=\sum_{m=0}^M\sum_{n=-N}^N f_{mn}\left(r\right) cos\left(m\theta+n\zeta\right) $$

or [math](math) f\\left(r,\\theta,\\zeta\\right)=\\sum
f\_{mn}\\left(r\\right) sin\\left(m\\theta+n\\zeta\\right). [math](math)
Here the angle per field period (zeta) has been utilized in place of the
total toroidal angle (phi). The two are related by the periodicity of
the toroidal domain [math](math) \\zeta=N\_{FP}\\phi. [math](math)
Finally, it is often more convenient to write the kernel of the
trigonometric functions in terms of normalized values instead of
radians. [math](math) u=\\frac{\\theta}{2\\pi} [math](math) [math](math)
v=\\frac{\\zeta}{2\\pi} [math](math)

The choice of trigonometric function for a given quantity is determined
by the symmetry of the problem. For systems with stellarator symmetry
(up-down in the phi=0 plane), the cylindrical radial coordinate (R) has
a even symmetry (cosine) while the vertical coordinate (Z) has a odd
symmetry (sine). In general, toroidal coordinates do not require this
symmetry and quantities are functions of a series of both odd and even
coefficients [math](math) f\\left(r,\\theta,\\zeta\\right)=\\sum
f\^C\_{mn}\\left(r\\right) cos\\left(m\\theta+n\\zeta\\right)+\\sum
f\^S\_{mn}\\left(r\\right) sin\\left(m\\theta+n\\zeta\\right).
[math](math) Here the superscripts denote that the coefficients are
different values.

There are three choices of kernel. Although they are simply a choice of
sign, they have been named after the codes which employ them. Here they
are \|\| Name \|\| Convention \|\| \|\| VMEC \|\| (mu-nv) \|\| \|\|
NESCOIL \|\| (mu+nv) \|\| \|\| PIES \|\| (nv-mu) \|\| Conversion between
these conventions is straightforward. Going from VMEC to NESCOIL simply
requires that an array be flipped about the toroidal mode index (n=0),
thus n=-n and -n=n (NESCOIL TO VMEC is the same conversion). The VMEC
convention is just the negative kernel of the PIES convention, so only
the odd (sine) coefficients need be multiplied by -1, remember:
[math](math) cos\\left(-x)\\right)=cos\\left(x\\right) [math](math)
[math](math) sin\\left(-x)\\right)=-sin\\left(x\\right) [math](math)
Thus conversion from one convention is simply a matter of flipping
arrays about the toroidal mode index (n) and negating odd coefficients
(sin).

------------------------------------------------------------------------

Curvilinear coordinates
-----------------------

In the toroidal domain a cylindrical coordinates system is used to
express the location of points in space. This position vector in this
coordinate system may be written [math](math)
\\vec{x}\\left(s,u,v\\right)=R\\left(s,u,v\\right)cos\\left(\\phi\\right)\\hat{x}+R\\left(s,u,v\\right)sin\\left(\\phi\\right)\\hat{y}+Z\\left(s,u,v\\right)s\\hat{z}
[math](math) where s is the normalized minor radial coordinate, u and v
are the normalized angular coordinates, and R and Z are functions of the
toroidal coordinates. The position function R and Z are written
[math](math) R\\left(s,u,v\\right)=\\sum R\^C\_{mn}\\left(s\\right)
cos\\left(mu+nv\\right)+\\sum R\^S\_{mn}\\left(s\\right)
sin\\left(mu+nv\\right) [math](math) \<span style=\"margin: 0px;
padding: 0px;\"\>and\</span\> [math](math) Z\\left(s,u,v\\right)=\\sum
Z\^S\_{mn}\\left(s\\right) sin\\left(mu+nv\\right)+\\sum
Z\^C\_{mn}\\left(s\\right) cos\\left(mu+nv\\right). [math](math) The
second term in each equation can be dropped if stellarator symmetry is
assumed. Vectors may be represented in terms of their contravariant
(sup, up, top) components or covariant (sub, dn, bottom) components
through the relation [math](math)
\\vec{A}=A\^s\\hat{e}\_s+A\^u\\hat{e}\_u+A\^v\\hat{e}\_v=A\_s\\hat{e}\^s+A\_u\\hat{e}\^u+A\_v\\hat{e}\^v.
[math](math) The covariant and contravariant basis vectors may be
written [math](math) \\hat{e}\_k=\\frac{\\partial \\vec{x}}{\\partial
x\_k} [math](math) and [math](math) \\hat{e}\^k=\\nabla x\_k.
[math](math)

### Contravariant Vector Components

This allows us to write the covariant basis vectors in terms of
cartesian unit vectors [math](math) \\hat{e}\_s=\\frac{\\partial
R}{\\partial s}cos\\left(\\phi\\right)\\hat{x}+\\frac{\\partial
R}{\\partial s}sin\\left(\\phi\\right)\\hat{y}+\\frac{\\partial
Z}{\\partial s}\\hat{z}, [math](math) [math](math)
\\hat{e}\_u=\\frac{\\partial R}{\\partial
u}cos\\left(\\phi\\right)\\hat{x}+\\frac{\\partial R}{\\partial
u}sin\\left(\\phi\\right)\\hat{y}+\\frac{\\partial Z}{\\partial
u}\\hat{z}, [math](math) [math](math)
\\hat{e}\_v=\\left(\\frac{\\partial R}{\\partial
v}cos\\left(\\phi\\right)-R\\frac{\\partial \\phi}{\\partial
v}sin\\left(\\phi\\right)\\right)\\hat{x}+\\left(\\frac{\\partial
R}{\\partial v}sin\\left(\\phi\\right)+R\\frac{\\partial
\\phi}{\\partial
v}cos\\left(\\phi\\right)\\right)\\hat{y}+\\frac{\\partial Z}{\\partial
v}\\hat{z}. [math](math) Here the derivative of phi with respect to the
normalized toroidal angle is kept general. This allows the cartesian
components of a vector to be written in terms of the contravariant
components: [math](math) A\_x=\\left(A\^s\\frac{\\partial R}{\\partial
s}+A\^u\\frac{\\partial R}{\\partial u}+A\^v\\frac{\\partial
R}{\\partial v}\\right)cos\\left(\\phi\\right)-A\^vR\\frac{\\partial
\\phi}{\\partial v}sin\\left(\\phi\\right), [math](math) [math](math)
A\_y=\\left(A\^s\\frac{\\partial R}{\\partial s}+A\^u\\frac{\\partial
R}{\\partial u}+A\^v\\frac{\\partial R}{\\partial
v}\\right)sin\\left(\\phi\\right)+A\^vR\\frac{\\partial \\phi}{\\partial
v}cos\\left(\\phi\\right), [math](math) and [math](math)
A\_z=A\^s\\frac{\\partial Z}{\\partial s}+A\^u\\frac{\\partial
Z}{\\partial u}+A\^v\\frac{\\partial Z}{\\partial v}. [math](math) The
components in cylindrical coordinates may also be written in terms of
the contravariant components: [math](math) A\_\\rho=A\^s\\frac{\\partial
R}{\\partial s}+A\^u\\frac{\\partial R}{\\partial
u}+A\^v\\frac{\\partial R}{\\partial v}, [math](math) [math](math)
A\_\\phi=A\^vR\\frac{\\partial \\phi}{\\partial v}, [math](math) and
[math](math) A\_z=A\^s\\frac{\\partial Z}{\\partial
s}+A\^u\\frac{\\partial Z}{\\partial u}+A\^v\\frac{\\partial
Z}{\\partial v}. [math](math) It is important to note that when working
with different coordinate systems a chain-rule can be used to convert
derivatives [math](math) \\frac{\\partial Z}{\\partial
v}=\\frac{\\partial Z}{\\partial v\_k}\\frac{\\partial v\_k}{\\partial
v\_l} [math](math) The surface normal vector (N, which does not have
unit length) can be written as the cross product of the covariant basis
vectors [math](math) \\vec{N}=\\frac{\\partial \\vec{x}}{\\partial
u}\\times\\frac{\\partial \\vec{x}}{\\partial
v}=\\hat{e}\_u\\times\\hat{e}\_v [math](math) allowing the cartesian
surface normal components to be written [math](math)
N\_x=-\\left(\\frac{\\partial R}{\\partial u}\\frac{\\partial
Z}{\\partial v}-\\frac{\\partial R}{\\partial v}\\frac{\\partial
Z}{\\partial u}\\right)sin\\left(\\phi\\right)+R\\frac{\\partial
\\phi}{\\partial v}\\frac{\\partial Z}{\\partial
u}cos\\left(\\phi\\right) [math](math) [math](math)
N\_y=\\left(\\frac{\\partial R}{\\partial u}\\frac{\\partial
Z}{\\partial v}-\\frac{\\partial R}{\\partial v}\\frac{\\partial
Z}{\\partial u}\\right)cos\\left(\\phi\\right)+R\\frac{\\partial
\\phi}{\\partial v}\\frac{\\partial Z}{\\partial
u}sin\\left(\\phi\\right) [math](math) [math](math)
N\_z=-R\\frac{\\partial \\phi}{\\partial v}\\frac{\\partial R}{\\partial
u}. [math](math) This vector integrated can be treated as the product of
the surface normal vector (unit length) and the differential surface
element [math](math) \\vec{N}=\\hat{n}\\cdot dA. [math](math)
