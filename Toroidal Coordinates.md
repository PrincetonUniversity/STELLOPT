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
or
\$$ f\left(r,\theta,\zeta\right)=\sum f\_{mn}\left(r\right) sin\left(m\theta+n\zeta\right). $$

Here the angle per field period (zeta) has been utilized in place of the
total toroidal angle (phi). The two are related by the periodicity of
the toroidal domain $$ \zeta=N\_{FP}\phi. $$
Finally, it is often more convenient to write the kernel of the
trigonometric functions in terms of normalized values instead of
radians. $$ u=\frac{\theta}{2\pi} $$ $$
v=\frac{\zeta}{2\pi} $$

The choice of trigonometric function for a given quantity is determined
by the symmetry of the problem. For systems with stellarator symmetry
(up-down in the phi=0 plane), the cylindrical radial coordinate (R) has
a even symmetry (cosine) while the vertical coordinate (Z) has a odd
symmetry (sine). In general, toroidal coordinates do not require this
symmetry and quantities are functions of a series of both odd and even
coefficients $$ f\left(r,\theta,\zeta\right)=\sum
f\^C\_{mn}\left(r\right) cos\left(m\theta+n\zeta\right)+\sum
f\^S\_{mn}\left(r\right) sin\left(m\theta+n\zeta\right).
$$ Here the superscripts denote that the coefficients are
different values.

There are three choices of kernel. Although they are simply a choice of
sign, they have been named after the codes which employ them. Here they
are 

| Name | Convention |
|------|----------|
| VMEC | (mu-nv) | 
|NESCOIL | (mu+nv) |
| PIES | (nv-mu) |

Conversion between these conventions is straightforward. Going from VMEC to NESCOIL simply
requires that an array be flipped about the toroidal mode index (n=0),
thus n=-n and -n=n (NESCOIL TO VMEC is the same conversion). The VMEC
convention is just the negative kernel of the PIES convention, so only
the odd (sine) coefficients need be multiplied by -1, remember:
$$ cos\left(-x)\right)=cos\left(x\right) $$
$$ sin\left(-x)\right)=-sin\left(x\right) $$
Thus conversion from one convention is simply a matter of flipping
arrays about the toroidal mode index (n) and negating odd coefficients
(sin).

------------------------------------------------------------------------

Curvilinear coordinates
-----------------------

In the toroidal domain a cylindrical coordinates system is used to
express the location of points in space. This position vector in this
coordinate system may be written 
\$$ \vec{x}\left(s,u,v\right)=R\left(s,u,v\right)cos\left(\phi\right)\hat{x}+R\left(s,u,v\right)sin\left(\phi\right)\hat{y}+Z\left(s,u,v\right)s\hat{z} $$ 
where s is the normalized minor radial coordinate, u and v
are the normalized angular coordinates, and R and Z are functions of the
toroidal coordinates. The position function R and Z are written

\$$ R\left(s,u,v\right)=\sum R\^C\_{mn}\left(s\right)
cos\left(mu+nv\right)+\sum R\^S\_{mn}\left(s\right)
sin\left(mu+nv\right) $$ 
\$$ Z\left(s,u,v\right)=\sum
Z\^S\_{mn}\left(s\right) sin\left(mu+nv\right)+\sum
Z\^C\_{mn}\left(s\right) cos\left(mu+nv\right). $$ 

Thesecond term in each equation can be dropped if stellarator symmetry is
assumed. Vectors may be represented in terms of their contravariant
(sup, up, top) components or covariant (sub, dn, bottom) components
through the relation 
$$ \vec{A}=A\^s\hat{e}\_s+A\^u\hat{e}\_u+A\^v\hat{e}\_v=A\_s\hat{e}\^s+A\_u\hat{e}\^u+A\_v\hat{e}\^v.$$ 
The covariant and contravariant basis vectors may be
written $$ \hat{e}\_k=\frac{\partial \vec{x}}{\partial
x\_k} $$ and $$ \hat{e}\^k=\nabla x\_k.$$

### Contravariant Vector Components

This allows us to write the covariant basis vectors in terms of
cartesian unit vectors \$$ \hat{e}\_s=\frac{\partial
R}{\partial s}cos\left(\phi\right)\hat{x}+\frac{\partial
R}{\partial s}sin\left(\phi\right)\hat{y}+\frac{\partial
Z}{\partial s}\hat{z}, $$ 
\$$
\hat{e}\_u=\frac{\partial R}{\partial
u}cos\left(\phi\right)\hat{x}+\frac{\partial R}{\partial
u}sin\left(\phi\right)\hat{y}+\frac{\partial Z}{\partial
u}\hat{z}, $$ 
\$$
\hat{e}\_v=\left(\frac{\partial R}{\partial
v}cos\left(\phi\right)-R\frac{\partial \phi}{\partial
v}sin\left(\phi\right)\right)\hat{x}+\left(\frac{\partial
R}{\partial v}sin\left(\phi\right)+R\frac{\partial
\phi}{\partial
v}cos\left(\phi\right)\right)\hat{y}+\frac{\partial Z}{\partial
v}\hat{z}. $$ Here the derivative of phi with respect to the
normalized toroidal angle is kept general. This allows the cartesian
components of a vector to be written in terms of the contravariant
components: $$ A\_x=\left(A\^s\frac{\partial R}{\partial
s}+A\^u\frac{\partial R}{\partial u}+A\^v\frac{\partial
R}{\partial v}\right)cos\left(\phi\right)-A\^vR\frac{\partial
\phi}{\partial v}sin\left(\phi\right), $$ 
\$$
A\_y=\left(A\^s\frac{\partial R}{\partial s}+A\^u\frac{\partial
R}{\partial u}+A\^v\frac{\partial R}{\partial
v}\right)sin\left(\phi\right)+A\^vR\frac{\partial \phi}{\partial
v}cos\left(\phi\right), $$ and $$
A\_z=A\^s\frac{\partial Z}{\partial s}+A\^u\frac{\partial
Z}{\partial u}+A\^v\frac{\partial Z}{\partial v}. $$ The
components in cylindrical coordinates may also be written in terms of
the contravariant components: $$ A\_\rho=A\^s\frac{\partial
R}{\partial s}+A\^u\frac{\partial R}{\partial
u}+A\^v\frac{\partial R}{\partial v}, $$ 
\$$
A\_\phi=A\^vR\frac{\partial \phi}{\partial v}, $$ and
$$ A\_z=A\^s\frac{\partial Z}{\partial
s}+A\^u\frac{\partial Z}{\partial u}+A\^v\frac{\partial
Z}{\partial v}. $$ It is important to note that when working
with different coordinate systems a chain-rule can be used to convert
derivatives $$ \frac{\partial Z}{\partial
v}=\frac{\partial Z}{\partial v\_k}\frac{\partial v\_k}{\partial
v\_l} $$ The surface normal vector (N, which does not have
unit length) can be written as the cross product of the covariant basis
vectors $$ \vec{N}=\frac{\partial \vec{x}}{\partial
u}\times\frac{\partial \vec{x}}{\partial
v}=\hat{e}\_u\times\hat{e}\_v $$ allowing the cartesian
surface normal components to be written $$
N\_x=-\left(\frac{\partial R}{\partial u}\frac{\partial
Z}{\partial v}-\frac{\partial R}{\partial v}\frac{\partial
Z}{\partial u}\right)sin\left(\phi\right)+R\frac{\partial
\phi}{\partial v}\frac{\partial Z}{\partial
u}cos\left(\phi\right) $$ 
\$$
N\_y=\left(\frac{\partial R}{\partial u}\frac{\partial
Z}{\partial v}-\frac{\partial R}{\partial v}\frac{\partial
Z}{\partial u}\right)cos\left(\phi\right)+R\frac{\partial
\phi}{\partial v}\frac{\partial Z}{\partial
u}sin\left(\phi\right) $$ 
\$$
N\_z=-R\frac{\partial \phi}{\partial v}\frac{\partial R}{\partial
u}. $$ This vector integrated can be treated as the product of
the surface normal vector (unit length) and the differential surface
element $$ \vec{N}=\hat{n}\cdot dA. $$
