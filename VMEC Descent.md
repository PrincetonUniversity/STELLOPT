VMEC Descent Algorithm
==========================

The page provide documentation of the algorithm used by the VMEC code 
to descend in the energy space. This description is based upon work
by Florian Hindenlang.

---

The published algorithm by (Hirshman and Betancourt, 1991) is a
Richardson scheme of the form:

$$\frac{\partial^2 X}{\partial t^2} + \frac{1}{\tau}\frac{\partial X}{\partial t} = F\left(X\right)$$

where the substitution is made

$$V=\frac{\partial X}{\partial t},~\frac{\partial V}{\partial t} + \frac{1}{\tau}V = F\left(X\right)$$

From the appendix of that paper the best choice for $$\tau$$ is found to be:

$$\frac{1}{\tau} = -\frac{d}{dt}\left(ln\left(|F|^2\right)\right)$$

In the paper they proposd and algorithm with $$P=V/\Delta t$$ giving,

$$P_n = \beta_nP_{n-1}+F_n$$

$$X_{n+1} = X_n + \left(\Delta t\right)^2P_n$$

where $$\beta_n=\lvertF\rvert^2_n/\left|F\right|^2_{n-1}$$

which resembles the momentum method ($$\beta\lt1$$ fixed), but with a 
dependence on $$F_n$$ and $$F_{n-1}$$ which are likely to vary a
great deal ($$\beta\ge 1$$, can occur).

---

The actual VMEC code itself uses the following formulation with $$P=V/\Delta t$$

$$P_n = \frac{1}{1+\bar\tau_n}\left(\left(1-\bar\tau_n\right)P_{n-1}+F_n\right)$$

$$X_{n+1} = X_n + \left(\Delta t\right)^2P_n$$

with $$\bar\tau_n$$ averaged over the last 10 iterations.

The timestep within a given iteration can be assumed to be constant which implies a specific time discretization of velocity:

$$\frac{\partial V}{\partial t} + \frac{1}{\tau}V = F\left(X\right)$$