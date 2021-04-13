VMEC Descent Algorithm
==========================

The page provide documentation of the algorithm used by the VMEC code 
to descend in the energy space. This description is based upon work
by Florian Hindenlang.

---

The published algorithm by ([Hirshman and Betancourt, 1991](https://doi.org/10.1016/0021-9991(91)90267-O)) is a
Richardson scheme of the form:

$$\frac{\partial^2 X}{\partial t^2} + \frac{1}{\tau}\frac{\partial X}{\partial t} = F\left(X\right)$$

where the substitution is made

$$V=\frac{\partial X}{\partial t},~\frac{\partial V}{\partial t} + \frac{1}{\tau}V = F\left(X\right)$$

From the appendix of that paper the best choice for $$\tau$$ is found to be:

$$\frac{1}{\tau} = -\frac{d}{dt}\left(ln\left(|F|^2\right)\right) \Rightarrow \frac{\Delta t}{\tau}=-ln\left(|F|^2_n/|F|^2_{n-1}\right)$$

In the paper they proposed an algorithm with $$P=V/\Delta t$$ giving,

$$P_n = \beta_nP_{n-1}+F_n$$

$$X_{n+1} = X_n + \left(\Delta t\right)^2P_n$$

where 

$$\beta_n=|F|^2_n/|F|^2_{n-1}$$

which resembles the momentum method ($$\beta\lt1$$ fixed), but with a 
dependence on $$F_n$$ and $$F_{n-1}$$ which are likely to vary a
great deal ($$\beta\ge 1$$, can occur).

---

The actual VMEC code itself uses the following formulation with $$P=V/\Delta t$$

$$P_n = \frac{1}{1+\bar\tau_n}\left(\left(1-\bar\tau_n\right)P_{n-1}+F_n\right)$$

$$X_{n+1} = X_n + \left(\Delta t\right)^2P_n$$

with $$\bar\tau_n$$ averaged over the last 10 iterations.

The timestep within a given iteration can be assumed to be constant. The velocity equation

$$\frac{\partial V}{\partial t} + \frac{1}{\tau}V = F\left(X\right)$$

with the specific discretization choice

$$\frac{\left(V_n-V_{n-1}\right)}{\Delta t}+\frac{1}{\tau_n}\frac{\left(V_n+V_{n-1}\right)}{2} = F_n$$

can be rewritten as

$$\left(1+\frac{\Delta t}{2\tau_n}\right)P_n - \left(1-\frac{\Delta t}{2\tau_n}\right)P_{n-1} = F_n $$

Thus the code defines $$\frac{\Delta t}{2\tau_n}=\bar\tau=<\hat\tau_n,...,\hat\tau_{n-9}>\Delta t_n/2$$ where

$$ \hat\tau_n = min\left(0.15, |ln\left(\frac{|F|^2_n}{|F|_{n-1}}\right)|\right)/\Delta t_n$$

this implies that $$\hat\tau_n\in\left[0.00,0.15\right]/\Delta t_n$$.  

We note that in this formulation $$F$$ is the vector of the preconditioned force.  If $$\Delta t$$ is reduced (because of the Jacobian becoming negative), the iteration is restarted form the initial state.