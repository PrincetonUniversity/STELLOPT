# Example of optimizing stellarator design using pySOT and STELLOPT's RPC
# server.

import numpy as np
from pySOT.surrogate import RBFInterpolant, CubicKernel, LinearTail, SurrogateUnitBox
from pySOT.experimental_design import SymmetricLatinHypercube
from poap.controller import SerialController
from pySOT.strategy import SRBFStrategy

from stellopt import StelloptProblem

# 2D STELLOPT problem (e.g.
# 2DOF_circularCrossSection_varyAxis_targetIotaAndQuasisymmetry)
prob = StelloptProblem(2, 'localhost:5923')
prob.ub = np.array([0.2, 0.2])

rbf = SurrogateUnitBox(RBFInterpolant(dim=prob.dim, kernel=CubicKernel(), 
                                      tail=LinearTail(prob.dim)),
                       lb=prob.lb, ub=prob.ub)
slhd = SymmetricLatinHypercube(dim=prob.dim, num_pts=2*(prob.dim+1))
controller = SerialController(prob.eval)
controller.strategy = SRBFStrategy(max_evals=50, opt_prob=prob,
                                   exp_design=slhd, surrogate=rbf)
result = controller.run()

print('params: %s' % result.params[0])
print('value: %s' % result.value)
