import numpy as np
from pySOT.optimization_problems import OptimizationProblem
import capnp
import eval_capnp

class StelloptProblem(OptimizationProblem):
    """STELLOPT weighted objective via RPC

    Connect to a Cap'n Proto RPC server and forward `eval` calls to an RPC
    `eval` method.  Intended for use with STELLOPT in "rpc-eval" mode.  Ideally
    each worker thread should connect to a different server, which will be the
    leader of its own MPI group.
    """
    def __init__(self, dim, server_addr):
        # TODO: Read dim, bounds, etc. from server
        self.dim = dim
        # For now, default bounds are [0,1] and may be safely overridden by the
        # user.
        self.lb = np.zeros(self.dim)
        self.ub = np.ones(self.dim)

        # All parameters are continuous.
        self.int_var = np.array([])
        self.cont_var = np.arange(0, self.dim)

        # Connect to RPC server and get handle to `eval` method.
        client = capnp.TwoPartyClient(server_addr)
        self.evaluator = client.bootstrap().cast_as(eval_capnp.Eval)

    def eval(self, x):
        """Evaluate objective function via RPC

        Call RPC `eval` method on configured server and wait for and return
        result.
        """
        self.__check_input__(x)
        # TODO: Handle errors and disconnects
        return self.evaluator.eval(x.tolist()).wait().objective

    def evalAll(self, x):
        self.__check_input__(x)
        return np.array(self.evaluator.evalAll(x.tolist()).wait().objectives)
