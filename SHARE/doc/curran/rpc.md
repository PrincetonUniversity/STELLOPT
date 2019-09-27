RPC server for external optimizers
==================================

Dependencies
------------

Server side: The RPC server depends on Cap'n Proto, as well as a C++ compiler and runtime library.  The Docker container and make rules provide these.

Client side: The pySOT example depends on the pySOT and Cap'n Proto Python libraries.  Cap'n Proto additionally requires some native packages; on Ubuntu, these can be installed with `sudo apt install python-dev libcapnp-dev`.  To install the Python packages, activate your Python environment (e.g. `conda activate`) and run `pip install pySOT && pip install pycapnp`.

Getting started
---------------

Server (Docker container):

1. You may need to re-build and re-run your container to pick up the new dependencies and port forwarding (see `docker-build.sh` and `docker-run.sh`).
2. In the scenario directory, edit `input.test` and change the value of `opt_type` to `rpc_server`.
3. Run `xstelloptv2` on your scenario's `input.test`.

Client (local host):

1. Activate Python environment and install client dependencies (see above).
2. Based on your scenario, edit `SHARE/pySOT/example.py` to adjust the dimension (first argument to `StelloptProblem` constructor) and parameter bounds (`prob.ub`) accordingly.  Defaults are appropriate for `2DOF_circularCrossSection_varyAxis_targetIotaAndQuasisymmetry`.
3. `Run SHARE/pySOT/example.py`.

Discussion: RPC backend
-----------------------

Several cross-language serializaton and RPC libraries were considered in order to add RPC capabilities to STELLOPT.  For the time being we have settled on Cap'n Proto.  Here are some brief comments on other options considered:

* gRPC: The go-to standard, but Ubuntu's packages appear to be broken (server example segfaults).
* Cap'n Proto: Provides fast serialization and RPC for both C++ and Python, and is available as an Ubuntu package.  Python supports runtime ICD parsing.
* Apache Thrift: gRPC is often preferred over Thrift due to support for streaming RPC (though we don't need this).  Possibly worth reconsidering.
* Apache Arrow: Targets serialization of arrays; worth considering if we start moving more data (but overkill for now).
* MessagePack: RPC implementations don't appear to be mature.
* JSON: Poor support for numeric types; serialization overhead.

TODO
----

* Allow server port to be configured
* Serve parameters from xstelloptv2 to client so that problem dimension and parameter bounds only need to be configured in the scenario.
* Add MPI support to server (idea is for each node to be an independent MPI group with a single server process, and for each pySOT worker thread to connect to a different node's server).
* Provide method that returns individual objective values, to support multi-objective optimizers.
* Commandable server shutdown
