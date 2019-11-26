#include "eval.capnp.h"
#include <capnp/ez-rpc.h>
#include <cassert>
#include <cmath>

/** Server implementation for Eval RPC methods */
class EvalImpl final: public Eval::Server {
public:
  /**
   * Construct a new server implementation for Eval methods.
   * 
   * @param fcn Callback to evaluate objective function
   * @param m Number of objectives
   * @param n Number of parameters
   * @param x Vector of length n to use when passing parameters to fcn
   * @param fvec Vector of length m to use for extracting objective values from
   *          fcn
   */
  EvalImpl(
    void(*fcn)(int const * m,
               int const * n,
               double x[],
               double fvec[],
               int * iflag,
               int const * ncnt),
    int const m,
    int const n,
    double x[],
    double fvec[]):
  _fcn(fcn),
  _m(m),
  _n(n),
  _x(x),
  _fvec(fvec),
  _ncnt() {}

  /**
   * Implementation of `eval` RPC method.  Returns Euclidian norm of objective
   * values at specified parameters.
   */
  kj::Promise<void> eval(EvalContext context) override {
    auto params = context.getParams().getParams();
    assert(params.size() == _n);

    ++_ncnt;

    // Evaluate `fcn` at provided parameters
    for (size_t i = 0; i < params.size(); ++i) {
      _x[i] = params[i];
    }
    int iflag = 0;
    _fcn(&_m, &_n, _x, _fvec, &iflag, &_ncnt);

    // Compute norm of objectives (weighting is performed by `fcn`)
    double obj = 0.0;
    for (size_t i = 0; i < _m; ++i) {
      obj += _fvec[i]*_fvec[i];
    }
    context.getResults().setObjective(std::sqrt(obj));

    return kj::READY_NOW;
  }

  kj::Promise<void> evalAll(EvalAllContext context) override {
    auto params = context.getParams().getParams();
    assert(params.size() == _n);

    ++_ncnt;

    // Evaluate `fcn` at provided parameters
    for (size_t i = 0; i < params.size(); ++i) {
      _x[i] = params[i];
    }
    int iflag = 0;
    _fcn(&_m, &_n, _x, _fvec, &iflag, &_ncnt);

    context.getResults().setObjectives(kj::ArrayPtr<double>(_fvec, _m));

    return kj::READY_NOW;
  }

private:
  /** Callback to evaluate objective function */
  void(*_fcn)(int const * m,
              int const * n,
              double x[],
              double fvec[],
              int * iflag,
              int const * ncnt);

  /** Number of objectives */
  int const _m;

  /** Number of parameters */
  int const _n;

  /** Vector of length n to use when passing parameters to fcn */
  double * _x;

  /** Vector of length m to use for extracting objective values from fcn */
  double * _fvec;

  /** Number of times objective function has been called */
  int _ncnt;
};


extern "C" {

/**
 * Start a Cap'n Proto RPC server for Eval methods.
 * 
 * @param fcn Callback to evaluate objective function
 * @param m Number of objectives
 * @param n Number of parameters
 * @param x Vector of length n to use when passing parameters to fcn
 * @param fvec Vector of length m to use for extracting objective values from
 *          fcn
 */
void serve_rpc(
    void(*fcn)(int const * m,
               int const * n,
               double x[],
               double fvec[],
               int * iflag,
               int const * ncnt),
    int * m,
    int * n,
    double x[],
    double fvec[]) {
  // TODO: Configurable port
  capnp::EzRpcServer server(kj::heap<EvalImpl>(fcn, *m, *n, x, fvec), "0.0.0.0", 5923);
  auto & waitScope = server.getWaitScope();
  kj::NEVER_DONE.wait(waitScope);
}

}
