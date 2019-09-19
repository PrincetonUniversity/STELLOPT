#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include <stdio.h>

struct fcn_data {
    void(*fcn)(int const * m,
               int const * n,
               double x[],
               double fvec[],
               int * iflag,
               int const * ncnt);
    double * fvec;
    int ncnt;
};

static int fcn_f(gsl_vector const * x, void * data, gsl_vector * f) {
    struct fcn_data * d = (struct fcn_data *)data;
    int const m = f->size;
    int const n = x->size;

    assert(x->stride == 1);
    double * ff = f->data;
    if (f->stride != 1) {
      ff = d->fvec;
    }

    int iflag = 0;
    ++d->ncnt;
    printf("ncnt = %d\n", d->ncnt);

    d->fcn(&m, &n, x->data, ff, &iflag, &d->ncnt);

    if (f->stride != 1) {
      gsl_vector_view fvec_view = gsl_vector_view_array(ff, m);
      gsl_vector_memcpy(f, &fvec_view.vector);
    }
    return GSL_SUCCESS;
}

static void
callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    gsl_vector *x = gsl_multifit_nlinear_position(w);
    printf("iter = %ld\n", iter);
    printf(" x =\n");
    gsl_vector_fprintf(stdout, x, "%g");
    printf(" f =\n");
    gsl_vector_fprintf(stdout, f, "%g");
}

void clmdif(
    void(*fcn)(int const * m,
               int const * n,
               double x[],
               double fvec[],
               int * iflag,
               int const * ncnt),
    int * m,
    int * n,
    double x[], // DIMENSION(n)
    double fvec[], // DIMENSION(m)
    double const * ftol,
    double const * xtol,
    double const * gtol,
    int * maxfev,
    double * epsfcn,
    double diag[], // DIMENSION(n)
    int * mode,
    double const * factor,
    int * nprint,
    int * info,
    int * nfev,
    double * fjac,
    int * ldfjac,
    int ipvt[],
    double qtf[],
    double wa1[],
    double wa2[],
    double wa3[],
    double wa4[],
    double const xvmin[],
    double const xvmax[]) {

    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params =
        gsl_multifit_nlinear_default_parameters();
    fdf_params.h_df = sqrt(*epsfcn);
    fdf_params.h_df *= 0.1;  // FIXME
    //fdf_params.trs = gsl_multifit_nlinear_trs_dogleg;
    //fdf_params.scale = gsl_multifit_nlinear_scale_levenberg;
    //fdf_params.fdtype = GSL_MULTIFIT_NLINEAR_CTRDIFF;

    struct fcn_data d = { fcn, fvec, 0 };
    gsl_vector_view x_init = gsl_vector_view_array (x, *n);
    int status, driver_info;

    /* define the function to be minimized */
    fdf.f = fcn_f;
    fdf.df = NULL;   /* set to NULL for finite-difference Jacobian */
    fdf.fvv = NULL;     /* not using geodesic acceleration */
    fdf.n = *m;
    fdf.p = *n;
    fdf.params = &d;

    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, *m, *n);

    /* initialize solver with starting point */
    gsl_multifit_nlinear_init (&x_init.vector, &fdf, w);

    /* solve the system with a maximum of maxfev iterations */
    status = gsl_multifit_nlinear_driver(*maxfev, *xtol, *gtol, *ftol,
                                         &callback, NULL, &driver_info, w);
    printf("Status = %s; Convergence = %d\n", gsl_strerror(status), driver_info);

    memcpy(x, w->x->data, w->x->size*sizeof(double));
    memcpy(fvec, w->f->data, w->f->size*sizeof(double));

    gsl_multifit_nlinear_free (w);
}
