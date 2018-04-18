#ifndef HAVE_GSLFUNCTIONS_H
#define HAVE_GSLFUNCTIONS_H

double rosenbrock(const gsl_vector *v, void *params);
void rosenbrock_df(const gsl_vector *v, void *params, gsl_vector *df);
void rosenbrock_fdf (const gsl_vector *v, void *params, double *f, gsl_vector *df);

double himmelblau(const gsl_vector *v, void *params);
void himmelblau_df(const gsl_vector *v, void *params, gsl_vector *df);
void himmelblau_fdf (const gsl_vector *v, void *params, double *f, gsl_vector *df);


int gsl_multiroot(double f(const gsl_vector *v, void *params), void df(const gsl_vector *v, void *params, gsl_vector *df), void fdf (const gsl_vector *v, void *params, double *f, gsl_vector *df), double x0, double y0 );

#endif
