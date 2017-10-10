/* monte/miser.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Michael Booth
 * Copyright (C) 2009 Brian Gough
 * Copyright 2017 Ilja Honkonen
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* MISER.  Based on the algorithm described in the following article,

   W.H. Press, G.R. Farrar, "Recursive Stratified Sampling for
   Multidimensional Monte Carlo Integration", Computers in Physics,
   v4 (1990), pp190-195.

*/

/* Author: MJB */
/* Modified by Brian Gough 12/2000 */
/* Modified by IH to return a suggested split dimension for hdintegrator */
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl_monte_miser2.h>

static int
estimate_corrmc (gsl_monte_function * f,
                 const double xl[], const double xu[],
                 size_t dim, size_t calls,
                 gsl_rng * r,
                 gsl_monte_miser_state * state,
                 double *result, double *abserr,
                 const double xmid[], double sigma_l[], double sigma_r[])
{
  size_t i, n;
  
  double *x = state->x;
  double *fsum_l = state->fsum_l;
  double *fsum_r = state->fsum_r;
  double *fsum2_l = state->fsum2_l;
  double *fsum2_r = state->fsum2_r;
  size_t *hits_l = state->hits_l;
  size_t *hits_r = state->hits_r;

  double m = 0.0, q = 0.0; 
  double vol = 1.0;

  for (i = 0; i < dim; i++)
    {
      vol *= xu[i] - xl[i];
      hits_l[i] = hits_r[i] = 0;
      fsum_l[i] = fsum_r[i] = 0.0;
      fsum2_l[i] = fsum2_r[i] = 0.0;
      sigma_l[i] = sigma_r[i] = -1;
    }

  for (n = 0; n < calls; n++)
    {
      double fval;
      
      unsigned int j = (n/2) % dim;
      unsigned int side = (n % 2);

      for (i = 0; i < dim; i++)
        {
          double z = gsl_rng_uniform_pos (r) ;

          if (i != j) 
            {
              x[i] = xl[i] + z * (xu[i] - xl[i]);
            }
          else
            {
              if (side == 0) 
                {
                  x[i] = xmid[i] + z * (xu[i] - xmid[i]);
                }
              else
                {
                  x[i] = xl[i] + z * (xmid[i] - xl[i]);
                }
            }
        }

      fval = GSL_MONTE_FN_EVAL (f, x);

      /* recurrence for mean and variance */
      {
        double d = fval - m;
        m += d / (n + 1.0);
        q += d * d * (n / (n + 1.0));
      }

      /* compute the variances on each side of the bisection */
      for (i = 0; i < dim; i++)
        {
          if (x[i] <= xmid[i])
            {
              fsum_l[i] += fval;
              fsum2_l[i] += fval * fval;
              hits_l[i]++;
            }
          else
            {
              fsum_r[i] += fval;
              fsum2_r[i] += fval * fval;
              hits_r[i]++;
            }
        }
    }

  for (i = 0; i < dim; i++)
    {
      double fraction_l = (xmid[i] - xl[i]) / (xu[i] - xl[i]);

      if (hits_l[i] > 0)
        {
          fsum_l[i] /= hits_l[i];
          sigma_l[i] = sqrt (fsum2_l[i] - fsum_l[i] * fsum_l[i] / hits_l[i]);
          sigma_l[i] *= fraction_l * vol / hits_l[i];
        }

      if (hits_r[i] > 0)
        {
          fsum_r[i] /= hits_r[i];
          sigma_r[i] = sqrt (fsum2_r[i] - fsum_r[i] * fsum_r[i] / hits_r[i]);
          sigma_r[i] *= (1 - fraction_l) * vol / hits_r[i];
        }
    }

  *result = vol * m;

  if (calls < 2)
    {
      *abserr = GSL_POSINF;
    }
  else
    {
      *abserr = vol * sqrt (q / (calls * (calls - 1.0)));
    }

  return GSL_SUCCESS;
}

int
gsl_monte_miser_integrate2 (gsl_monte_function * f,
                           const double xl[], const double xu[],
                           size_t dim, size_t calls,
                           gsl_rng * r,
                           gsl_monte_miser_state * state,
                           double *result, double *abserr, int* split_dims)
{
  size_t n, estimate_calls, calls_l, calls_r;
  const size_t min_calls = state->min_calls;
  size_t i;
  size_t i_bisect;
  int found_best;

  double res_est = 0, err_est = 0;
  double res_r = 0, err_r = 0, res_l = 0, err_l = 0;
  double xbi_l, xbi_m, xbi_r, s;

  double vol;
  double weight_l, weight_r;

  double *x = state->x;
  double *xmid = state->xmid;
  double *sigma_l = state->sigma_l, *sigma_r = state->sigma_r;

  if (dim != state->dim)
    {
      GSL_ERROR ("number of dimensions must match allocated size", GSL_EINVAL);
    }

  for (i = 0; i < dim; i++)
    {
      if (xu[i] <= xl[i])
        {
          GSL_ERROR ("xu must be greater than xl", GSL_EINVAL);
        }

      if (xu[i] - xl[i] > GSL_DBL_MAX)
        {
          GSL_ERROR ("Range of integration is too large, please rescale",
                     GSL_EINVAL);
        }
    }

  if (state->alpha < 0)
    {
      GSL_ERROR ("alpha must be non-negative", GSL_EINVAL);
    }

  /* Compute volume */

  vol = 1;

  for (i = 0; i < dim; i++)
    {
      vol *= xu[i] - xl[i];
    }

  if (calls < state->min_calls_per_bisection)
    {
      double m = 0.0, q = 0.0;

      if (calls < 2)
        {
          GSL_ERROR ("insufficient calls for subvolume", GSL_EFAILED);
        }

      for (n = 0; n < calls; n++)
        {
          /* Choose a random point in the integration region */

          for (i = 0; i < dim; i++)
            {
              x[i] = xl[i] + gsl_rng_uniform_pos (r) * (xu[i] - xl[i]);
            }

          {
            double fval = GSL_MONTE_FN_EVAL (f, x);

            /* recurrence for mean and variance */

            double d = fval - m;
            m += d / (n + 1.0);
            q += d * d * (n / (n + 1.0));
          }
        }

      *result = vol * m;

      *abserr = vol * sqrt (q / (calls * (calls - 1.0)));

      return GSL_SUCCESS;
    }

  estimate_calls = GSL_MAX (min_calls, calls * (state->estimate_frac));

  if (estimate_calls < 4 * dim)
    {
      GSL_ERROR ("insufficient calls to sample all halfspaces", GSL_ESANITY);
    }

  /* Flip coins to bisect the integration region with some fuzz */

  for (i = 0; i < dim; i++)
    {
      s = (gsl_rng_uniform (r) - 0.5) >= 0.0 ? state->dither : -state->dither;
      state->xmid[i] = (0.5 + s) * xl[i] + (0.5 - s) * xu[i];
    }

  /* The idea is to chose the direction to bisect based on which will
     give the smallest total variance.  We could (and may do so later)
     use MC to compute these variances.  But the NR guys simply estimate
     the variances by finding the min and max function values 
     for each half-region for each bisection. */

  estimate_corrmc (f, xl, xu, dim, estimate_calls,
                   r, state, &res_est, &err_est, xmid, sigma_l, sigma_r);

  /* We have now used up some calls for the estimation */

  calls -= estimate_calls;

  /* Now find direction with the smallest total "variance" */

  {
    double best_var = GSL_DBL_MAX;
    double beta = 2.0 / (1.0 + state->alpha);
    found_best = 0;
    i_bisect = 0;
    weight_l = weight_r = 1.0;

    for (i = 0; i < dim; i++)
      {
        if (sigma_l[i] >= 0 && sigma_r[i] >= 0)
          {
            /* estimates are okay */
            double var = pow (sigma_l[i], beta) + pow (sigma_r[i], beta);

            if (var <= best_var)
              {
                found_best = 1;
                best_var = var;
                i_bisect = i;
                weight_l = pow (sigma_l[i], beta);
                weight_r = pow (sigma_r[i], beta);

                if (weight_l == 0 && weight_r == 0)
                  {
                    weight_l = 1;
                    weight_r = 1;
                  }
              }
          }
        else
          {
            if (sigma_l[i] < 0)
              {
                GSL_ERROR ("no points in left-half space!", GSL_ESANITY);
              }
            if (sigma_r[i] < 0)
              {
                GSL_ERROR ("no points in right-half space!", GSL_ESANITY);
              }
          }
      }
  }

  if (!found_best)
    {
      /* All estimates were the same, so chose a direction at random */

      i_bisect = gsl_rng_uniform_int (r, dim);
    }

  (*(split_dims + i_bisect))++;

  xbi_l = xl[i_bisect];
  xbi_m = xmid[i_bisect];
  xbi_r = xu[i_bisect];

  /* Get the actual fractional sizes of the two "halves", and
     distribute the remaining calls among them */

  {
    double fraction_l = fabs ((xbi_m - xbi_l) / (xbi_r - xbi_l));
    double fraction_r = 1 - fraction_l;

    double a = fraction_l * weight_l;
    double b = fraction_r * weight_r;

    calls_l = min_calls + (calls - 2 * min_calls) * a / (a + b);
    calls_r = min_calls + (calls - 2 * min_calls) * b / (a + b);
  }

  /* Compute the integral for the left hand side of the bisection */

  /* Due to the recursive nature of the algorithm we must allocate
     some new memory for each recursive call */

  {
    int status;

    double *xu_tmp = (double *) malloc (dim * sizeof (double));

    if (xu_tmp == 0)
      {
        GSL_ERROR_VAL ("out of memory for left workspace", GSL_ENOMEM, 0);
      }

    for (i = 0; i < dim; i++)
      {
        xu_tmp[i] = xu[i];
      }

    xu_tmp[i_bisect] = xbi_m;

    status = gsl_monte_miser_integrate2 (f, xl, xu_tmp,
                                        dim, calls_l, r, state,
                                        &res_l, &err_l, split_dims);
    free (xu_tmp);

    if (status != GSL_SUCCESS)
      {
        return status;
      }
  }

  /* Compute the integral for the right hand side of the bisection */

  {
    int status;

    double *xl_tmp = (double *) malloc (dim * sizeof (double));

    if (xl_tmp == 0)
      {
        GSL_ERROR_VAL ("out of memory for right workspace", GSL_ENOMEM, 0);
      }

    for (i = 0; i < dim; i++)
      {
        xl_tmp[i] = xl[i];
      }

    xl_tmp[i_bisect] = xbi_m;

    status = gsl_monte_miser_integrate2 (f, xl_tmp, xu,
                                        dim, calls_r, r, state,
                                        &res_r, &err_r, split_dims);
    free (xl_tmp);

    if (status != GSL_SUCCESS)
      {
        return status;
      }
  }

  *result = res_l + res_r;
  *abserr = sqrt (err_l * err_l + err_r * err_r);

  return GSL_SUCCESS;
}
