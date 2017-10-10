/* monte/plain.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2009 Michael Booth
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

/* Plain Monte-Carlo. */

/* Author: MJB */
/* Modified by IH to return a suggested split dimension for hdintegrator */
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl_monte_plain2.h>

int
gsl_monte_plain_integrate2 (const gsl_monte_function * f,
                           const double xl[], const double xu[],
                           const size_t dim,
                           const size_t calls,
                           gsl_rng * r,
                           gsl_monte_plain_state * state,
                           double *result, double *abserr, int* split_dims)
{
  double vol, m = 0, q = 0;
  double *x = state->x;
  size_t n, i;

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

  /* Compute the volume of the region */

  vol = 1;

  for (i = 0; i < dim; i++)
    {
      vol *= xu[i] - xl[i];
    }

  double* quad_avgs = (double*) calloc(2*dim, sizeof(double));
  int* quad_nr = (int*) calloc(2*dim, sizeof(int));
  if (quad_avgs == NULL || quad_nr == NULL) {
    return GSL_FAILURE;
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

        for (unsigned int d = 0; d < dim; d++) {
          if (x[d] - xl[d] < xu[d] - x[d]) {
            quad_avgs[d] += fval;
            quad_nr[d]++;
          } else {
            quad_avgs[d+1] += fval;
            quad_nr[d+1]++;
          }
        }

        /* recurrence for mean and variance */

        double d = fval - m;
        m += d / (n + 1.0);
        q += d * d * (n / (n + 1.0));
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

  double max_diff = -1;
  int max_diff_d = 0;
  for (size_t d = 0; d < dim; d++) {
    quad_avgs[d] /= quad_nr[d];
    quad_avgs[d+1] /= quad_nr[d+1];
    const double diff = abs(quad_avgs[d] - quad_avgs[d+1]);
    if (max_diff < diff) {
      max_diff = diff;
      max_diff_d = d;
    }
  }
  (*(split_dims + max_diff_d))++;

  return GSL_SUCCESS;
}
