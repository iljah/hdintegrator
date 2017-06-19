/* monte/gsl_monte_vegas.h
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

/* header for the gsl "vegas" routines.  Mike Booth, May 1998 */
/* Modified by IH to return a suggested split dimension for hdintegrator */
#ifndef __GSL_MONTE_VEGAS2_H__
#define __GSL_MONTE_VEGAS2_H__

#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

int gsl_monte_vegas_integrate2(gsl_monte_function * f, 
                              double xl[], double xu[], 
                              size_t dim, size_t calls,
                              gsl_rng * r,
                              gsl_monte_vegas_state *state,
                              double* result, double* abserr, int* split_dims);
__END_DECLS

#endif /* __GSL_MONTE_VEGAS2_H__ */

