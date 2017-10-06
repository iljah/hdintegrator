---
title: 'High-Dimensional Integrator'
tags:
  - numerical integration
  - message passing interface
  - arbitrary number of dimensions
authors:
 - name: Ilja Honkonen
   orcid: 0000-0002-9542-5866
   affiliation: 1
affiliations:
 - name: Finnish Meteorological Institute
   index: 1
date: 2017-10-06
bibliography: paper.bib
---

# Summary

HDIntegrator [@hdint] is a program for parallel numerical integration of functions in arbitrary number of dimensions.
It is implemented in Python and parallelized using the Message Passing Interface [@mpi] via mpi4py [@mpi4py].
HDIntegrator is a wrapper program that divides the integration region into rectangular subvolumes and uses a separate serial program (the integrand) for evaluating the integral within those subvolumes.
Integration of a subvolume is stopped when one of the user-specified criteria for convergence is fulfilled, otherwise the subvolume is divided into two smaller subvolumes and integration is continued.
Communication between HDIntegrator and integrand is handled via standard input and output in ASCII format.


# Statement of need

This program was developed due to the need for a free and open source software (FOSS) for adaptive and parallel calculation of numerical integrals.
Several FOSS packages exist for numerical integration [e.g. @gsl, @cubature] but none seem to support parallel integration in arbitrary number of dimensions using a supercomputer.
HDIntegrator is implemented as a parallel wrapper program that can take advantage of existing serial software for numerical integration.
This also allows one to easily develop and use efficient integrands separately from the parallel wrapper, that use specialized hardware or environments such as GPUs and even cloud computing.
This program is currently used for studying turbulence via numerical integration of functional integrals.


# References
