Integrands to be used by HDIntegrator [main_mpi.py](../main_mpi.py).

# I/O format

Integrands read integration parameters from stdin and output the results to
stdout. For every line of input one line of output is produced. Some integrands
take additional arguments on the command line, use `--help` as an argument to
get help on the rest. Input format for every line is:

    number_of_points x_min x_max y_min y_max z_min ...

Where mins and maxs are the minimum and maximum extent of integration over the
respective dimensions. Every line of output is:

    value error dim

where dim is the suggested dimension to split given integration
volume for more accurate calculation.

# Examples

Command:

    echo 1e6 -1 1 -2 2 -3 3 | ./turbulence2_miser --corr1 -1 --corr2 -1

Expected output:

    1.386146698197008e+00 1.502145041181887e-03 1

Command:

    echo 1e6 -1 1 -2 2 -3 3 | ./turbulence2_vegas --corr1 1 --corr2 2

Expected output:

    2.769778993537326e-01 3.291349703842326e-05 2

Command:

    echo 1e6 -1 1 -2 2 -3 3 | ./turbulence2_miser --corr1 1 --corr2 2

Expected output:

    2.774590113466751e-01 2.751452266693708e-04 2
