# Legendre_polynomials

Code Information
================

For a maximum degree $N$, the following code computes the fully normalized Legendre polynomials $\bar{P_n}(\cos(\theta))$ as a function of co-latitude $\theta$.

There are two source codes provided for the computation: one in C and one in Fortran.

The source codes have the following names:

- **geolp.c**
- **geolp.f90**

Code in C
---------

Before running the *geolp.c* code, you must insert the values for the maximum degree $N$ and the co-latitude $\theta$ in lines **33** and **34**, respectively.
For example:

```c
N = 10; /* The maximum degree */
theta = 30.0; /* The co-latitude in degrees */
```

Enter the following command at the Linux command line to compile the code:

```bash
gcc -Wall geolp.c -o geolpc -lm
```

This command will build a single executable with the name *geolpc*. Enter the following command in the Linux command line to execute the program:

```bash
./geolpc
```

Below are the results:

```
N = 10
theta =  3.00000000e+01 deg

n       Pn(theta)       
0        1.0000000000000000e+00
1        1.5000000000000000e+00
2        1.3975424859373691e+00
3        8.5923294280422069e-01
4        7.0312500000000708e-02
5       -7.4051002865529225e-01
6       -1.3485606821315501e+00
7       -1.5886127476866183e+00
8       -1.3968077255243430e+00
9       -8.2633914801326269e-01
10      -3.2252701405311168e-02
```

---

Code in Fortran
---------------

Similarly, you need to enter the values for the co-latitude $\theta$ and maximum degree $N$ in lines **32** and **33**, respectively, in order to run the *geolp.f90* code.

For instance:

```fortran
NX=10 ! The maximum degree, less than NMAX
theta=30.0d0 ! The co-latitude in degrees
```

In the Linux command line, enter the following command to compile this code:

```bash
gfortran -Wall geolp.f90 -o geolpf90 -lm
```

One executable with the name *geolpf90* will be produced when you type this command.

In the Linux command line, enter the following command to execute the program:

```bash
./geolpf90
```

The outcomes are displayed below:

```
                 NX=                  10
         theta(deg)=      3.00000000E+01
         n                Pn(theta)
         0   1.0000000000000000E+00
         1   1.5000000000000000E+00
         2   1.3975424859373691E+00
         3   8.5923294280422069E-01
         4   7.0312500000000708E-02
         5  -7.4051002865529225E-01
         6  -1.3485606821315501E+00
         7  -1.5886127476866183E+00
         8  -1.3968077255243430E+00
         9  -8.2633914801326269E-01
        10  -3.2252701405311231E-02
```

*For individuals who are unfamiliar with the Linux command line, using an Integrated Development Environment (IDE) can facilitate the compilation and execution process.
