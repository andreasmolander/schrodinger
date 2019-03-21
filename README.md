# Numerical solver of the Schrödinger equation

A simple Python program for numerically solving the [Schrödinger equation](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation) in one dimension.

This was originally developed as part of my Bachelor's thesis. For the original code used for the thesis, see branch [bachelor's-thesis](/andreasmolander/schrodinger/tree/bachelor's-thesis).

## About

When the program is run it will open an animated plot that visualizes the evolvement of the chosen initial wave function under the influence of a chosen potential according to the Schrödinger equation.

The currently available algorithms for solving the Schrödinger equation are the [finite difference methods](https://en.wikipedia.org/wiki/Finite_difference_method):

* [Euler method](https://en.wikipedia.org/wiki/Euler_method)
* [Crank-Nicolson method](https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method)

## <a name="usage"></a> Usage

To run the program, in the src/ directory simply run `simulate.py`, e.g.

```bash
./simulate.py
```

The default setup is an arbitrary harmonic oscillator with a gaussian wave packet used as the wave function.

## Future plans

* More dimensions
* More solving algorithms
* More potentials