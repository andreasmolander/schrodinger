# Numerical solver for the Schrödinger equation

A simple Python program for numerically solving the <a href="https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation" target="_blank">Schrödinger equation</a> in one dimension.

This was originally developed as part of my Bachelor's thesis. For the original code used for the thesis, see branch <a href="https://github.com/andreasmolander/schrodinger/tree/bachelor's-thesis" target="_blank">bachelor's-thesis</a>.

## About

When the program is run it will open an animated plot that visualizes the evolvement of the chosen initial wave function under the influence of a chosen potential according to the Schrödinger equation.

The currently available algorithms for solving the Schrödinger equation are the <a href="https://en.wikipedia.org/wiki/Finite_difference_method" target="_blank">finite difference methods</a>:

* <a href="https://en.wikipedia.org/wiki/Euler_method" target="_blank">Euler method</a>
* <a href="https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method" target="_blank">Crank-Nicolson method</a>

## <a name="usage"></a> Usage

To run the program, in the src/ directory simply run `simulate.py`, e.g.

```shell
./simulate.py
```

The default setup is an arbitrary harmonic oscillator with a gaussian wave packet used as the wave function.

## Future plans

* More dimensions
* More solving algorithms
* More potentials