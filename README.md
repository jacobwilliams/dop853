# dop853

This is a modern Fortran (2003/2008) implementation of Hairer's DOP853 ODE solver. The original FORTRAN 77 code has been extensively refactored, and is now object-oriented and thread-safe, with an easy-to-use class interface.  DOP853 is an explicit Runge-Kutta method of order 8(5,3) due to Dormand & Prince (with stepsize control and dense output).

This project is hosted on [GitHub](https://github.com/jacobwilliams/dop853).

## References

1. E. Hairer, S.P. Norsett and G. Wanner, "[Solving ordinary
   Differential Equations I. Nonstiff Problems](http://www.unige.ch/~hairer/books.html)", 2nd edition.
   Springer Series in Computational Mathematics,
   Springer-Verlag (1993).
2. Ernst Hairer's website: [Fortran and Matlab Codes](http://www.unige.ch/~hairer/software.html)

## License

* [Original license for Hairer's codes](http://www.unige.ch/~hairer/prog/licence.txt).
* The updates are released under a [similar BSD-style license](https://raw.githubusercontent.com/jacobwilliams/dop853/master/LICENSE?token=AGLts7b7U0xmmh6R8kZCsRNyVG7m71KZks5WkIeIwA%3D%3D).
