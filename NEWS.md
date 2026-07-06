# Changes in version 1.1.2

* `Poisson_RDP()` now returns `splitvec`, which encodes the final hereditary dyadic partition. The function's output is now a list of two vectors: the estimated intensity and the final hereditary dyadic partition.
* Improved documentation for the `Poisson_RDP()` estimator and its tree representation.
* Corrected handling of zero-valued fitted Poisson means in selected terminal intervals.
* All calls within the package to the function `Poisson_RDP()` have been corrected accordingly.
