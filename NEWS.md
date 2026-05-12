# cdom 0.1.1

- `cdom_fit_exponential()` has been renamed to `cdom_exponential()`.

- All fitting functions are now using `purrr:safely()` to safely fit non-linear models.

- Implemented S3 funnctions `coef()`, `predict()` and `plot()` for object returned by `cdom_exponential()`.

- Internal changes to fix CRAN check issues.

# cdom 0.1.0

- First version of cdom
