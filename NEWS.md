# cdom 0.1.1 (Unreleased)

- `cdom_fit_exponential()` has been renamed to `cdom_exponential()`.

- All fitting functions are now using `purrr:safely()` to safely fit non-linear models.

- Implemented S3 funnctions `coef()`, `predict()` and `plot()` for object returned by `cdom_exponential()`. 
 
# cdom 0.1.0

- First version of cdom
