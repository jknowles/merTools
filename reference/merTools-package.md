# merTools: Provides methods for extracting and exploring results from merMod objects in the lme4 package.

The merTools package contains convenience tools for extracting useful
information from and exploring the implications of merMod objects
created by the lme4 package. These convenience functions are especially
useful for merMod objects that take a long time to estimate due to their
complexity or because they are estimated on very large samples.

## Details

See the vignettes for usage examples

## merMod extraction/utility functions

- [`fastdisp`](https://jknowles.github.io/merTools/reference/fastdisp.md)

- [`superFactor`](https://jknowles.github.io/merTools/reference/superFactor.md)

- [`REextract`](https://jknowles.github.io/merTools/reference/REextract.md)

- [`REsim`](https://jknowles.github.io/merTools/reference/REsim.md)

- [`FEsim`](https://jknowles.github.io/merTools/reference/FEsim.md)

- [`RMSE.merMod`](https://jknowles.github.io/merTools/reference/RMSE.merMod.md)

- [`thetaExtract`](https://jknowles.github.io/merTools/reference/thetaExtract.md)

- [`REquantile`](https://jknowles.github.io/merTools/reference/REquantile.md)

## merMod exploration functions

- [`plotREsim`](https://jknowles.github.io/merTools/reference/plotREsim.md)

- [`plotFEsim`](https://jknowles.github.io/merTools/reference/plotFEsim.md)

- [`plotREimpact`](https://jknowles.github.io/merTools/reference/plotREimpact.md)

- [`draw`](https://jknowles.github.io/merTools/reference/draw.md)

- [`wiggle`](https://jknowles.github.io/merTools/reference/wiggle.md)

- [`subBoot`](https://jknowles.github.io/merTools/reference/subBoot.md)

- [`predictInterval`](https://jknowles.github.io/merTools/reference/predictInterval.md)

- [`expectedRank`](https://jknowles.github.io/merTools/reference/expectedRank.md)

- [`REimpact`](https://jknowles.github.io/merTools/reference/REimpact.md)

- [`shinyMer`](https://jknowles.github.io/merTools/reference/shinyMer.md)

## Influences and acknowledgements

The simulation approach used by
[`predictInterval`](https://jknowles.github.io/merTools/reference/predictInterval.md)
to propagate uncertainty from both the fixed and random effects draws on
the framework described by Gelman and Hill (2007) and implemented in the
arm package's [`sim()`](https://rdrr.io/pkg/arm/man/sim.html) function.
The package is built on top of the model objects fit by lme4 (Bates et
al. 2015) and is intended to complement
[`lme4::bootMer()`](https://rdrr.io/pkg/lme4/man/bootMer.html) by
offering a faster, simulation-based alternative for large models. We are
grateful to the authors and maintainers of those packages.

## References

Gelman, A. and Hill, J. (2007). *Data Analysis Using Regression and
Multilevel/Hierarchical Models*. Cambridge University Press.

Bates, D., Maechler, M., Bolker, B., and Walker, S. (2015). Fitting
Linear Mixed-Effects Models Using lme4. *Journal of Statistical
Software*, 67(1), 1-48.
[doi:10.18637/jss.v067.i01](https://doi.org/10.18637/jss.v067.i01)

## See also

Useful links:

- <https://jknowles.github.io/merTools/>

- <https://github.com/jknowles/merTools>

- Report bugs at <https://github.com/jknowles/merTools/issues>

## Author

**Maintainer**: Jared E. Knowles <jared@civilytics.com>

Authors:

- Jared E. Knowles <jared@civilytics.com>

- Carl Frederick <carlbfrederick@gmail.com>

Other contributors:

- Alex Whitworth <whitworth.alex@gmail.com> \[contributor\]
