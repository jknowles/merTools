# Pin RNG algorithm and version across R releases so that results generated
# with set.seed() are reproducible regardless of which R version runs the
# tests. Without this, R-release vs R-devel vs R-oldrel-1 on CI can drift on
# the same seed because default RNG defaults have shifted across R versions
# (e.g., sample.kind changed in R 3.6).
RNGversion("4.1.0")
RNGkind(kind = "Mersenne-Twister",
        normal.kind = "Inversion",
        sample.kind = "Rejection")
