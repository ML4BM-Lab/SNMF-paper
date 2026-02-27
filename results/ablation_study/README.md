## When does the Gaussian kernel become the identity matrix?

We construct the spatial similarity matrix using a Gaussian (RBF) kernel:

\[
S_{ij} = \exp(-\gamma d_{ij}^2)
\]

where:
- \( d_{ij} \) is the Euclidean distance between spots \( i \) and \( j \)
- \( \gamma \) controls how fast similarity decays with distance

### Diagonal entries

Since \( d_{ii} = 0 \),

\[
S_{ii} = \exp(0) = 1
\]

This holds for any value of \( \gamma \).

### Off-diagonal entries

For \( i \ne j \),

\[
S_{ij} = \exp(-\gamma d_{ij}^2)
\]

The exponential function never reaches exactly zero for finite \( \gamma \), so the matrix can only become the identity matrix in the limit:

\[
\gamma \to \infty
\]

### Practical identity (with thresholding)

In practice, we threshold small values:

```r
S[S < 1e-3] <- 0
```

To force all off-diagonal entries below this threshold, we require:

\[
\exp(-\gamma d_{\min}^2) < 10^{-3}
\]

Taking logs:

\[
\gamma > \frac{\log(1000)}{d_{\min}^2}
\]

Since:

\[
\log(1000) \approx 6.91
\]

we obtain the practical condition:

\[
\gamma > \frac{6.91}{d_{\min}^2}
\]

where \( d_{\min} \) is the smallest non-zero pairwise distance between spots.

### Interpretation

If

\[
\gamma > \frac{6.91}{d_{\min}^2}
\]

then all off-diagonal entries fall below `1e-3`, are thresholded to zero, and the matrix becomes the identity matrix after row normalization.