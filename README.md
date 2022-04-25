# Curvogram

Envisioned by Ryan Schenck with development of method from Jeffrey West.

Method for plotting histograms on a curve using a West Transform. Originally, the term 'histo' means mast, as in a sailboat mast that holds the luff of a sail. Here we call these a curvogram, which is a histogram on a curve that is given by a user defined function. The transformed coordinates are a conformal mapping to the curve/function defined by the user, termed the West Coordinate Transform.

Together, the West Coordinate Transform and the Curvogram represent a fun way to visualize high density information that is often plotted separately for sake of accuracy. However, this separation, may not be necessary if the take away for a plot is easily understood using a curvogram or the curvogram is supplemented with the original histogram elsewhere.

## Examples
Here are two versions with matlab and python code.

### Matlab

Here we show the functions for a normal distribution whose mean is 1.5 and variance of 0.5. The curve functions are f(x)=sin(x), f(x)=cos(x), and f(x)=x^2.
![]("./MATLAB/MATLAB_example.png")


### Python

Here we show the same curve functions are f(x)=sin(x), f(x)=cos(x), and f(x)=x^2. However, we show a beta distriution whose alpha and beta shape parameters are 2 and 5, respectively.



### Important considerations:
1. The curve cannot be that complex. So simple x/y relationships work best.
2. The distributions that you plot on the curve needs to be clear, but how clear depends entirely on the steepness of the curve.
3. The ***plotting area must be square***.
