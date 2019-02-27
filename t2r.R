## In a linear regression model, 
## the conversion between t-statistic and partial R2
t2r = function(t, df) {
  sign(t)*sqrt(t^2/(t^2+df))
}

## computation of standar error from r-squared and n
rsq2se = function(r, n) {
  sqrt((1-r^2)/(n-2))
}