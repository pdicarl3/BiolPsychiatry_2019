aprime <- function(hit, fa)
{ 
  if(fa <= .5 & hit >= .5)
  {
    a <- .75 + (hit -fa)/4 - fa*(1-hit)
  } else if(fa <= hit & hit <= .5)
  {
    a <- .75 + (hit -fa)/4 - fa/(4*hit)
  } else {
    a <- .75 + (hit -fa)/4 - (1-hit)/(4 * (1-fa))
  } 
  return(a)
}
#######################################################

beta <- function(hit, fa)
{
  if(fa <= .5 & hit >= .5)
  {
    b <-(5-4*hit)/(1+4*fa)
  } else if(fa <= hit & hit <= .5)
  {
    b <-(hit^2+hit)/(hit^2+fa) 
  } else { 
    b <- ((1-fa)^2 + (1-hit))/((1-fa)^2 + (1-fa))
  } 
  return(b)
}
#######################################################

dprime <- function(hit, fa) 
{
  qnorm(hit) - qnorm(fa)
}
#####################################################



