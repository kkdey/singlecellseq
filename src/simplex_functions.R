
###  Transformation and Reverse transform functions


## the reverse transform function takes a simplex vector and un-simplexes it on (-infty,infty)

reverse_transform=function(x) 
{
  out=array(0,length(x)-1);
  for(k in 2:length(x))
  {
    out[k-1]=log((x[k] + 1e-7)/(x[1] + 1e-7));
  }
  return(out)
}

# the transform function simplexes a vector

transform <- function(y) 
{
  temp =c(1,exp(y));
  out=temp/sum(temp);
  return(out)
}
