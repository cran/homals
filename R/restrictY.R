`restrictY` <-
function(d,y,r,level) {
  if (sum(y^2) == 0) return(y)
  switch(level,
  "nominal"=return(nominalY(d,y,r)),
  "ordinal"=return(ordinalY(d,y,r)),
  "numerical"=return(numericalY(d,y,r)),
  "polynomial"=return(polynomialY(d,y,r)))
}

