`computeY` <-
function(g,x) apply(x, 2, function(z) tapply(z,g,mean))

