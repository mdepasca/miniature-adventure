CFoM.Ia <- function(N.Ia.tot, N.Ia.true, N.Ia.false, W.Ia.false=3){
  fom <- (1/N.Ia.tot)*(N.Ia.true**2/(N.Ia.true+W.Ia.false*N.Ia.false))
  return(fom)
}



