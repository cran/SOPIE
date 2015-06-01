von_mises_density <-
function(x_grid,c,kappa,mu=pi)
{
fv_x<-((1/(2*pi*besselI(kappa,0)))*exp(kappa*cos(x_grid-pi))-(exp(+1*kappa*cos(c-mu))/(2*pi*besselI(kappa,0))))/(ppvonmises((2*pi)-c,pi,kappa) - ppvonmises(c,pi,kappa) - (((pi-c)/(pi*besselI(kappa,0)))*exp(+1*kappa*cos(c-mu))))
return (fv_x)
}
