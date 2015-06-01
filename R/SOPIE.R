SOPIE <-
function(data = NULL ,h=1,to=1,alpha=0.05,g=20,r=10,m=1,grid=512)
{
if (is.null(data)) 
        stop("\n", "A data vector is needed")

if ((trunc(g) != g) || (trunc(r) != r) || (trunc(m) != m) || (trunc(grid) != grid) || is.na(r) || is.na(g))
        stop("\n", "g, r, m and grid must be integers")

if (length(data) < (g*r)) 
{
  g<-1 
  r<-1
}

cl<-match.call()

sp<-findh(data,h,to)
min_points<-circ.kernel(data,sp,to=1,grid,m)

cl<-match.call()

a_hat<-a.estimate(data,to,min_points$minimum,alpha,g,r)
b_hat<-b.estimate(data,to,min_points$minimum,alpha,g,r)



table<-matrix(c(a_hat$summary,b_hat$summary),ncol=5,byrow=T,dimnames=list(c("a-hat","b-hat"),c("Cramer von Mises","Kolmogorov-Smirnoff",
"Anderson-Darling","Rayleigh","MEDIAN")))

a_hat$General$call<-cl
a_hat$General$kernel_bandwidth<-sp
a_hat$General$kernel_grid<-grid

result<-list(Summary=table,General=a_hat$General)
graphics.off()
#tiff("C:\\crab.tiff", width = 16, height = 15, units = 'cm', res = 300)
hist(data, nclass=100,freq = FALSE, col = "grey",main="Histogram, kernel density estimator and est. off-pulse interval",xlab="",ylab="")
height<-0.5*(max(hist(data, nclass=20,plot=FALSE)$density)+min(hist(data, nclass=20,plot=FALSE)$density))

lines(c((table[1,5]+0.05),table[1,5],table[1,5],table[1,5]+0.05),c(height,height,-0.015,-0.015),col="blue",lwd=2)
lines(c((table[2,5]-0.05),table[2,5],table[2,5],table[2,5]-0.05),c(height,height,-0.015,-0.015),col="blue",lwd=2)


#lines(c((table[1,1]+0.05),table[1,1],table[1,1],table[1,1]+0.05),c(height,height,-0.15,-0.15),col="blue")
#lines(c(table[1,1],table[1,1]),c(height,-0.15),col="blue")
#lines(c((table[1,1]+0.05),table[1,1]),c(-0.15,-0.15),col="blue")

#lines(c((table[1,1]+0.05),table[1,1]),c(0.5*(max(hist(data, nclass=20,plot=FALSE)$density)),0.5*(max(hist(data, nclass=20,plot=FALSE)$density))),col="blue")
#lines(c(table[1,1],table[1,1]),c(0.5*(max(hist(data, nclass=20,plot=FALSE)$density)),-0.15),col="blue")
#lines(c((table[1,1]+0.05),table[1,1]),c(-0.15,-0.15),col="blue")

#lines(c((table[2,1]-0.05),table[2,1],table[2,1],table[2,1]-0.05),c(height,height,-0.15,-0.15),col="blue")
#lines(c((table[2,1]-0.05),table[2,1]),c(height,height),col="blue")
#lines(c(table[2,1],table[2,1]),c(height,-0.15),col="blue")
#lines(c((table[2,1]-0.05),table[2,1]),c(-0.15,-0.15),col="blue")

#lines(c((table[1,2]+0.05),table[1,2],table[1,2],table[1,2]+0.05),c(height,height,-0.15,-0.15),col="green")
#lines(c((table[1,2]+0.05),table[1,2]),c(height,height),col="green")
#lines(c(table[1,2],table[1,2]),c(0.5*(max(hist(data, nclass=20,plot=FALSE)$density)),-0.15),col="green")
#lines(c((table[1,2]+0.05),table[1,2]),c(-0.15,-0.15),col="green")


#lines(c((table[2,2]-0.05),table[2,2],table[2,2],table[2,2]-0.05),c(height,height,-0.15,-0.15),col="green")
#lines(c((table[2,2]-0.05),table[2,2]),c(height,height),col="green")
#lines(c(table[2,2],table[2,2]),c(height,-0.15),col="green")
#lines(c((table[2,2]-0.05),table[2,2]),c(-0.15,-0.15),col="green")

#lines(c((table[1,3]+0.05),table[1,3],table[1,3],table[1,3]+0.05),c(height,height,-0.15,-0.15),col="yellow")
#lines(c((table[1,3]+0.05),table[1,3]),c(height,height),col="yellow")
#lines(c(table[1,3],table[1,3]),c(height,-0.15),col="yellow")
#lines(c((table[1,3]+0.05),table[1,3]),c(-0.15,-0.15),col="yellow")

#lines(c((table[2,3]-0.05),table[2,3],table[2,3],table[2,3]-0.05),c(height,height,-0.15,-0.15),col="yellow")
#lines(c((table[2,3]-0.05),table[2,3]),c(height,height),col="yellow")
#lines(c(table[2,3],table[2,3]),c(height,-0.15),col="yellow")
#lines(c((table[2,3]-0.05),table[2,3]),c(-0.15,-0.15),col="yellow")

#lines(c((table[1,4]+0.05),table[1,4],table[1,4],table[1,4]+0.05),c(height,height,-0.15,-0.15),col="purple")
#lines(c((table[1,4]+0.05),table[1,4]),c(height,height),col="purple")
#lines(c(table[1,4],table[1,4]),c(height,-0.15),col="purple")
#lines(c((table[1,4]+0.05),table[1,4]),c(-0.15,-0.15),col="purple")

#lines(c((table[2,4]-0.05),table[2,4],table[2,4],table[2,4]-0.05),c(height,height,-0.15,-0.15),col="purple")
#lines(c((table[2,4]-0.05),table[2,4]),c(height,height),col="purple")
#lines(c(table[2,4],table[2,4]),c(height,-0.15),col="purple")
#lines(c((table[2,4]-0.05),table[2,4]),c(-0.15,-0.15),col="purple")




lines(min_points$x,min_points$y/(sum(min_points$y)*(1/grid)),col="red",lwd=2,lty=1)
#dev.off()
result

}
