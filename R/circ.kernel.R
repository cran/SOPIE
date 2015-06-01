circ.kernel <-
function(data,sp,to=1,grid=512,m=1)
{

if (round(to,2)==round(2*pi,2))
{
data<-unique(data/(2*pi))
}

feta<-seq(0,1,length=(grid+1))
f_feta_calc<-matrix(nrow=(grid+1),ncol=length(data))
f_feta<-vector()

for (i in 1:(grid+1))
{ 
      f_feta_calc[i,]<-pmin(abs(feta[i]-data),(1)-(abs(feta[i]-data)))/sp 



f_feta[i]<- sum(1-((1-cos(f_feta_calc[i,][f_feta_calc[i,]<=1])))^2)/(sp*length(data))
}


minimum_punte<-feta[order(f_feta)[1:m]]
results<-list(x=feta,y=f_feta,minimum=minimum_punte)
return(results)
}
