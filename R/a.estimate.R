a.estimate <-
function(data,to=1,min_points,alpha=0.05,g=1,r=1)
{

cl<-match.call()

if (round(max(data),1)==round(2*pi,1))
{
data<-unique(data/(2*pi))
}

data<-sort(unique(data))

data<-c(data-1,data,data+1)
close_to_min_index<-unique(findInterval(min_points,data)+1)

CVM_rejectionpnt<-numeric(length(close_to_min_index))
AD_rejectionpnt<-numeric(length(close_to_min_index))
KS_rejectionpnt<-numeric(length(close_to_min_index))
rayleigh_rejectionpnt<-numeric(length(close_to_min_index))

for (k in 1:length(close_to_min_index))
{
end<-close_to_min_index[k]
begin<-close_to_min_index[k]-(length(data)/3)



#if (begin == 0) {start<-1} else if (begin<=1) {start<- findInterval(begin,data)} else {start<-begin}
if (begin == 0) {start<-1} else {start<-begin}

if (end<=1) endpoint<-min(findInterval(end,data),length(data)) else endpoint<-end

lengte<-endpoint-start-1

ks<-numeric(trunc(lengte/g))
ks_p<-numeric(trunc(lengte/g))

b<-0
c<-numeric(trunc(lengte/g))
d<-numeric(trunc(lengte/g))
e<-0
f<-0
tel<-0
asq<-numeric(trunc(lengte/g))
asq_p<-numeric(trunc(lengte/g))
CVM_pvalue<-numeric(trunc(lengte/g))
rayleigh_p<-numeric(trunc(lengte/g))
tol<-0.000


j<-1
#j is die versameling nommer

CVM_rejection<-FALSE
KS_rejection<-FALSE
AD_rejection<-FALSE
rayleigh_rejection<-FALSE
rejected<-(CVM_rejection & KS_rejection & AD_rejection & rayleigh_rejection)






while (!rejected & (j<=trunc(lengte/g)))
{
sn<-0
n<-0
sn_star<-0
star_pval<-0

start<-endpoint-(j*g)-1
if ((start < begin) & (begin > 1 )) 
{
start<-begin
} else if (start<1)
{
 start<-1
}

teller<-1
n<-endpoint-start-1

#Bereken W square wat gelyk is aan Cramer Von Mises
#Bereken W square soos in Stephens p101 1ste deel sum_i{ (z_(i) - (2i-1)/(2n)}^2

sn<-((data[(start+1):(endpoint-1)]-data[start])/(data[endpoint]-data[start])-(1:n-0.5)/(n))^2

#cat("start",start, " endpoint",endpoint," sn" ,sn, "\n")

if (n<8)
{
asq[j] <- 0
asq_p[j] <-0.99
}
else 
{
asq[j]<-adgofteststat(data[(start+1):(endpoint-1)],punif,data[start],data[endpoint],tol)
asq_p[j]<-1-adgoftestpval(asq[j],n)
}


#Bereken W square wat gelyk is aan Cramer Von Mises
#Bereken W square soos in Stephens p101 2de deel sum_i{ (z_(i) - (2i-1)/(2n)}^2 + 1/12n

c[j]<-sum(sn)+(1/(12*(n)))

#Bereken KS met ingeboude funksie

ks[j]<-ks.test(data[(start+1):(endpoint-1)],punif,data[start],data[endpoint])$statistic
ks_p[j]<-ks.test(data[(start+1):(endpoint-1)],punif,data[start],data[endpoint])$p.value

d[j]<-n

#Bereken CVM pvalues met Stephens se W^2 ster op p105 Stephens
#cat("j=",j," c[j]=",c[j]," n=",n, "\n")


CVM_pvalue[j]<-CVM_pvall(c[j],n)
#cat("j=",j," CVM_pal=",CVM_pvalue[j],"\n")

#Bereken Rayleigh pvalue folgens ingebeboude toets in CircStat
rayleigh_p[j]<-r.test(((data[(start+1):(endpoint-1)])-data[(start+1)])/(data[(endpoint-1)]-data[(start+1)])*2*pi)$p.value
if (rayleigh_p[j]=="NaN")
{rayleigh_p[j]<-1}



if (j>=r)
{
#rejection point calc

if (CVM_rejection == FALSE)
{
 if (max(CVM_pvalue[(j-r+1):j]) < alpha)
{
  CVM_rejection <- TRUE
  if (endpoint-1-((j-r+1)*g) >= ((length(data)/3)+1))
  {
   CVM_rejectionpnt[k]<-(data[endpoint-1-((j-r+1)*g)])/(1)
  }
  else
  {
  CVM_rejectionpnt[k]<-(data[endpoint-1-((j-r+1)*g)+(length(data)/3)])/(1)

  }

}
}

if (KS_rejection == FALSE)
{
 if (max(ks_p[(j-r+1):j]) < alpha)
{
  KS_rejection <- TRUE
    if (endpoint-1-((j-r+1)*g) >= ((length(data)/3)+1))
  {
   KS_rejectionpnt[k]<-(data[endpoint-1-((j-r+1)*g)])/(1)
  }
  else
  {
  KS_rejectionpnt[k]<-(data[endpoint-1-((j-r+1)*g)+(length(data)/3)])/(1)

  }


}
}

if (AD_rejection == FALSE)
{
 if (max(asq_p[(j-r+1):j]) < alpha)
{
  AD_rejection <- TRUE
  if (endpoint-1-((j-r+1)*g) >= ((length(data)/3)+1))
  {
   AD_rejectionpnt[k]<-(data[endpoint-1-((j-r+1)*g)])/(1)
  }
  else
  {
  AD_rejectionpnt[k]<-(data[endpoint-1-((j-r+1)*g)+(length(data)/3)])/(1)

  }
}
}

if (rayleigh_rejection == FALSE)
{
 if (max(rayleigh_p[(j-r+1):j]) < alpha)
{
  rayleigh_rejection <- TRUE
    if (endpoint-1-((j-r+1)*g) >= ((length(data)/3)+1))
  {
   rayleigh_rejectionpnt[k]<-(data[endpoint-1-((j-r+1)*g)])/(1)
  }
  else
  {
  rayleigh_rejectionpnt[k]<-(data[endpoint-1-((j-r+1)*g)+(length(data)/3)])/(1)

  }


}
}



}
# if j > rejection loop end



rejected<-(CVM_rejection & KS_rejection & AD_rejection & rayleigh_rejection)
j<-j+1
}
# while rejection loop end

if (!(is.na(which(c==0)[1])))
{
CVM_pvalue<-CVM_pvalue[1:((which(CVM_pvalue==0)[1])-1)]
}
if (!(is.na(which(d==0)[1])))
{
d<-d[1:((which(d==0)[1]-1))]
}
if (!(is.na(which(ks==0)[1])))
{
ks<-ks[1:((which(ks==0)[1]-1))]
}
if (!(is.na(which(ks_p==0)[1])))
{
ks_p<-ks_p[1:((which(ks_p==0)[1]-1))]
}
if (!(is.na(which(asq==0)[1])))
{
asq<-asq[1:((which(asq==0)[1]-1))]
}
if (!(is.na(which(asq_p==0)[1])))
{
asq_p<-asq_p[1:((which(asq_p==0)[1]-1))]
}
if (!(is.na(which(rayleigh_p==0)[1])))
{
rayleigh_p<-rayleigh_p[1:((which(rayleigh_p==0)[1]-1))]
}



}
#for k loop end (aantal minimum punte)
vec<-c(mean(CVM_rejectionpnt),mean(KS_rejectionpnt),mean(AD_rejectionpnt),mean(rayleigh_rejectionpnt))
sum<-matrix(c(mean(CVM_rejectionpnt),mean(KS_rejectionpnt),mean(AD_rejectionpnt),mean(rayleigh_rejectionpnt),median(vec)),ncol=5,dimnames=list(c("a-hat"),
c("Cramer von Mises","Kolmogorov-Smirnoff","Anderson-Darling","Rayleigh","MEDIAN")))

#CVM<-list(rejection=CVM_rejectionpnt,mean=mean(CVM_rejectionpnt),stdev=sqrt(var(CVM_rejectionpnt)),minimums=data#[close_to_min_index],grow=g,nr_reject=r,alpha=alpha,kernel="1-cos di")


#KS<-list(rejection=KS_rejectionpnt,mean=mean(KS_rejectionpnt),stdev=sqrt(var(KS_rejectionpnt)),minimums=data#[close_to_min_index],grow=g,nr_reject=r,alpha=alpha,kernel="1-cos di")

#AD<-list(rejection=AD_rejectionpnt,mean=mean(AD_rejectionpnt),stdev=sqrt(var(AD_rejectionpnt)),minimums=data#[close_to_min_index],grow=g,nr_reject=r,alpha=alpha,kernel="1-cos di")

#rayleigh<-list(rejection=rayleigh_rejectionpnt,mean=mean(rayleigh_rejectionpnt),stdev=sqrt(var(rayleigh_rejectionpnt)),minimums=data#[close_to_min_index],grow=g,nr_reject=r,alpha=alpha,kernel="1-cos di")

CVM<-list(rejection=CVM_rejectionpnt,mean=mean(CVM_rejectionpnt))

KS<-list(rejection=KS_rejectionpnt,mean=mean(KS_rejectionpnt))

AD<-list(rejection=AD_rejectionpnt,mean=mean(AD_rejectionpnt))

rayleigh<-list(rejection=rayleigh_rejectionpnt,mean=mean(rayleigh_rejectionpnt))

general<-list(call=cl,Minimums=data[close_to_min_index],alpha=alpha,grow=g,nr_reject=r,kernel_function="Epanechnikov")


comb<-list(summary=sum,General=general)


return(comb)

}
