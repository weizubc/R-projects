# rectangle(country,year,population,gdp,longitude,latitude,position1,position2,mark,count)

# position is for locating the position of labels by controlling adj in the text function.
# position1 determines the horizontal direction,0.5 denotes the middle, lower than 0.5 denotes "to the right",greater than 0.5 denotes "to the left"
# position2 determines the vertical direction,0.5 denotes the middle, lower than 0.5 denotes the top,greater than 0.5 denotes the bottom
# mark,count are used to solve the overlapping problem,just needed to keep in the function


rm(list=ls())

plot(c(-10,36),c(28,64),axes=FALSE, xlab=NA,ylab=NA,type='n')
par(mar=c(0,0,0,0)+0.1)
box(lty='solid')

move=function(mark,x,y,a,b)

{ time=length(mark)/4

for(i in 1:3)

{for ( i in 1:time )

{ x1=mark[i*4-3]

x2=mark[i*4-2]

y1=mark[i*4-1]

y2=mark[i*4]

newx1=x-a/2

newx2=x+a/2

newy1=y-b/2

newy2=y+b/2

distance=1


if ((x1<=(newx2+0.1))&(newx1<=(x2+0.1))&(y1<=(newy2+0.1))&(newy1<=(y2+0.1)))

{ if (( x1<=newx1)&(newx2<=x2)&(y1<=newy1)&(newy2<=y2))

{if ((x2-x1)<=(y2-y1)){x=x1-a/2-distance}else{y=y2+b/2+distance}}


  if ((y1>newy1)&( x1<=newx1)&(newx2<=x2)) {y=y1-b/2-distance}

  if ((newy2>y2)&(x1<=newx1)&(newx2<=x2)) {y=y2+b/2+distance}

  if (newx2>x2) {x=x2+a/2+distance}

  if (x1>newx1) {x=x1-a/2-distance}


}

  

} # 2 for ends

} # 1 for ends


combine=c(x,y)

combine

} # move function ends




rectangle=function(country,year,population,gdp,longitude,latitude,position1,position2,mark,count)

{ x=longitude
  y=latitude
  a=population/(10)
  b=gdp/(10)^4

 change=move(mark,x,y,a,b)
  x=change[1]
  y=change[2]
  
  mark[count]=x-a/2
  mark[count+1]=x+a/2
 mark[count+2]=y-b/2
 mark[count+3]=y+b/2
  

  lty=1
  lwd=1
  col=0

  

   if(is.na(year)){col=0} else
{  if(year==1952)
  {col='grey20'}

  if((year<=1995)&(year>=1973))
  {col='grey50'}

  if(year==2004)
  {col='grey90'}

  
  if(year>2007)
  {lty=2
   lwd=2}
}
   
  
 
  polygon(c((x-a/2),(x-a/2),(x+a/2),(x+a/2)),c((y-b/2),(y+b/2),(y+b/2),(y-b/2)),col=col,lty=lty,lwd=lwd)
  text(x,y,labels=country,cex=1,adj=c(position1,position2))

  mark 

} # rectangle function ends



mark=c(0,0,0,0)

count=1


library(foreign)

data=read.dta("eu.dta")

attach(data)


for(i in 1:length(country))

{ mark=rectangle(country[i],year[i],population[i],gdp[i],longitude[i],latitude[i],position1[i],position2[i],mark,count)

count=count+4

}





# for legend


mark=rectangle("GDP per capita",2007,30,10000,-3,63,position1=-0.3,position2=0.5,mark ,count)

count=count+4

text(-3,63,"Population",adj=c(0.5,2))

mark=rectangle("1950s",1952,10,10000,32.5,63,-0.415,0.5,mark,count)

count=count+4

mark=rectangle("1973-1995",1973,10,10000,32.5,61.5,-0.24,0.5,mark,count)

count=count+4

mark=rectangle("2004",2004,10,10000,32.5,60,-0.55,0.5,mark,count)

count=count+4

mark=rectangle("2007",2007,10,10000,32.5,58.5,-0.55,0.5,mark,count)

count=count+4

mark=rectangle("20??",2009,10,10000,32.5,57,-0.55,0.5,mark,count)

count=count+4




