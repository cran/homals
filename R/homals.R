sum.set<-function(g,n,p,y,set,active) {
z<-array(0.0,c(n,p))
for (j in set)
	if (active[j]){	    
	    gg<-g[,j]; ii<-which(!is.na(gg)); 
		z[ii,]<-z[ii,]+y[[j]][gg[ii],]}
return(z)
}

pca.update.y<-function(data,x,y,totsum,active,rank,level,nset){
for (j in 1:nset){
	if (active[j]){	    
		gg<-data[,j]; ycen<-compute.y(gg,x); d<-as.vector(table(gg)); ii<-which(!is.na(gg))
		yy<-restrict.y(d,ycen,rank[j],level[j])$y
		y[[j]]<-yy
		totsum[ii,]<-totsum[ii,]+yy[gg[ii],]
		}
	}
return(list(y=y,totsum=totsum))
}

cca.update.y<-function(data,x,y,totsum,active,rank,level,sets){
nobj<-dim(x)[1]; ndim<-dim(x)[2]; nset<-length(sets)
for (l in 1:nset) {
	indi<-sets[[l]]
	ss<-sum.set(data,nobj,ndim,y,indi,active)
	for (j in indi) {
		if (active[j]){	
			gg<-data[,j]; yy<-y[[j]]; ii<-which(!is.na(gg)); d<-as.vector(table(gg))
			ss<-ss-yy[gg[ii],]
			ycen<-compute.y(data[,j],x-ss)
			yy<-restrict.y(d,ycen,rank[j],level[j])$y
			ss<-ss+yy[gg[ii],]
			y[[j]]<-yy
			}
		}
	totsum<-totsum+ss
	}
return(list(y=y,totsum=totsum))
}

compute.y<-function(g,x) apply(x, 2, function(z) tapply(z,g,mean))

restrict.y<-function(d,y,r,level) {
if (sum(y^2) == 0) return(y)
switch(level,
"NO"=return(nominal.y(d,y,r)),
"OR"=return(ordinal.y(d,y,r)),
"NU"=return(numerical.y(d,y,r)),
"PO"=return(polynomial.y(d,y,r)))
}

nominal.y<-function(d,y,r) {
qq<-La.svd(sqrt(d)*y,r,r)
zz<-(1/sqrt(d))*qq$u
aa<-qq$d[1:r]*qq$vt
list(y=zz%*%aa,z=zz,a=aa)
}

ordinal.y<-function(d,y,r,itermax=100,eps=1e-6) {
qq<-La.svd(sqrt(d)*y,r,r)
a1=cbind(qq$vt[1,])
z1=cbind(orthogonal.polynomials(d,1:length(d),1))
if (r > 1) {
	a2<-cbind(qq$vt[2:r,])
	z2<-cbind(qq$u[,2:r])}
iter<-1; sold<-Inf
repeat{
	if (r == 1) ytilde<-y else ytilde<-y-(z2%*%t(a2))
	z1<-ytilde%*%a1/sum(a1*2)
	zo<-pava(z1,d)
	qo<-sum(d*zo^2)
	if (qo > 0) ao<-crossprod(ytilde,d*zo)/qo else ao<-a1
	so<-sum(d*(ytilde-zo%*%t(ao))^2)
	zp<-pava(-z1,d)
	qp<-sum(d*zp^2)
    if (qp > 0) ap<-crossprod(ytilde,d*zp)/qp else ap<-a1
	sp<-sum(d*(ytilde-zp%*%t(ap))^2)
	if (so < sp) {a1<-ao; z1<-zo; snew<-so}
		else {a1<-ap; z1<-zp; snew<-sp}
	if (r > 1) {
    	ytilde=y-(z1%*%t(a1))
    	qq<-La.svd(sqrt(d)*ytilde,r-1,r-1)
    	z2<-(1/sqrt(d))*qq$u
    	a2<-t(qq$d[1:(r-1)]*qq$vt)
    	snew<-sum(d*(ytilde-z2%*%t(a2))^2)}
	if ((iter == itermax) || ((sold - snew) < eps)) break()
		else {iter<-iter+1; sold<-snew}
}
if (r > 1) z<-cbind(z1,z2) else z<-cbind(z1)
z<-weighted.gram.schmidt(z,d)$pol
a<-crossprod(y,d*z)
list(yhat=z%*%t(a),z=z,a=a)
}

numerical.y<-function(d,y,r) {
z0<-orthogonal.polynomials(d,1:length(d),1)
a0<-as.vector(crossprod(z0*d,y))
if (r == 1)
	return(list(y=z0%o%a0,z=cbind(z0),a=rbind(a0)))
else {
	yy<-y-z0%o%a0
	qq<-La.svd(sqrt(d)*yy,r,r)
	zz<-cbind(z0,array((1/sqrt(d))*qq$u[,1:(r-1)],c(dim(y)[1],r-1)))
	aa<-rbind(a0,array((qq$d[1:r]*qq$vt)[1:(r-1),],c(r-1,dim(y)[2])))
	return(list(y=zz%*%aa,z=zz,a=aa))}
}

polynomial.y<-function(d,y,r) {
k<-length(d)
zz<-orthogonal.polynomials(d,1:k,min(r,k-1))
aa<-crossprod(zz,d*y)
list(y=zz%*%aa,z=zz,a=aa)
}

catplot<-function(name,vname,g,y,xx,yy,s,t) {
plot(y,type="l",xlim=xx,ylim=yy,main=paste("Category plot for",name,":",vname),
xlab=paste("dimension",s),ylab=paste("dimension",t))
text(y,levels(g))
abline(h=0)
abline(v=0)
}

trfplot<-function(name,vname,g,y) {
yy<-c(min(y),max(y))
p<-dim(y)[2]
r<-rainbow(p-1)
plot(y[,1],ylim=yy,type="l",
main=paste("Transformation plot for",name,":",vname),
xlab="original",ylab="transformed")
text(y[,1],levels(g))
if (p > 1)
	for (s in 2:p) {lines(y[,s],col=r[s-1]); text(y[,s],levels(g))}
abline(h=0)
}

starplot<-function(name,vname,g,y,x,s,t) {
plot(x,col="GREEN",pch=8,main=paste("Starplot for",name,":",vname),
xlab=paste("dimension",s),ylab=paste("dimension",t))
z<-compute.y(g, x)
points(z,type="n")
text(z,levels(g),col="RED")
for (i in 1:length(g)) lines(rbind(x[i,],z[g[i],]),col="BLUE")
}

lossplot<-function(name,vname,g,y,x,s,t) {
z<-compute.y(g, x); k<-dim(z)[1]
xx<-c(min(c(z[,1],y[,1])),max(c(z[,1],y[,1])))
yy<-c(min(c(z[,2],y[,2])),max(c(z[,2],y[,2])))
plot(y,type="n",main=paste("Lossplot for",name,":",vname),
xlab=paste("dimension",s),ylab=paste("dimension",t),xlim=xx,ylim=yy)
text(y,levels(g),col="RED");lines(y,col="RED")
points(z,type="n")
text(z,levels(g),col="GREEN");lines(z,col="GREEN")
for (i in 1:k) lines(rbind(y[i,],z[i,]),col="BLUE")
abline(h=0)
abline(v=0)
}

graphplot<-function(name,g,x,s,t) {
plot(x,col="GREEN",pch=8,main=paste("Graphplot for",name),
xlab=paste("dimension",s),ylab=paste("dimension",t))
for (j in 1:dim(g)[2]){
	y<-compute.y(g[,j], x)
	points(y,col="RED",pch=16)
	for (i in 1:dim(g)[1]) lines(rbind(x[i,],y[g[i,j],]))}
}

objplot<-function(name,rname,objlabel,offset,x,s,t) {
xx<-c(min(x[,1]),max(x[,1]))
yy<-c(min(x[,2]),max(x[,2]))
if (objlabel) {
plot(x,type="n",main=paste("Object score plot for",name),
xlab=paste("dimension",s),ylab=paste("dimension",t),xlim=offset*xx,ylim=offset*yy)
text(x,rname)}
else
plot(x,col="GREEN",pch=8,main=paste("Object score plot for",name),
xlab=paste("dimension",s),ylab=paste("dimension",t))
}

write.head<-function(name,vname,p,a,r,c,ofile) {
s<-NULL; sl<-35+10*p; for (i in 1:sl) s<-paste(s,"*",sep="")
cat("\n",formatC(s,format="s"),"\n",file=ofile,sep="")
cat(formatC(vname,format="s"),
if (a) "(Active)," else "(Passive),",
"Rank =", formatC(r,width=1),"and",
"Level is", 
switch(c,"NO"="nominal","OR"="ordinal","NU"="numerical","PO"="polynomial"),
"\n",file=ofile)
cat(formatC(s,format="s"),"\n",file=ofile)
}

write.y<-function(g,y,t,ofile) {
d<-as.vector(table(g));l<-levels(g)
s<-NULL; sl<-35+10*dim(y)[2]; for (i in 1:sl) s<-paste(s,"*",sep="")
u<-NULL; sl<-35+10*dim(y)[2]; for (i in 1:sl) u<-paste(u,"-",sep="")
switch(t,"Z"=cat("Lower Rank Quantifications\n",file=ofile),"C"=cat("Category Centroids\n",file=ofile),
"Y"=cat("Rank-restricted Category Quantifications\n",file=ofile))
cat(formatC(u,format="s"),"\n",file=ofile)
for (k in 1:length(d))
cat(formatC(l[k],width=10,format="s"),
formatC(d[k],width=5)," *** ",
formatC(y[k,],digits=6,width=10,format="f")," *** ",
formatC(sum(d[k]*y[k,]^2),digits=6,width=10,format="f"),"\n",sep="",file=ofile)
cat(formatC(u,format="s"),"\n",file=ofile)
cat(formatC(" ",format="s",width=10),
formatC(sum(d),width=5)," *** ",
formatC(d%*%(y^2),digits=6,width=10,format="f")," *** ",
formatC(sum(d%*%(y^2)),digits=6,width=10,format="f"),"\n",sep="",file=ofile)
cat(formatC(s,format="s"),"\n",file=ofile)
}

write.a<-function(a,ofile) {
s<-NULL; sl<-35+10*dim(a)[2]; for (i in 1:sl) s<-paste(s,"*",sep="")
u<-NULL; sl<-35+10*dim(a)[2]; for (i in 1:sl) u<-paste(u,"-",sep="")
cat("Category Loadings\n",file=ofile)
cat(formatC(u,format="s"),"\n",file=ofile)
for (k in 1:dim(a)[1])
cat(formatC(" ",format="s",width=10),formatC(k,width=5)," *** ",
formatC(a[k,],digits=6,width=10,format="f")," *** ",
formatC(sum(a[k,]^2),digits=6,width=10,format="f"),"\n",sep="",file=ofile)
cat(formatC(u,format="s"),"\n",file=ofile)
cat(formatC(" ",format="s",width=15)," *** ",
formatC(apply(a^2,2,sum),digits=6,width=10,format="f")," *** ",
formatC(sum(a^2),digits=6,width=10,format="f"),"\n",sep="",file=ofile)
cat(formatC(s,format="s"),"\n",file=ofile)
}

weighted.gram.schmidt<-function(x,w) {
ss<-NULL; 
for (j in 1:dim(x)[2]) {
if (j > 1) {xx<-x[,1:(j-1)]; x[,j]<-x[,j]-xx%*%(crossprod(xx,(w*x[,j])))}
s<-sqrt(sum(w*x[,j]^2)); ss<-c(ss,s); x[,j]<-x[,j]/s;}
list(pol=x,fac=ss)
}

orthogonal.polynomials<-function(w,x,p) {
z<-weighted.gram.schmidt(outer(x,0:p,"^"),w)$pol[,2:(p+1)]
}

center.x<-function(x,w) apply(x,2,function(z) z-weighted.mean(z,w))

norm.x<-function(x,w) {
qq<-La.svd(sqrt(w)*x); list(q=(1/sqrt(w))*(qq$u),r=qq$d)}

pava<-function(x,w=rep(1,length(x)),block=weighted.mean){
nblock<-n<-length(x); blocklist<-array(1:n,c(n,2)); blockvalues<-x; active<-1
repeat{
	if (!is.up.satisfied(blockvalues,active)) {
		blockmerge<-mergeBlockup(blocklist,blockvalues,x,w,active,block)
		blockvalues<-blockmerge$v; blocklist<-blockmerge$l
		nblock<-nblock-1
		while (!is.down.satisfied(blockvalues,active)) {
			blockmerge<-mergeBlockup(blocklist,blockvalues,x,w,active-1,block)
			blockvalues<-blockmerge$v; blocklist<-blockmerge$l; 
			nblock<-nblock-1; active<-active-1;
			}
		}
	else if (active == nblock) break() else active<-active+1
	}	
put.back(n,blocklist,blockvalues)
}

mergeBlockup<-function(blocklist,blockvalues,x,w,i,block){
n<-length(blockvalues); nn<-1:n; ii<-which(i+1!=nn)
blocklist[i,]<-c(blocklist[i,1],blocklist[i+1,2])
indi<-blocklist[i,1]:blocklist[i+1,2]
blockvalues[i]<-block(x[indi],w[indi])
blocklist<-blocklist[ii,]
if (length(ii) == 1) dim(blocklist)<-c(1,2)
blockvalues<-blockvalues[ii]
list(v=blockvalues,l=blocklist)
}

put.back<-function(n,blocklist,blockvalues){
x<-rep(0,n);nb<-length(blockvalues)
for (i in 1:nb) {
		x[blocklist[i,1]:blocklist[i,2]]<-blockvalues[i]}
return(x)
}

is.up.satisfied<-function(x,i) (i == length(x))||(x[i]<=x[i+1])

is.down.satisfied<-function(x,i) (i == 1)||(x[i-1]<=x[i])

tkhomals<-function(data)
	{
	require(tcltk,quietly=T); require(tkrplot,quietly=T)
	ncat <- length(txtlabel <- c("active","rank","level","starplot","catplot","trfplot","lossplot","set"))
	tkvar<-allvars<-""
	#need to have global window objects
	output.Win <<- graphplot.Win <<- objplot.Win  <<- voronoi.Win <<- FALSE
	starplots <- catplots <-trfplots <-lossplots <- FALSE
	active <- TRUE; sets <- 0; rank <- ndim <- 2; level <- "NO"
	plotgraph <- tclVar(0)
	plotobj <- tclVar(0)
	objscores <- tclVar(0)
	objlabel <- tclVar(0)
	voronoi <- tclVar(0)
	timer <- tclVar(0)
	saveme <- tclVar(0)
	tclvalue(ndim) <- ndim
		
	menu.Win <<- tktoplevel()
	tktitle(menu.Win) <- "Menu"
	main.frm <- tkframe(menu.Win)
	tframe <- tkframe(menu.Win)
	top.frm <- tkframe(main.frm, relief="groove", borderwidth=2)
	left.frm <- tkframe(main.frm)
	right.frm <- tkframe(top.frm)
	bottom.frm <- tkframe(main.frm)
	bfile<- tkmenubutton(tframe, text="File")
	bmenu <- tkmenu(bfile, tearoff=0, relief="raised")
	tkpack(tframe, fill="x", side = "top")
	tkpack(bfile, side="left")
	tkconfigure(bfile, menu=bmenu)
	tkadd (bmenu, "command", label="Quit", command=tkHom.quit)
	tkpack(tklabel(top.frm, text="HomalsTk 1.0"))
	graphplot.Label <- tklabel(right.frm, text="Graphplot")
	graphplot.Button <- tkcheckbutton(right.frm, variable=plotgraph)
	objplot.Label <- tklabel(right.frm, text="Objplot")
	objplot.Button <- tkcheckbutton(right.frm, variable=plotobj)
	objscores.Label <- tklabel(right.frm, text="Objscores")
	objscores.Button <- tkcheckbutton(right.frm, variable=objscores)
	objlabel.Label <- tklabel(right.frm, text="Objlabel")
	objlabel.Button <- tkcheckbutton(right.frm, variable=objlabel)
	voronoi.Label <- tklabel(right.frm, text="Voronoi")
	voronoi.Button <- tkcheckbutton(right.frm, variable=voronoi)
	timer.Label <- tklabel(right.frm, text="Timer")
	timer.Button <- tkcheckbutton(right.frm, variable=timer)
	save.Label <- tklabel(right.frm, text="Save")
	save.Button <- tkcheckbutton(right.frm, variable=saveme)
	ndim.Label <- tklabel(right.frm, text="Ndim")
	ndim.Entry <- tkentry(right.frm, textvariable=ndim, width=3)
	
	tkgrid(graphplot.Label, graphplot.Button, objplot.Label, objplot.Button, objscores.Label, objscores.Button, objlabel.Label, objlabel.Button)
	tkgrid(voronoi.Label, voronoi.Button, timer.Label, timer.Button, save.Label, save.Button, ndim.Label, ndim.Entry )

	# Create variable info section
	name <- deparse(substitute(data))
	message <- paste("Data selected: ", name)
	nvar<-dim(data)[2];
	vname<-attr(data,"names") 
	varlen <- max(nchar(vname))
	
	for(i in 1:nvar) {
		label <- vname[i]
		for (j in 1:ncat) {
			allvars[j + (ncat * (i-1))] <- tkvar[j] <- paste(i,txtlabel[j],sep=".")
			}
		tclvalue(tkvar[1]) <- 1
		tclvalue(tkvar[2]) <- 2
		tclvalue(tkvar[3]) <- "NO"
		tclvalue(tkvar[4]) <- 0
		tclvalue(tkvar[5]) <- 0
		tclvalue(tkvar[6]) <- 0
		tclvalue(tkvar[7]) <- 0
		tclvalue(tkvar[8]) <- i
		}	

	headerf <- tkframe(menu.Win, borderwidth=2)
	froptions <- tkframe(menu.Win, borderwidth=2, relief="groove")
	scr <- tkscrollbar(froptions, orient="vertical", command=function(...)tkyview(can,...))
	top.can <- tkcanvas(headerf)
	tkpack(top.can, fill="both", expand=TRUE, side="left")
	can <- tkcanvas(froptions)
	tkconfigure(can, yscrollcommand=function(...)tkset(scr,...))
	tkpack(scr, fill="y", side="right")
	tkpack(can, fill="both", expand=TRUE, side="left")
	tframe <- tkframe(menu.Win)
	varlen <- varlen * 12
	
	# Positioning checkboxes hack
	x <- c(2,varlen,(varlen+50),(varlen+90),(varlen+130),(varlen+180),(varlen+230),(varlen+275),(varlen+330),(varlen+335))
	y <- 2
	variable.Label <- tklabel(top.can, text="Variable"); active.Label <- tklabel(top.can, text="Active"); rank.Label <- tklabel(top.can, text="Rank")
	level.Label <- tklabel(top.can, text="Level"); starplot.Label <- tklabel(top.can, text="Starplot"); catplot.Label <- tklabel(top.can, text="Catplot")
	trfplot.Label <- tklabel(top.can, text="Trfplot"); lossplot.Label <- tklabel(top.can, text="Lossplot"); sets.Label <- tklabel(top.can, text="Sets")
	tkcreate(top.can, "window", x[1], y, anchor = "nw", window = variable.Label)
	tkcreate(top.can, "window", x[2], y, anchor = "nw", window = active.Label)
	tkcreate(top.can, "window", x[3], y, anchor = "nw", window = rank.Label)
	tkcreate(top.can, "window", x[4], y, anchor = "nw", window = level.Label)
	tkcreate(top.can, "window", x[5], y, anchor = "nw", window = starplot.Label)
	tkcreate(top.can, "window", x[6], y, anchor = "nw", window = catplot.Label)
	tkcreate(top.can, "window", x[7], y, anchor = "nw", window = trfplot.Label)
	tkcreate(top.can, "window", x[8], y, anchor = "nw", window = lossplot.Label)
	tkcreate(top.can, "window", x[9], y, anchor = "nw", window = sets.Label)
	
	for(i in 1:nvar) {
		label <- vname[i]
		for (j in 1:ncat) {
			allvars[j + (ncat * (i-1))] <<- tkvar[j] <- paste(i,txtlabel[j],sep=".")
			}
		variable.Header <- tklabel(can, text=label); active.Button <- tkcheckbutton(can, variable=tkvar[1]); rank.Entry <- tkentry(can, textvariable=tkvar[2], width=3)
		level.Entry <- tkentry(can, textvariable=tkvar[3], width=3); starplot.Button <- tkcheckbutton(can, variable=tkvar[4]); catplot.Button <- tkcheckbutton(can, variable=tkvar[5])
		trfplot.Button <- tkcheckbutton(can, variable=tkvar[6]); lossplot.Button <- tkcheckbutton(can, variable=tkvar[7]); sets.Entry <- tkentry(can, textvariable=tkvar[8], width=3)	
		tkcreate(can, "window", x[1], y, anchor = "nw", window = variable.Header)
		tkcreate(can, "window", x[2], y, anchor = "nw", window = active.Button)
		tkcreate(can, "window", x[3], y, anchor = "nw", window = rank.Entry)
		tkcreate(can, "window", x[4], y, anchor = "nw", window = level.Entry)
		tkcreate(can, "window", x[5], y, anchor = "nw", window = starplot.Button)
		tkcreate(can, "window", x[6], y, anchor = "nw", window = catplot.Button)
		tkcreate(can, "window", x[7], y, anchor = "nw", window = trfplot.Button)
		tkcreate(can, "window", x[8], y, anchor = "nw", window = lossplot.Button)
		tkcreate(can, "window", x[9], y, anchor = "nw", window = sets.Entry)
		y <- y + 20
		}
	
	rwidth = x[9] + 30
	top.rwidth = rwidth + 22
	coord <- sprintf("0 0 %s %s", as.character(y), as.character(y))

	fileinfo <- tklabel(main.frm, text=message, width=25, height=1)   
	bsubmit <- tkbutton(bottom.frm, text="SUBMIT", command=function()tkHom.submit(data,name,plotgraph,plotobj,objlabel,objscores,voronoi,timer,saveme,ndim,nvar,ncat,allvars,tkvar,active,rank,level,starplots,catplots,trfplots,lossplots,sets))
	tkpack(right.frm, fill="x", side = "left")
	tkpack(top.frm, fill="x", side = "top")
	tkpack(bottom.frm, fill="x", side="bottom")
	tkpack(main.frm, fileinfo)

	tkconfigure(top.can, width=top.rwidth, height=20)
	tkconfigure(can, width=rwidth, height=250,"-scrollregion",coord)
	tkpack(tframe, fill="x", side = "top")
	tkpack(headerf)
	tkpack(froptions)
	tkpack(bsubmit)
	tkfocus(menu.Win)
	}

tkHom.submit <- function(data,name,plotgraph,plotobj,objlabel,objscores,voronoi,timer,saveme,ndim,nvar,ncat,allvars,tkvar,active,rank,level,starplots,catplots,trfplots,lossplots,sets) {
    if (data!="") { 
    	graphplot.Image <- objplot.Image <- voronoi.Image <- ""
    	graphplot <- as.numeric(tclvalue(plotgraph)); objplot <- as.numeric(tclvalue(plotobj)); objlabel <- as.numeric(tclvalue(objlabel))
		objscores <- as.numeric(tclvalue(objscores)); voronoi <- as.numeric(tclvalue(voronoi)); timer <- as.numeric(tclvalue(timer)); saveme <- as.numeric(tclvalue(saveme))
		ndim <- as.numeric(tclvalue(ndim))
        for(i in 1:nvar) {
			arr1 <- as.numeric(tclvalue(allvars[1 + (ncat * (i-1))]))
			if (arr1==1) { active[i] <- TRUE } else { active[i] <- FALSE }
			rank[i] <- as.numeric(tclvalue(allvars[2 + (ncat * (i-1))]))
			level[i] <- tclvalue(allvars[3 + (ncat * (i-1))])
			arr2 <- as.numeric(tclvalue(allvars[4 + (ncat * (i-1))]))
			if (arr2==1) { starplots[i] <- TRUE } else { starplots[i] <- FALSE }
			arr3 <- as.numeric(tclvalue(allvars[5 + (ncat * (i-1))]))
			if (arr3==1) { catplots[i] <- TRUE } else { catplots[i] <- FALSE }
			arr4 <- as.numeric(tclvalue(allvars[6 + (ncat * (i-1))]))
			if (arr4==1) { trfplots[i] <- TRUE } else { trfplots[i] <- FALSE }
			arr5 <- as.numeric(tclvalue(allvars[7 + (ncat * (i-1))]))
			if (arr5==1) { lossplots[i] <- TRUE } else { lossplots[i] <- FALSE }
            sets[i] <- as.numeric(tclvalue(allvars[8 + (ncat * (i-1))]))
		 	}
		if (graphplot) {
			if (is.tkwin(graphplot.Win)) tkdestroy(graphplot.Win)
			graphplot.Win <<- tktoplevel()
			tktitle(graphplot.Win) <- "Homals Graphplot"
			graphplot.Image <- tkrplot(graphplot.Win, function() plot(0,0))
			tkgrid(graphplot.Image)
			}
		if (objplot) {
			if (is.tkwin(objplot.Win)) tkdestroy(objplot.Win)
			objplot.Win <<- tktoplevel()
			tktitle(objplot.Win) <- "Homals Objplot"
			objplot.Image <- tkrplot(objplot.Win, function() plot(0,0))
			tkgrid(objplot.Image)
            }
        if (voronoi) {
			if (is.tkwin(voronoi.Win)) tkdestroy(voronoi.Win)
			voronoi.Win <<- tktoplevel()
			tktitle(voronoi.Win) <- "Voronoi"
			voronoi.Image <- tkrplot(voronoi.Win, function() plot(0,0))
			tkgrid(voronoi.Image)
            }
        homals(data, sets = sets, ndim = ndim, active = active, rank = rank, level = level, starplots = starplots, catplots = catplots, trfplots = trfplots, lossplots = lossplots, graphplot = graphplot, objplot = objplot, objscores = objscores, objlabel = objlabel, voronoi = voronoi, timer = timer, save.me = saveme, tk = TRUE, img1 = graphplot.Image, img2 = objplot.Image, img3 = voronoi.Image, name = name)
		if (is.tkwin(output.Win)) tkdestroy(output.Win)
		output.Win <<- tktoplevel()
		tktitle(output.Win) <- "Results"
		outputFrame <- tkframe(output.Win)
		view <- tklistbox(outputFrame, width = 70, height = 40)
		scr <- tkscrollbar(outputFrame, orient = "vertical", command = function(...) tkyview(view,...))
		tkconfigure(view, yscrollcommand = function(...) tkset(scr, ...))
		tkpack(scr, side = "right", fill = "y")
		tkpack(view, side = "left", fill = "both", expand = TRUE)
		tkpack(outputFrame)
		filename <- paste(name,"out",sep=".")
		out <- scan(filename, what="", sep="\n")
		tkdelete(view, 0, "end")
    	tkinsert(view, "end", out)
		}
    }
    
tkHom.quit <- function() {
	tkdestroy(menu.Win)
	if (is.tkwin(output.Win)) tkdestroy(output.Win) 
	if (is.tkwin(graphplot.Win)) tkdestroy(graphplot.Win)
	if (is.tkwin(objplot.Win)) tkdestroy(objplot.Win) 
   	if (is.tkwin(voronoi.Win)) tkdestroy(voronoi.Win) 
	}

homals<-function(data,                 # data (in data-frame)
                  sets=0,              # list of vectors of indices
                  ndim=2,              # dimensionality (default 2)
                  active=T,            # which variables are active (single T means all)
                  rank=ndim,           # which quantification ranks (default all p)
                  level="NO",          # which quantification levels (default all nominal)            
                  starplots=F,         # which starplots (default none)
                  catplots=F,          # which category plots (default none)
                  trfplots=F,          # which transformation plots (default none)
                  lossplots=F,         # which loss plots (default none)
                  graphplot=F,         # graphplot (default no)
                  objplot=F,           # object score plot (default no)
                  objscores=F,         # object scores written to file (default no)
                  objlabel=F,          # object score plot labeled (default no)
                  offset=1.20,         # offset for labeled plots (default 1.20)
                  eps1=-Inf,           # iteration precision eigenvalues (default 1e-6)
                  eps2=1e-6,           # iteration precision eigenvectors (default 1e-6)
                  itermax=100,         # maximum number of iterations (default 100)
                  voronoi=F,           # voronoi diagram
                  save.me=F,           # do we return the results
                  demo=F,              # animated iteration demo (default no)
                  timer=F,             # time the steps of program (default no)
                  tk=F,				   # create tk homals output (default no)
                  img1,				   # tkrplot image1 placeholder
                  img2,                # tkrplot image2 placeholder
                  img3,				   # tkrplot image3 placeholder
                  name				   # dataframe name from tkhomals
                  )               
                  {            
if (timer) stime<-proc.time()
if (!tk) name<-deparse(substitute(data))
nobj<-dim(data)[1]; nvar<-dim(data)[2]; iter<-0; pops<-0; 
for (j in 1:nvar) data[,j]<-as.factor(data[,j]) 
if (length(sets) == 1) sets<-lapply(1:nvar,"c")
pca<-max(sapply(sets,length)) <= 1
nset<-length(sets)
if (pca) mis<-apply(data,1,function (x) sum(ifelse(is.na(x),0,1)))
	else mis<-apply(sapply(1:nset,
		function(s) apply(cbind(data[,sets[[s]]]),1,
		function (x) prod(ifelse(is.na(x),0,1)))),1,sum)
vname<-attr(data,"names") 
rname<-attr(data,"row.names")
if (ndim == 1) starplots<-catplots<-lossplots<-graphplot<-objplot<-F
if (length(active)==1) active<-rep(active,nvar)
if (length(starplots)==1) starplots<-rep(starplots,nvar)
if (length(catplots)==1) catplots<-rep(catplots,nvar)
if (length(trfplots)==1) trfplots<-rep(trfplots,nvar)
if (length(lossplots)==1) lossplots<-rep(lossplots,nvar)
do.pdf<-any(starplots,catplots,trfplots,lossplots,graphplot,objplot,voronoi)
if (length(rank)==1) rank<-rep(rank,nvar)
if (length(level)==1) level<-rep(level,nvar)
for (j in 1:nvar) {
    k<-length(levels(data[,j]))
	if (rank[j] > min(ndim,k-1)) rank[j]<-min(ndim,k-1)
	}
outfile<-file(paste(name,"out",sep="."),"w")
if (objscores)
	objfile<-file(paste(name,"obj",sep="."),"w")
x<-cbind(orthogonal.polynomials(mis,1:nobj,ndim))
x<-norm.x(center.x(x,mis),mis)$q
y<-lapply(1:nvar, function(j) compute.y(data[,j],x))
if (demo) plot(x, col="GREEN", pch=8)
if (timer) itime<-proc.time()
cat("Iterations:\n",file=outfile)
repeat {
	iter<-iter+1
	totsum<-array(0.0,dim(x))
	if (pca) up.y<-pca.update.y(data,x,y,totsum,active,rank,level,nset)
		else up.y<-cca.update.y(data,x,y,totsum,active,rank,level,sets)
	y<-up.y$y; totsum<-up.y$totsum
	qv<-norm.x(center.x((1/mis)*totsum,mis),mis)
	z<-qv$q;r<-qv$r;ops=sum(r);aps<-sum(La.svd(crossprod(x,mis*z),0,0)$d)/ndim	
	cat("Iteration: ",formatC(iter,digits=3,width=3)," Eigenvalues: ",
		formatC(r,digits=6,width=9,format="f"),
		" Gain: ",formatC(c(ops,aps),digits=6,width=9, format="f"),"\n",file=outfile)
	if (!tk) cat("Iteration: ",formatC(iter,digits=3,width=3)," Eigenvalues: ",
				formatC(r,digits=6,width=9,format="f")," Gain: ",
				formatC(c(ops,aps),digits=6,width=9, format="f"),"\n")
		if (demo) {
			points(z, col="RED",pch=8)
			for (i in 1:nobj) lines(rbind(x[i,],z[i,]))
		}
	if (((ops - pops) < eps1) || ((1.0 - aps) < eps2) || (iter == itermax)) break 
		else {x<-z; pops<-ops}}
if (timer) ctime<-proc.time()
if (do.pdf)
	pdf(file=paste(name,"pdf",sep="."),encoding="MacRoman")
xlim<-c(min(z[,1]),max(z[,1])); if (ndim > 1) ylim<-c(min(z[,2]),max(z[,2]))
if (graphplot) if (ndim > 1) for (s in 1:(ndim-1)) for (t in (s+1):ndim) {
	graphplot(name,data,z[,c(s,t)],s,t); 
	if (tk) tkrreplot(img1, function() graphplot(name,data,z[,c(s,t)],s,t))
	}
if (objplot) if (ndim > 1) for (s in 1:(ndim-1)) for (t in (s+1):ndim) {
	objplot(name,rname,objlabel,offset,z[,c(s,t)],s,t); 
	if (tk) tkrreplot(img2, function() objplot(name,rname,objlabel,offset,z[,c(s,t)],s,t))
	}
if (voronoi) {
	require("deldir",quietly=T)
	if (ndim > 1) for (s in 1:(ndim-1)) for (t in (s+1):ndim) {
		plot(deldir(z[,s],z[,t]),wlines="tess",wpoints="real",col=c(1,2,1),lty=c(1,1),pch=8)
		if (tk) tkrreplot(img3, function() plot(deldir(z[,s],z[,t]),wlines="tess",wpoints="real",col=c(1,2,1),lty=c(1,1),pch=8))
		}
	}
ylist<-alist<-clist<-ulist<-NULL
for (j in 1:nvar) {
    gg<-data[,j]; c<-compute.y(gg,z); d<-as.vector(table(gg))
    lst<-restrict.y(d,c,rank[j],level[j])
    y<-lst$y; a<-lst$a; u<-lst$z
    ylist<-c(ylist,list(y)); alist<-c(alist,list(a)); clist<-c(clist,list(c)); ulist<-c(ulist,list(u))
    write.head(name,vname[j],ndim,active[j],rank[j],level[j],outfile)
    write.y(data[,j],y,"Y",outfile)
    write.y(data[,j],c,"C",outfile)
    write.y(data[,j],u,"Z",outfile)
    write.a(t(a),outfile)
	if (starplots[j]) {
		if (ndim > 1) for (s in 1:(ndim-1)) for (t in (s+1):ndim) {
			starplot(name,vname[j],data[,j],y[,c(s,t)],z[,c(s,t)],s,t)}}
	if (lossplots[j]) {
		if (ndim > 1) for (s in 1:(ndim-1)) for (t in (s+1):ndim) {
			lossplot(name,vname[j],data[,j],y[,c(s,t)],z[,c(s,t)],s,t)}}
	if (catplots[j]) {
		if (ndim > 1) for (s in 1:(ndim-1)) for (t in (s+1):ndim) {
			catplot(name,vname[j],data[,j],y[,c(s,t)],xlim,ylim,s,t)}}
	if (trfplots[j]) trfplot(name,vname[j],data[,j],y)}
if (do.pdf) dev.off()
if (objscores){
	write(x,file=objfile,ncolumns=p)
	close(objfile)}
if (timer) {otime<-proc.time()
	cat("Input time:     ",formatC(itime-stime,digits=6,width=10,format="f"),"\n",file=outfile)
	cat("Compute time:   ",formatC(ctime-itime,digits=6,width=10,format="f"),"\n",file=outfile)
	cat("Iteration time: ",formatC((ctime-itime)/iter,digits=6,width=10,format="f"),"\n",file=outfile)
	cat("Output time:    ",formatC(otime-ctime,digits=6,width=10,format="f"),"\n",file=outfile)}
close(outfile)
if (save.me) return(x=z,y=ylist,c=clist,a=alist,z=ulist)
}

