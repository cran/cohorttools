#'
#' Plots cumulative incidence rates
#'
#' @param ftime failure time variable
#' @param fstatus variable with distinct codes for different causes of failure and also a distinct code for censored observations
#' @param group plots will be made for each group. If missing then treated as all one group
#' @param cencode value of fstatus variable which indicates the failure time is censored.
#' @param pop.length number of population sizes shown
#' @param ... additional parameters
#' @seealso \code{\link{survival}}  \code{\link{pyears}}
#' @note package cmprsk and ggplot2 are utilized
#' @keywords survival
#' @author Jari Haukka \email{jari.haukka@@helsinki.fi}
#' @return if missing group ggplot2 object or if group given named list of ggplot2 objects
#' @export
#' @examples
#' set.seed(2)
#' ss <- rexp(100)
#' gg <- factor(sample(1:3,100,replace=TRUE),1:3,c('a','b','c'))
#' cc <- sample(0:2,100,replace=TRUE)
#' print(plotcuminc(ftime=ss,fstatus=cc,cencode=0))
#' print(plotcuminc(ftime=ss,fstatus=cc,cencode=0,group=gg))
plotcuminc<-function(ftime,fstatus,cencode,pop.length=50,group,...){
  if(!missing(group)){
    lv.gr<-sort(unique(group))
    lv.ulos<-lapply(lv.gr,function(x){
      plotcuminc(ftime[group==x],fstatus[group==x],cencode=cencode,pop.length=pop.length)})
    names(lv.ulos)<-lv.gr
    return(lv.ulos)
  }
  lv.cinc<-cmprsk::cuminc(ftime=ftime,fstatus=fstatus,cencode=cencode)
  lv.df<-data.frame(lv.time=unlist(lapply(lv.cinc,function(x)x$time)),
                    lv.est=unlist(lapply(lv.cinc,function(x)x$est)))
  lv.R<-rep(names(lv.cinc),sapply(lv.cinc,function(x)length(x$time)))
  lv.df$lv.Type<-substr(lv.R,3,nchar(lv.R))
  lv.pop.time<-seq(min(lv.df$lv.time),max(lv.df$lv.time),length=pop.length)
  lv.pop<-sapply(lv.pop.time,function(x)sum(ftime>=x))
  lv.pop<-data.frame(lv.time=lv.pop.time,lv.pop=lv.pop)
  lv.scale<-max(lv.pop$lv.pop)/max(lv.df$lv.est)
  lv.time<-NA;lv.est<-NA;lv.Type<-NA # This for check
  lv.p1<-ggplot(lv.df,aes(x=lv.time))+geom_line(aes(y=lv.est,linetype=lv.Type))+
    xlab("Time") + ylab("Cumulative incidence")
  lv.p2<- lv.p1 + geom_line(data=lv.pop,aes(x=lv.time,y = lv.pop/(lv.scale)),size=1)
  lv.p3<- lv.p2 + scale_y_continuous(sec.axis = sec_axis(~.*(lv.scale), name = "Population size"))
  lv.p3
}
#'
#' Function makes rate table with confidence intervals for crude incidences (rates)
#'
#' @param formula  where Surv object is on lhs and marginal variable(s)
#' on rhs. Marginal variables should usually be factors
#' @param data data.frame to be used
#' @param alpha confidence level, default is 0.05
#' @param add.RR should rate ratio (RR) be added
#' @param lowest.N lowest frequency to be shown
#' @param ... additional parameter for function survival::pyears
#' @seealso \code{\link{survival}}  \code{\link{pyears}}
#' @note packages survival is utilized.
#' Frequencies lower than lowest.N replaced by 999999
#' Person-years scaled by default with 365.25
#' @keywords survival
#' @author Jari Haukka \email{jari.haukka@@helsinki.fi}
#' @return table with columns named after marginal variables and n, event, incidence, se, exact.lower95ci
#' and  exact.upper95ci variables
#' @export
#' @examples
#' library(survival)
#' tmp.lt1<-mkratetable(Surv(time,status)~ sex,data=lung)
#' tmp.lt2<-mkratetable(Surv(time,status)~ sex+ph.ecog,data=lung,add.RR=TRUE,lowest.N=10)
mkratetable<-function(formula,data,alpha=0.05,add.RR=FALSE,lowest.N=0,...)
{
  # Exact CI for Poisson distribution (originally from epitools package)
  lv.pois.exact<-function (x, pt = 1, conf.level = 0.95)
  {
    xc <- cbind(x, conf.level, pt)
    pt2 <- xc[, 3]
    results <- matrix(NA, nrow(xc), 6)
    f1 <- function(x, ans, alpha = alp) {
      ppois(x, ans) - alpha/2
    }
    f2 <- function(x, ans, alpha = alp) 1 - ppois(x, ans) + dpois(x,
                                                                  ans) - alpha/2
    for (i in 1:nrow(xc)) {
      alp <- 1 - xc[i, 2]
      interval <- c(0, xc[i, 1] * 5 + 4)
      uci <- uniroot(f1, interval = interval, x = xc[i, 1])$root/pt2[i]
      if (xc[i, 1] == 0) {
        lci <- 0
      }
      else {
        lci <- uniroot(f2, interval = interval, x = xc[i,
                                                       1])$root/pt2[i]
      }
      results[i, ] <- c(xc[i, 1], pt2[i], xc[i, 1]/pt2[i],
                        lci, uci, xc[i, 2])
    }
    coln <- c("x", "pt", "rate", "lower", "upper", "conf.level")
    colnames(results) <- coln
    data.frame(results)
  }

  lv.kutsu<-deparse(substitute(formula))
  lv.data<-deparse(substitute(data))

  lv.00<-as.character(formula)
  lv.0<-strsplit(x = lv.00[3],split = "\\+")[[1]]
  if(length(lv.0)>1){
    lv.fun.1<-function(x){
      lv.1<-data.frame(Var=c(colnames(x)[1],rep("",nrow(x)-1)),
                       Group=as.character(x[,1]),x[,-1],
                       stringsAsFactors =FALSE,row.names = NULL)
      lv.1
    }
    lv.ulos<-lapply(lv.0,function(x){
      mkratetable(formula(paste0(paste(rev(lv.00[-3]),collapse = ""),x)),data=data,add.RR=add.RR,lowest.N=lowest.N,...)
    })
    lv.len<-sapply(lv.ulos,nrow)
    lv.ulos<-do.call(rbind,lapply(lv.ulos,lv.fun.1))

    # Add attributes
    attributes(lv.ulos)$ratetable.len<-lv.len
    attributes(lv.ulos)$ratetable.call<-lv.kutsu
    attributes(lv.ulos)$ratetable.data<-lv.data
    attributes(lv.ulos)$ratetable.time<-date()
    # class(lv.ulos)<-"ratetable"
    return(lv.ulos)
  }
  lv.1<-survival::pyears(formula,data=data,data.frame=TRUE,...)$data
  lv.2<-lv.pois.exact(x=lv.1$event,pt=lv.1$pyears,conf.level =1-alpha)[,-6]
  lv.2<-cbind(lv.1,lv.2[,-c(1,2)])
  lv.3<-rev(names(lv.2))
  lv.3[1:2]<-paste(c("high","low"),sep=".")
  names(lv.2)<-rev(lv.3)
  lv.2<-lv.2[,is.na(match(names(lv.2),"n"))]
  colnames(lv.2)<-gsub("incidence","inci.",colnames(lv.2))
  if(add.RR)
  {
    lv.m<-glm(lv.1$event~offset(log(pyears))+factor(seq(nrow(lv.1))),data=lv.1,family=poisson)
    lv.RR<-rbind(c(1,1,1),Epi::ci.exp(lv.m,subset=-1))
    colnames(lv.RR)<-gsub("%","p",colnames(lv.RR))
    colnames(lv.RR)[1]<-"RR"
    lv.2<-cbind(lv.2,lv.RR)
  }

  # Event cell frequencies, if too low replace it with high number
  if(min(lv.2$event)<lowest.N) warning(paste("lowest cell frequency under",
                                             lowest.N,", replaced by 999999"))
  lv.2$event<-ifelse(lv.2$event<lowest.N,9999999,lv.2$event)
  # Add attributes
  attributes(lv.2)$ratetable.call<-lv.kutsu
  attributes(lv.2)$ratetable.data<-lv.data
  attributes(lv.2)$ratetable.time<-date()
  rownames(lv.2)<-NULL
  # class(lv.2)<-"ratetable"
  lv.2
}


#' Boxes plot summarizing Lexis object
#'
#' Creates boxes graph describing Lexis
#'
#' @param x Lexis object
#' @param layout Graphviz layout "circo", "dot", "twopi" or, "neato".
#' It determines general layout of graph.
#' @param prop.penwidth use line width relative to incidence. If TRUE linewidths of
#' showing transition rates beween states are relative to log of rate.
#' @param scale.Y scale for incidence. Scale factor rates, default is 1.
#' @param rankdir for graph, default is TB. NOTE! this works best with layout "dot"
#' @param node.attr general node attributers.
#' Attributes like shape, color, fillcolor, etc. for nodes.
#' Consult Graphviz documentation for details
#' \url{https://www.graphviz.org/doc/info/attrs.html}.
#' @param edge.attr general edge (line) attributers.
#' Attributes like  color, arrowhead, fontcolor  etc. for edges.
#' Consult Graphviz documentation for details
#' \url{https://www.graphviz.org/doc/info/attrs.html}
#' @param show.loop ,  should loop (staying in same state be shown), default FALSE
#' @param show.persons , should number of persons be shown (entry->exit), default FALSE
#' @param fontsizeN font size for nodes
#' @param fontsizeL font size for edges
#' @param show.gr  should graph be shown. If TRUE,
#' function DiagrammeR::grViz is used to show graph.
#' @author Jari Haukka <jari.haukka@@helsinki.fi>
#' @return Character vector containing Graphviz script. This may used
#' to create graph by  DiagrammeR::grViz function.
#' @seealso \code{\link{grViz}}
#' @export
#' @examples
#' library(DiagrammeR)
#' library(survival)
#' library(Epi)
#' library(mstate)
#' data(ebmt3)
#' bmt <- Lexis(exit = list(tft = rfstime/365.25),
#'              exit.status = factor(rfsstat, labels = c("Tx", "RD")),
#'                           data = ebmt3)
#' bmtr <- cutLexis(bmt, cut = bmt$prtime/365.25, precursor.states = "Tx",
#'                                            new.state = "PR")
#'
#' summary(bmtr)
#' kk<-boxesLx(bmtr)
#' \dontrun{
#' # Graph to file
#' gv2image(kk, file="k1", type="pdf")
#' }
#' boxesLx(bmtr,layout="dot",rankdir = "LR",show.loop = FALSE,show.persons = TRUE)
#' boxesLx(bmtr,node.attr='shape=hexagon color=navy style=filled fillcolor=lightblue',
#' edge.attr = ' color=steelblue arrowhead=vee fontcolor="#8801d7" ',
#' layout="circo",prop.penwidth=TRUE)
boxesLx<-function(x,
                  layout="circo",
                  prop.penwidth=FALSE,
                  scale.Y = 1,
                  rankdir="TB",
                  node.attr="shape=box",
                  edge.attr="minlen=1",
                  show.loop=FALSE,
                  show.persons=FALSE,
                  fontsizeN=14,fontsizeL=8,
                  show.gr=TRUE)
{
  # layout<-c("dot","neato","twopi","circo")
  if(!("Lexis"%in%class(x))) stop("x must be Lexis object")
  lv.1<-summary(x)
  lv.1$Transitions[,"Risk time:"]<-lv.1$Transitions[,"Risk time:"]/scale.Y
  lv.len<-nrow(lv.1$Transitions)-1
  lv.N<-lv.1$Transitions[1:lv.len,!(colnames(lv.1$Transitions)%in%c(" Records:", " Events:","Risk time:", " Persons:"))]

  lv.N.Persons<-lv.1$Transitions[1:lv.len," Persons:" ]
  lv.rate<-lv.N/lv.1$Transitions[1:lv.len,"Risk time:"]
  lv.nd.names<-paste0('"',colnames(lv.N),'"')
  lv.ent.N<-table(status(x, at="entry", by.id=TRUE))
  lv.ex.N<-table(status(x, at="exit" , by.id=TRUE))

  lv.RT<-lv.1$Transitions[1:lv.len,"Risk time:"]
  lv.RT<-lv.RT[match(names(lv.ent.N),names(lv.RT))]
  lv.RT<-paste(round(lv.RT,2))
  lv.RT<-ifelse(grepl("NA",lv.RT),"",lv.RT)

  lv.head.gv<-c("digraph G {\n" ,
                "layout=",layout,";\n",
                "rankdir=",rankdir,";\n",
                "node [",node.attr,"];\n",
                "edge [",edge.attr,"];\n")
  lv.main.gv<-NULL

  if(show.persons)  lv.main.gv<-c(lv.main.gv,paste0('"',
                                                    names(lv.ent.N),
                                                    '"',' [label="',
                                                    paste(names(lv.ent.N),'\\n',
                                                          lv.RT,'\\n',
                                                          lv.ent.N,'&rarr;',
                                                          lv.ex.N),
                                                    '" , fontsize=',fontsizeN,'];\n'))
  else{
    lv.main.gv<-c(lv.main.gv,paste0('"',
                                    names(lv.ent.N),
                                    '"',' [label="',
                                    names(lv.ent.N),'\\n',lv.RT,
                                    '" , fontsize=',fontsizeN,'];\n'))

  }
  lv.pen1<-log(lv.rate/min(lv.rate[lv.rate>0]))
  lv.pen2<-3*(lv.pen1/max(lv.pen1))+0.2

  for(lv.from in seq(nrow(lv.N))) for(lv.to  in seq(ncol(lv.N))){
    if(lv.N[lv.from,lv.to]>0){
      lv.pen3<-ifelse(prop.penwidth,paste0(" penwidth=",lv.pen2[lv.from,lv.to])," ")
      if((!show.loop & (rownames(lv.N)[lv.from]==colnames(lv.N)[lv.to]))){lv.main.gv<-lv.main.gv}
      else {
        lv.main.gv<-c(lv.main.gv,paste0('"',rownames(lv.N)[lv.from],'"'),"->",
                      paste0('"',colnames(lv.N)[lv.to],'"'),
                      paste0('[label="',lv.N[lv.from,lv.to],
                             '\\n(',
                             round(lv.rate[lv.from,lv.to],2),
                             ')",fontsize=',fontsizeL," ",lv.pen3,']'),";\n")
      }
    }
  }
  lv.gv<-c(lv.head.gv,lv.main.gv,"}\n")
  if(show.gr){
    print(DiagrammeR::grViz(lv.gv,engine=layout))
  }
  invisible(lv.gv)
}

#'
#' Function makes image from graphviz code
#'
#' @param gv  character string containing graphviz code
#' @param file file name for image, character string
#' @param type type of ('pdf', 'png', 'ps', 'raw','svg','webp') as character string
#' @param engine grViz engine, defaults is 'dot'
#' @param ... parameters for rsvg_
#' @author Jari Haukka \email{jari.haukka@@helsinki.fi}
#' @return Invisible name of file created.
#' @export
gv2image<-function(gv,file="gv",type="png",engine="dot",...){
  lv.svg<-DiagrammeRsvg::export_svg(DiagrammeR::grViz(gv,engine=engine))
  lv.file<-paste0(file,".",type)
  lv.txt<-paste0('rsvg::rsvg_',type,'(svg=charToRaw(lv.svg),file=lv.file,...)')
  eval(parse(text=lv.txt))
  invisible(lv.file)
}
#'
#' Estimates hazard function using Poisson model
#' @param formula formula with Surv in LHS, NOTE! only one variable in RHS
#' @param data data used by formula
#' @param time  time variables
#' @param status status indicator Lowest value used as sensoring.
#' If only one unique value detected, all are assumed events
#' @param breaks time is splitted with these values
#' @param knots knots for natural splines used in estimation of hazard function
#' @param time.eval in which time points hazard function is evaluate.
#' @param alpha significance level for confidence intervals
#' @param use.GAM logical determining if generalized additive model (GAM) is used
#' @param print.GAM.summary logical determining if summary of GAM is printed
#' @param ... parameters for glm
#' @author Jari Haukka \email{jari.haukka@@helsinki.fi}
#' @return Returns data frame with time and hazard function values with attribute 'estim.hazard.param'
#' containing estimation parameters (breaks and knots)
#' @export
#' @examples
#' library(survival)
#' tmp.hz<-estim.hazard(time=lung$time,status=lung$status)
#' head(tmp.hz,2)
#' attributes(tmp.hz)$estim.hazard.param # estimation parameters
#' tmp.hz2<-estim.hazard(formula=Surv(time,status)~sex,data=lung)
#' head(tmp.hz2,2)
estim.hazard<-function(formula,data,time,status,breaks,knots,
                           time.eval=breaks,alpha=0.05,use.GAM=FALSE,print.GAM.summary=FALSE,...){

  # require(splines)
  #------------------------------
  # Function to calculate hazard estimate using Poisson regression
  #------------------------------
  lv.fun1<-function(time,status,breaks=breaks,knots=knots,time.eval=breaks,alpha=0.05){
    # Check if there is only one unique status value
    if(length(unique(status))==1)lv.status<-1
    else {
      lv.min.status<-min(status)
      lv.status<-ifelse(status==lv.min.status,0,1)
    }
    # Make Lexis
    lv.Lx <- Epi::Lexis( duration=time,
                         entry=list(fu.time=0),
                         exit.status=lv.status,
                         entry.status = 0)
    # Split Lexis, default is 100 breaks
    if(missing(breaks)) breaks<-with(lv.Lx,c(0,seq(from=min(lex.dur),to=max(lex.dur),
                                                   length=100)))
    lv.Lx.s <- Epi::splitLexis( lv.Lx, "fu.time", breaks=breaks )
    lv.Lx.s.agg<-stats::aggregate(list(lex.dur=lv.Lx.s$lex.dur,lex.Xst=lv.Lx.s$lex.Xst),list(fu.time=lv.Lx.s$fu.time),sum)

    # Poisson model with splines
    # Define knots, if missing
    if(missing(knots)) knots<-with(lv.Lx,seq(from=min(lex.dur),to=max(lex.dur),
                                             length=7))[-c(1,7)]
    lv.Xmat<-splines::ns(lv.Lx.s.agg$fu.time,knots=knots,intercept = TRUE)
    lv.m<-glm(lex.Xst~lv.Xmat+offset(log(lex.dur))-1,family = poisson,data=lv.Lx.s.agg,maxit=50,...)

    if(missing(time.eval))time.eval<-breaks
    lv.Xmat.eval<-splines::ns(time.eval,knots=knots,intercept = TRUE)
    lv.lambda <- Epi::ci.lin( lv.m, ctr.mat=lv.Xmat.eval, Exp=TRUE,alpha=alpha )
    lv.ulos<-data.frame(time.eval=time.eval,lv.lambda)
    names(lv.ulos)[6:8]<-c("haz","haz.lo","haz.hi")
    lv.ulos<-lv.ulos[,-c(4,5)]
    attributes(lv.ulos)$estim.hazard.param<-list(breaks=breaks,knots=knots)
    lv.ulos
  }

  #------------------------------
  # If missing formula, old style
  #------------------------------
  if(missing(formula)){
    return(lv.fun1(time=time,status=status,breaks=breaks,knots=knots,time.eval=breaks,alpha=alpha))
  }

  #------------------------------
  # With formula, new style
  #------------------------------

  # Create survival object and make working data frame
  lv.Surv<-eval(parse(text=paste0("with(data,",as.character(formula)[2],")")))

  # If no grouping in formula RHS and no GAM
  if(as.character(formula)[3]=="1"&!use.GAM){
    lv.dt<-data.frame(time=lv.Surv[,1],status=lv.Surv[,2])
    return(lv.fun1(time=lv.dt$time,status=lv.dt$status,breaks=breaks,knots=knots,time.eval=breaks,alpha=alpha))
  }

  # Make data.frame to use in analyses
  if(as.character(formula)[3]=="1"&use.GAM){
    lv.dt<-data.frame(time=lv.Surv[,1],status=lv.Surv[,2])
  }
  else{
    lv.cvar<-unlist(strsplit(as.character(formula)[3],fixed = TRUE,split = "+"))[1]
    lv.var2<-eval(parse(text=paste0("with(data,",lv.cvar,")")))
    lv.dt<-data.frame(time=lv.Surv[,1],status=lv.Surv[,2],Indeksi=lv.var2)
  }
  #------------------------------
  # Using GAM
  #------------------------------
  if(use.GAM){
    # Make Lexis
    lv.Lx <- with(lv.dt,Epi::Lexis( duration=time,
                                    entry=list(fu.time=0),
                                    exit.status=status,
                                    entry.status = 0,data=lv.dt))

    # Split Lexis, default is 100 breaks
    if(missing(breaks)) breaks<-with(lv.Lx,c(0,seq(from=min(lex.dur),to=max(lex.dur),length=100)))
    lv.Lx.s <- Epi::splitLexis( lv.Lx, "fu.time", breaks=breaks )

    if(as.character(formula)[3]=="1"){
      lv.Lx.s.agg<-with(lv.Lx,stats::aggregate(list(lex.dur=lv.Lx.s$lex.dur,lex.Xst=lv.Lx.s$lex.Xst),list(fu.time=lv.Lx.s$fu.time),sum))
      lv.m<-mgcv::gam(lex.Xst ~ s(fu.time)+offset(log(lex.dur)),data=lv.Lx.s.agg,family=poisson)
      lv.dt.new<-data.frame(fu.time=lv.Lx.s.agg$fu.time,lex.dur=1)
      lv.apu.p<-mgcv::predict.gam(lv.m,type="link",se.fit = TRUE,newdata=lv.dt.new)
      if(print.GAM.summary){print(summary(lv.m))}
      lv.apu.ulos<-with(lv.Lx.s.agg,
                        data.frame(time.eval=fu.time,
                                   haz=exp(lv.apu.p$fit),
                                   haz.lo=exp(lv.apu.p$fit+lv.apu.p$se.fit*qnorm(alpha/2)),
                                   haz.hi=exp(lv.apu.p$fit+lv.apu.p$se.fit*qnorm(1-alpha/2))
                        ))
      return(lv.apu.ulos)
    }
    else{

      lv.Lx.s.agg<-with(lv.Lx,stats::aggregate(list(lex.dur=lv.Lx.s$lex.dur,lex.Xst=lv.Lx.s$lex.Xst),
                                               list(fu.time=lv.Lx.s$fu.time,Indeksi=lv.Lx.s$Indeksi),sum))
      lv.m<-by(data = lv.Lx.s.agg,INDICES = lv.Lx.s.agg$Indeksi,
               function(x){
                 lv.m<-mgcv::gam(lex.Xst ~ s(fu.time)+offset(log(lex.dur)),data=x,family=poisson)
                 lv.dt.new<-data.frame(fu.time=x$fu.time,lex.dur=1)
                 lv.apu.p<-mgcv::predict.gam(lv.m,type="link",se.fit = TRUE,newdata=lv.dt.new)

                 if(print.GAM.summary){
                   cat("\n--------------------\nGAM summary for ",
                       lv.cvar," with level ",x$Indeksi[1],
                       "\n--------------------\n")
                   print(summary(lv.m))
                   cat("\n--------------------\n")
                 }
                 lv.apu.ulos<-with(x,
                                   data.frame(time.eval=fu.time,
                                              haz=exp(lv.apu.p$fit),
                                              haz.lo=exp(lv.apu.p$fit+lv.apu.p$se.fit*qnorm(alpha/2)),
                                              haz.hi=exp(lv.apu.p$fit+lv.apu.p$se.fit*qnorm(1-alpha/2))
                                   ))
                 return(lv.apu.ulos)
               })
      lv.apu.ulos<-data.frame(do.call("rbind",lv.m),Indeksi=rep(names(lv.m),sapply(lv.m,nrow)))
      names(lv.apu.ulos)[ncol(lv.apu.ulos)]<-lv.cvar
      return(lv.apu.ulos)
    }
  }

  #------------------------------
  # No GAM used
  #------------------------------

  if (missing(breaks)) breaks <- with(lv.dt, c(0, seq(from = min(time), to = max(time),length = 100)))
  if (missing(knots)) knots <- with(lv.dt, seq(from = min(time), to = max(time),length = 7))[-c(1, 7)]
  lv.2<-by(data=lv.dt,INDICES =lv.dt$Indeksi,
           function(x){
             lv.fun1(time=x$time,status=x$status,breaks=breaks,knots=knots,time.eval=breaks,alpha=alpha)
           })

  lv.ulos<-do.call("rbind",lv.2)
  lv.ulos$lv.NimiPitka<-factor(rep(names(lv.2),sapply(lv.2,nrow)))
  names(lv.ulos)[ncol(lv.ulos)]<-lv.cvar
  lv.ulos
}
#'
#' Function makes flowchart in graphviz
#'
#' @param N  Population sizes
#' @param text.M Text for exclusions, length one less than N
#' @param text.P Text for main boxes, must be same length with N
#' @param type flowchart type (1 or 2)
#' @author Jari Haukka \email{jari.haukka@@helsinki.fi}
#' @return Character string, graphviz language
#' @export
#' @export
#' @examples
#' DiagrammeR::grViz(mkflowchart(N=c(743,32,20),
#' text.M=c("Excluded","Excluded \n other with reasons"),
#' text.P=c("Studies","Relevant studies","Included in final review"),type=1))
mkflowchart<-function(N,text.M,text.P,type=1){
  lv.alku<-"digraph graphname { \n rankdir=TB \n style=filled; \n color=white; \n len=0.001; \n "
  lv.loppu<-" \n }"
  lv.ero<- (-diff(N))
  lv.len<-length(lv.ero)

  lv.4<-paste0("P",seq(lv.len+1),' [label= "',text.P,'\n N=',prettyNum(N,big.mark = ","),'",fontname=Calibri,fontsize=12,shape=box] \n')
  lv.6<-paste0("M",seq(lv.len),'[label="N=',prettyNum(lv.ero,big.mark = ","),'",fontname=Calibri,fontsize=10,shape=box] \n')

  if(type==1){
    lv.rank<-paste0("{ rank=same; V",seq(lv.len),"; M",seq(lv.len)," }\n")
    lv.1<-paste0("P",seq(lv.len),"->V",seq(lv.len),' [arrowhead="none"] \n')
    lv.2<-paste0("V",seq(lv.len),"->M",seq(lv.len),'[label="',text.M,'",fontname=Calibri,fontsize=10] \n')
    lv.3<-paste0("V",seq(lv.len),' [label="",shape="point",width=0.001] \n')
    lv.5<-paste0("V",seq(lv.len),"->P",seq(lv.len)+1,' \n')
    return(c(lv.alku,lv.1,lv.2,lv.3,lv.4,lv.5,lv.6,lv.rank,lv.loppu))
  }
  else{
    lv.rank<-paste0("{ rank=same; P",seq(lv.len),"; M",seq(lv.len)," }\n")
    lv.2<-paste0("P",seq(lv.len),"->M",seq(lv.len),'[label="',text.M,'",fontname=Calibri,fontsize=10,arrowsize=.7] \n')
    lv.5<-paste0("P",seq(lv.len),"->P",seq(lv.len)+1,' \n')
    return(c(lv.alku,lv.2,lv.4,lv.5,lv.6,lv.rank,lv.loppu))
  }
}
#'
#' Function makes plot(s) from ratetable
#'
#' @param rt  Rate table produced by function mkratetable
#' @param RR Boolean, if TRUE rate ratios plotted
#' @return ggplot object, or list if multiple variables in rate table
#' @export
#' @examples
#' library(ggplot2)
#' library(survival)
#' tmp.lt1<-mkratetable(Surv(time,status)~ ph.ecog,data=lung,add.RR = FALSE)
#' plotratetable(tmp.lt1)
#' tmp.lt2<-mkratetable(Surv(time,status)~ sex+ph.ecog+cut(age,4),data=lung,add.RR=TRUE,lowest.N=1)
#' plotratetable(tmp.lt2,TRUE)
plotratetable<-function(rt,RR=FALSE){
  if(RR & (ncol(rt)<9)) stop("No RR in ratetable")

  if(names(rt)[1]=="Var"){
    # lv.1<-rep(seq(attributes(rt)$ratetable.len),attributes(rt)$ratetable.len)
    lv.0<-unique(rt$Var)
    lv.1<-rep(lv.0[lv.0!=""],attributes(rt)$ratetable.len)
    lv.2<-by(rt,lv.1,function(x){
      if(RR) {
        names(x)[9:10]<-c("RR.low","RR.high")
        eval(parse(text=paste("lv.p<-ggplot(x,aes(x=factor(Group),y=RR))")))
        with(x,lv.p+geom_point()+geom_pointrange(aes(ymin=RR.low,ymax=RR.high))+xlab(x$Var[1])+coord_flip()+
          geom_hline(yintercept = 1, linetype="dotted"))
      }
      else{
      eval(parse(text=paste("lv.p<-ggplot(x,aes(x=factor(Group),y=rate))")))
      with(x,lv.p+geom_point()+geom_pointrange(aes(ymin=low,ymax=high))+xlab(x$Var[1])+coord_flip())
      }
    })
    return(lv.2)
  }
  if(RR) {
    names(rt)[8:9]<-c("RR.low","RR.high")
    eval(parse(text=paste("lv.p<-ggplot(rt,aes(x=factor(",names(rt)[1],"),y=RR))")))
    with(rt,lv.p+geom_point()+geom_pointrange(aes(ymin=RR.low,ymax=RR.high))+xlab(names(rt)[1])+coord_flip()+
    geom_hline(yintercept = 1, linetype="dotted"))
  }
  else{
    eval(parse(text=paste("lv.p<-ggplot(rt,aes(x=factor(",names(rt)[1],"),y=rate))")))
    with(rt,lv.p+geom_point()+geom_pointrange(aes(ymin=low,ymax=high))+xlab(names(rt)[1])+coord_flip())
  }
}
