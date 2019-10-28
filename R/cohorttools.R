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
#' @note packages survival and epitools are utilized.
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
#' tmp.lt2<-mkratetable(Surv(time,status)~ sex,data=lung,add.RR=TRUE,lowest.N=60)
mkratetable<-function(formula,data,alpha=0.05,add.RR=FALSE,lowest.N=0,...)
{
  lv.1<-survival::pyears(formula,data=data,data.frame=TRUE,...)$data
  lv.2<-epitools::pois.exact(x=lv.1$event,pt=lv.1$pyears,conf.level =1-alpha)[,-6]
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

  lv.kutsu<-deparse(substitute(formula))
  lv.data<-deparse(substitute(data))
  # Event cell frequencies, if too low replace it with high number
  if(min(lv.2$event)<lowest.N) warning(paste("lowest cell frequency under",
                                             lowest.N,", replaced by 999999"))
  lv.2$event<-ifelse(lv.2$event<lowest.N,9999999,lv.2$event)
  # Add attributes
  attributes(lv.2)$ratetable.call<-lv.kutsu
  attributes(lv.2)$ratetable.data<-lv.data
  attributes(lv.2)$ratetable.time<-date()
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
#' @param node.attr general node attributers.
#' Attributes like shape, color, fillcolor, etc. for nodes.
#' Consult Graphviz documentation for details
#' \url{https://www.graphviz.org/doc/info/attrs.html}.
#' @param edge.attr general edge (line) attributers.
#' Attributes like  color, arrowhead, fontcolor  etc. for edges.
#' Consult Graphviz documentation for details
#' \url{https://www.graphviz.org/doc/info/attrs.html}
#' @param fontsizeN font size for nodes
#' @param fontsizeL font size for edges
#' @param show.gr  should graph be shown. If TRUE,
#' function DiagrammeR::grViz is used to show graph.
#' @author Jari Haukka <jari.haukka@@helsinki.fi>
#' @return Character vector containing Graphviz script. This may used
#' to create graph by  DiagrammeR::grViz function.
#' @seealso \code{\link{grViz}}
#' @seealso
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
#' summary(bmtr)
#' boxesLx(bmtr)
#' boxesLx(bmtr,layout="dot")
#' boxesLx(bmtr,node.attr='shape=hexagon color=navy style=filled fillcolor=lightblue',
#' edge.attr = ' color=steelblue arrowhead=vee fontcolor="#8801d7" ',
#' layout="circo",prop.penwidth=TRUE)
boxesLx<-function(x,
                  layout="circo",
                  prop.penwidth=FALSE,
                  scale.Y = 1,
                  node.attr="shape=box",
                  edge.attr="minlen=1",
                  fontsizeN=14,fontsizeL=8,
                  show.gr=TRUE)
{
  # layout<-c("dot","neato","twopi","circo")
  if(!("Lexis"%in%class(x))) stop("x must be Lexis object")
  lv.1<-summary(x)
  lv.1$Transitions[,"Risk time:"]<-lv.1$Transitions[,"Risk time:"]/scale.Y
  lv.len<-nrow(lv.1$Transitions)-1
  lv.N<-lv.1$Transitions[1:lv.len,!(colnames(lv.1$Transitions)%in%c(" Records:", " Events:","Risk time:", " Persons:"))]
  lv.RT<-lv.1$Transitions[1:lv.len,"Risk time:"]
  lv.rate<-lv.N/lv.1$Transitions[1:lv.len,"Risk time:"]
  lv.nd.names<-paste0('"',colnames(lv.N),'"')
  lv.head.gv<-c("digraph G {\n" ,"layout=",layout,";\n",
                "node [",node.attr,"];\n","edge [",edge.attr,"];\n")
  lv.main.gv<-NULL
  lv.main.gv<-c(lv.main.gv,paste0('"',rownames(lv.N),'"'," [label='",
                                  paste(rownames(lv.N),'\\n',
                                        round(lv.RT,2)),
                                  '" , fontsize=',fontsizeN=10,'];\n'))
  lv.pen1<-log(lv.rate/min(lv.rate[lv.rate>0]))
  lv.pen2<-3*(lv.pen1/max(lv.pen1))+0.2

  for(lv.from in seq(nrow(lv.N))) for(lv.to  in seq(ncol(lv.N))){
    if(lv.N[lv.from,lv.to]>0){
      lv.pen3<-ifelse(prop.penwidth,paste0(" penwidth=",lv.pen2[lv.from,lv.to])," ")
      lv.main.gv<-c(lv.main.gv,paste0('"',rownames(lv.N)[lv.from],'"'),"->",paste0('"',colnames(lv.N)[lv.to],'"'),
                    paste0('[label="',lv.N[lv.from,lv.to],'\\n(',
                           round(lv.rate[lv.from,lv.to],2),')",fontsize=',fontsizeL," ",lv.pen3,']'),";\n")
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
#'
#' @param time  time variables
#' @param status status indicator Lowest value used as sensoring.
#' If only one unique value detected, all are assumed events
#' @param breaks time is splitted with these values
#' @param knots knots for natural splines used in estimation of hazard function
#' @param time.eval in which time points hazard function is evaluate.
#' @param alpha significance level for confidence intervals
#' @param ... parameters for glm
#' @author Jari Haukka \email{jari.haukka@@helsinki.fi}
#' @return Returns data frame with time and hazard function values with attribute 'estim.hazard.param'
#' containing estmation parameters (breaks and knots)
#' @export
#' @examples
#' library(survival)
#' tmp.hz<-estim.hazard(time=lung$time,status=lung$status)
#' head(tmp.hz,2)
#' attributes(tmp.hz)$estim.hazard.param # estimation parameters
estim.hazard<-function(time,status,breaks,knots,time.eval=breaks,alpha=0.05,...){
  # require(survival)
  # require(splines)
  # require(Epi)
  # Check if there is only one unique status value
  if(length(unique(status))==1)lv.status<-1
  else {
    lv.min.status<-min(status)
    lv.status<-ifelse(status==lv.min.status,0,1)
  }
  # Make Lexis
  lv.Lx <- Epi::Lexis( exit=list(fu.time=time),
                  exit.status=lv.status)
  # Split Lexis, default is 100 breaks
  if(missing(breaks)) breaks<-with(lv.Lx,c(0,seq(from=min(lex.dur),to=max(lex.dur),
                                             length=100)))
  lv.Lx.s <- Epi::splitLexis( lv.Lx, "fu.time", breaks=breaks )
  # Poisson model with splines
  # Define knots, if missing
  if(missing(knots)) knots<-with(lv.Lx,seq(from=min(lex.dur),to=max(lex.dur),
                                             length=7))[-c(1,7)]
  lv.y<-Epi::status(lv.Lx.s, "exit")
  lv.offs<-Epi::dur(lv.Lx.s)
  lv.tb<-Epi::timeBand( lv.Lx.s, "fu.time", "left" )
  lv.Xmat<-splines::ns(lv.tb,knots=knots,intercept = TRUE)
  lv.m<-glm(lv.y~lv.Xmat+offset(log(lv.offs))-1,family = poisson,maxit=50,...)
  if(missing(time.eval))time.eval<-breaks
  lv.Xmat.eval<-splines::ns(time.eval,knots=knots,intercept = TRUE)
  lv.lambda <- Epi::ci.lin( lv.m, ctr.mat=lv.Xmat.eval, Exp=TRUE,alpha=alpha )
  lv.ulos<-data.frame(time.eval=time.eval,lv.lambda)
  names(lv.ulos)[6:8]<-c("haz","haz.lo","haz.hi")
  lv.ulos<-lv.ulos[,-c(4,5)]
  attributes(lv.ulos)$estim.hazard.param<-list(breaks=breaks,knots=knots)
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
