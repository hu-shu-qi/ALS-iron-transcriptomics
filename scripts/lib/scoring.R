adj <- function(p){z<-p.adjust(ifelse(is.finite(p),p,1),'BH');z[!is.finite(p)]<-NA_real_;z}
score3<-function(x,seed){
 set.seed(seed); usable<-lapply(sets,intersect,y=rownames(x))
 u<-as.matrix(ScoreSignatures_UCell(x,features=usable,maxRank=1500,missing_genes='skip',BPPARAM=SerialParam(),ncores=1))
 colnames(u)<-sub('_UCell$','',colnames(u))
 set.seed(seed)
 a<-t(as.matrix(getAUC(AUCell_calcAUC(usable,AUCell_buildRankings(as.matrix(x),plotStats=FALSE,verbose=FALSE),aucMaxRank=ceiling(.05*nrow(x)),nCores=1,verbose=FALSE))))
 rk<-rankGenes(as.matrix(x));si<-sapply(usable,function(gs)simpleScore(rk,upSet=gs,centerScore=TRUE)$TotalScore);rownames(si)<-colnames(x)
 ans<-list(UCell=u[,names(sets),drop=FALSE],AUCell=a[,names(sets),drop=FALSE],singscore=si[,names(sets),drop=FALSE])
 stopifnot(all(vapply(ans,function(z)identical(rownames(z),colnames(x)) && all(is.finite(z)),TRUE)))
 ans
}
fit_score<-function(md,form){
 md<-copy(md);md[,group:=factor(group,levels=c('Control','ALS'))]
 X<-model.matrix(form,md);stopifnot(qr(X)$rank==ncol(X),nrow(X)==nrow(md))
 df<-nrow(X)-ncol(X);j<-match('groupALS',colnames(X));stopifnot(!is.na(j),df>0)
 if(sd(md$score)<1e-14)return(data.table(effect=NA_real_,SE=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,P_HC3=NA_real_,LODO_same_direction=NA_real_,status='constant_score',n=nrow(md),df=df))
 f<-lm.fit(X,md$score);bread<-solve(crossprod(X));h<-rowSums(qr.Q(qr(X))^2)
 se<-sqrt(bread[j,j]*sum(f$residuals^2)/df);ef<-unname(f$coefficients[j]);se3<-sqrt((bread%*%crossprod(X,X*as.numeric(f$residuals^2/(1-h)^2))%*%bread)[j,j])
 li<-sapply(seq_len(nrow(X)),function(i)if(qr(X[-i,,drop=FALSE])$rank==ncol(X))lm.fit(X[-i,,drop=FALSE],md$score[-i])$coefficients[j] else NA_real_)
 data.table(effect=ef,SE=se,CI_low=ef-qt(.975,df)*se,CI_high=ef+qt(.975,df)*se,P=2*pt(-abs(ef/se),df),P_HC3=2*pt(-abs(ef/se3),df),LODO_same_direction=mean(sign(li)==sign(ef),na.rm=TRUE),status='tested',n=nrow(md),df=df)
}

collect<-function(x,md,context,layer,forms,seed){
 scores<-score3(x,seed)
 coverage[[length(coverage)+1]]<<-data.table(layer=layer,context=context,panel=names(sets),requested=lengths(sets),detected=lengths(lapply(sets,intersect,y=rownames(x))))
 for(meth in names(scores))for(pan in names(sets)){
  d<-copy(md);d[,score:=scores[[meth]][,pan]];d[,`:=`(layer=layer,context=context,method=meth,panel=pan)]
  allsc[[length(allsc)+1]]<<-d
  for(mod in names(forms)){
   z<-copy(d)
   if(mod=='sex_complete')z<-z[!is.na(sex)]
   if(mod=='age_sex_complete')z<-z[complete.cases(z[,.(age,sex)])]
   if(mod=='exclude_comorbidity')z<-z[sensitivity_eligible==TRUE]
   t<-fit_score(z,forms[[mod]]);t[,`:=`(layer=layer,context=context,method=meth,panel=pan,model=mod)]
   tests[[length(tests)+1]]<<-t
  }
 }
}
