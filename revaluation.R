evalcg_10000_0.05 = function(graph){
  tgdata <- rnorm.cg(10000,graph,get.normal.dist(graph))
  #lcd方法
  startlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.05)
  tg.jtree <- ug.to.jtree(tgug)
  #startlcd2 = proc.time()
  tg.patlcd <- learn.mec.norm(tg.jtree,cov(tgdata),10000,0.05,"CG")
  learntLCGlcd<-studeny_rules(tg.patlcd)
  endlcd1 = proc.time() -startlcd1
  #endlcd2 = proc.time() -startlcd2
  time1_lcd = endlcd1[[3]][[1]]
  #time2_lcd = endlcd2[[3]][[1]]
  #mmlcd方法(道德图出发)
  #startmmlcd1 = proc.time()
 # tgmor <- mor(as.data.frame(tgdata))
  #tgminitri <- minimal_triangMAT(tgmor)
  #tg.ctree <- ug.to.jtree1(tgminitri)
  #startmmlcd2 = proc.time()
  #tg.patmmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),10000,0.05,"CG")
  #learntLCGmmlcd<-studeny_rules(tg.patmmlcd)
  #endmmlcd1 = proc.time() -startmmlcd1
  #endmmlcd2 = proc.time() -startmmlcd2
  #time1_mmlcd = endmmlcd1[[3]][[1]]
  #time2_mmlcd = endmmlcd2[[3]][[1]]  
  #mlcd方法(无向独立图出发)
  startmlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.05)
  tgminitri <- minimal_triangMAT(tgug)
  tg.ctree <- ug.to.jtree1(tgminitri)
  #startmlcd2 = proc.time()
  tg.patmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),10000,0.05,"CG")
  learntLCGmlcd<-studeny_rules(tg.patmlcd)
  endmlcd1 = proc.time() -startmlcd1
  #endmlcd2 = proc.time() -startmlcd2
  time1_mlcd = endmlcd1[[3]][[1]]
  #time2_mlcd = endmlcd2[[3]][[1]]  
  #SPC方法
  startspc = proc.time()
  tg.patspc <- learn.lwf.norm(tgdata,0.05,method="stable",LCG=TRUE)
  endspc = proc.time() -startspc
  time_spc = endspc[[3]][[1]]
  #mbcsp方法
  startmbcsp = proc.time()
  mb=mbcsp(data.frame(tgdata),alpha = 0.05)
  tg.mpscp = learn.vstruct(data.frame(tgdata),mb,alpha = 0.05)
  endmbcsp = proc.time() -startmbcsp
  time_mbcsp = endmbcsp[[3]][[1]]
  #指标
  fnlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fn
  fplcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fp
  tprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$TPR
  fprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$FPR
  SHDlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$SHD
  prelcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$PRE
  f1lcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$F1
  acclcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$ACC
  #indexlcd <-c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd)  
  #fnmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fn
  #fpmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fp
  #tprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$TPR
  #fprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$FPR
  #SHDmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$SHD
  #premmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$PRE
  #f1mmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$F1
  #indexmmlcd <- c(fnmmlcd,fpmmlcd,tprmmlcd,fprmmlcd,SHDmmlcd,time_mmlcd)  
  fnmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fn
  fpmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fp
  tprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$TPR
  fprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$FPR
  SHDmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$SHD
  premlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$PRE
  f1mlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$F1
  accmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$ACC
  #indexmlcd <- c(fnmlcd,fpmlcd,tprmlcd,fprmlcd,SHDmlcd,time_mlcd)  
  fnspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fn
  fpspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fp
  tprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$TPR
  fprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$FPR
  SHDspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$SHD
  prespc <- comp.cgs(pattern(graph),tg.patspc$matrix)$PRE
  f1spc <- comp.cgs(pattern(graph),tg.patspc$matrix)$F1
  accspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$ACC
  #indexspc <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc) 
  fnmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fn
  fpmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fp
  tprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$TPR
  fprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$FPR
  SHDmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$SHD
  prembcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$PRE
  f1mbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$F1
  accmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$ACC
  #indexmbcsp <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc)
  index <- c(fnlcd,fplcd,tprlcd,fprlcd,prelcd,f1lcd,acclcd,SHDlcd,time1_lcd,
             fnmlcd,fpmlcd,tprmlcd,fprmlcd,premlcd,f1mlcd,accmlcd,SHDmlcd,time1_mlcd,
             fnspc,fpspc,tprspc,fprspc,prespc,f1spc,accspc,SHDspc,time_spc,
             fnmbcsp,fpmbcsp,tprmbcsp,fprmbcsp,prembcsp,f1mbcsp,accmbcsp,SHDmbcsp,time_mbcsp)
  #index1 <- c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd,fnmlcd1,fpmlcd1,tprmlcd1,fprmlcd1,SHDmlcd1,time_mlcd1)
  #return(index)
  return(index)
}

evalcg_10000_0.01 = function(graph){
  tgdata <- rnorm.cg(10000,graph,get.normal.dist(graph))
  #lcd方法
  startlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.01)
  tg.jtree <- ug.to.jtree(tgug)
  #startlcd2 = proc.time()
  tg.patlcd <- learn.mec.norm(tg.jtree,cov(tgdata),10000,0.01,"CG")
  learntLCGlcd<-studeny_rules(tg.patlcd)
  endlcd1 = proc.time() -startlcd1
  #endlcd2 = proc.time() -startlcd2
  time1_lcd = endlcd1[[3]][[1]]
  #time2_lcd = endlcd2[[3]][[1]]
  #mmlcd方法(道德图出发)
  #startmmlcd1 = proc.time()
  # tgmor <- mor(as.data.frame(tgdata))
  #tgminitri <- minimal_triangMAT(tgmor)
  #tg.ctree <- ug.to.jtree1(tgminitri)
  #startmmlcd2 = proc.time()
  #tg.patmmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),10000,0.05,"CG")
  #learntLCGmmlcd<-studeny_rules(tg.patmmlcd)
  #endmmlcd1 = proc.time() -startmmlcd1
  #endmmlcd2 = proc.time() -startmmlcd2
  #time1_mmlcd = endmmlcd1[[3]][[1]]
  #time2_mmlcd = endmmlcd2[[3]][[1]]  
  #mlcd方法(无向独立图出发)
  startmlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.01)
  tgminitri <- minimal_triangMAT(tgug)
  tg.ctree <- ug.to.jtree1(tgminitri)
  #startmlcd2 = proc.time()
  tg.patmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),10000,0.01,"CG")
  learntLCGmlcd<-studeny_rules(tg.patmlcd)
  endmlcd1 = proc.time() -startmlcd1
  #endmlcd2 = proc.time() -startmlcd2
  time1_mlcd = endmlcd1[[3]][[1]]
  #time2_mlcd = endmlcd2[[3]][[1]]  
  #SPC方法
  startspc = proc.time()
  tg.patspc <- learn.lwf.norm(tgdata,0.01,method="stable",LCG=TRUE)
  endspc = proc.time() -startspc
  time_spc = endspc[[3]][[1]]
  #mbcsp方法
  startmbcsp = proc.time()
  mb=mbcsp(data.frame(tgdata),alpha = 0.01)
  tg.mpscp = learn.vstruct(data.frame(tgdata),mb,alpha = 0.01)
  endmbcsp = proc.time() -startmbcsp
  time_mbcsp = endmbcsp[[3]][[1]]
  #指标
  fnlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fn
  fplcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fp
  tprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$TPR
  fprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$FPR
  SHDlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$SHD
  prelcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$PRE
  f1lcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$F1
  acclcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$ACC
  #indexlcd <-c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd)  
  #fnmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fn
  #fpmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fp
  #tprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$TPR
  #fprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$FPR
  #SHDmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$SHD
  #premmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$PRE
  #f1mmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$F1
  #indexmmlcd <- c(fnmmlcd,fpmmlcd,tprmmlcd,fprmmlcd,SHDmmlcd,time_mmlcd)  
  fnmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fn
  fpmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fp
  tprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$TPR
  fprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$FPR
  SHDmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$SHD
  premlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$PRE
  f1mlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$F1
  accmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$ACC
  #indexmlcd <- c(fnmlcd,fpmlcd,tprmlcd,fprmlcd,SHDmlcd,time_mlcd)  
  fnspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fn
  fpspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fp
  tprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$TPR
  fprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$FPR
  SHDspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$SHD
  prespc <- comp.cgs(pattern(graph),tg.patspc$matrix)$PRE
  f1spc <- comp.cgs(pattern(graph),tg.patspc$matrix)$F1
  accspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$ACC
  #indexspc <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc) 
  fnmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fn
  fpmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fp
  tprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$TPR
  fprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$FPR
  SHDmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$SHD
  prembcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$PRE
  f1mbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$F1
  accmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$ACC
  #indexmbcsp <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc)
  index <- c(fnlcd,fplcd,tprlcd,fprlcd,prelcd,f1lcd,acclcd,SHDlcd,time1_lcd,
             fnmlcd,fpmlcd,tprmlcd,fprmlcd,premlcd,f1mlcd,accmlcd,SHDmlcd,time1_mlcd,
             fnspc,fpspc,tprspc,fprspc,prespc,f1spc,accspc,SHDspc,time_spc,
             fnmbcsp,fpmbcsp,tprmbcsp,fprmbcsp,prembcsp,f1mbcsp,accmbcsp,SHDmbcsp,time_mbcsp)
  #index1 <- c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd,fnmlcd1,fpmlcd1,tprmlcd1,fprmlcd1,SHDmlcd1,time_mlcd1)
  #return(index)
  return(index)
}

evalcg_5000_0.05 = function(graph){
  tgdata <- rnorm.cg(5000,graph,get.normal.dist(graph))
  #lcd方法
  startlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.05)
  tg.jtree <- ug.to.jtree(tgug)
  #startlcd2 = proc.time()
  tg.patlcd <- learn.mec.norm(tg.jtree,cov(tgdata),5000,0.05,"CG")
  learntLCGlcd<-studeny_rules(tg.patlcd)
  endlcd1 = proc.time() -startlcd1
  #endlcd2 = proc.time() -startlcd2
  time1_lcd = endlcd1[[3]][[1]]
  #time2_lcd = endlcd2[[3]][[1]]
  #mmlcd方法(道德图出发)
  #startmmlcd1 = proc.time()
  # tgmor <- mor(as.data.frame(tgdata))
  #tgminitri <- minimal_triangMAT(tgmor)
  #tg.ctree <- ug.to.jtree1(tgminitri)
  #startmmlcd2 = proc.time()
  #tg.patmmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),10000,0.05,"CG")
  #learntLCGmmlcd<-studeny_rules(tg.patmmlcd)
  #endmmlcd1 = proc.time() -startmmlcd1
  #endmmlcd2 = proc.time() -startmmlcd2
  #time1_mmlcd = endmmlcd1[[3]][[1]]
  #time2_mmlcd = endmmlcd2[[3]][[1]]  
  #mlcd方法(无向独立图出发)
  startmlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.05)
  tgminitri <- minimal_triangMAT(tgug)
  tg.ctree <- ug.to.jtree1(tgminitri)
  #startmlcd2 = proc.time()
  tg.patmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),5000,0.05,"CG")
  learntLCGmlcd<-studeny_rules(tg.patmlcd)
  endmlcd1 = proc.time() -startmlcd1
  #endmlcd2 = proc.time() -startmlcd2
  time1_mlcd = endmlcd1[[3]][[1]]
  #time2_mlcd = endmlcd2[[3]][[1]]  
  #SPC方法
  startspc = proc.time()
  tg.patspc <- learn.lwf.norm(tgdata,0.05,method="stable",LCG=TRUE)
  endspc = proc.time() -startspc
  time_spc = endspc[[3]][[1]]
  #mbcsp方法
  startmbcsp = proc.time()
  mb=mbcsp(data.frame(tgdata),alpha = 0.05)
  tg.mpscp = learn.vstruct(data.frame(tgdata),mb,alpha = 0.05)
  endmbcsp = proc.time() -startmbcsp
  time_mbcsp = endmbcsp[[3]][[1]]
  #指标
  fnlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fn
  fplcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fp
  tprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$TPR
  fprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$FPR
  SHDlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$SHD
  prelcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$PRE
  f1lcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$F1
  acclcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$ACC
  #indexlcd <-c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd)  
  #fnmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fn
  #fpmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fp
  #tprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$TPR
  #fprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$FPR
  #SHDmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$SHD
  #premmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$PRE
  #f1mmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$F1
  #indexmmlcd <- c(fnmmlcd,fpmmlcd,tprmmlcd,fprmmlcd,SHDmmlcd,time_mmlcd)  
  fnmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fn
  fpmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fp
  tprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$TPR
  fprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$FPR
  SHDmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$SHD
  premlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$PRE
  f1mlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$F1
  accmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$ACC
  #indexmlcd <- c(fnmlcd,fpmlcd,tprmlcd,fprmlcd,SHDmlcd,time_mlcd)  
  fnspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fn
  fpspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fp
  tprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$TPR
  fprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$FPR
  SHDspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$SHD
  prespc <- comp.cgs(pattern(graph),tg.patspc$matrix)$PRE
  f1spc <- comp.cgs(pattern(graph),tg.patspc$matrix)$F1
  accspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$ACC
  #indexspc <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc) 
  fnmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fn
  fpmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fp
  tprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$TPR
  fprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$FPR
  SHDmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$SHD
  prembcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$PRE
  f1mbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$F1
  accmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$ACC
  #indexmbcsp <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc)
  index <- c(fnlcd,fplcd,tprlcd,fprlcd,prelcd,f1lcd,acclcd,SHDlcd,time1_lcd,
             fnmlcd,fpmlcd,tprmlcd,fprmlcd,premlcd,f1mlcd,accmlcd,SHDmlcd,time1_mlcd,
             fnspc,fpspc,tprspc,fprspc,prespc,f1spc,accspc,SHDspc,time_spc,
             fnmbcsp,fpmbcsp,tprmbcsp,fprmbcsp,prembcsp,f1mbcsp,accmbcsp,SHDmbcsp,time_mbcsp)
  #index1 <- c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd,fnmlcd1,fpmlcd1,tprmlcd1,fprmlcd1,SHDmlcd1,time_mlcd1)
  #return(index)
  return(index)
}

evalcg_5000_0.01 = function(graph){
  tgdata <- rnorm.cg(5000,graph,get.normal.dist(graph))
  #lcd方法
  startlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.01)
  tg.jtree <- ug.to.jtree(tgug)
  #startlcd2 = proc.time()
  tg.patlcd <- learn.mec.norm(tg.jtree,cov(tgdata),5000,0.01,"CG")
  learntLCGlcd<-studeny_rules(tg.patlcd)
  endlcd1 = proc.time() -startlcd1
  #endlcd2 = proc.time() -startlcd2
  time1_lcd = endlcd1[[3]][[1]]
  #time2_lcd = endlcd2[[3]][[1]]
  #mmlcd方法(道德图出发)
  #startmmlcd1 = proc.time()
  # tgmor <- mor(as.data.frame(tgdata))
  #tgminitri <- minimal_triangMAT(tgmor)
  #tg.ctree <- ug.to.jtree1(tgminitri)
  #startmmlcd2 = proc.time()
  #tg.patmmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),10000,0.05,"CG")
  #learntLCGmmlcd<-studeny_rules(tg.patmmlcd)
  #endmmlcd1 = proc.time() -startmmlcd1
  #endmmlcd2 = proc.time() -startmmlcd2
  #time1_mmlcd = endmmlcd1[[3]][[1]]
  #time2_mmlcd = endmmlcd2[[3]][[1]]  
  #mlcd方法(无向独立图出发)
  startmlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.01)
  tgminitri <- minimal_triangMAT(tgug)
  tg.ctree <- ug.to.jtree1(tgminitri)
  #startmlcd2 = proc.time()
  tg.patmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),5000,0.01,"CG")
  learntLCGmlcd<-studeny_rules(tg.patmlcd)
  endmlcd1 = proc.time() -startmlcd1
  #endmlcd2 = proc.time() -startmlcd2
  time1_mlcd = endmlcd1[[3]][[1]]
  #time2_mlcd = endmlcd2[[3]][[1]]  
  #SPC方法
  startspc = proc.time()
  tg.patspc <- learn.lwf.norm(tgdata,0.01,method="stable",LCG=TRUE)
  endspc = proc.time() -startspc
  time_spc = endspc[[3]][[1]]
  #mbcsp方法
  startmbcsp = proc.time()
  mb=mbcsp(data.frame(tgdata),alpha = 0.01)
  tg.mpscp = learn.vstruct(data.frame(tgdata),mb,alpha = 0.01)
  endmbcsp = proc.time() -startmbcsp
  time_mbcsp = endmbcsp[[3]][[1]]
  #指标
  fnlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fn
  fplcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fp
  tprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$TPR
  fprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$FPR
  SHDlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$SHD
  prelcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$PRE
  f1lcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$F1
  acclcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$ACC
  #indexlcd <-c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd)  
  #fnmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fn
  #fpmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fp
  #tprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$TPR
  #fprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$FPR
  #SHDmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$SHD
  #premmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$PRE
  #f1mmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$F1
  #indexmmlcd <- c(fnmmlcd,fpmmlcd,tprmmlcd,fprmmlcd,SHDmmlcd,time_mmlcd)  
  fnmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fn
  fpmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fp
  tprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$TPR
  fprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$FPR
  SHDmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$SHD
  premlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$PRE
  f1mlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$F1
  accmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$ACC
  #indexmlcd <- c(fnmlcd,fpmlcd,tprmlcd,fprmlcd,SHDmlcd,time_mlcd)  
  fnspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fn
  fpspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fp
  tprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$TPR
  fprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$FPR
  SHDspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$SHD
  prespc <- comp.cgs(pattern(graph),tg.patspc$matrix)$PRE
  f1spc <- comp.cgs(pattern(graph),tg.patspc$matrix)$F1
  accspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$ACC
  #indexspc <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc) 
  fnmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fn
  fpmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fp
  tprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$TPR
  fprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$FPR
  SHDmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$SHD
  prembcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$PRE
  f1mbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$F1
  accmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$ACC
  #indexmbcsp <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc)
  index <- c(fnlcd,fplcd,tprlcd,fprlcd,prelcd,f1lcd,acclcd,SHDlcd,time1_lcd,
             fnmlcd,fpmlcd,tprmlcd,fprmlcd,premlcd,f1mlcd,accmlcd,SHDmlcd,time1_mlcd,
             fnspc,fpspc,tprspc,fprspc,prespc,f1spc,accspc,SHDspc,time_spc,
             fnmbcsp,fpmbcsp,tprmbcsp,fprmbcsp,prembcsp,f1mbcsp,accmbcsp,SHDmbcsp,time_mbcsp)
  #index1 <- c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd,fnmlcd1,fpmlcd1,tprmlcd1,fprmlcd1,SHDmlcd1,time_mlcd1)
  #return(index)
  return(index)
}

evalcg_1000_0.05 = function(graph){
  tgdata <- rnorm.cg(1000,graph,get.normal.dist(graph))
  #lcd方法
  startlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.05)
  tg.jtree <- ug.to.jtree(tgug)
  #startlcd2 = proc.time()
  tg.patlcd <- learn.mec.norm(tg.jtree,cov(tgdata),1000,0.05,"CG")
  learntLCGlcd<-studeny_rules(tg.patlcd)
  endlcd1 = proc.time() -startlcd1
  #endlcd2 = proc.time() -startlcd2
  time1_lcd = endlcd1[[3]][[1]]
  #time2_lcd = endlcd2[[3]][[1]]
  #mmlcd方法(道德图出发)
  #startmmlcd1 = proc.time()
  # tgmor <- mor(as.data.frame(tgdata))
  #tgminitri <- minimal_triangMAT(tgmor)
  #tg.ctree <- ug.to.jtree1(tgminitri)
  #startmmlcd2 = proc.time()
  #tg.patmmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),10000,0.05,"CG")
  #learntLCGmmlcd<-studeny_rules(tg.patmmlcd)
  #endmmlcd1 = proc.time() -startmmlcd1
  #endmmlcd2 = proc.time() -startmmlcd2
  #time1_mmlcd = endmmlcd1[[3]][[1]]
  #time2_mmlcd = endmmlcd2[[3]][[1]]  
  #mlcd方法(无向独立图出发)
  startmlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.05)
  tgminitri <- minimal_triangMAT(tgug)
  tg.ctree <- ug.to.jtree1(tgminitri)
  #startmlcd2 = proc.time()
  tg.patmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),1000,0.05,"CG")
  learntLCGmlcd<-studeny_rules(tg.patmlcd)
  endmlcd1 = proc.time() -startmlcd1
  #endmlcd2 = proc.time() -startmlcd2
  time1_mlcd = endmlcd1[[3]][[1]]
  #time2_mlcd = endmlcd2[[3]][[1]]  
  #SPC方法
  startspc = proc.time()
  tg.patspc <- learn.lwf.norm(tgdata,0.05,method="stable",LCG=TRUE)
  endspc = proc.time() -startspc
  time_spc = endspc[[3]][[1]]
  #mbcsp方法
  startmbcsp = proc.time()
  mb=mbcsp(data.frame(tgdata),alpha = 0.05)
  tg.mpscp = learn.vstruct(data.frame(tgdata),mb,alpha = 0.05)
  endmbcsp = proc.time() -startmbcsp
  time_mbcsp = endmbcsp[[3]][[1]]
  #指标
  fnlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fn
  fplcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fp
  tprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$TPR
  fprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$FPR
  SHDlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$SHD
  prelcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$PRE
  f1lcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$F1
  acclcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$ACC
  #indexlcd <-c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd)  
  #fnmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fn
  #fpmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fp
  #tprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$TPR
  #fprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$FPR
  #SHDmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$SHD
  #premmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$PRE
  #f1mmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$F1
  #indexmmlcd <- c(fnmmlcd,fpmmlcd,tprmmlcd,fprmmlcd,SHDmmlcd,time_mmlcd)  
  fnmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fn
  fpmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fp
  tprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$TPR
  fprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$FPR
  SHDmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$SHD
  premlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$PRE
  f1mlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$F1
  accmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$ACC
  #indexmlcd <- c(fnmlcd,fpmlcd,tprmlcd,fprmlcd,SHDmlcd,time_mlcd)  
  fnspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fn
  fpspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fp
  tprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$TPR
  fprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$FPR
  SHDspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$SHD
  prespc <- comp.cgs(pattern(graph),tg.patspc$matrix)$PRE
  f1spc <- comp.cgs(pattern(graph),tg.patspc$matrix)$F1
  accspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$ACC
  #indexspc <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc) 
  fnmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fn
  fpmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fp
  tprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$TPR
  fprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$FPR
  SHDmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$SHD
  prembcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$PRE
  f1mbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$F1
  accmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$ACC
  #indexmbcsp <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc)
  index <- c(fnlcd,fplcd,tprlcd,fprlcd,prelcd,f1lcd,acclcd,SHDlcd,time1_lcd,
             fnmlcd,fpmlcd,tprmlcd,fprmlcd,premlcd,f1mlcd,accmlcd,SHDmlcd,time1_mlcd,
             fnspc,fpspc,tprspc,fprspc,prespc,f1spc,accspc,SHDspc,time_spc,
             fnmbcsp,fpmbcsp,tprmbcsp,fprmbcsp,prembcsp,f1mbcsp,accmbcsp,SHDmbcsp,time_mbcsp)
  #index1 <- c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd,fnmlcd1,fpmlcd1,tprmlcd1,fprmlcd1,SHDmlcd1,time_mlcd1)
  #return(index)
  return(index)
}

evalcg_1000_0.01 = function(graph){
  tgdata <- rnorm.cg(1000,graph,get.normal.dist(graph))
  #lcd方法
  startlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.01)
  tg.jtree <- ug.to.jtree(tgug)
  #startlcd2 = proc.time()
  tg.patlcd <- learn.mec.norm(tg.jtree,cov(tgdata),1000,0.01,"CG")
  learntLCGlcd<-studeny_rules(tg.patlcd)
  endlcd1 = proc.time() -startlcd1
  #endlcd2 = proc.time() -startlcd2
  time1_lcd = endlcd1[[3]][[1]]
  #time2_lcd = endlcd2[[3]][[1]]
  #mmlcd方法(道德图出发)
  #startmmlcd1 = proc.time()
  # tgmor <- mor(as.data.frame(tgdata))
  #tgminitri <- minimal_triangMAT(tgmor)
  #tg.ctree <- ug.to.jtree1(tgminitri)
  #startmmlcd2 = proc.time()
  #tg.patmmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),10000,0.05,"CG")
  #learntLCGmmlcd<-studeny_rules(tg.patmmlcd)
  #endmmlcd1 = proc.time() -startmmlcd1
  #endmmlcd2 = proc.time() -startmmlcd2
  #time1_mmlcd = endmmlcd1[[3]][[1]]
  #time2_mmlcd = endmmlcd2[[3]][[1]]  
  #mlcd方法(无向独立图出发)
  startmlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.01)
  tgminitri <- minimal_triangMAT(tgug)
  tg.ctree <- ug.to.jtree1(tgminitri)
  #startmlcd2 = proc.time()
  tg.patmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),1000,0.01,"CG")
  learntLCGmlcd<-studeny_rules(tg.patmlcd)
  endmlcd1 = proc.time() -startmlcd1
  #endmlcd2 = proc.time() -startmlcd2
  time1_mlcd = endmlcd1[[3]][[1]]
  #time2_mlcd = endmlcd2[[3]][[1]]  
  #SPC方法
  startspc = proc.time()
  tg.patspc <- learn.lwf.norm(tgdata,0.01,method="stable",LCG=TRUE)
  endspc = proc.time() -startspc
  time_spc = endspc[[3]][[1]]
  #mbcsp方法
  startmbcsp = proc.time()
  mb=mbcsp(data.frame(tgdata),alpha = 0.01)
  tg.mpscp = learn.vstruct(data.frame(tgdata),mb,alpha = 0.01)
  endmbcsp = proc.time() -startmbcsp
  time_mbcsp = endmbcsp[[3]][[1]]
  #指标
  fnlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fn
  fplcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fp
  tprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$TPR
  fprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$FPR
  SHDlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$SHD
  prelcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$PRE
  f1lcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$F1
  acclcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$ACC
  #indexlcd <-c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd)  
  #fnmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fn
  #fpmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fp
  #tprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$TPR
  #fprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$FPR
  #SHDmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$SHD
  #premmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$PRE
  #f1mmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$F1
  #indexmmlcd <- c(fnmmlcd,fpmmlcd,tprmmlcd,fprmmlcd,SHDmmlcd,time_mmlcd)  
  fnmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fn
  fpmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fp
  tprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$TPR
  fprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$FPR
  SHDmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$SHD
  premlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$PRE
  f1mlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$F1
  accmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$ACC
  #indexmlcd <- c(fnmlcd,fpmlcd,tprmlcd,fprmlcd,SHDmlcd,time_mlcd)  
  fnspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fn
  fpspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fp
  tprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$TPR
  fprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$FPR
  SHDspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$SHD
  prespc <- comp.cgs(pattern(graph),tg.patspc$matrix)$PRE
  f1spc <- comp.cgs(pattern(graph),tg.patspc$matrix)$F1
  accspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$ACC
  #indexspc <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc) 
  fnmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fn
  fpmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fp
  tprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$TPR
  fprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$FPR
  SHDmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$SHD
  prembcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$PRE
  f1mbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$F1
  accmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$ACC
  #indexmbcsp <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc)
  index <- c(fnlcd,fplcd,tprlcd,fprlcd,prelcd,f1lcd,acclcd,SHDlcd,time1_lcd,
             fnmlcd,fpmlcd,tprmlcd,fprmlcd,premlcd,f1mlcd,accmlcd,SHDmlcd,time1_mlcd,
             fnspc,fpspc,tprspc,fprspc,prespc,f1spc,accspc,SHDspc,time_spc,
             fnmbcsp,fpmbcsp,tprmbcsp,fprmbcsp,prembcsp,f1mbcsp,accmbcsp,SHDmbcsp,time_mbcsp)
  #index1 <- c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd,fnmlcd1,fpmlcd1,tprmlcd1,fprmlcd1,SHDmlcd1,time_mlcd1)
  #return(index)
  return(index)
}

evalcg_500_0.05 = function(graph){
  tgdata <- rnorm.cg(500,graph,get.normal.dist(graph))
  #lcd方法
  startlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.05)
  tg.jtree <- ug.to.jtree(tgug)
  #startlcd2 = proc.time()
  tg.patlcd <- learn.mec.norm(tg.jtree,cov(tgdata),500,0.05,"CG")
  learntLCGlcd<-studeny_rules(tg.patlcd)
  endlcd1 = proc.time() -startlcd1
  #endlcd2 = proc.time() -startlcd2
  time1_lcd = endlcd1[[3]][[1]]
  #time2_lcd = endlcd2[[3]][[1]]
  #mmlcd方法(道德图出发)
  #startmmlcd1 = proc.time()
  # tgmor <- mor(as.data.frame(tgdata))
  #tgminitri <- minimal_triangMAT(tgmor)
  #tg.ctree <- ug.to.jtree1(tgminitri)
  #startmmlcd2 = proc.time()
  #tg.patmmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),10000,0.05,"CG")
  #learntLCGmmlcd<-studeny_rules(tg.patmmlcd)
  #endmmlcd1 = proc.time() -startmmlcd1
  #endmmlcd2 = proc.time() -startmmlcd2
  #time1_mmlcd = endmmlcd1[[3]][[1]]
  #time2_mmlcd = endmmlcd2[[3]][[1]]  
  #mlcd方法(无向独立图出发)
  startmlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.05)
  tgminitri <- minimal_triangMAT(tgug)
  tg.ctree <- ug.to.jtree1(tgminitri)
  #startmlcd2 = proc.time()
  tg.patmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),500,0.05,"CG")
  learntLCGmlcd<-studeny_rules(tg.patmlcd)
  endmlcd1 = proc.time() -startmlcd1
  #endmlcd2 = proc.time() -startmlcd2
  time1_mlcd = endmlcd1[[3]][[1]]
  #time2_mlcd = endmlcd2[[3]][[1]]  
  #SPC方法
  startspc = proc.time()
  tg.patspc <- learn.lwf.norm(tgdata,0.05,method="stable",LCG=TRUE)
  endspc = proc.time() -startspc
  time_spc = endspc[[3]][[1]]
  #mbcsp方法
  startmbcsp = proc.time()
  mb=mbcsp(data.frame(tgdata),alpha = 0.05)
  tg.mpscp = learn.vstruct(data.frame(tgdata),mb,alpha = 0.05)
  endmbcsp = proc.time() -startmbcsp
  time_mbcsp = endmbcsp[[3]][[1]]
  #指标
  fnlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fn
  fplcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fp
  tprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$TPR
  fprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$FPR
  SHDlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$SHD
  prelcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$PRE
  f1lcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$F1
  acclcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$ACC
  #indexlcd <-c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd)  
  #fnmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fn
  #fpmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fp
  #tprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$TPR
  #fprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$FPR
  #SHDmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$SHD
  #premmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$PRE
  #f1mmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$F1
  #indexmmlcd <- c(fnmmlcd,fpmmlcd,tprmmlcd,fprmmlcd,SHDmmlcd,time_mmlcd)  
  fnmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fn
  fpmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fp
  tprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$TPR
  fprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$FPR
  SHDmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$SHD
  premlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$PRE
  f1mlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$F1
  accmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$ACC
  #indexmlcd <- c(fnmlcd,fpmlcd,tprmlcd,fprmlcd,SHDmlcd,time_mlcd)  
  fnspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fn
  fpspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fp
  tprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$TPR
  fprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$FPR
  SHDspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$SHD
  prespc <- comp.cgs(pattern(graph),tg.patspc$matrix)$PRE
  f1spc <- comp.cgs(pattern(graph),tg.patspc$matrix)$F1
  accspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$ACC
  #indexspc <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc) 
  fnmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fn
  fpmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fp
  tprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$TPR
  fprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$FPR
  SHDmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$SHD
  prembcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$PRE
  f1mbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$F1
  accmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$ACC
  #indexmbcsp <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc)
  index <- c(fnlcd,fplcd,tprlcd,fprlcd,prelcd,f1lcd,acclcd,SHDlcd,time1_lcd,
             fnmlcd,fpmlcd,tprmlcd,fprmlcd,premlcd,f1mlcd,accmlcd,SHDmlcd,time1_mlcd,
             fnspc,fpspc,tprspc,fprspc,prespc,f1spc,accspc,SHDspc,time_spc,
             fnmbcsp,fpmbcsp,tprmbcsp,fprmbcsp,prembcsp,f1mbcsp,accmbcsp,SHDmbcsp,time_mbcsp)
  #index1 <- c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd,fnmlcd1,fpmlcd1,tprmlcd1,fprmlcd1,SHDmlcd1,time_mlcd1)
  #return(index)
  return(index)
}

evalcg_500_0.01 = function(graph){
  tgdata <- rnorm.cg(500,graph,get.normal.dist(graph))
  #lcd方法
  startlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.01)
  tg.jtree <- ug.to.jtree(tgug)
  #startlcd2 = proc.time()
  tg.patlcd <- learn.mec.norm(tg.jtree,cov(tgdata),500,0.01,"CG")
  learntLCGlcd<-studeny_rules(tg.patlcd)
  endlcd1 = proc.time() -startlcd1
  #endlcd2 = proc.time() -startlcd2
  time1_lcd = endlcd1[[3]][[1]]
  #time2_lcd = endlcd2[[3]][[1]]
  #mmlcd方法(道德图出发)
  #startmmlcd1 = proc.time()
  # tgmor <- mor(as.data.frame(tgdata))
  #tgminitri <- minimal_triangMAT(tgmor)
  #tg.ctree <- ug.to.jtree1(tgminitri)
  #startmmlcd2 = proc.time()
  #tg.patmmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),10000,0.05,"CG")
  #learntLCGmmlcd<-studeny_rules(tg.patmmlcd)
  #endmmlcd1 = proc.time() -startmmlcd1
  #endmmlcd2 = proc.time() -startmmlcd2
  #time1_mmlcd = endmmlcd1[[3]][[1]]
  #time2_mmlcd = endmmlcd2[[3]][[1]]  
  #mlcd方法(无向独立图出发)
  startmlcd1 = proc.time()
  tgug <- naive.getug.norm(tgdata,0.01)
  tgminitri <- minimal_triangMAT(tgug)
  tg.ctree <- ug.to.jtree1(tgminitri)
  #startmlcd2 = proc.time()
  tg.patmlcd <- learn.mec.norm(tg.ctree,cov(tgdata),500,0.01,"CG")
  learntLCGmlcd<-studeny_rules(tg.patmlcd)
  endmlcd1 = proc.time() -startmlcd1
  #endmlcd2 = proc.time() -startmlcd2
  time1_mlcd = endmlcd1[[3]][[1]]
  #time2_mlcd = endmlcd2[[3]][[1]]  
  #SPC方法
  startspc = proc.time()
  tg.patspc <- learn.lwf.norm(tgdata,0.01,method="stable",LCG=TRUE)
  endspc = proc.time() -startspc
  time_spc = endspc[[3]][[1]]
  #mbcsp方法
  startmbcsp = proc.time()
  mb=mbcsp(data.frame(tgdata),alpha = 0.01)
  tg.mpscp = learn.vstruct(data.frame(tgdata),mb,alpha = 0.01)
  endmbcsp = proc.time() -startmbcsp
  time_mbcsp = endmbcsp[[3]][[1]]
  #指标
  fnlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fn
  fplcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$fp
  tprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$TPR
  fprlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$FPR
  SHDlcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$SHD
  prelcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$PRE
  f1lcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$F1
  acclcd <- comp.cgs(pattern(graph),learntLCGlcd$matrix)$ACC
  #indexlcd <-c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd)  
  #fnmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fn
  #fpmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$fp
  #tprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$TPR
  #fprmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$FPR
  #SHDmmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$SHD
  #premmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$PRE
  #f1mmlcd <- comp.cgs(pattern(graph),learntLCGmmlcd$matrix)$F1
  #indexmmlcd <- c(fnmmlcd,fpmmlcd,tprmmlcd,fprmmlcd,SHDmmlcd,time_mmlcd)  
  fnmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fn
  fpmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$fp
  tprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$TPR
  fprmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$FPR
  SHDmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$SHD
  premlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$PRE
  f1mlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$F1
  accmlcd <- comp.cgs(pattern(graph),learntLCGmlcd$matrix)$ACC
  #indexmlcd <- c(fnmlcd,fpmlcd,tprmlcd,fprmlcd,SHDmlcd,time_mlcd)  
  fnspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fn
  fpspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$fp
  tprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$TPR
  fprspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$FPR
  SHDspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$SHD
  prespc <- comp.cgs(pattern(graph),tg.patspc$matrix)$PRE
  f1spc <- comp.cgs(pattern(graph),tg.patspc$matrix)$F1
  accspc <- comp.cgs(pattern(graph),tg.patspc$matrix)$ACC
  #indexspc <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc) 
  fnmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fn
  fpmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$fp
  tprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$TPR
  fprmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$FPR
  SHDmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$SHD
  prembcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$PRE
  f1mbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$F1
  accmbcsp <- comp.cgs(pattern(graph),tg.mpscp$pattern)$ACC
  #indexmbcsp <- c(fnspc,fpspc,tprspc,fprspc,SHDspc,time_spc)
  index <- c(fnlcd,fplcd,tprlcd,fprlcd,prelcd,f1lcd,acclcd,SHDlcd,time1_lcd,
             fnmlcd,fpmlcd,tprmlcd,fprmlcd,premlcd,f1mlcd,accmlcd,SHDmlcd,time1_mlcd,
             fnspc,fpspc,tprspc,fprspc,prespc,f1spc,accspc,SHDspc,time_spc,
             fnmbcsp,fpmbcsp,tprmbcsp,fprmbcsp,prembcsp,f1mbcsp,accmbcsp,SHDmbcsp,time_mbcsp)
  #index1 <- c(fnlcd,fplcd,tprlcd,fprlcd,SHDlcd,time_lcd,fnmlcd1,fpmlcd1,tprmlcd1,fprmlcd1,SHDmlcd1,time_mlcd1)
  #return(index)
  return(index)
}

