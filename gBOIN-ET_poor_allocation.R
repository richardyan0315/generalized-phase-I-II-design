rm(list=ls())

scenario_gBOINET_fix <- readRDS("./Outputs/scenario_gBOINET_fix.rds")
scenario_gBOINET_all <- readRDS("./Outputs/scenario_gBOINET_all.rds")

OBD_true_random_HE <- readRDS("./Outputs/OBD_true_random_HE.rds") # true OBD for all 729 random scenarios (HE) via scenario generator (line 595)
OBD_true_random_LT <- readRDS("./Outputs/OBD_true_random_LT.rds") # true OBD for all 729 random scenarios (LT) via scenario generator (line 595)

# the manually updated gboinet function added poor allocation rate calculation
gboinet_with_poor_allocation <- function(
    n.dose, start.dose, size.cohort, n.cohort,
    toxprob, effprob, sev.weight, res.weight,
    phi, phi1=phi*0.1, phi2=phi*1.4, delta, delta1=delta*0.6,
    alpha.T1=0.5, alpha.E1=0.5, tau.T, tau.E,
    te.corr=0.2, gen.event.time="weibull",
    accrual, gen.enroll.time="uniform",
    stopping.npts=size.cohort*n.cohort,
    stopping.prob.T=0.95, stopping.prob.E=0.99,
    estpt.method, obd.method,
    w1=0.33, w2=1.09,
    plow.ast=phi1, pupp.ast=phi2, qlow.ast=delta1/2, qupp.ast=delta,
    psi00=40, psi11=60,
    n.sim=1000, seed.sim=100)
{
  if(ncol(toxprob)!=n.dose){
    stop("Number of dose must be the same as the length of true toxicity probability.")
  }else if(ncol(effprob)!=n.dose){
    stop("Number of dose must be the same as the length of true efficacy probability.")
  }else if(nrow(toxprob)!=length(sev.weight)){
    stop("Number of toxicity category must be the same as the length of weight.")
  }else if(nrow(effprob)!=length(res.weight)){
    stop("Number of efficacy category must be the same as the length of weight.")
  }else if(!((phi1<phi)&(phi<phi2))){
    stop("Design parameters must satisfy a condition of phi1 < phi < phi2.")
  }else if(!(delta1<delta)){
    stop("Design parameters must satisfy a condition of delta1 < delta.")
  }else if((!is.matrix(toxprob))|(!is.matrix(effprob))){
    stop("True toxicity and efficacy probability must be specified as matrix.")
  }else{
    
    dosen <- 1:n.dose
    dose  <- paste("Dose", dosen, sep = "")
    toxp <- data.frame(toxprob)
    colnames(toxp) <- dose
    effp <- data.frame(effprob)
    colnames(effp) <- dose
    ncat.T <- nrow(toxp)-1
    ncat.E <- nrow(effp)-1
    ncoh <- size.cohort
    nesc <- n.cohort
    nmax <- ncoh*nesc
    design.par <- gridoptim(phi=phi,phi1=phi1,phi2=phi2,delta=delta,delta1=delta1)
    lambda1 <- design.par$lambda1
    lambda2 <- design.par$lambda2
    eta1    <- design.par$eta1
    pr.alpha <- 1
    pr.beta  <- 1
    alpha.T2 <- 0.5
    alpha.E2 <- 0.5
    efftoxp <- list(toxp=toxp,effp=effp)
    ncop    <- copula::normalCopula(te.corr,dim=2,dispstr="ex")
    mv.ncop <- NULL
    
    if(gen.event.time=="weibull"){
      for(i in 1:n.dose){
        psi.T    <- sum(efftoxp$toxp[-1,i])
        zetta.T1 <- log(log(1-psi.T)/log(1-psi.T+alpha.T1*psi.T))/log(1/(1-alpha.T2))
        zetta.T2 <- tau.T/(-log(1-psi.T))^(1/zetta.T1)
        psi.E    <- sum(efftoxp$effp[-1,i])
        zetta.E1 <- log(log(1-psi.E)/log(1-psi.E+alpha.E1*psi.E))/log(1/(1-alpha.E2))
        zetta.E2 <- tau.E/(-log(1-psi.E))^(1/zetta.E1)
        mv.ncop <- append(mv.ncop, copula::mvdc(copula = ncop,
                                                margins = c("weibull","weibull"),
                                                paramMargins = list(list(shape=zetta.T1,scale=zetta.T2),
                                                                    list(shape=zetta.E1,scale=zetta.E2))))
      }
    }else if(gen.event.time=="uniform"){
      for(i in 1:n.dose){
        psi.T <- sum(efftoxp$toxp[-1,i])
        psi.E <- sum(efftoxp$effp[-1,i])
        mv.ncop <- append(mv.ncop, copula::mvdc(copula = ncop,
                                                margins = c("unif","unif"),
                                                paramMargins = list(list(min=0,max=tau.T*(1/psi.T)),
                                                                    list(min=0,max=tau.E*(1/psi.E)))))
      }
    }
    
    data.obs.n <- array(0, dim=c(n.sim, n.dose))
    data.dur   <- array(0, dim=c(n.sim))
    obd <- array(0, dim=c(n.sim))
    poor_allocation <- numeric(n.sim)
    
    set.seed(seed.sim)
    
    for(ss in 1:n.sim){
      obs.n   <- numeric(n.dose)
      obs.tox <- numeric(n.dose)
      obs.eff <- numeric(n.dose)
      pe      <- numeric(n.dose)
      pt      <- numeric(n.dose)
      t.enter <- NULL
      t.decision <- 0
      curdose <- start.dose
      early.stop <- 0
      
      for(i in 1:nesc){
        dlab <- paste("Dose",curdose,sep="")
        obs.n[curdose] <- obs.n[curdose] + ncoh
        
        for(j in 1:ncoh){
          if(j==1){
            t.enter <- c(t.enter, t.decision)
          }else{
            if(gen.enroll.time=="uniform"){
              t.enter <- c(t.enter, t.enter[length(t.enter)]+runif(1, 0, 2*accrual))
            }else if(gen.enroll.time=="exponential"){
              t.enter <- c(t.enter, t.enter[length(t.enter)]+rexp(1,1/accrual))
            }
          }
        }
        t.decision <- t.enter[length(t.enter)]+max(tau.T, tau.E)
        
        time.te <- copula::rMvdc(ncoh, mv.ncop[[curdose]])
        event.T  <- as.numeric(time.te[,1] <= tau.T)
        grade    <- event.T * ((1:ncat.T) %*% rmultinom(ncat.T, 1, efftoxp$toxp[-1, dlab])) + 1
        ETS      <- apply(grade, 2, function(x) {return(sev.weight[x])})
        nETS     <- ETS / max(sev.weight)
        event.E  <- as.numeric(time.te[,2] <= tau.E)
        response <- event.E * ((1:ncat.E) %*% rmultinom(ncat.E, 1, efftoxp$effp[-1, dlab])) + 1
        EES      <- apply(response, 2, function(x) {return(res.weight[x])})
        nEES     <- EES / max(res.weight)
        obs.tox[curdose] <- obs.tox[curdose] + sum(nETS)
        pt[curdose] <- obs.tox[curdose] / obs.n[curdose]
        obs.eff[curdose] <- obs.eff[curdose] + sum(nEES)
        pe[curdose] <- obs.eff[curdose] / obs.n[curdose]
        
        if((pt[curdose]<=lambda1)&(pe[curdose]<=eta1)){
          nxtdose <- curdose+1
        }else if((pt[curdose]<lambda2)&(pe[curdose]>eta1)){
          nxtdose <- curdose
        }else if(pt[curdose]>=lambda2){
          nxtdose <- curdose-1
        }else if((pt[curdose]>lambda1)&(pt[curdose]<lambda2)&(pe[curdose]<=eta1)){
          if(curdose==n.dose){
            three   <- c(curdose-1,curdose)
            maxpe   <- max(pe[three])
            nxtdose <- sample(dosen[which((pe==maxpe)&(is.element(dosen,three)))],1)
          }else if(obs.n[curdose+1]==0){
            nxtdose <- curdose+1
          }else if(curdose==1){
            three   <- c(curdose,curdose+1)
            maxpe   <- max(pe[three])
            nxtdose <- sample(dosen[which((pe==maxpe)&(is.element(dosen,three)))],1)
          }else{
            three   <- c(curdose-1,curdose,curdose+1)
            maxpe   <- max(pe[three])
            nxtdose <- sample(dosen[which((pe==maxpe)&(is.element(dosen,three)))],1)
          }
        }
        
        po.shape1 <- pr.alpha + obs.tox
        po.shape2 <- pr.beta  + (obs.n - obs.tox)
        tterm     <- pbeta(phi, po.shape1, po.shape2)
        po.shape1 <- pr.alpha + obs.eff
        po.shape2 <- pr.beta  + (obs.n - obs.eff)
        eterm     <- 1 - pbeta(delta1, po.shape1, po.shape2)
        admflg  <- !((eterm < (1 - stopping.prob.E)) | (tterm < (1 - stopping.prob.T)))
        admdose <- dosen[admflg]
        
        if(sum(admflg)==0){
          early.stop <- 1
          break
        }else if(sum(obs.n >= stopping.npts) > 0){
          break
        }else{
          if(nxtdose == 0){
            if(admflg[1]){
              curdose <- 1
            }else{
              early.stop <- 1
              break
            }
          }else if(nxtdose == (n.dose + 1)){
            curdose <- n.dose
          }else if(is.element(nxtdose, admdose)){
            curdose <- nxtdose
          }else if(curdose < nxtdose){
            if(sum(admdose >= nxtdose) != 0){
              curdose <- min(admdose[admdose >= nxtdose])
            }
          }else if(curdose >= nxtdose){
            if(sum(admdose <= nxtdose) != 0){
              curdose <- max(admdose[admdose <= nxtdose])
            }else{
              early.stop <- 1
              break
            }
          }
        }
      }
      
      data.obs.n[ss, ] <- obs.n
      data.dur[ss] <- t.decision
      evadose <- dosen[obs.n != 0]
      obspt <- obs.tox[evadose] / obs.n[evadose]
      obspe <- obs.eff[evadose] / obs.n[evadose]
      tterm.obd <- numeric(n.dose)
      eterm.obd <- numeric(n.dose)
      
      for(i in evadose){
        po.shape1 <- pr.alpha + obs.tox[i]
        po.shape2 <- pr.beta + (obs.n[i] - obs.tox[i])
        tterm.obd[i] <- pbeta(phi, po.shape1, po.shape2)
        po.shape1 <- pr.alpha + obs.eff[i]
        po.shape2 <- pr.beta + (obs.n[i] - obs.eff[i])
        eterm.obd[i] <- 1 - pbeta(delta1, po.shape1, po.shape2)
      }
      
      if(early.stop == 1){
        obd[ss] <- 0
      }else if(length(evadose) == 1){
        if((tterm.obd[evadose] >= (1 - stopping.prob.T)) & (eterm.obd[evadose] >= (1 - stopping.prob.E))){
          obd[ss] <- evadose
        }
      }else if(sum((tterm.obd[evadose] >= (1 - stopping.prob.T)) & (eterm.obd[evadose] >= (1 - stopping.prob.E))) >= 1){
        estpt <- Iso::pava(obspt)
        if(estpt.method == "multi.iso"){
          estpe <- multi.iso(obs = obs.eff[evadose], n = obs.n[evadose])
        }else if(estpt.method == "fp.logistic"){
          estpe <- fp.logit(obs = obs.eff[evadose], n = obs.n[evadose], dose = evadose)
        }else if(estpt.method == "obs.prob"){
          estpe <- obspe
        }
        obd[ss] <- obd.select(
          probt = estpt, probe = estpe, method = obd.method,
          phi = phi, phi1 = phi1, phi2 = phi2, delta = delta, delta1 = delta1,
          tterm = tterm.obd[evadose], eterm = eterm.obd[evadose],
          stopT = stopping.prob.T, stopE = stopping.prob.E,
          w1 = w1, w2 = w2,
          plow.ast = plow.ast, pupp.ast = pupp.ast, qlow.ast = qlow.ast, qupp.ast = qupp.ast,
          psi00 = psi00, psi11 = psi11)
      }
      
      if (obd[ss] > 0) {
        allocated_patients <- obs.n[obd[ss]]
        total_patients <- sum(obs.n)
        allocation_rate <- allocated_patients / total_patients
        poor_allocation[ss] <- ifelse(allocation_rate < 0.2, 1, 0)
      }
    }
    
    prop.select <- array(0, dim=c(n.dose))
    for(i in 1:n.dose){
      prop.select[i] <- round(mean(obd == i) * 100, digits = 1)
    }
    names(prop.select) <- dose
    
    prop.stop <- round(mean(obd == 0) * 100, digits = 1)
    names(prop.stop) <- "Stop %"
    
    n.patient <- round(apply(data.obs.n, 2, mean), digits = 1)
    names(n.patient) <- dose
    
    duration  <- round(mean(data.dur), digits = 1)
    names(duration) <- "Trial duration (days)"
    
    t.nets <- round(apply(apply(toxprob, 2, function(x) x * sev.weight / max(sev.weight)), 2, sum), digits = 2)
    t.nees <- round(apply(apply(effprob, 2, function(x) x * res.weight / max(res.weight)), 2, sum), digits = 2)
    
    dimnames(toxprob) <- list(paste("Tox.cat", 1:(ncat.T + 1), sep = ""), dose)
    dimnames(effprob) <- list(paste("Eff.cat", 1:(ncat.E + 1), sep = ""), dose)
    names(phi) <- "Target toxicity prob."
    names(delta) <- "Target efficacy prob."
    names(lambda1) <- "Lower toxicity boundary"
    names(lambda2) <- "Upper toxicity boundary"
    names(eta1) <- "Lower efficacy boundary"
    names(tau.T) <- "Tox. assessment window (days)"
    names(tau.E) <- "Eff. assessment window (days)"
    names(accrual) <- "Accrual rate (days)"
    names(ncat.T) <- "Number of toxicity category"
    names(ncat.E) <- "Number of efficacy category"
    names(estpt.method) <- "Efficacy prob. estimation"
    names(obd.method) <- "OBD selection"
    
    poor_allocation_rate <- mean(poor_allocation)
    
    result <- list(toxprob = toxprob,
                   effprob = effprob,
                   nETS = t.nets,
                   nEES = t.nees,
                   phi = phi,
                   delta = delta,
                   lambda1 = lambda1,
                   lambda2 = lambda2,
                   eta1 = eta1,
                   tau.T = tau.T,
                   tau.E = tau.E,
                   accrual = accrual,
                   ncat.T = ncat.T + 1,
                   ncat.E = ncat.E + 1,
                   estpt.method = estpt.method,
                   obd.method = obd.method,
                   n.patient = n.patient,
                   prop.select = prop.select,
                   prop.stop = prop.stop,
                   duration = duration,
                   poor_allocation_rate = poor_allocation_rate)
    
    class(result) <- "gboinet"
    result
  }
}

# parameter settings
{
  n.dose      <- 6                 # total dose levels
  start.dose  <- 1                 # initial dose level
  size.cohort <- 3                 # cohort size
  n.cohort    <- 9                 # cohort number
  n.sim       <- 100               # simulation replicates
  
  estpt.method <- "obs.prob"    
  obd.method   <- "utility.scoring"  # OBD selection method
  
  phi   <- 0.30  
  delta <- 0.50  
  
  tau.T   <- 30     
  tau.E   <- 45     
  accrual <- 10    
  
  psi00_HE <- 30  # w_t (30 or 55)
  psi11_HE <- 70  # w_e (70 or 45)
  
  psi00_LT <- 55  # w_t (30 or 55)
  psi11_LT <- 45  # w_e (70 or 45)
  
  sev.weight.A <- c(0.90, 1.30, 1.80, 2.30)  # toxicity A # weight sum up for grade 0 and 1
  sev.weight.B <- c(0.57, 1.00, 1.50, 2.00)  # toxicity B 
  sev.weight.C <- c(0.07, 0.25, 0.50, 1.00)  # toxicity C 
  
  res.weight.1 <- c(0.05, 0.20, 0.335, 0.47)  # efficacy 1, the third as the average of second and fourth
  res.weight.2 <- c(0.07, 0.14, 0.275, 0.41)  # efficacy 2 
  res.weight.3 <- c(0.03, 0.15, 0.300, 0.45)  # efficacy 3 
  
  sev.weight.merged <- sev.weight.A + sev.weight.B + sev.weight.C
  
  res.weight.merged <- res.weight.1 + res.weight.2 + res.weight.3
  
  merge_rows <- function(mat) { new_mat <- rbind(colSums(mat[1:2, ]), mat[-c(1, 2), ]); return(new_mat)}
  split_second_row <- function(mat) { new_mat <- rbind(mat[1, ], mat[2, ] / 2, mat[2, ] / 2, mat[-c(1, 2), ]); return(new_mat)}
}





# ---- Fix 
results_fix_HE <- vector("list", 9)
results_fix_LT <- vector("list", 9)

set.seed(315) # set the same seed to ensure same outputs as boinet::gboinet
for (i in 1:9) {
  toxprobA <- merge_rows(t(scenario_gBOINET_fix[[i]]$TOX_A))
  toxprobB <- merge_rows(t(scenario_gBOINET_fix[[i]]$TOX_B))
  toxprobC <- merge_rows(t(scenario_gBOINET_fix[[i]]$TOX_C))
  
  effprob1 <- split_second_row(t(scenario_gBOINET_fix[[i]]$EFF_1))
  effprob2 <- split_second_row(t(scenario_gBOINET_fix[[i]]$EFF_2))
  effprob3 <- split_second_row(t(scenario_gBOINET_fix[[i]]$EFF_3))
  
  # merged toxicity (Integrate information on multiple toxicities without discrimination)
  toxprob_merged <- round((toxprobA + toxprobB + toxprobC) / 3, 3)
  
  # merged efficacy (Integrate information on multiple efficacies without discrimination)
  effprob_merged <- round((effprob1 + effprob2 + effprob3) / 3, 3)
  
  results_fix_HE[[i]] <- gboinet_with_poor_allocation(
    n.dose = n.dose, 
    start.dose = start.dose,
    size.cohort = size.cohort, 
    n.cohort = n.cohort,
    
    toxprob = toxprob_merged, 
    effprob = effprob_merged,
    
    sev.weight = sev.weight.merged, 
    res.weight = res.weight.merged,
    
    phi = phi, 
    delta = delta,
    
    tau.T = tau.T,
    tau.E = tau.E,
    accrual = accrual,
    
    estpt.method = estpt.method, 
    obd.method = obd.method,
    
    n.sim = n.sim,
    psi00 = psi00_HE, # w_t
    psi11 = psi11_HE  # w_e
  )
  
  results_fix_LT[[i]] <- gboinet_with_poor_allocation(
    n.dose = n.dose, 
    start.dose = start.dose,
    size.cohort = size.cohort, 
    n.cohort = n.cohort,
    
    toxprob = toxprob_merged, 
    effprob = effprob_merged,
    
    sev.weight = sev.weight.merged, 
    res.weight = res.weight.merged,
    
    phi = phi, 
    delta = delta,
    
    tau.T = tau.T,
    tau.E = tau.E,
    accrual = accrual,
    
    estpt.method = estpt.method, 
    obd.method = obd.method,
    
    n.sim = n.sim,
    psi00 = psi00_LT, # w_t
    psi11 = psi11_LT  # w_e
  )
  
  print(paste0(" --- Results for Scenario ", i, " at HE --- "))
  print(results_fix_HE[[i]]$poor_allocation_rate)

  print(paste0(" --- Results for Scenario ", i, " at LT --- "))
  print(results_fix_LT[[i]]$poor_allocation_rate)

}

saveRDS(results_fix_HE, file = "./Outputs/gBOIN-ET_poorAllocation_fix_HE.rds")
saveRDS(results_fix_LT, file = "./Outputs/gBOIN-ET_poorAllocation_fix_LT.rds")





# ---- Random (All)
results_all_HE <- vector("list", 729)
results_all_LT <- vector("list", 729)

set.seed(315)
for (i in 1:729) {
  toxprobA <- merge_rows(t(scenario_gBOINET_all[[i]]$TOX_A))
  toxprobB <- merge_rows(t(scenario_gBOINET_all[[i]]$TOX_B))
  toxprobC <- merge_rows(t(scenario_gBOINET_all[[i]]$TOX_C))
  
  effprob1 <- split_second_row(t(scenario_gBOINET_all[[i]]$EFF_1))
  effprob2 <- split_second_row(t(scenario_gBOINET_all[[i]]$EFF_2))
  effprob3 <- split_second_row(t(scenario_gBOINET_all[[i]]$EFF_3))
  
  # merged toxicity (Integrate information on multiple toxicities without discrimination)
  toxprob_merged <- round((toxprobA + toxprobB + toxprobC) / 3, 3)
  
  # merged efficacy (Integrate information on multiple efficacies without discrimination)
  effprob_merged <- round((effprob1 + effprob2 + effprob3) / 3, 3)
  
  results_all_HE[[i]] <- gboinet_with_poor_allocation(
    n.dose = n.dose, 
    start.dose = start.dose,
    size.cohort = size.cohort, 
    n.cohort = n.cohort,
    
    toxprob = toxprob_merged, 
    effprob = effprob_merged,
    
    sev.weight = sev.weight.merged, 
    res.weight = res.weight.merged,
    
    phi = phi, 
    delta = delta,
    
    tau.T = tau.T,
    tau.E = tau.E,
    accrual = accrual,
    
    estpt.method = estpt.method, 
    obd.method = obd.method,
    
    n.sim = n.sim,
    psi00 = psi00_HE, # w_t
    psi11 = psi11_HE  # w_e
  )
  
  results_all_LT[[i]] <- gboinet_with_poor_allocation(
    n.dose = n.dose, 
    start.dose = start.dose,
    size.cohort = size.cohort, 
    n.cohort = n.cohort,
    
    toxprob = toxprob_merged, 
    effprob = effprob_merged,
    
    sev.weight = sev.weight.merged, 
    res.weight = res.weight.merged,
    
    phi = phi, 
    delta = delta,
    
    tau.T = tau.T,
    tau.E = tau.E,
    accrual = accrual,
    
    estpt.method = estpt.method, 
    obd.method = obd.method,
    
    n.sim = n.sim,
    psi00 = psi00_LT, # w_t
    psi11 = psi11_LT  # w_e
  )
  
  print(paste0(" --- Results for Scenario ", i, " at HE --- "))
  print(results_all_HE[[i]]$poor_allocation_rate)
  
  print(paste0(" --- Results for Scenario ", i, " at LT --- "))
  print(results_all_LT[[i]]$poor_allocation_rate)

}


saveRDS(results_all_HE, file = "./Outputs/gBOIN-ET_poorAllocation_all_HE.rds")
saveRDS(results_all_LT, file = "./Outputs/gBOIN-ET_poorAllocation_all_LT.rds")


