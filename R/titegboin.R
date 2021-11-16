
#'@title Obtain the operating characteristics of the Time-to-event gBOIN design for late-onset toxicity grades by simulating trials.
#' @param target the target ETS
#' @param prob a matrix containing the true toxicity probabilities of each grade x the investigational dose levels.
#' @param ncohort the total number of cohorts
#' @param cohortsize the cohort size
#' @param maxt the maximum follow-up time 
#' @param prior.p a vector of length 3, which specifies the prior probability that the time to toxicity lies
#'              inside the time interval (0,\code{maxt}/3), (\code{maxt}/3,\code{2maxt}/3), (\code{2maxt}/3,1).
#'              The default value is \code{prior.p=c(1/3,1/3,1/3)}. 
#' @param accrual the accrual rate, i.e., the number of patients accrued in 1 unit of time,
#' @param maxpen the maximum proportion of pending patients
#' @param alpha a number from (0,1) that controls alpha*100% events in (0, 1/2T). 
#'              The default is \code{alpha=0.5}.             
#' @param startdose the starting dose level for the trial
#' @param cutoff.eli the cutoff to eliminate an overly toxic dose for safety.
#'                  We recommend the default value of (\code{cutoff.eli=0.95}) for general use.
#' @param ntrial the total number of trials to be simulated.
#' @param seed set seed
#' @details This function generates he operating characteristics of the Time-To-Event Bayesian Optimal 
#'          Interval design for toxicity grades(TITE-gBOIN)
#'          with delayed toxicities grades by simulating trials under the prespecified true ETS of the investigational doses. 
#' @return \code{get.oc.tite.gBOIN()} returns the operating characteristics of the TITE-gBOIN design as a data frame,
#'         including: (1) selection percentage at each dose level (\code{selpercent}),
#'         (2) the number of patients treated at each dose level (\code{nptsdose}),
#'         (3) the number of toxicities observed at each dose level (\code{ntoxdose}),
#'         (4) the average number of toxicities (\code{totaltox}),
#'         (5) the average number of patients (\code{totaln}),
#'         (6) the percentage of early stopping without selecting the MTD (\code{pctearlystop}).   
#'         (7) the average trial duration needed for the trial based on the TITE-gBOIN design (\code{duration})  
#'         (8)the standard deviation of trial duration needed for the trial based on the TITE-gBOIN design (\code{sdduration})
#' @import stat  MASS 
#' @export 

get.oc.tite.gBOIN <- function(target, prob, ncohort, cohortsize, maxt=1, prior.p=NA, accrual=3, maxpen=0.5, 
                              alpha=0.5, startdose = 1,  cutoff.eli = 0.95, ntrial = 1000, seed=seed)
{
  ### simple error checking
  if(!is.na(prior.p[1])){if(length(prior.p)!=3){cat("Error: The length of the prior probabilities should be 3! \n"); return();}}
  if(is.na(maxpen)){maxpen=0.5;}
  if(maxpen<0 || maxpen>0.65) {cat("Error: the value of maxpen should lie within (0,0.65)!  \n"); return();}
  
  select.mtd <- function(target, y, n, cutoff.eli=0.95)
  {
    ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
    pava <- function (x, wt = rep(1, length(x))) 
    {
      n <- length(x)
      if (n <= 1) 
        return(x)
      if (any(is.na(x)) || any(is.na(wt))) {
        stop("Missing values in 'x' or 'wt' not allowed")
      }
      lvlsets <- (1:n)
      repeat {
        viol <- (as.vector(diff(x)) < 0)
        if (!(any(viol))) 
          break
        i <- min((1:(n - 1))[viol])
        lvl1 <- lvlsets[i]
        lvl2 <- lvlsets[i + 1]
        ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
        x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
        lvlsets[ilvl] <- lvl1
      }
      x
    }
    ## determine whether the dose has been eliminated during the trial
    ndose=length(n);
    elimi=rep(0, ndose);
    for(i in 1:ndose)
    {
      if(n[i]>2) {if(y[i]>n[i]|1-pbeta(target, y[i]+1, n[i]-y[i]+1)>cutoff.eli) {elimi[i:ndose]=1; break;}}
    }
    
    if(elimi[1]==1) { selectdose=99; } ## no dose should be selected if the first dose is already very toxic
    else
    {
      nadmis = min(max(which(elimi==0)), max(which(n!=0))); ## the highest admissble (or un-eliminated) dose level
      ## poster mean and variance of toxicity probabilities using beta(0.005, 0.005) as the prior 
      phat = (y[1:nadmis]+0.005)/(n[1:nadmis]+0.01); 
      phat.var = (y[1:nadmis]+0.005)*(n[1:nadmis]-y[1:nadmis]+0.005)/((n[1:nadmis]+0.01)^2*(n[1:nadmis]+0.01+1))
      
      ## perform the isotonic transformation using PAVA
      phat = pava(phat, wt=1/phat.var) 
      phat = phat + (1:nadmis)*1E-10 ## break ties by adding an increasingly small number 
      selectdose = sort(abs(phat-target), index.return=T)$ix[1]  ## select dose closest to the target as the MTD
    }
    return(selectdose);  
  }
  
  gen.tite<-function(n, prob, alpha=0.5, Tobs=1)
  {
    ############ subroutines ############
    weib<-function(n, pi, pihalft)
    {
      ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
      alpha = log(log(1-pi)/log(1-pihalft))/log(2);
      lambda = -log(1-pi)/(Tobs^alpha);
      t = (-log(runif(n))/lambda)^(1/alpha);
      return(t);
    }
    
    tox = rep(0, n);
    t.tox = rep(0, n);
    grade <- rep(0, n)
    pi <- sum(prob[c(3,4)])
    
    pihalft = alpha*pi;  # alpha*100% event in (0, 1/2T)
    t.tox = weib(n, pi, pihalft);
    tox[t.tox<=Tobs]=1;    
    t.tox[tox==0]=Tobs;
    
    for (i in 1:n){
      if (tox[i]==0){
        grade[i] <- (1-rbinom(1,1,prob[1]/sum(prob[1:2])))*0.5
      }else {
        grade[i] <- (1-rbinom(1,1,prob[3]/sum(prob[3:4])))*0.5+1
        
      }
    }
    y = sum(grade);
    
    return(list(t.tox=t.tox,y=y,tox=tox,grade=grade))
  }
  
  
  set.seed(seed);
  if(is.na(prior.p[1])){prior.p = rep(1/3,3)}
  prior.p = prior.p/sum(prior.p)
  nprior=1; # the effective sample size for the Beta prior.
  priortox=target/2; # the prior mean for the Beta prior.
  ndose=ncol(prob);	
  p.true <- apply(prob*c(0,0,1,1),2,sum)
  Y = matrix(rep(0, ndose * ntrial), ncol = ndose);
  N = matrix(rep(0, ndose * ntrial), ncol = ndose);
  dselect = rep(0, ntrial);
  durationV = rep(0, ntrial);
  npendV = rep(0, ntrial);
  npts = ncohort*cohortsize;
  ntarget <- target/1.5 #normalized target ETS
  np.saf = ntarget*0.6; 
  np.tox = ntarget*1.4;
  lambda1 = log((1-np.saf)/(1-ntarget))/log(ntarget*(1-np.saf)/(np.saf*(1-ntarget)));
  lambda2 = log((1-ntarget)/(1-np.tox))/log(np.tox*(1-ntarget)/(ntarget*(1-np.tox)));
  signal.exists <- file.exists('sig.cmd');
  
  
  for(trial in 1:ntrial)
  {
    if (signal.exists & trial %% 500 == 1) shell('sig.cmd')
    y=NULL;  #toxicity indicator for each subject
    dv=NULL;  #dose for each subject
    n.d = rep(0, ndose);  # number of patient at each dose
    y.d = rep(0, ndose);  # number of toxicity at each dose
    t.enter=NULL; # time enter the study
    t.event=NULL; # time to event
    t.decision = 0; # decision making time
    d = startdose;  # current dose level
    earlystop = 0; #indicate if trial stops early
    elimi = rep(0, ndose)
    npend = 0;
    
    for(i in 1:ncohort)
    {
      # generate data for the new patient
      for(j in 1:cohortsize)
      {
        if(j==1) { t.enter = c(t.enter, t.decision); }
        else { 
          t.enter = c(t.enter, t.enter[length(t.enter)] + rexp(1, rate=accrual))
        }
      }
      obscohort = gen.tite(cohortsize, prob[,d], alpha=alpha,Tobs=maxt);
      t.event = c(t.event, obscohort$t.tox);
      y = c(y, obscohort$grade)
      dv = c(dv, rep(d, cohortsize));
      t.decision = t.enter[length(t.enter)];
      nobs=-1; pending=1;
      d.curr=d;
      npend = npend-1;
      
      while(pending==1)
      {
        npend = npend+1;
        pending = 0;
        if(i==ncohort) { t.decision = t.decision + maxt; }
        else {
          t.decision = t.decision + rexp(1, rate=accrual)
        }
        
        # determine which observation are observed
        delta = ((t.enter+t.event)<=t.decision);
        t = pmin(t.event, t.decision-t.enter, maxt);  ## used for recording potential censoring time
        cset = (dv==d);
        delta.curr = delta[cset];
        t.curr = t[cset];
        ntox.curr = sum((y[cset])[delta.curr==1]); 
        totalt = t.curr[delta.curr==0]
        totalt = 3*prior.p[1]*totalt*(totalt<=maxt/3)+
          ((prior.p[1]-prior.p[2])*maxt+3*prior.p[2]*totalt)*(maxt/3<totalt & totalt<=2*maxt/3)+
          ((prior.p[1]+prior.p[2]-2*prior.p[3])*maxt+3*prior.p[3]*totalt)*(2*maxt/3<totalt & totalt<=maxt)
        totalt = sum(totalt)/maxt
        n.curr = sum(cset);
        n.pend = sum(delta[cset]==0);
        nobs = sum(delta[cset]);
        m <- sum((y[cset])[delta.curr==1]==0) # #of patients finish assessment but with no DLT
        # determine which dose level should be eliminated
        for(dd in 1:ndose){
          cset1 = dv==dd;
          delta.curr1 = delta[cset1];
          ntox.curr1 = sum((y[cset1])[delta.curr1==1]) #
          n.curr1 = sum(cset1);
          if (n.curr1>0){
            nobs.curr = sum(delta.curr1)
            if (1-pbeta(ntarget, ntox.curr1/1.5+1, n.curr1-ntox.curr1/1.5+1)>cutoff.eli && nobs.curr>=3){ #################Fix here
              elimi[dd:ndose]=1;
              break;
            }
          }
        }
        #check whether extra safey rule should be applied
        
        if(d==1){
          if(1-pbeta(ntarget, ntox.curr/1.5+1, n.curr-ntox.curr/1.5+1)>cutoff.eli && nobs>=3) {
            earlystop = 1;   break;}
        }
        
        
        # check whether the current dose level should be eliminated
        if(elimi[d.curr]==1) {
          d=which(elimi==1)[1]-1
          if(d==0){earlystop = 1; break;}
          next;
        }
        
        #check whether the trial should be early terminated
        # if(n.curr>=n.earlystop){break;}
        #if(n.curr>=n.earlystop && sum(y[dv==d.curr])/sum(dv==d.curr) >lambda1 && sum(y[dv==d.curr])/sum(dv==d.curr)<lambda2){break;}
        
        #check whether the current dose is toxic based on observed ata
        if((ntox.curr/n.curr/1.5)>=lambda2){
          if(d==1){d=d; if(n.pend>0){pending=1};} else{d=d-1};
          next;
        }
        
        
        #check whether the trial should be suspended
        if(n.pend>n.curr*maxpen) {pending=1;}
        else
        {
          #compute the estimate of toxicity rate
          if(n.pend==0){
            phat=ntox.curr/n.curr/1.5
          } else {
            #compute the estimated toxicity rate based on the imputed data
            phat = ntox.curr/1.5/((ntox.curr/1.5)+m+totalt)
          }
          
          
          
          #make dose assignment decisions
          if(phat<=lambda1 & (ntox.curr/n.curr/1.5)<=ntarget){
            if(n.pend==n.curr){
              pending=1;
            } else {
              #check whether the current dose is the highest
              if(d==ndose){d=d;} else{
                if(elimi[d+1]==1){d=d} else{ d=d+1}
              }}} else if (phat>=lambda2 & (ntox.curr/n.curr/1.5)>=ntarget){
                if(d==1){d=d; if(n.pend>0){pending=1}} else{d=d-1}
              } else {d=d}
          if(elimi[d]==1){d=which(elimi==1)[1]-1}
        }
      }
      if(earlystop==1){t.decision = max(t.enter)+maxt; break;}
      if(earlystop==2){t.decision = max(t.enter)+maxt; break;}
    }
    
    for(k in 1:ndose){
      y.d[k] = sum(y[dv==k]);
      n.d[k] = sum(dv==k);
    }
    
    npendV[trial]= npend;
    Y[trial, ] = y.d
    N[trial, ] = n.d
    durationV[trial] = t.decision
    if (earlystop == 1) {
      dselect[trial] = 99
    }
    else { dselect[trial] = select.mtd(ntarget, y.d/1.5, n.d, cutoff.eli)}
  }
  
  selpercent = rep(0, ndose)
  selpercent=rep(0, ndose);
  nptsdose = apply(N,2,mean);
  ntoxdose = apply(Y,2,mean);
  
  for(i in 1:ndose) { selpercent[i]=sum(dselect==i)/ntrial*100; }
  
  if(length(which(p.true==target))>0) # if MTD exists, calculate risk of overdosing
  {
    if (which(p.true==target) == ndose-1) {
      overdosing60=mean(N[,p.true>target]>0.6*npts)*100;
      overdosing80=mean(N[,p.true>target]>0.8*npts)*100;
    } else {
      if(ntrial>1){
        overdosing60=mean(rowSums(N[,p.true>target])>0.6*npts)*100;
        overdosing80=mean(rowSums(N[,p.true>target])>0.8*npts)*100;
      } else {
        overdosing60=(sum(N[p.true>target])>0.6*npts)*100;
        overdosing80=(sum(N[p.true>target])>0.8*npts)*100;
      }
    }
    
    out=list(selpercent=selpercent, npatients=nptsdose, ntox=ntoxdose, totaltox=sum(Y)/ntrial, totaln=sum(N)/ntrial,
             percentstop=sum(dselect== 99)/ntrial*100, # poorallocation=mean(N[, p.true==target]<npts/ndose)*100,
             overdose60=overdosing60, overdose80=overdosing80, duration=mean(durationV),sdduration=sqrt(var(durationV)));
  }else {
    out=list(selpercent=selpercent, npatients=nptsdose, ntox=ntoxdose, totaltox=sum(Y)/ntrial, totaln=sum(N)/ntrial,
             percentstop=sum(dselect== 99)/ntrial*100, duration=mean(durationV),sdduration=sqrt(var(durationV)))
  }
  return(out);
  
}
