##################################################################### SIMULATION STUDY TO COMPARE DIFFERENT ANALYSIS METHODS TO ASSESS THE TREATMENT EFFECT ###############################################################
library(survival);library(MASS);library(snow);library(MonteCarlo);library(Matrix);library(geeM);library(sandwich);library(lmtest);library(pracma)
##################################################################### AUXILIARY FUNCTION FOR THE DATA GENERATING MECHANISM: CREATING EVENT TIMES ##########################################################################
create_event_times<-function(n,A,L1,max_fup,bas_haz_out,coef_trt_out,shape){
  frailty<-rgamma(n,shape=shape,scale=1/shape) ## GENERATE FRAILTY SO THAT THE NUMBER OF PAST INFECTIONS IS PREDICTIVE FOR FUTURE INFECTIONS (GAMMA DISTRIBUTION WITH 1 PARAMETER).
  event_times<-lapply(1:n,function(sigma){
    vector_of_event_times_per_individual<-c();cumulative_event_time<-0 ## THE EVENT TIMES WILL BE STORED IN "vector_of_event_times_per_individual". Whilst the generated event times themselves are not of cumulative nature, the event times stored in this vector are cumulative!
    while(cumulative_event_time<max_fup){ ## DRAW EVENT TIMES AS LONG AS THE CUMULATIVE EVENT TIME IS SMALLER THAN THE MAXIMUM FOLLOW-UP TIME.
      y<-rexp(1,rate=exp(bas_haz_out+coef_trt_out*A[sigma]+0.5*L1[sigma])*frailty[sigma]) ## Y_t|X FOLLOWS A NEGATIVE BINOMIAL DISTRIBUTION, AT EACH TIME POINT t. (AT TIME POINT t, WE EXPECT THAT Y_t IS POISSON DISTRIBUTED WITH PARAMETER exp{\beta'X+log(t)}*frailty)
      aux_var<-length(vector_of_event_times_per_individual) ## Auxiliary covariate
      if(y<(max_fup-ifelse(aux_var!=0,vector_of_event_times_per_individual[aux_var],0))){ ## CHECK WHETHER THERE IS STILL "SPACE" FOR AN EVENT.
        vector_of_event_times_per_individual<-c(vector_of_event_times_per_individual,y+ifelse(aux_var!=0,vector_of_event_times_per_individual[aux_var],0))}
      cumulative_event_time<-cumulative_event_time+y}
    return(vector_of_event_times_per_individual)})
  return(event_times)}
##################################################################### transformation FUNCTION FOR CENSORING MODEL ##############################################################################################################
transformation<-function(x,delta,sigmoid_type){if(sigmoid_type==0){return(erf(delta*x*sqrt(pi)/2))}else{return(delta*x/(1+delta*x))}}
##################################################################### SIMULATE CENSORING TIMES BASED ON THE TREATMENT, BASELINE COVARIATES IN OUTCOME MODEL AND n.infect{t-}! #############################################
simulate_censoring_times<-function(n,A,L1,max_fup,event_times,bas_haz_cens,censoring_strength_trt,censoring_strength_bas_cov,censoring_strength_time_varying_covariate,cens_interaction,delta,sigmoid_type){
  ## Let t_{i,1},t_{i,2},...,t_{i,j} be the event times for the i-th individual.
  ## Complete event_times by the maximum follow-up time.
  event_times<-lapply(1:n,function(sigma){if(length(event_times[[sigma]])>0){return(c(event_times[[sigma]],max_fup))}else{return(max_fup)}})
  ## As we wish to include y_{t-} in the model for ct, we must update the rate when an individual crosses a time point at which s/he experiences an event.
  return(sapply(1:n,function(sigma){
    j<-0;cumulative_time<-0
    while(j<=(length(event_times[[sigma]])-1)){
      current_time<-rexp(1,rate=exp(bas_haz_cens+censoring_strength_trt*A[sigma]+censoring_strength_bas_cov*L1[sigma]+censoring_strength_time_varying_covariate*transformation(j,delta,sigmoid_type)+cens_interaction*A[sigma]*transformation(j,delta,sigmoid_type))) ## Probability of getting censored raises with time as the time-varying rate also raises with time.
      cumulative_time<-cumulative_time+current_time
      if(cumulative_time>event_times[[sigma]][j+1]){ ## Individual sigma is uncensored at the (j+1)-th event time.
        cumulative_time<-event_times[[sigma]][j+1]
      }else{ ## Individual is censored before the j-th event time.
        return(cumulative_time)}
      j<-j+1}
    return(cumulative_time)}))} ## The censoring time of an individual who passes the entire while statement equals the follow-up time; and should still be returned at the end.
##################################################################### RESTRICT THE EVENT TIMES TO THOSE WHICH LIE BEFORE THE CENSORING TIMES ##############################################################################
restrict_event_times_to_events_before_censoring_times<-function(n,event_times,censoring_times){
  event_times<-lapply(1:n,function(sigma){
    ## Restrict event times to those before the censoring time of each individual.
    aux_var<-findInterval(censoring_times[sigma],event_times[[sigma]])
    if(aux_var>0){return(event_times[[sigma]][1:aux_var])}else{return(NULL)}})
  return(list(event_times,censoring_times))}
##################################################################### CREATE A DATASET WHICH STORES FOR EACH INDIVIDUAL THE FIRST, SECOND, THIRD,... EVENT TIMES ##########################################################
create_dataset_with_covariates_representing_at_which_time_point_the_q_th_event_took_place<-function(n,list,max_fup){
  ## INDIVIDUALS WHO REACH THE PRESPECIFIED FOLLOW-UP TIME ARE ADMINISTRATIVELY CENSORED; THEIR TRUE CENSORING TIME IS HENCE UNKNOWN.
  list_event_times<-list[[1]];fup<-list[[2]];cens.time<-sapply(1:n,function(sigma){if(fup[sigma]!=max_fup){return(fup[sigma])}else{return(NA)}})
  number_of_covariates<-max(sapply(1:n,function(sigma){return(length(list_event_times[[sigma]]))})) ## DETERMINE MAXIMUM NUMBER OF EVENT TIMES ACROSS ALL INDIVIDUALS.
  covariate_names<-sapply(1:number_of_covariates,function(rho){paste("event.dt.",rho,sep="")}) ## CREATE COLUMN NAMES FOR EACH EVENT NUMBER.
  list_with_event_times<-lapply(1:number_of_covariates,function(kappa){
    return(sapply(1:n,function(lambda){if(kappa<=length(list_event_times[[lambda]])){return(list_event_times[[lambda]][kappa])}else{return(NA)}}))})
  data<-data.frame(cbind(1:n,data.frame(do.call("cbind",list_with_event_times)),cens.time,fup));colnames(data)<-c("id",covariate_names,"cens.time","fup")
  return(data)}
##################################################################### AUXILIARY FUNCTION TO DETERMINE AT WHICH ROWS THE IDENTIFIER VARIABLE CHANGES #######################################################################
determine_row_indices_at_which_identifier_variable_changes<-function(data.long){
  return(c(which(data.long$tstart==0),nrow(data.long)+1))}
##################################################################### TURN DATA INTO THE DESIRED FORMAT SO THAT SURVSPLIT CAN BE APPLIED TO IT ############################################################################
turn_dataset_in_format_so_that_survsplit_can_be_used<-function(n,data.baseline,data.event.time.covariates){
  number.of.event.time.covariates<-ncol(data.event.time.covariates)-3 ## data.event.time.covariates has also an id, censoring time and follow-up time column.
  new.data<-tmerge(data.baseline,data.event.time.covariates,id=id,tstop=fup) ## Add tstart and tstop. (tstop corresponds to the follow-up time).
  ## Split each individual's time window in intervals [tstart,t1], [t1,t2],...[tk,tstop], where each endpoint corresponds to a point in time at which an individual experiences an event.
  for(l in 2:(number.of.event.time.covariates+1)){
    cov<-names(data.event.time.covariates)[l]
    new.data<-tmerge(new.data,data.event.time.covariates,id=id,infect=event(data.event.time.covariates[,cov]))}
  new.data<-tmerge(new.data,data.event.time.covariates,id=id,cens=event(cens.time))
  ## The stopping times are the time points at which an individual experiences an event, or at which s/he gets censored, or at which s/he reaches the prespecified follow-up time.
  cutting_times<-sort(unique(new.data$tstop)) ## Time points at which an individual gets censored, or at which his/her time-varying covariates change, or at which she reaches the pre-specified follow-up time.
  ## Convert new.data to a long data frame form.
  data.long<-survSplit(new.data,cut=cutting_times,end="tstop",start="tstart",event="infect") ## Create an indicator infect representing whether an individual had an event at "cutting_times".
  ## The censoring model includes n.infect-.
  n.infect<-ave(data.long$infect,data.long$id,FUN=cumsum);indices<-determine_row_indices_at_which_identifier_variable_changes(data.long)
  data.long$n.infect<-unlist(lapply(1:n,function(sigma){
    if(length((indices[sigma]):(indices[sigma+1]))==2){ ## The individual under consideration gets immediately censored at the end of the first interval and therefore cannot have any COPD attack.
      return(0)
    }else{ ## The individual under consideration has at least two time intervals (not necessarily corresponding to intervals at which his/her n.infect value changes)
      return(c(0,n.infect[(indices[sigma]):(indices[sigma+1]-2)]))}}))
  data.long.cens<-survSplit(new.data,cut=cutting_times,end="tstop",start="tstart",event="cens") ## Create an indicator cens representing at which cutting time an individual was censored.
  data.long$cens<-data.long.cens$cens
  return(list("data.long"=data.long,"cutting_times"=cutting_times,"indices"=indices))}
##################################################################### DATA-GENERATING MECHANISM ###########################################################################################################################
## n (number of individuals), mu (mean of a normal distributed continuous baseline covariate), sigma (corresponding variance)
## coef_trt_out is the coefficient corresponding to the exposure A in the outcome model
## strength of censoring; coefficient which determines the censoring strength depending on the treatment assignment
## fup (follow up time of the study)
## shape (The shape parameter of a one-parameter gamma distribution which represents a frailty)
generate_data<-function(n,mu,sigma,max_fup,bas_haz_out,coef_trt_out,bas_haz_cens,censoring_strength_trt,censoring_strength_bas_cov,censoring_strength_time_varying_covariate,cens_interaction,shape,delta,sigmoid_type){
  ################################################################### DATA GENERATION PART ################################################################################################################################
  L1<-rnorm(n,mu,sigma) ## BASELINE COVARIATE (STRONGLY PREDICITIVE FOR CENSORING)
  A<-rbinom(n,size=1,prob=1/2) ## DETERMINE TREATMENT INDICATOR, WITH PROBABILITY 1/2 SINCE WE WORK IN A RANDOMISED STUDY.
  event_times<-create_event_times(n,A,L1,max_fup,bas_haz_out,coef_trt_out,shape) ## CREATE EVENT TIMES
  cens_times<-simulate_censoring_times(n,A,L1,max_fup,event_times,bas_haz_cens,censoring_strength_trt,censoring_strength_bas_cov,censoring_strength_time_varying_covariate,cens_interaction,delta,sigmoid_type) ## CREATE CENSORING TIMES
  list<-restrict_event_times_to_events_before_censoring_times(n,event_times,cens_times) ## RESTRICT EVENT TIMES TO THOSE BEFORE CENSORING TIMES
  data.event.time.covariates<-create_dataset_with_covariates_representing_at_which_time_point_the_q_th_event_took_place(n,list,max_fup) ## CREATE COVARIATES WHICH REPRESENT FOR THE J-TH EVENT ITS OCCURRENCE TIME.
  data.baseline<-data.frame("id"=1:n,"A"=A,"L1"=L1,"fup"=data.event.time.covariates[,ncol(data.event.time.covariates)]) ## CREATE DATASET SO THAT TMERGE CAN BE APPLIED.
  final_list<-turn_dataset_in_format_so_that_survsplit_can_be_used(n,data.baseline,data.event.time.covariates) ## TURN THE DATASET IN LONG FORMAT SO THAT A COX PROPORTIONAL HAZARDS MODEL CAN BE FITTED/ STORE ALSO THE CUTTING TIMES.
  cens.times<-unique(sort(list[[2]]));final_list[["cens.times"]]<-cens.times[!(cens.times %in% max_fup)] ## ADD THE CENSORING TIMES FOR EACH INDIVIDUAL TO REDUCE TIME COMPLEXITY LATER.
  return(final_list)}
##################################################################### RESTRICT DATA LONG TO CENSORING TIME POINTS #########################################################################################################
restrict_data_long_to_censoring_time_points<-function(n,data,indices,cens.times){
  ## DETERMINE THE DATA FRAME IN LONG FORMAT WHERE THE TIME WINDOW OF EACH INDIVIDUAL IS ONLY SPLITTED ACROSS CENSORING TIMES.
  data.long.cens.times<-do.call("rbind",lapply(1:n,function(sigma){
    current_data<-data[(indices[sigma]):(indices[sigma+1]-1),]
    current_data_at_cens_times<-current_data[which(current_data$tstop %in% cens.times),];number_of_rows<-length(current_data_at_cens_times$tstop)
    current_data_at_cens_times$tstart<-if(number_of_rows>1){c(0,current_data_at_cens_times$tstop[1:(length(current_data_at_cens_times$tstop)-1)])}else{0}
    return(current_data_at_cens_times)}))
  ## DETERMINE THE NUMBER OF TIMES A PARTICULAR WEIGHT WILL REMAIN UNCHANGED WHEN THE TIME WINDOW OF EACH INDIVIDUAL IS ALSO SPLITTED ACROSS TIME POINTS
  ## AT WHICH TIME-VARYING COVARIATES CHANGE. DO THIS FOR EACH POSSIBLE WEIGHT OF EACH INDIVIDUAL AND FOR EACH INDIVIDUAL
  number.of.repetitions<-list();cens.indices<-strtoi(rownames(subset(data.long.cens.times,id==1)))
  if(cens.indices[1]!=1){ ## THERE ARE TIME POINTS AT WHICH TIME-VARYING COVARIATES CHANGE BEFORE THE FIRST CENSORING TIME POINT OCCURS.
    number.of.repetitions<-lapply(1:n,function(sigma){
      current_data<-subset(data.long.cens.times,id==sigma);cens.indices<-strtoi(rownames(current_data))
      if(length(cens.indices)>1){ ## INDIVIDUAL SIGMA CROSSES AT LEAST ONE CENSORING TIME POINT BEFORE S/HE GETS HIM/HERSELF CENSORED
        if(data[(indices[sigma+1]-1),"cens"]!=1){ ## INDIVIDUAL SIGMA REACHES THE PRE-SPECIFIED FOLLOW-UP TIME.
          return(c(cens.indices[1]-indices[sigma],diff(cens.indices),indices[sigma+1]-1-cens.indices[length(cens.indices)]+1))
        }else{ ## INDIVIDUAL SIGMA DOES NOT REACH THE PRE-SPECIFIED FOLLOW-UP TIME
          return(c(cens.indices[1]-indices[sigma],diff(cens.indices),1))
        }
      }else{## INDIVIDUAL SIGMA GETS IMMEDIATELY CENSORED.
        ## BUT THERE IS AT LEAST ONE TIME POINT AT WHICH THE TIME VARYING-COVARIATE OF AN INDIVIDUAL CHANGES BEFORE S/HE GETS CENSORED.
        return(c(cens.indices[1]-indices[sigma],1))
      }
    })
  }else{ ## THE FIRST TIME POINT CORRESPONDS TO A POINT IN TIME AT WHICH CENSORING OCCURS.
    number.of.repetitions<-lapply(1:n,function(sigma){
      current_data<-subset(data.long.cens.times,id==sigma);cens.indices<-strtoi(rownames(current_data))
      if(length(cens.indices)>1){ ## INDIVIDUAL SIGMA CROSSES AT LEAST ONE CENSORING TIME POINT BEFORE S/HE GETS HIM/HERSELF CENSORED
        if(data[(indices[sigma+1]-1),"cens"]!=1){ ## INDIVIDUAL SIGMA REACHES THE PRE-SPECIFIED FOLLOW-UP TIME
          return(c(diff(cens.indices),indices[sigma+1]-1-cens.indices[length(cens.indices)]+1))
        }else{ ## INDIVIDUAL SIGMA DOES NOT REACH THE PRE-SPECIFIED FOLLOW-UP TIME
          return(c(diff(cens.indices),1))
        }
      }else{ ## INDIVIDUAL SIGMA GETS IMMEDIATELY CENSORED.
        return(1)
      }
    })
  }
  indices.data.long.cens<-determine_row_indices_at_which_identifier_variable_changes(data.long.cens.times)
  return(list("data.long.cens.times"=data.long.cens.times,"number.of.repetitions"=number.of.repetitions,"indices.data.long.cens"=indices.data.long.cens))}
##################################################################### DETERMINE CENSORING MODELS ##########################################################################################################################
determine_censoring_model_baseline_covariates<-function(data){
  return(coxph(Surv(tstart,tstop,cens)~A+L1,id=id,data=data,control=coxph.control(timefix=FALSE)))
} ## The censoring model is generated in function of A and L1.
determine_censoring_model_time_varying_covariates<-function(data,delta,sigmoid_type){
  return(coxph(Surv(tstart,tstop,cens)~A+L1+transformation(n.infect,delta,sigmoid_type)+A*transformation(n.infect,delta,sigmoid_type),id=id,data=data,control=coxph.control(timefix=FALSE)))
} ## The censoring model is generated in function of A,L1 and n.infect-! (n.infect is encoded as n.infect-)
##################################################################### DETERMINE IPCW ######################################################################################################################################
determine_IPCW<-function(n,data,cens_mod_time_varying_covariates,indices,data.long.cens.times,number.of.repetitions,indices.data.long.cens){
  if(data[1,"tstop"]!=data.long.cens.times[1,"tstop"]){ ## There is at least one point in time at which a time-varying covariate of an individual changes which precedes the first censoring time point.
    return(do.call("rbind",lapply(1:n,function(sigma){ ## The path of the time-varying covariate must be specified in the newdata argument of survfit to obtain correct survival probabilities!
      current_data<-data[(indices[sigma]):(indices[sigma+1]-1),] ## Restrict data.long to the long data frame for individual sigma.
      current_data_at_cens_times<-data.long.cens.times[(indices.data.long.cens[sigma]):(indices.data.long.cens[sigma+1]-1),] ## Restrict data.long.cens.times to the long data frame for individual sigma which only includes censoring time points.
      current_data$Unstab_IPCW<-rep(c(1,1/summary(survfit(cens_mod_time_varying_covariates,newdata=current_data_at_cens_times,id=id,timefix=FALSE))$surv),number.of.repetitions[[sigma]])
      return(current_data)})))
  }else{ ## The first point in time corresponds to a point in time at which censoring occurs.
    return(do.call("rbind",lapply(1:n,function(sigma){ ## The path of the time-varying covariate must be specified in the newdata argument of survfit to obtain correct survival probabilities!
      current_data<-data[(indices[sigma]):(indices[sigma+1]-1),]
      current_data_at_cens_times<-data.long.cens.times[(indices.data.long.cens[sigma]):(indices.data.long.cens[sigma+1]-1),] ## Restrict data.long.cens.times to the long data frame for individual sigma which only includes censoring time points.
      current_data$Unstab_IPCW<-rep(1/summary(survfit(cens_mod_time_varying_covariates,newdata=current_data_at_cens_times,id=id,timefix=FALSE))$surv,number.of.repetitions[[sigma]])
      return(current_data)})))}}
##################################################################### DETERMINE IPCW STABILIZED ###########################################################################################################################
determine_IPCW_stabilized<-function(n,data,cens_mod_baseline_covariates,indices,data.long.cens.times,number.of.repetitions,indices.data.long.cens){
  if(data[1,"tstop"]!=data.long.cens.times[1,"tstop"]){ ## There is at least one point in time at which a time-varying covariate of an individual changes which precedes the first censoring time point.
    return(do.call("rbind",lapply(1:n,function(sigma){ ## The path of the time-varying covariate must be specified in the newdata argument of survfit to obtain correct survival probabilities!
      current_data<-data[(indices[sigma]):(indices[sigma+1]-1),] ## Restrict data.long to the long data frame for individual sigma.
      current_data_at_cens_times<-data.long.cens.times[(indices.data.long.cens[sigma]):(indices.data.long.cens[sigma+1]-1),] ## Restrict data.long.cens.times to the long data frame for individual sigma which only includes censoring time points.
      current_data$Stab_IPCW<-current_data$Unstab_IPCW*rep(c(1,summary(survfit(cens_mod_baseline_covariates,newdata=current_data_at_cens_times,id=id,timefix=FALSE))$surv),number.of.repetitions[[sigma]])
      return(current_data)})))
  }else{ ## The first point in time corresponds to a point in time at which censoring occurs.
    return(do.call("rbind",lapply(1:n,function(sigma){ ## The path of the time-varying covariate must be specified in the newdata argument of survfit to obtain correct survival probabilities!
      current_data<-data[(indices[sigma]):(indices[sigma+1]-1),] ## Restrict data.long to the long data frame for individual sigma.
      current_data_at_cens_times<-data.long.cens.times[(indices.data.long.cens[sigma]):(indices.data.long.cens[sigma+1]-1),] ## Restrict data.long.cens.times to the long data frame for individual sigma which only includes censoring time points.
      current_data$Stab_IPCW<-current_data$Unstab_IPCW*rep(summary(survfit(cens_mod_baseline_covariates,newdata=current_data_at_cens_times,id=id,timefix=FALSE))$surv,number.of.repetitions[[sigma]])
      return(current_data)})))}}
##################################################################### AUXILIARY FUNCTION FOR FIXED TIME POINT #############################################################################################################
create_uncensored_wide_data_frame_at_a_given_time_point<-function(data,cutting_times,tau,max_fup,column_name_weights){
  row_index_smallest_time_point_bigger_than_tau<-findInterval(tau,cutting_times)+if(cutting_times[findInterval(tau,cutting_times)]==max_fup){0}else{1} ## Row index at which we find the smallest cutting time (in the column tstop) bigger than tau.
  if(cutting_times[row_index_smallest_time_point_bigger_than_tau]!=max_fup){ ## The fixed time point tau is located between two cutting times, where the right endpoint does not correspond to the pre-specified follow-up time.
    current_data<-subset(data,id %in% unique((data[which(data$fup>tau),])$id)) ## Restrict data to subjects whose censoring time lies beyond time point tau.
    indices_current_data<-determine_row_indices_at_which_identifier_variable_changes(current_data) ## Determine the row indices at which the identifier variable changes for current data.
    return(do.call("rbind",lapply(1:(length(indices_current_data)-1),function(sigma){
      aux_data<-current_data[(indices_current_data[sigma]):(indices_current_data[sigma+1]-1),]
      ## We need the weight for individual sigma at the cutting time which preceeds time point tau as stopping time,
      ## but the records for that individual at the cutting time which exceeds time point tau as stopping time.
      final_data<-aux_data[row_index_smallest_time_point_bigger_than_tau,]
      final_data[,column_name_weights]<-if(row_index_smallest_time_point_bigger_than_tau>1){aux_data[row_index_smallest_time_point_bigger_than_tau-1,column_name_weights]}else{1}
      return(final_data)
    })))
  }else{ ## The fixed time point tau equals EITHER the maximum follow-up time OR is located between the penultimate and ultimate cutting time. Note that the ultimate cutting time corresponds to the pre-specified follow-up time.
    ## For the latter case, individuals who are uncensored at time point tau, are also uncensored at the pre-specified follow-up time,
    ## as no censoring can occur at this pre-specified follow-up time. We hence only include individuals who reach the pre-specified follow-up time. 
    return(subset(data,data$tstop==max_fup))}}
##################################################################### AUXILIARY FUNCTION FOR LAST TIME POINT ##############################################################################################################
create_uncensored_wide_data_frame_at_the_last_time_point<-function(n,data,column_name_weights){
  return(do.call("rbind",lapply(1:n,function(rho){
    current_dataset<-subset(data,id==rho);n_of_rows<-nrow(current_dataset)
    if((current_dataset[n_of_rows,"cens"])==1){ ## Individual rho gets censored at last time point.
      if(n_of_rows>1){ ## Individual rho has previous uncensored time points. Which data has to be returned follows the same reasoning as outlined in the previous auxiliary function.
        data.wide<-current_dataset[n_of_rows,];data.wide[,column_name_weights]<-current_dataset[n_of_rows-1,column_name_weights]
        return(data.wide)
      }else{ ## Individual rho is censored at the first possible time point, meaning that his/her weight should be adapted to 1.
        current_dataset[,column_name_weights]<-1; return(current_dataset)}
    }else{ ## Individual rho is uncensored at the last point in time (i.e. s/he reaches the pre-specified follow-up time.)
      return(current_dataset[n_of_rows,])}})))}
##################################################################### AUXILIARY FUNCTION FOR LAST TIME POINT WITH ADDITIONAL HAZARD #######################################################################################
adapt_weights_with_additional_hazard<-function(n,data,data_last_time_point,cens_mod_baseline_covariates,cens_mod_time_varying_covariates,indices,max_fup,delta,sigmoid_type){
  aux_list_time_varying_covariates<-cens_mod_time_varying_covariates[["coefficients"]]
  aux_list_fixed_covariates<-cens_mod_baseline_covariates[["coefficients"]]
  for(rho in 1:n){
    if(data_last_time_point[rho,"fup"]!=max_fup){ ## Individual rho did not reach the pre-specified follow-up time.
      current_data<-data[indices[rho+1]-1,]
      data_last_time_point[rho,"Stab_IPCW"]<-data_last_time_point[rho,"Stab_IPCW"]*(exp(aux_list_fixed_covariates[["A"]]*current_data$A+aux_list_fixed_covariates[["L1"]]*current_data$L1))/(exp(aux_list_time_varying_covariates[["A"]]*current_data$A+aux_list_time_varying_covariates[["L1"]]*current_data$L1+aux_list_time_varying_covariates[[3]]*transformation(current_data$n.infect,delta,sigmoid_type)+aux_list_time_varying_covariates[[4]]*(current_data$A)*transformation(current_data$n.infect,delta,sigmoid_type)))
    }## There is no else branche, as the weight for an individual who reaches the pre-specified follow-up time should not be adapted.
  }
  return(data_last_time_point)}
##################################################################### AUXILIARY FUNCTION FOR MULTIPLE TIME POINTS USING GEE ROUTINES ######################################################################################
determine_long_data_frame_in_balanced_form_for_the_vector_of_outcomes<-function(n,data.long,set_of_partition_points,indices,cutting_times,max_fup,increment){ ## Censoring time points excluded.
  data.long.balanced<-lapply(1:n,function(sigma){
    if(data.long[indices[sigma],"fup"]!=max_fup){ ## Individual sigma gets censored before the pre-specified follow-up time.
      aux_var<-findInterval(data.long[indices[sigma],"fup"],set_of_partition_points) ## Determine for individual sigma the largest visit time at which one is still uncensored.
      if(aux_var>0){ ## Individual sigma does not get censored before the first visit time. (The first visit time is assumed to be strictly larger than 0).
        set_of_partition_points_for_ind_sigma<-set_of_partition_points[1:aux_var] ## Partition points at which records of individual sigma are observed.
        current_data<-data.long[(indices[sigma]):(indices[sigma+1]-1),c("id","A","L1","n.infect","cens")] ## The relevant data in long format for individual sigma.
        list_of_data_at_visit_times<-lapply(1:(length(set_of_partition_points_for_ind_sigma)),function(rho){ ## Go through all visit times for individual sigma.
          data_at_visit_times<-current_data[findInterval(set_of_partition_points_for_ind_sigma[rho],cutting_times)+1,] ## Determine the largest cutting time smaller than the rho-th visit time; go one line deeper to get the number of events one experienced by the rho-th visit time.
          return(data_at_visit_times)
        })
        balanced_data_for_ind_sigma<-do.call("rbind",list_of_data_at_visit_times) ## Convert list_of_data_at_visit_times to a data frame via row binding.
        balanced_data_for_ind_sigma$tstart<-if(length(set_of_partition_points_for_ind_sigma)>1){c(0,set_of_partition_points_for_ind_sigma[1:(length(set_of_partition_points_for_ind_sigma)-1)])}else{0} ## Add a column tstart representing the starting time points for each interval.
        balanced_data_for_ind_sigma$tstop<-set_of_partition_points_for_ind_sigma ## Add a column tstop representing the stopping time points for each interval.
        balanced_data_for_ind_sigma$cens<-0 ## Individual sigma is uncensored at each of these visit times, but could have a cens value equal to 1 at the last visit time. (The values of the last visit time are determined by the cutting time following this last visit time, at which the individual might get censored.)
        return(balanced_data_for_ind_sigma)
      } ## The individual gets censored before the first visit time, and should therefore not be included in the analysis.
    }else{ ## Individual sigma reaches the pre-specified follow-up time.
      current_data<-data.long[(indices[sigma]):(indices[sigma+1]-1),c("id","A","L1","n.infect","cens")] ## The relevant data in long format for individual sigma.
      list_of_data_at_visit_times<-lapply(1:(length(set_of_partition_points)),function(rho){ ## Go through all visit times for individual sigma.
        aux_var<-findInterval(set_of_partition_points[rho],cutting_times) ## Determine largest cutting time smaller than the rho-th visit time.
        data_at_visit_times<-current_data[if(aux_var>0){aux_var+if(cutting_times[aux_var]!=max_fup){1}else{0}}else{1},] ## Records of individual sigma corresponding to visit time rho.
        return(data_at_visit_times)
      })
      balanced_data_for_ind_sigma<-do.call("rbind",list_of_data_at_visit_times) ## Convert list_of_data_at_visit_times to a data frame via row binding.
      balanced_data_for_ind_sigma$tstart<-c(0,set_of_partition_points[1:(length(set_of_partition_points)-1)]) ## Add a column tstart representing the starting time points for each interval.
      balanced_data_for_ind_sigma$tstop<-set_of_partition_points ## Add a column tstop representing the stopping time points for each interval.
      return(balanced_data_for_ind_sigma)
    }
  })
  return(do.call("rbind",data.long.balanced))}
##################################################################### AUXILIARY FUNCTION TO RESTRICT A LONG DATA FRAME TO UNCENSORED OBSERVATIONS #########################################################################
restrict_data_to_uncensored_observations<-function(data.long){return(subset(data.long,cens==0))}
##################################################################### AUXILIARY FUNCTION WHICH DETERMINES WHICH PERCENTAGE OF INDIVIDUALS WAS CENSORED ####################################################################
determine_percentage_censored<-function(n,max_fup,data.long){return(length(unique(data.long[which(data.long[,"fup"]<max_fup),]$id))/n)}
##################################################################### AUXILIARY FUNCTION WHICH DETERMINES THE MAXIMUM NUMBER OF INFECTIONS OVER [0,max_fup] ###############################################################
determine_max_number_of_infections<-function(data.long){return(max(data.long[,"n.infect"]))}
##################################################################### AUXILIARY FUNCTION WHICH WILL BE USED TO COMPUTE THE RUNTIME FOR BALANCED DATA MORE ACCURATELY ######################################################
###########################################################################################################################################################################################################################
##################################################################### ASSESSING THE TREATMENT EFFECT ######################################################################################################################
###########################################################################################################################################################################################################################
## Actually, we need to fit a negative binomial regression (i.e. we should use glm.nb() in each method below.)
## But, when the shape parameter becomes extremely large, glm.nb() returns an error (caused by the fact that the variance is not significantly larger than the mean).
## Indeed, when the shape parameter \theta becomes large, the variance (mu+mu^2/\theta) tends to \mu.
## To overcome this problem, one can reason that the estimating equation corresponding to poisson regression will yield consistent estimates of the regression coefficients too.
## This therefore justifies the use of glm(formula,family="poisson",...) instead of glm.nb(formula,...) to avoid errors due to extremely large shape parameters.
##################################################################### APPROACH I, III, V: ASSESS TREATMENT EFFECT AT A FIXED THEORETICAL TIME POINT, NO WEIGHTS ###########################################################
Method1<-function(data,cutting_times,tau,column_name_weights,max_fup){
  data.wide<-create_uncensored_wide_data_frame_at_a_given_time_point(data,cutting_times,tau,max_fup,column_name_weights)
  pr<-glm(n.infect~A+L1+offset(log(fup)),family="poisson",data.wide,control=glm.control(maxit=10000))
  return(list("pr"=pr[["coefficients"]][["A"]],"min.weight"=1,"max.weight"=1,"Sandwich.SE"=coeftest(pr,vcov=sandwich)[2,2]))}
##################################################################### APPROACH II, IV, VI: ASSESS THE TREATMENT EFFECT AT A FIXED THEORETICAL TIME POINT, STABILIZED WEIGHTS ##############################################
Method2<-function(data,cutting_times,tau,column_name_weights,max_fup){
  data.wide<-create_uncensored_wide_data_frame_at_a_given_time_point(data,cutting_times,tau,max_fup,column_name_weights)
  pr.IPCW<-glm(n.infect~A+L1+offset(log(fup)),family="poisson",data.wide,weights=data.wide[,column_name_weights],control=glm.control(maxit=10000))
  min.weight<-min(data.wide[,column_name_weights]);max.weight<-max(data.wide[,column_name_weights]) ## Compute minimal and maximal weights
  return(list("pr.IPCW"=pr.IPCW[["coefficients"]][["A"]],"min.weight"=min.weight,"max.weight"=max.weight,"Sandwich.SE"=coeftest(pr.IPCW,vcov=sandwich)[2,2]))}
##################################################################### APPROACH VII, VIII, IX, X: ASSESS THE TREATMENT EFFECT BY SPLITTING THE TIME INTERVAL IN A FIXED NUMBER OF TIME WINDOWS, STABILIZED WEIGHTS #########
Method3<-function(data.long.weights.balanced,column_name_weights,corstr){
  pr.IPCW.fup.splitted<-geem(n.infect~A+L1+offset(log(tstop)),family="poisson",data=data.long.weights.balanced,weights=data.long.weights.balanced[,column_name_weights],id=id,corstr=corstr)
  pr.fup.splitted<-geem(n.infect~A+L1+offset(log(tstop)),family="poisson",data=data.long.weights.balanced,id=id,corstr=corstr)
  return(list("pr.fup.splitted"=summary(pr.fup.splitted)[[1]][2],"Sandwich.SE"=summary(pr.fup.splitted)[[3]][2],
              "pr.IPCW.fup.splitted"=summary(pr.IPCW.fup.splitted)[[1]][2],"Sandwich.SE"=summary(pr.IPCW.fup.splitted)[[3]][2]))}
##################################################################### AUXILIARY FUNCTION TO DETERMINE HOW INFLUENTIAL INDIVIDUALS ARE FOR THE LONG DATA FRAME FORMAT ######################################################
determine_ratio_of_weights_for_long_format<-function(data.long.balanced.Stab,indices.data.long.balanced,increment){
  vector_of_averaged_weights<-sapply(1:(length(indices.data.long.balanced)-1),function(rho){
    current.data<-data.long.balanced.Stab[(indices.data.long.balanced[rho]):(indices.data.long.balanced[rho+1]-1),]
    return((sum(current.data[,"Stab_IPCW"])*increment)/(current.data[nrow(current.data),"tstop"]))})
  return(list("min.weight"=min(vector_of_averaged_weights),"max.weight"=max(vector_of_averaged_weights)))}
##################################################################### AUXILIARY FUNCTION TO DETERMINE THE ESTIMATED TREATMENT EFFECT WHEN SPLITTING THE FOLLOW-UP TIME FOR A GIVEN INCREMENT ##############################
determine_treatment_coefficient_by_splitting_the_follow_up_time<-function(n,data.long,set_of_partition_points,indices,cutting_times,max_fup,increment,data.long.Stab){
  ## DETERMINE THE LONG DATA FRAME FORMAT WHICH ONE SHOULD PASS TO ANY GEE ROUTINE
  data.long.balanced<-determine_long_data_frame_in_balanced_form_for_the_vector_of_outcomes(n,data.long,set_of_partition_points,indices,cutting_times,max_fup,increment)
  indices.data.long.balanced<-determine_row_indices_at_which_identifier_variable_changes(data.long.balanced)
  ## DETERMINE THE STABILISED WEIGHTS FOR THE GEE PROCEDURE
  Stab.IPCW<-unlist(sapply(1:(length(indices.data.long.balanced)-1),function(sigma){
    current.data.balanced<-data.long.balanced[(indices.data.long.balanced[sigma]):(indices.data.long.balanced[sigma+1]-1),]
    sigma_bis<-unique(current.data.balanced[,"id"]) ## ID variable of individual sigma in data.long.Stab
    current.data.long.Stab<-data.long.Stab[(indices[sigma_bis]):(indices[sigma_bis+1]-1),]
    return(sapply(1:(nrow(current.data.balanced)),function(rho){
      aux_var<-findInterval(current.data.balanced[rho,"tstop"],current.data.long.Stab[,"tstop"])
      if(aux_var>0){return(current.data.long.Stab[aux_var,"Stab_IPCW"])}else{return(1)}}))}))
  ## CREATE THE DATA FRAME OF INTEREST WITH STABILISED WEIGHTS
  data.long.balanced.Stab<-data.long.balanced;data.long.balanced.Stab$Stab_IPCW<-Stab.IPCW;list_of_min_and_max_weights<-determine_ratio_of_weights_for_long_format(data.long.balanced.Stab,indices.data.long.balanced,increment)
  res.object.balanced.stabilized<-Method3(data.long.balanced.Stab,"Stab_IPCW","independence"); res.object.balanced.stabilized.ar1<-Method3(data.long.balanced.Stab,"Stab_IPCW","ar1")
  ## DETERMINE THE TREATMENT EFFECT BY MAKING USE OF THE OUTCOME VALUES AT SEVERAL PRE-SPECIFIED TIME POINTS, STABILIZED WEIGHTS
  return(list("Naive.trt.coef.est"=res.object.balanced.stabilized[[1]],"min.weight.naive"=1,"max.weight.naive"=1,"Sandwich.SE.Naive"=res.object.balanced.stabilized[[2]],
              "Stab.trt.coef.est"=res.object.balanced.stabilized[[3]],"min.weight.stab"=list_of_min_and_max_weights[[1]],"max.weight.stab"=list_of_min_and_max_weights[[2]],"Sandwich.SE.Stabilized"=res.object.balanced.stabilized[[4]],
              "Naive.trt.coef.est.ar1"=res.object.balanced.stabilized.ar1[[1]],"min.weight.naive.ar1"=1,"max.weight.naive.ar1"=1,"Sandwich.SE.Naive.ar1"=res.object.balanced.stabilized.ar1[[2]],
              "Stab.trt.coef.est.ar1"=res.object.balanced.stabilized.ar1[[3]],"min.weight.stab.ar1"=list_of_min_and_max_weights[[1]],"max.weight.stab.ar1"=list_of_min_and_max_weights[[2]],"Sandwich.SE.Stabilized.ar1"=res.object.balanced.stabilized.ar1[[4]]))}
##################################################################### APPROACH XI: ASSESS THE TREATMENT EFFECT BY MAKING USE OF THE LAST UNCENSORED TIME POINT, NO WEIGHTS ################################################
Method4<-function(data.wide.weights){
  pr.last.time.point<-glm(n.infect~A+L1+offset(log(fup)),family="poisson",data.wide.weights,control=glm.control(maxit=10000))
  return(list("pr.last.time.point"=pr.last.time.point[["coefficients"]][["A"]],"min.weight"=1,"max.weight"=1,"Sandwich.SE"=coeftest(pr.last.time.point,vcov=sandwich)[2,2]))}
##################################################################### APPROACH XII: ASSESS THE TREATMENT EFFECT BY MAKING USE OF THE LAST UNCENSORED TIME POINT, STABILIZED WEIGHTS, NO ADDITIONAL HAZARD #################
Method5<-function(data.wide.weights,column_name_weights){
  pr.IPCW.last.time.point<-glm(n.infect~A+L1+offset(log(fup)),family="poisson",data.wide.weights,weights=data.wide.weights[,column_name_weights],control=glm.control(maxit=10000))
  min.weight<-min(data.wide.weights[,column_name_weights]);max.weight<-max(data.wide.weights[,column_name_weights]) ## Compute minimal and maximal weights
  return(list("pr.IPCW.last.time.point"=pr.IPCW.last.time.point[["coefficients"]][["A"]],"min.weight"=min.weight,"max.weight"=max.weight,"Sandwich.SE"=coeftest(pr.IPCW.last.time.point,vcov=sandwich)[2,2]))}
##################################################################### APPROACH XIII: ASSESS THE TREATMENT EFFECT BY MAKING USE OF THE LAST UNCENSORED TIME POINT, STABILIZED WEIGHTS, ADDITIONAL HAZARD ####################
Method6<-function(n,data.long.weights,data_last_time_point,cens_mod_baseline_covariates,cens_mod_time_varying_covariates,indices,column_name_weights,max_fup,truncated,delta,sigmoid_type){
  data.wide<-adapt_weights_with_additional_hazard(n,data.long.weights,data_last_time_point,cens_mod_baseline_covariates,cens_mod_time_varying_covariates,indices,max_fup,delta,sigmoid_type)
  if(truncated){
    lower<-quantile(data.wide[,column_name_weights],0.025)[[1]]; upper<-quantile(data.wide[,column_name_weights],0.975)[[1]]
    data.wide[,column_name_weights]<-sapply(1:n,function(sigma){if(data.wide[sigma,column_name_weights]<lower){return(lower)}else if(data.wide[sigma,column_name_weights]>upper){return(upper)}else{return(data.wide[sigma,column_name_weights])}})
  }
  pr.IPCW.last.time.point.additional.hazard<-glm(n.infect~A+L1+offset(log(fup)),family="poisson",data.wide,weights=data.wide[,column_name_weights],control=glm.control(maxit=10000))
  min.weight<-min(data.wide[,column_name_weights]);max.weight<-max(data.wide[,column_name_weights]) ## Compute minimal and maximal weights
  return(list("pr.IPCW.last.time.point.additional.hazard"=pr.IPCW.last.time.point.additional.hazard[["coefficients"]][["A"]],"min.weight"=min.weight,"max.weight"=max.weight,"Sandwich.SE"=coeftest(pr.IPCW.last.time.point.additional.hazard,vcov=sandwich)[2,2]))}
##################################################################### APPROACH XIV: NAIVE ANDERSEN-GILL APPROACH ##########################################################################################################
Method7<-function(data){
  ## Fit the Andersen-Gill model.
  AG<-coxph(Surv(tstart,tstop,infect)~A+L1+cluster(id),data=data,timefix=FALSE)
  ## Return the point estimate corresponding to treatment; a conservative standard error for the coefficient corresponding to treatment ###################################################################################
  return(list("coef_AG"=AG[["coefficients"]][["A"]],"min.weight"=1,"max.weight"=1,"Sandwich.SE"=summary(AG)[["coefficients"]][1,4]))}
##################################################################### APPROACH XV: ANDERSEN-GILL APPROACH WITH IPCW WEIGHTS ###############################################################################################
Method8<-function(data,indices){
  ## Fit the Andersen-Gill model with weights the stabilised IPCW weights.
  data<-do.call("rbind",lapply(1:(length(indices)-1),function(sigma){
    current_data<-data[(indices[sigma]):(indices[sigma+1]-1),]
    current_data[,"Stab_IPCW"]<-if(length(current_data[,"Stab_IPCW"])>1){c(1,current_data[,"Stab_IPCW"][1:(length(current_data[,"Stab_IPCW"])-1)])}else{1}
    return(current_data)
  }))
  AG_IPCW<-coxph(Surv(tstart,tstop,infect)~A+L1+cluster(id),data=data,weights=data[,"Stab_IPCW"],timefix=FALSE)
  ## Determine the average weight for each individual
  vector_of_averaged_weights<-sapply(1:(length(indices)-1),function(rho){
    ## The degree to which a weight counts in computing the average weight of an individual depends on how long each weight is observed.
    current.data<-data[(indices[rho]):(indices[rho+1]-1),]
    increment<-c(current.data[1,"tstop"],diff(current.data[,"tstop"]))
    return((sum(current.data[,"Stab_IPCW"]*increment))/(current.data[nrow(current.data),"tstop"]))})
  ## Determine the minimum and the maximum average weight.
  min_weight<-min(vector_of_averaged_weights); max_weight<-max(vector_of_averaged_weights)
  ## Return the point estimate corresponding to treatment; the minimum and maximum average weight; a conservative standard error for the coefficient corresponding to treatment #######################################
  return(list("coef_AG_IPCW"=AG_IPCW[["coefficients"]][["A"]],"min.weight"=min_weight,"max.weight"=max_weight,"Sandwich.SE"=summary(AG_IPCW)[["coefficients"]][1,4]))}
###########################################################################################################################################################################################################################
##################################################################### SIMULATION PART #####################################################################################################################################
MonteCarlo_sim<-function(n,max_fup,bas_haz_out,coef_trt_out,bas_haz_cens,censoring_strength_trt,censoring_strength_bas_cov,censoring_strength_time_varying_covariate,cens_interaction,shape,delta,sigmoid_type){
  return(tryCatch(expr={ ## TRY THE EXPRESSION PART. IF SUCCESFULL, THE LAST LINE OF THE TRY PART WILL BE RETURNED
    ################################################################# PART I ##############################################################################################################################################
    #######################################################################################################################################################################################################################
    ## GENERATE THE DATA AND STORE THE LONG DATA FRAME FORMAT, THE CUTTING TIMES, THE INDICES AND THE CENSORING TIMES IN THE CORRESPONDING OBJECTS.
    object<-generate_data(n,mu,sigma,max_fup,bas_haz_out,coef_trt_out,bas_haz_cens,censoring_strength_trt,censoring_strength_bas_cov,censoring_strength_time_varying_covariate,cens_interaction,shape,delta,sigmoid_type); data.long<-object[[1]];cutting_times<-object[[2]];indices<-object[[3]];cens.times<-object[[4]]
    ## AUXILIARY OBJECT (IN VIEW OF REDUCING TIME COMPLEXITY)
    aux_object_part_1<-restrict_data_long_to_censoring_time_points(n,data.long,indices,cens.times);data.long.cens.times<-aux_object_part_1[[1]];number.of.repetitions<-aux_object_part_1[[2]];indices.data.long.cens<-aux_object_part_1[[3]]
    ## CENSORING MODEL NECESSARY TO COMPUTE THE STABILIZED WEIGHTS
    cens_mod_baseline_covariates<-determine_censoring_model_baseline_covariates(data.long.cens.times)
    ## CENSORING MODEL NECESSARY TO COMPUTE THE (UN)STABILIZED WEIGHTS.
    cens_mod_time_varying_covariates<-determine_censoring_model_time_varying_covariates(data.long.cens.times,delta,sigmoid_type)
    ## DETERMINE THE UNSTABILIZED WEIGHTS FOR THE LONG DATA FRAME FORMAT.
    data.long.Unstab<-determine_IPCW(n,data.long,cens_mod_time_varying_covariates,indices,data.long.cens.times,number.of.repetitions,indices.data.long.cens)
    ## DETERMINE THE STABILIZED WEIGHTS FOR THE LONG DATA FRAME FORMAT.
    data.long.Stab<-determine_IPCW_stabilized(n,data.long.Unstab,cens_mod_baseline_covariates,indices,data.long.cens.times,number.of.repetitions,indices.data.long.cens)
    ## DETERMINE THE TREATMENT EFFECT AT HALF OFTHE PRE-SPECIFIED FOLLOW-UP TIME, UNSTABILIZED WEIGHTS.
    aux_Method_1<-Method1(data.long.Unstab,cutting_times,1/2*max_fup,"Unstab_IPCW",max_fup)
    res1<-aux_Method_1[[1]];min.weight.Appr.1<-aux_Method_1[[2]];max.weight.Appr.1<-aux_Method_1[[3]];Sandwich.SE.A.1<-aux_Method_1[[4]]
    ## DETERMINE THE TREATMENT EFFECT AT HALF OF THE PRE-SPECIFIED FOLLOW-UP TIME, STABILIZED WEIGHTS.
    aux_Method_2<-Method2(data.long.Stab,cutting_times,1/2*max_fup,"Stab_IPCW",max_fup)
    res2<-aux_Method_2[[1]];min.weight.Appr.2<-aux_Method_2[[2]];max.weight.Appr.2<-aux_Method_2[[3]];Sandwich.SE.A.2<-aux_Method_2[[4]]
    ## DETERMINE THE TREATMENT EFFECT AT THREE QUARTERS OF THE PRE-SPECIFIED FOLLOW-UP TIME, UNSTABILIZED WEIGHTS.
    aux_Method_3<-Method1(data.long.Unstab,cutting_times,3/4*max_fup,"Unstab_IPCW",max_fup)
    res3<-aux_Method_3[[1]];min.weight.Appr.3<-aux_Method_3[[2]];max.weight.Appr.3<-aux_Method_3[[3]];Sandwich.SE.A.3<-aux_Method_3[[4]]
    ## DETERMINE THE TREATMENT EFFECT AT THREE QUARTERS OF THE PRE-SPECIFIED FOLLOW-UP TIME, STABILIZED WEIGHTS.
    aux_Method_4<-Method2(data.long.Stab,cutting_times,3/4*max_fup,"Stab_IPCW",max_fup)
    res4<-aux_Method_4[[1]];min.weight.Appr.4<-aux_Method_4[[2]];max.weight.Appr.4<-aux_Method_4[[3]];Sandwich.SE.A.4<-aux_Method_4[[4]]
    ## DETERMINE THE TREATMENT EFFECT AT THE PRE-SPECIFIED FOLLOW-UP TIME, UNSTABILIZED WEIGHTS.
    aux_Method_5<-Method1(data.long.Unstab,cutting_times,max_fup,"Unstab_IPCW",max_fup)
    res5<-aux_Method_5[[1]];min.weight.Appr.5<-aux_Method_5[[2]];max.weight.Appr.5<-aux_Method_5[[3]];Sandwich.SE.A.5<-aux_Method_5[[4]]
    ## DETERMINE THE TREATMENT EFFECT AT THE PRE-SPECIFIED FOLLOW-UP TIME, STABILIZED WEIGHTS.
    aux_Method_6<-Method2(data.long.Stab,cutting_times,max_fup,"Stab_IPCW",max_fup)
    res6<-aux_Method_6[[1]];min.weight.Appr.6<-aux_Method_6[[2]];max.weight.Appr.6<-aux_Method_6[[3]];Sandwich.SE.A.6<-aux_Method_6[[4]]
    ################################################################# PART II #############################################################################################################################################
    #######################################################################################################################################################################################################################
    ## DETERMINE THE TREATMENT EFFECT BY MAKING USE OF THE OUTCOME VALUES AT SEVERAL PRE-SPECIFIED TIME POINTS (SPLIT THE FOLLOW-UP TIME OVER MULTIPLE TIME WINDOWS FOR EACH INDIVIDUAL)
    increment_grid<-c(50,10) ## SET OF PARTITION POINTS WHICH REPRESENT THE VISIT TIMES
    res_for_each_increment<-lapply(1:(length(increment_grid)),function(psi){
      set_of_partition_points<-seq(increment_grid[psi],max_fup,by=increment_grid[psi])
      return(determine_treatment_coefficient_by_splitting_the_follow_up_time(n,data.long,set_of_partition_points,indices,cutting_times,max_fup,increment_grid[psi],data.long.Stab))})
    ################################################################# PART III ############################################################################################################################################
    #######################################################################################################################################################################################################################
    ## DETERMINE THE TREATMENT EFFECT BY MAKING USE OF THE RECORDS AT THE LAST UNCENSORED TIME POINT FOR EACH INDIVIDUAL NO WEIGHTS
    data.wide.Stab<-create_uncensored_wide_data_frame_at_the_last_time_point(n,data.long.Stab,"Stab_IPCW"); aux_object<-Method4(data.wide.Stab)
    res11<-aux_object[[1]];min.weight.Appr.11<-aux_object[[2]];max.weight.Appr.11<-aux_object[[3]];Sandwich.SE.A.11<-aux_object[[4]]
    ## DETERMINE THE TREATMENT EFFECT BY MAKING USE OF THE RECORDS AT THE LAST UNCENSORED TIME POINT FOR EACH INDIVIDUAL, STABILIZED WEIGHTS (APPROACH ANDREA)
    res.object.stabilized<-Method5(data.wide.Stab,"Stab_IPCW")
    res12<-res.object.stabilized[[1]];min.weight.Appr.12<-res.object.stabilized[[2]];max.weight.Appr.12<-res.object.stabilized[[3]];Sandwich.SE.A.12<-res.object.stabilized[[4]]
    ## DETERMINE THE TREATMENT EFFECT BY MAKING USE OF THE RECORDS AT THE LAST UNCENSORED TIME POINT FOR EACH INDIVIDUAL, STABILIZED WEIGHTS, STABILIZED ADDITIONAL HAZARD.
    aux_object_stabilized<-Method6(n,data.long.Stab,data.wide.Stab,cens_mod_baseline_covariates,cens_mod_time_varying_covariates,indices,"Stab_IPCW",max_fup,FALSE,delta,sigmoid_type)
    res13<-aux_object_stabilized[[1]];min.weight.Appr.13<-aux_object_stabilized[[2]];max.weight.Appr.13<-aux_object_stabilized[[3]];Sandwich.SE.A.13<-aux_object_stabilized[[4]]
    ## DETERMINE THE TREATMENT EFFECT BY MAKING USE OF THE RECORDS AT THE LAST UNCENSORED TIME POINT FOR EACH INDIVIDUAL, STABILIZED WEIGHTS, STABILIZED ADDITIONAL HAZARD, TRUNCATED WEIGHTS
    aux_object_stabilized_truncated<-Method6(n,data.long.Stab,data.wide.Stab,cens_mod_baseline_covariates,cens_mod_time_varying_covariates,indices,"Stab_IPCW",max_fup,TRUE,delta,sigmoid_type)
    res14<-aux_object_stabilized_truncated[[1]];min.weight.Appr.14<-aux_object_stabilized_truncated[[2]];max.weight.Appr.14<-aux_object_stabilized_truncated[[3]];Sandwich.SE.A.14<-aux_object_stabilized_truncated[[4]]
    ################################################################# PART IV #############################################################################################################################################
    #######################################################################################################################################################################################################################
    ## Determine the treatment effect via (i) the Andersen-Gill model and (ii) the Andersen-Gill model with IPCW weights
    aux_AG<-Method7(data.long.Stab); aux_AG_IPCW<-Method8(data.long.Stab,indices)
    ################################################################# PART IV #############################################################################################################################################
    #######################################################################################################################################################################################################################
    ## WHAT THE FUNCTION WILL RETURN
    list("A.1"=res1,"A.2"=res2,"A.3"=res3,"A.4"=res4,"A.5"=res5,"A.6"=res6,"A.7"=res_for_each_increment[[1]][[1]],"A.8"=res_for_each_increment[[1]][[5]],"A.9"=res_for_each_increment[[2]][[1]],"A.10"=res_for_each_increment[[2]][[5]],"A.11"=res_for_each_increment[[1]][[9]],"A.12"=res_for_each_increment[[1]][[13]],"A.13"=res_for_each_increment[[2]][[9]],"A.14"=res_for_each_increment[[2]][[13]],"A.15"=res11,"A.16"=res12,"A.17"=res13,"A.18"=res14,
         "min.weight.Appr.1"=min.weight.Appr.1,"max.weight.Appr.1"=max.weight.Appr.1,"min.weight.Appr.2"=min.weight.Appr.2,"max.weight.Appr.2"=max.weight.Appr.2,"min.weight.Appr.3"=min.weight.Appr.3,"max.weight.Appr.3"=max.weight.Appr.3,"min.weight.Appr.4"=min.weight.Appr.4,"max.weight.Appr.4"=max.weight.Appr.4,"min.weight.Appr.5"=min.weight.Appr.5,"max.weight.Appr.5"=max.weight.Appr.5,"min.weight.Appr.6"=min.weight.Appr.6,"max.weight.Appr.6"=max.weight.Appr.6,
         "min.weight.Appr.7"=res_for_each_increment[[1]][[2]],"max.weight.Appr.7"=res_for_each_increment[[1]][[3]],"min.weight.Appr.8"=res_for_each_increment[[1]][[6]],"max.weight.Appr.8"=res_for_each_increment[[1]][[7]],"min.weight.Appr.9"=res_for_each_increment[[2]][[2]],"max.weight.Appr.9"=res_for_each_increment[[2]][[3]],"min.weight.Appr.10"=res_for_each_increment[[2]][[6]],"max.weight.Appr.10"=res_for_each_increment[[2]][[7]],
         "min.weight.Appr.11"=res_for_each_increment[[1]][[10]],"max.weight.Appr.11"=res_for_each_increment[[1]][[11]],"min.weight.Appr.12"=res_for_each_increment[[1]][[14]],"max.weight.Appr.12"=res_for_each_increment[[1]][[15]],"min.weight.Appr.13"=res_for_each_increment[[2]][[10]],"max.weight.Appr.13"=res_for_each_increment[[2]][[11]],"min.weight.Appr.14"=res_for_each_increment[[2]][[14]],"max.weight.Appr.14"=res_for_each_increment[[2]][[15]],"min.weight.Appr.15"=min.weight.Appr.11,"max.weight.Appr.15"=max.weight.Appr.11,"min.weight.Appr.16"=min.weight.Appr.12,"max.weight.Appr.16"=max.weight.Appr.12,"min.weight.Appr.17"=min.weight.Appr.13,"max.weight.Appr.17"=max.weight.Appr.13,"min.weight.Appr.18"=min.weight.Appr.14,"max.weight.Appr.18"=max.weight.Appr.14,
         "Sandwich.SE.A.1"=Sandwich.SE.A.1,"Sandwich.SE.A.2"=Sandwich.SE.A.2,"Sandwich.SE.A.3"=Sandwich.SE.A.3,"Sandwich.SE.A.4"=Sandwich.SE.A.4,"Sandwich.SE.A.5"=Sandwich.SE.A.5,"Sandwich.SE.A.6"=Sandwich.SE.A.6,
         "Sandwich.SE.A.7"=res_for_each_increment[[1]][[4]],"Sandwich.SE.A.8"=res_for_each_increment[[1]][[8]],"Sandwich.SE.A.9"=res_for_each_increment[[2]][[4]],"Sandwich.SE.A.10"=res_for_each_increment[[2]][[8]],
         "Sandwich.SE.A.11"=res_for_each_increment[[1]][[12]],"Sandwich.SE.A.12"=res_for_each_increment[[1]][[16]],"Sandwich.SE.A.13"=res_for_each_increment[[2]][[12]],"Sandwich.SE.A.14"=res_for_each_increment[[2]][[16]],
         "Sandwich.SE.A.15"=Sandwich.SE.A.11,"Sandwich.SE.A.16"=Sandwich.SE.A.12,"Sandwich.SE.A.17"=Sandwich.SE.A.13,"Sandwich.SE.A.18"=Sandwich.SE.A.14,
         "percentage.censored"=determine_percentage_censored(n,max_fup,data.long),"max.number.of.infections"=determine_max_number_of_infections(data.long),
         "A.19"=aux_AG[[1]],"min.weight.Appr.19"=aux_AG[[2]],"max.weight.Appr.19"=aux_AG[[3]],"Sandwich.SE.A.19"=aux_AG[[4]],"A.20"=aux_AG_IPCW[[1]],"min.weight.Appr.20"=aux_AG_IPCW[[2]],"max.weight.Appr.20"=aux_AG_IPCW[[3]],"Sandwich.SE.A.20"=aux_AG_IPCW[[4]])
  },error=function(e){ ## IF AN ERROR OCCURS, THE OUTPUT FOR EACH APPROACH WILL BE SET TO NOT AVAILABLE. (NA).
    list("A.1"=NA,"A.2"=NA,"A.3"=NA,"A.4"=NA,"A.5"=NA,"A.6"=NA,"A.7"=NA,"A.8"=NA,"A.9"=NA,"A.10"=NA,"A.11"=NA,"A.12"=NA,"A.13"=NA,"A.14"=NA,"A.15"=NA,"A.16"=NA,"A.17"=NA,"A.18"=NA,
         "min.weight.Appr.1"=NA,"max.weight.Appr.1"=NA,"min.weight.Appr.2"=NA,"max.weight.Appr.2"=NA,"min.weight.Appr.3"=NA,"max.weight.Appr.3"=NA,"min.weight.Appr.4"=NA,"max.weight.Appr.4"=NA,"min.weight.Appr.5"=NA,"max.weight.Appr.5"=NA,
         "min.weight.Appr.6"=NA,"max.weight.Appr.6"=NA,"min.weight.Appr.7"=NA,"max.weight.Appr.7"=NA,"min.weight.Appr.8"=NA,"max.weight.Appr.8"=NA,"min.weight.Appr.9"=NA,"max.weight.Appr.9"=NA,"min.weight.Appr.10"=NA,"max.weight.Appr.10"=NA,
         "min.weight.Appr.11"=NA,"max.weight.Appr.11"=NA,"min.weight.Appr.12"=NA,"max.weight.Appr.12"=NA,"min.weight.Appr.13"=NA,"max.weight.Appr.13"=NA,"min.weight.Appr.14"=NA,"max.weight.Appr.14"=NA,
         "min.weight.Appr.15"=NA,"max.weight.Appr.15"=NA,"min.weight.Appr.16"=NA,"max.weight.Appr.16"=NA,"min.weight.Appr.17"=NA,"max.weight.Appr.17"=NA,"min.weight.Appr.18"=NA,"max.weight.Appr.18"=NA,
         "Sandwich.SE.A.1"=NA,"Sandwich.SE.A.2"=NA,"Sandwich.SE.A.3"=NA,"Sandwich.SE.A.4"=NA,"Sandwich.SE.A.5"=NA,"Sandwich.SE.A.6"=NA,"Sandwich.SE.A.7"=NA,"Sandwich.SE.A.8"=NA,"Sandwich.SE.A.9"=NA,"Sandwich.SE.A.10"=NA,"Sandwich.SE.A.11"=NA,"Sandwich.SE.A.12"=NA,"Sandwich.SE.A.13"=NA,"Sandwich.SE.A.14"=NA,"Sandwich.SE.A.15"=NA,"Sandwich.SE.A.16"=NA,"Sandwich.SE.A.17"=NA,"Sandwich.SE.A.18"=NA,
         "percentage.censored"=NA,"max.number.of.infections"=NA,"A.19"=NA,"min.weight.Appr.19"=NA,"max.weight.Appr.19"=NA,"Sandwich.SE.A.19"=NA,"A.20"=NA,"min.weight.Appr.20"=NA,"max.weight.Appr.20"=NA,"Sandwich.SE.A.20"=NA)}))}
###########################################################################################################################################################################################################################
##################################################################### SET OF PARAMETERS TO BE PASSED TO MONTECARLO() ######################################################################################################
mu<--5;sigma<-1/4 ## SHAPE PARAMETER OF THE ONE PARAMETER GAMMA DISTRIBUTION;MEAN OF BASELINE COVARIATE L;VARIANCE OF BASELINE COVARIATE L
sample_grid<-c(250,500,750) ## GRID OF SAMPLE SIZES 
coef_trt_out_grid<-c(-0.4,-0.1,0.2) ## STRONG BENFICIAL VACCINE EFFECT; MILD BENEFICIAL VACCINE EFFECT; MILD ADVERSE VACCINE EFFECT
###########################################################################################################################################################################################################################
##################################################################### PART TO EVALUATE FOR THE SUPERCOMPUTER ##############################################################################################################
args<-commandArgs(trailingOnly = TRUE) ## GET PARAMETER VALUES FROM THE COMMAND LINE

MC_result<-MonteCarlo(MonteCarlo_sim,nrep=1000,ncpus=strtoi(Sys.getenv('PBS_NUM_PPN')),
                      param_list=list("n"=sample_grid,
                                      "max_fup"=as.double(args[1]),
                                      "bas_haz_out"=as.double(args[2]),
                                      "coef_trt_out"=coef_trt_out_grid,
                                      "bas_haz_cens"=as.double(args[3]),
                                      "censoring_strength_trt"=as.double(args[4]),
                                      "censoring_strength_bas_cov"=as.double(args[5]),
                                      "censoring_strength_time_varying_covariate"=as.double(args[6]),
                                      "cens_interaction"=as.double(args[7]),
                                      "shape"=as.double(args[8]),
                                      "delta"=as.double(args[9]),
                                      "sigmoid_type"=as.double(args[10])))

output_file<-paste(paste("Max.fup=",args[1],sep=""),
                   paste("Bas.haz.out=",args[2],sep=""),
                   paste("Bas.haz.cens=",args[3],sep=""),
                   paste("Cens.trt=",args[4],sep=""),
                   paste("Cens.bas.cov=",args[5],sep=""),
                   paste("Cens.time.var=",args[6],sep=""),
                   paste("Cens.interaction=",args[7],sep=""),
                   paste("Shape=",args[8],sep=""),
                   paste("delta=",args[9],sep=""),
                   paste("sigmoid_type=",args[10],sep=""),
                   sep=",") ## STORE THE RESULTS IN output_file.

saveRDS(MC_result,file=output_file) ## SAVE THE RESULTS
