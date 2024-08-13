function [LogPosterior] = logposterior_MCMC(Log_proposed_Likehood,Proposed_Nlayers)

    Log_Proposed_Likehood  = Log_proposed_Likehood;
     
     RJMCMC_prior           = RJ_MCMC_THBI_Prior(Proposed_Nlayers);
     LogPosterior           = Log_Proposed_Likehood + RJMCMC_prior;

end