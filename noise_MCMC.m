function [RJ_MCMC_THBI_state,accepted] = noise_MCMC(RJMCMC_Chain_Array_Updates,Markov_Chain)
%************************************************************************
%* Change Noise
%************************************************************************
RJ_MCMC_THBI_state              = RJMCMC_Chain_Array_Updates;
RJ_MCMC_THBI_state_old          = RJMCMC_Chain_Array_Updates;
RJ_MCMC_THBI_state.Iteration    = Markov_Chain.Iteration;
Log_Likelihood_Current          = RJ_MCMC_THBI_state.LogLikelihood;
Uncertainty_Current             = RJ_MCMC_THBI_state.Uncertainty_Current;
% Reset acceptance and rejection counts
RJ_MCMC_THBI_state.ABirth_count = [];
RJ_MCMC_THBI_state.RBirth_count = [];
RJ_MCMC_THBI_state.ADeath_count = [];
RJ_MCMC_THBI_state.RDeath_count = [];
RJ_MCMC_THBI_state.AValue_count = [];
RJ_MCMC_THBI_state.RValue_count = [];
RJ_MCMC_THBI_state.AMove_count  = [];
RJ_MCMC_THBI_state.RMove_count  = [];
RJ_MCMC_THBI_state.ANoise_count = [];
RJ_MCMC_THBI_state.RNoise_count = [];
RJ_MCMC_THBI_state.Rejection    = [];
RJ_MCMC_THBI_state.Acceptance   = [];
RJ_MCMC_THBI_state_old.ABirth_count = [];
RJ_MCMC_THBI_state_old.RBirth_count = [];
RJ_MCMC_THBI_state_old.ADeath_count = [];
RJ_MCMC_THBI_state_old.RDeath_count = [];
RJ_MCMC_THBI_state_old.AValue_count = [];
RJ_MCMC_THBI_state_old.RValue_count = [];
RJ_MCMC_THBI_state_old.AMove_count  = [];
RJ_MCMC_THBI_state_old.RMove_count  = [];
RJ_MCMC_THBI_state_old.ANoise_count = [];
RJ_MCMC_THBI_state_old.RNoise_count = [];
RJ_MCMC_THBI_state_old.Rejection    = [];
RJ_MCMC_THBI_state_old.Acceptance   = [];
%**************************************************************************
%  Current distribution parameters
%**************************************************************************
    Current_Layers      = [RJ_MCMC_THBI_state.Current_Layers];
    Current_Nlayers     = length([RJ_MCMC_THBI_state.Current_Layers]);  
    Current_Viscous     = [RJ_MCMC_THBI_state.Current_Viscosity];
    Fixed_Visc_depth    = RJ_MCMC_THBI_state.Fixed_Visc_depth;
    Fixed_Visc_Contrast = RJ_MCMC_THBI_state.Fixed_Visc_Contrast;
    NoiseType           = RJ_MCMC_THBI_state.NoiseType;
    Proposed_Geoid      = RJ_MCMC_THBI_state.Proposed_Geoid;
    Var_change          = Markov_Chain.Var_change;
%**************************************************************************
% Change Noise / Sigma randomly using uniform or gaussian distribution.
%************************************************************************** 
        switch NoiseType
            case 1
                disp({'Single error parameter estimation'})
                sigma_current = RJ_MCMC_THBI_state.SigmaCurrent;
                sigma_proposed = RJ_MCMC_THBI_state.SigmaProposed;
                Var_change = 0.3;
                Var_min  = 1e-5;
                Var_max = 10;
                dvar = Var_change*randn;
                %Varfakt = Var2/Var1;
                %L2 = sqrt(Var2)*L
                if (10^(log10(sigma_current)+dvar)>Var_max || 10^(log10(sigma_proposed)+dvar)<Var_min)
                    indx = 0;
                else
                    sigma_proposed = 10^(log10(sigma_current) + dvar);
                    RJ_MCMC_THBI_state.SigmaProposed = sigma_proposed;
                    indx = 1;
                end
           case 2
                disp({'Multiple error parameter estimation'})
                sigma_prior = [0.1:0.1:10]; % Apply flat prior noise
                idx = randperm(length(sigma_prior),1);
                sigma_proposed= sigma_prior(idx);
                RJ_MCMC_THBI_state.SigmaProposed = sigma_proposed;
                RJ_MCMC_THBI_state.NoiseType = 2;
                indx = 1;
        end

if (indx == 1)        
        
    Proposed_visco = Current_Viscous;
    Proposed_Layers = Current_Layers;
    
    % Update proposed structure
    RJ_MCMC_THBI_state.Proposed_Viscosity    =  Proposed_visco;
    RJ_MCMC_THBI_state.Proposed_Layers       =  Proposed_Layers;
    RJ_MCMC_THBI_state.NProposed_Layers      =  length(Proposed_Layers);  
%%
%**************************************************************************
% Compute alpha prior probability P(m) = k+1/k'+1
%**************************************************************************
    Proposed_Nlayers    = length(Proposed_Layers);
    LogAlpha_Prior      = log((Current_Nlayers+1)/(Proposed_Nlayers+1));
%if (indx == 1)
if Markov_Chain.Proposal_sample
    %**************************************************************************
    %Compute proposed geoid from birth layer and visocity
    %**************************************************************************
    if isempty(Proposed_Geoid)
       [Proposed_Layers_,ind]   = sort([Fixed_Visc_depth Proposed_Layers]);
        Proposed_visco_         = [Fixed_Visc_Contrast Proposed_visco];
        Proposed_visco_         = exp(Proposed_visco_(ind));
        %[~,state] for regional
        %[state,~] for global
        [~,RJ_MCMC_THBI_state.Proposed_Geoid] = geoidkernel_MCMC(RJ_MCMC_THBI_state.Maximum_Degree,RJ_MCMC_THBI_state.Noden,[Proposed_Layers_; Proposed_visco_]',RJ_MCMC_THBI_state.Region,RJ_MCMC_THBI_state.Tomography,RJ_MCMC_THBI_state.Scaling,RJ_MCMC_THBI_state.Step);
    end
    %**************************************************************************
    % Compute log-likelihood, Uncertainty, and Misfit
    %**************************************************************************
    [Log_Likehood_Proposed,Prefactor, Uncertainty_Proposed,Misfit,Residual,RMS_Error] = likelihood_MCMC(RJ_MCMC_THBI_state);
    %**************************************************************************
    % Calculate ratio term needed to evaluate the acceptance ratio and check 
    %for acceptance or rejection
    %**************************************************************************
        Log_Alpha_Ratio =  Prefactor + Log_Likehood_Proposed - Log_Likelihood_Current + LogAlpha_Prior;
    %**************************************************************************
    % calculate the posterior distrobution
    %**************************************************************************
    RJ_MCMC_THBI_state.LogPosterior = logposterior_MCMC(Log_Likehood_Proposed,Proposed_Nlayers);
    %
else
	Log_Alpha_Ratio = 1;
	Log_Likehood_Proposed = Log_Likelihood_Current;
	Uncertainty_Proposed = Uncertainty_Current;
	Misfit = [];
	Residual = [];
	RMS_Error = [];
end
% Metropolis Hasting method:
% =========================
    Noise = log(rand);
if (isfinite(Log_Likelihood_Current) && isfinite(Log_Likehood_Proposed))    
    if (Log_Alpha_Ratio > 0 || Log_Alpha_Ratio > Noise)
        accepted = 1;
        RJ_MCMC_THBI_state.LogLikelihood        = Log_Likehood_Proposed;
        RJ_MCMC_THBI_state.Current_Viscosity    = Proposed_visco;
        RJ_MCMC_THBI_state.Current_Layers       = Proposed_Layers;
        RJ_MCMC_THBI_state.NCurrent_Layers      = length(Proposed_Layers);
        RJ_MCMC_THBI_state.ANoise_count         = 1;
        RJ_MCMC_THBI_state.RNoise_count         = 0;
        RJ_MCMC_THBI_state.Acceptance           = 1;
        RJ_MCMC_THBI_state.Rejection            = 0;
        RJ_MCMC_THBI_state.Log_Alpha_Ratio      = Log_Alpha_Ratio;
        RJ_MCMC_THBI_state.LogPrior             = LogAlpha_Prior;
        RJ_MCMC_THBI_state.Misfit               = Misfit;
        RJ_MCMC_THBI_state.Residual             = Residual;        
        RJ_MCMC_THBI_state.Kernel_Geoid         = RJ_MCMC_THBI_state.Proposed_Geoid;
        RJ_MCMC_THBI_state.Uncertainty_Current  = Uncertainty_Proposed;
        RJ_MCMC_THBI_state.Uncertainty_Proposed = Uncertainty_Proposed;
        RJ_MCMC_THBI_state.RMS_Error            = RMS_Error;
        RJ_MCMC_THBI_state.SigmaCurrent         = RJ_MCMC_THBI_state.SigmaProposed;
    else
        accepted = 0;
        RJ_MCMC_THBI_state_old.ANoise_count     = 0;
        RJ_MCMC_THBI_state_old.RNoise_count     = 1;
        RJ_MCMC_THBI_state_old.Acceptance       = 0;
        RJ_MCMC_THBI_state_old.Rejection        = 1;
        RJ_MCMC_THBI_state_old.Iteration        = RJ_MCMC_THBI_state.Iteration;
        RJ_MCMC_THBI_state                      = RJ_MCMC_THBI_state_old;

    end
end
else
        accepted = -1;
        RJ_MCMC_THBI_state                      = RJ_MCMC_THBI_state_old;

end
end
