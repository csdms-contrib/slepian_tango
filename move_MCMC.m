function [RJ_MCMC_THBI_state,accepted]=move_MCMC(RJMCMC_Chain_Array_Updates,Markov_Chain)
%**************************************************************************
% 3). Location Move or change perturbs the location of a current layer
%**************************************************************************
RJ_MCMC_THBI_state = RJMCMC_Chain_Array_Updates;
RJ_MCMC_THBI_state_old          = RJMCMC_Chain_Array_Updates;
RJ_MCMC_THBI_state.Iteration = Markov_Chain.Iteration;
Log_Likelihood_Current = RJ_MCMC_THBI_state.LogLikelihood;
Uncertainty_Current     = RJ_MCMC_THBI_state.Uncertainty_Current;
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
    Zmax                = RJ_MCMC_THBI_state.Zmax;
    Zmin                = RJ_MCMC_THBI_state.Zmin;
    Fixed_Visc_depth    = RJ_MCMC_THBI_state.Fixed_Visc_depth;
    Fixed_Visc_Contrast = RJ_MCMC_THBI_state.Fixed_Visc_Contrast;
    NoiseType           = RJ_MCMC_THBI_state.NoiseType;
    Layer_interface_std = Markov_Chain.Layer_interface_std;
    Move_type           = Markov_Chain.Move_type;
    drmin               = Markov_Chain.drmin;
%**************************************************************************
%   Moves a Kth layer in the current model to a new location - Select a  
%   location a random and purturb it.
%**************************************************************************
    current_layer       = Current_Layers;
    idx                 = randperm(length(current_layer),1);
    Kth_layer           = current_layer(idx);
%     [Kth_layer, idx]    = datasample(current_layer,1);
    
if Move_type == 1 % Move step type I using normal distibution
    K_Move_Layer = Kth_layer+Layer_interface_std*randn;
elseif Move_type == 2 % Move step type II using uniform distibution
    K_Move_Layer = Kth_layer+unifrnd(-1,1);
end
         Proposed_layers = current_layer;
         Proposed_layers(idx)=K_Move_Layer;
         Vics_proposed = Current_Viscous(idx);

% Check for the minimum allowable layer
% checks to see that layers have not gotten too thin (which can
% cause convergence problems
if length(Proposed_layers)>1
    tmp_nrad = Proposed_layers;
    tmp_nrad = sort(tmp_nrad);
    tmp_nrad = min(diff(tmp_nrad));
    indx2 = tmp_nrad > drmin;
elseif length(Proposed_layers)==1
    indx2 = 1;
end

indx1 = ((K_Move_Layer >= Zmin) & (K_Move_Layer <= Zmax));
if indx1 == 1 && indx2 == 1
    Proposed_visco = Current_Viscous;
    Proposed_visco(idx) = Vics_proposed;
    Proposed_Layers = Proposed_layers;   

    % Update proposed structure
    RJ_MCMC_THBI_state.Proposed_Viscosity    =  Proposed_visco;
    RJ_MCMC_THBI_state.Proposed_Layers       =  Proposed_Layers;
    RJ_MCMC_THBI_state.NProposed_Layers      =  length(Proposed_Layers);   
%%
%**************************************************************************
%Compute alpha prior probability P(m) = k+1/k'+1
%**************************************************************************
    Proposed_Nlayers    = length(Proposed_Layers);
    LogAlpha_Prior      = log((Current_Nlayers+1)/(Proposed_Nlayers+1)); % multiplying factor k+1/k'+1
if Markov_Chain.Proposal_sample
%**************************************************************************
%Compute proposed geoid from birth layer and visocity
%**************************************************************************
   [Proposed_Layers_,ind]   = sort([Fixed_Visc_depth Proposed_Layers]);
    Proposed_visco_         = [Fixed_Visc_Contrast Proposed_visco];
    
    % For pushing top viscosity to surface
    if max(Proposed_Layers_)==RJ_MCMC_THBI_state.Zmax
       Proposed_visco_         = exp(Proposed_visco_(ind));
    elseif max(Proposed_Layers_) < RJ_MCMC_THBI_state.Zmax
      Proposed_visco_          = Proposed_visco_(ind);
      Proposed_Layers_         = [Proposed_Layers_ 1];
      Proposed_visco_          = exp([Proposed_visco_ Proposed_visco_(end)]);
    end

%[~,state] for regional
%[state,~] for global
   [~,RJ_MCMC_THBI_state.Proposed_Geoid] = geoidkernel_MCMC(RJ_MCMC_THBI_state.Maximum_Degree,RJ_MCMC_THBI_state.Noden,[Proposed_Layers_; Proposed_visco_]',RJ_MCMC_THBI_state.Region,RJ_MCMC_THBI_state.Tomography,RJ_MCMC_THBI_state.Scaling,RJ_MCMC_THBI_state.Step);

%**************************************************************************
% Compute log-likelihood, Uncertainty and Misfit
%**************************************************************************
[Log_Likehood_Proposed,Prefactor,Uncertainty_Proposed,Misfit,Residual,RMS_Error] = likelihood_MCMC(RJ_MCMC_THBI_state);
%**************************************************************************
% Calculate ratio term needed to evaluate the acceptance ratio and check for
%acceptance or rejection
%**************************************************************************
    Log_Alpha_Ratio = Prefactor + Log_Likehood_Proposed - Log_Likelihood_Current + LogAlpha_Prior;
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
    Move = log(rand);
if (isfinite(Log_Likelihood_Current) && isfinite(Log_Likehood_Proposed))    
    if (Log_Alpha_Ratio >0 || Log_Alpha_Ratio > Move)
        accepted = 1;
        RJ_MCMC_THBI_state.LogLikelihood        = Log_Likehood_Proposed;
        RJ_MCMC_THBI_state.Current_Viscosity    = Proposed_visco;
        RJ_MCMC_THBI_state.Current_Layers       = Proposed_Layers;
        RJ_MCMC_THBI_state.NCurrent_Layers      = length(Proposed_Layers);
        RJ_MCMC_THBI_state.AMove_count          = 1;
        RJ_MCMC_THBI_state.RMove_count          = 0;
        RJ_MCMC_THBI_state.Acceptance           = 1;
        RJ_MCMC_THBI_state.Rejection            = 0;
        RJ_MCMC_THBI_state.Log_Alpha_Ratio      = Log_Alpha_Ratio;
        RJ_MCMC_THBI_state.LogPrior             = LogAlpha_Prior;
        RJ_MCMC_THBI_state.Misfit               = Misfit;
        RJ_MCMC_THBI_state.Residual             = Residual;        
        RJ_MCMC_THBI_state.Kernel_Geoid         = RJ_MCMC_THBI_state.Proposed_Geoid;
        RJ_MCMC_THBI_state.Uncertainty_Current  = Uncertainty_Proposed;
        RJ_MCMC_THBI_state.RMS_Error            = RMS_Error;
  else
        accepted = 0;
        RJ_MCMC_THBI_state_old.AMove_count      = 0;
        RJ_MCMC_THBI_state_old.RMove_count      = 1;
        RJ_MCMC_THBI_state_old.Acceptance       = 0;
        RJ_MCMC_THBI_state_old.Rejection        = 1;
        RJ_MCMC_THBI_state_old.Iteration        = RJ_MCMC_THBI_state.Iteration;
        RJ_MCMC_THBI_state                      = RJ_MCMC_THBI_state_old;
    end
end    
else
        accepted = -1;
        RJ_MCMC_THBI_state = RJ_MCMC_THBI_state_old;
end
end

