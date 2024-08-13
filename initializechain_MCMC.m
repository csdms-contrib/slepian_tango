function [RJ_MCMC_THBI_state] = initializechain_MCMC(Markov_Chain)
% [RJ_MCMC_THBI_state]=INITIALIZECHAIN_MCMC(Markov_Chain)
%
% Primary function that runs the Markov chain Monte Carlo routine. Provide
% the total number of iterations you would like the computation to run for
% and existing variable files if you are starting from a stopped
% computation. You can start the chain from the beginning or from a
% previous set of saved files. The algorithm automatically saves files
% every 1000 iterations. You will need to change the first block of code to
% conform to the MCMC input variables. You will also need to change the
% block of code that preallocates PROPOSALSGLOBAL if you are running the
% chain for a number of iterations other than 1,000,000 (the default).
%
% INPUT:
%
% Markov_Chain      The current set of Markov chain parameters.
% 
% OUTPUT:
%
% RJ_MCMC_THBI_state   The new set of Markov chain parameters.
%
% Last modified by kengourley-at-arizona.edu, 8-13-2024

    RJ_MCMC_THBI_state = emptystate_MCMC();

% -------------------------------------------------------------------------
%       Return initial state as defined in settings
% -------------------------------------------------------------------------
    RJ_MCMC_THBI_state.LogPrior                 = [];
    RJ_MCMC_THBI_state.Zmax                     = Markov_Chain.InitialState.Zmax;
    RJ_MCMC_THBI_state.Zmin                     = Markov_Chain.InitialState.Zmin;
    RJ_MCMC_THBI_state.SigmaInit                = Markov_Chain.InitialState.SigmaInit;
    RJ_MCMC_THBI_state.SigmaProposed            = [];
    RJ_MCMC_THBI_state.SigmaCurrent             = [];
    RJ_MCMC_THBI_state.Proposed_Viscosity       = [];
    RJ_MCMC_THBI_state.Proposed_Layers          = [];    
    RJ_MCMC_THBI_state.NProposed_Layers         = [];
    RJ_MCMC_THBI_state.Accpetance               = [];
    RJ_MCMC_THBI_state.Rejection                = [];
    RJ_MCMC_THBI_state.Viscosity_Increase       = [];
    RJ_MCMC_THBI_state.Log_Alpha_Ratio          = [];
    RJ_MCMC_THBI_state.Misfit                   = [];
    RJ_MCMC_THBI_state.Kernel_Geoid             = [];
    RJ_MCMC_THBI_state.Observed_Geoid           = Markov_Chain.Observed_Geoid;
    RJ_MCMC_THBI_state.Proposed_Geoid           = Markov_Chain.Proposed_Geoid;
    RJ_MCMC_THBI_state.NoiseType                = Markov_Chain.NoiseType;
    RJ_MCMC_THBI_state.Current_Layers           = Markov_Chain.InitialState.Current_Layers;
    RJ_MCMC_THBI_state.Current_Viscosity        = Markov_Chain.InitialState.Current_Viscosity;
    RJ_MCMC_THBI_state.Maximum_Degree           = Markov_Chain.Maximum_Degree;
    RJ_MCMC_THBI_state.Noden                    = Markov_Chain.Noden;
    RJ_MCMC_THBI_state.Tomography               = Markov_Chain.Tomography;
    RJ_MCMC_THBI_state.Scaling                  = Markov_Chain.Scaling;
    RJ_MCMC_THBI_state.Step                     = Markov_Chain.Step;
    RJ_MCMC_THBI_state.NCurrent_Layers          = length(Markov_Chain.InitialState.Current_Layers);
    RJ_MCMC_THBI_state.Iteration                = Markov_Chain.InitialState.Iteration;
    RJ_MCMC_THBI_state.Max_Visc_value           = Markov_Chain.InitialState.Max_Visc_value;
    RJ_MCMC_THBI_state.Min_Visc_value           = Markov_Chain.InitialState.Min_Visc_value;
    RJ_MCMC_THBI_state.Fixed_Visc_Contrast      = Markov_Chain.InitialState.Fixed_Visc_Contrast;
    RJ_MCMC_THBI_state.Fixed_Visc_depth         = Markov_Chain.InitialState.Fixed_Visc_depth;
    RJ_MCMC_THBI_state.Uncertainty_Current      = [];
    RJ_MCMC_THBI_state.Uncertainty_Proposed     = [];
    RJ_MCMC_THBI_state.RMS_Error                = [];
% -------------------------------------------------------------------------  
% Compute and verify initial predicted geoid to be used in subsequent
% iterations
% -------------------------------------------------------------------------
disp('Compute and verify initial predicted geoid to be used in subsequent iterations')
    
    [Proposed_Layers_,ind]   = sort([RJ_MCMC_THBI_state.Fixed_Visc_depth RJ_MCMC_THBI_state.Current_Layers]);
    Proposed_visco_         = [RJ_MCMC_THBI_state.Fixed_Visc_Contrast RJ_MCMC_THBI_state.Current_Viscosity];
    
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

% -------------------------------------------------------------------------  
% Update current state with hyperparameters to be used in subsequent
% iterations
% -------------------------------------------------------------------------
    RJ_MCMC_THBI_state.HyperSigma_I         = Markov_Chain.HyperSigma_I;
    RJ_MCMC_THBI_state.HyperSigma_II        = Markov_Chain.HyperSigma_II;
    
        switch RJ_MCMC_THBI_state.NoiseType
            case 1
                disp({'Single error parameter estimation selected'})
                sigma_current = RJ_MCMC_THBI_state.SigmaInit;
                RJ_MCMC_THBI_state.SigmaCurrent    =  sigma_current;
                RJ_MCMC_THBI_state.NoiseType        = 1;
           case 2
                disp({'Error parameter estimation using min and max sigma selected'})
                sigma_current  = RJ_MCMC_THBI_state.SigmaInit;
                RJ_MCMC_THBI_state.SigmaCurrent    = sigma_current;
                RJ_MCMC_THBI_state.NoiseType        = 2;
        end
        
   [LogLikelihood,~,Uncertainty,Misfit,Residual,RMS_Error] = likelihood_MCMC(RJ_MCMC_THBI_state);
    RJ_MCMC_THBI_state.LogLikelihood        = LogLikelihood;
    RJ_MCMC_THBI_state.Uncertainty_Current  = Uncertainty;
    RJ_MCMC_THBI_state.Misfit               = Misfit;
    RJ_MCMC_THBI_state.Residual             = Residual;
    RJ_MCMC_THBI_state.RMS_Error            = RMS_Error;
    NCurrent_layers = RJ_MCMC_THBI_state.NCurrent_Layers;
    RJ_MCMC_THBI_state.SigmaProposed        = RJ_MCMC_THBI_state.SigmaInit;
   [RJ_MCMC_THBI_state.LogPosterior]  = logposterior_MCMC(LogLikelihood,NCurrent_layers);
end
