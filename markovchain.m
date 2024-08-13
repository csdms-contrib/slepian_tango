function [PROPOSALSGLOBAL] = markovchain(Step,Continuing_Step,cont_Proposal)
% [PROPOSALSGLOBAL]=MARKOVCHAIN(Step,Continuing_Step,cont_Proposal)
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
% Depending on whether or not you're running a global or regional
% inversion, you will also have to change your call of 'geoidkernel_MCMC'
% in 'initializechain_MCMC', 'birth_MCMC', 'death_MCMC', 'move_MCMC',
% 'noise_MCMC', and value_MCMC'.
%
% INPUT:
%
% Step              Total number of iterations for which the MCMC will run.
% Continuing_Step   The previous state of the MCMC model if you are
%                   starting from a saved step. This is saved by the 
%                   algorithm as 'RJMCMC_Chain_Array_UpdatesXXXX.m'.
%                   If starting a new chain, this should be '[]'.
% cont_Proposal     All of the previous ensemble information to be used
%                   when computing statistics when chain is complete. This
%                   is saved by the algorithm as 'PROPOSALS_GLOBALXXXX.m'.
%                   If starting a new chain, this should be '[]'.
% 
% OUTPUT:
%
% PROPOSALSGLOBAL   An object containing all of the ensemble variables to
%                   be used for statistical calculations.
%
% Last modified by kengourley-at-arizona.edu, 8-13-2024

% Set total number of iterations
Markov_Chain.Nsample = Step;
% Set burn in period
Markov_Chain.BurnIn             = 5000;
% Load observed geoid
Markov_Chain.Observed_Geoid = Observed_Geoid.Observed20;
% Initialize predicted geoid
Markov_Chain.Proposed_Geoid = [];
% Initialize iteration number
Markov_Chain.Iteration = [];
R=6371;         % Radius of Earth [km]
CMB=3476;       % Radius of Earth's CMB [km]
Markov_Chain.Maximum_Degree = 0;  % Maximum spherical harmonic degree
Markov_Chain.Noden = 0;           % Depth above which there is zero density contribution [km]
Markov_Chain.Tomography = '';     % Name of tomography model to use 
Markov_Chain.Scaling = 0;         % Scaling profile to use
Markov_Chain.Step = 10;           % Depth step to use [km]
Markov_Chain.Region = '';         % Name of Slepian computation region ('' for global)

min_visc=log([1e-2]);     % value of minimum viscosity for proposed layers
max_visc=log([1e2]);      % value of maximum viscosity for proposed layers
fixed_visc=[log([1e0])];  % value of fixed viscosity layer for all iterations
fixed_depth=[(R-180)/R];  % depth of fixed viscosity layer

%-------------------------------------------------------------------------%
           % set variables for the five types of Markov chain steps
%-------------------------------------------------------------------------%
% 1.Birth
Max_depth = R/R;
Min_depth = (CMB+50)/R;
Max_num_layers = 9; % maximum total number of layers (excluding fixed viscosity layer)
Max_num_layer  = [];
for Iteration = 2:Max_num_layers
    Max_num_layer = [Max_num_layer; Iteration*ones(Markov_Chain.BurnIn*Iteration,1)];
end
if(length(Max_num_layer)>Markov_Chain.Nsample)
    disp('Not enough steps to get to Max_num_vor');
else
    Max_num_layer(length(Max_num_layer):Markov_Chain.Nsample) = Max_num_layers;
end
Markov_Chain.Max_num_layer=Max_num_layer;
%-------------------------------------------------------------------------%
% 2.Death
% NB: This step is applied when number of layers are greater than 1
%-------------------------------------------------------------------------%
% 3.Value
% Can vary to see effect of imposed error on different geoids
Markov_Chain.Vicosity_change_std = [0.1];
%-------------------------------------------------------------------------%
% 4.Move --> Set the standard deviations of the perturbations to the coefficients
Markov_Chain.Layer_interface_std    = [0.1];
Markov_Chain.Move_type              = 1; % move with normal=1 or uniform=2 distribution
% How small can a layer become?
Markov_Chain.drmin                  = 1e-5;
%-------------------------------------------------------------------------%
% 5.Noise --> Set the standard deviations of the perturbations to the variance
% (uncertainty hyperparameter) in log units
Markov_Chain.Var_change = [1e-1];

Markov_Chain.NoiseType          = 1;   % whether to use correlated or uncorrelated noise
Markov_Chain.HyperSigma_I       = [];  % Noise type 1
Markov_Chain.HyperSigma_II      = [];  % Noise type 2

%--------------------------------------------------------------------------
% Collect all output for analysis                                         %
%--------------------------------------------------------------------------

% preallocate space for variables
PROPOSALSGLOBAL(1e6).LogLikelihood      = [];
PROPOSALSGLOBAL(1e6).Misfit             = [];
PROPOSALSGLOBAL(1e6).Proposed_Layers    = [];
PROPOSALSGLOBAL(1e6).Proposed_Viscosity = [];
PROPOSALSGLOBAL(1e6).Acceptance         = [];
PROPOSALSGLOBAL(1e6).RMS_Error          = [];
PROPOSALSGLOBAL(1e6).SigmaProposed      = [];
PROPOSALSGLOBAL(1e6).NProposed_Layers   = [];

% preallocate space for final variables
PROPOSALSGLOBALF(1e4).LogLikelihood      = [];
PROPOSALSGLOBALF(1e4).Misfit             = [];
PROPOSALSGLOBALF(1e4).Residual           = [];
PROPOSALSGLOBALF(1e4).Proposed_Layers    = [];
PROPOSALSGLOBALF(1e4).Proposed_Viscosity = [];
PROPOSALSGLOBALF(1e4).Acceptance         = [];
PROPOSALSGLOBALF(1e4).RMS_Error          = [];
PROPOSALSGLOBALF(1e4).SigmaProposed      = [];
PROPOSALSGLOBALF(1e4).NProposed_Layers   = [];
PROPOSALSGLOBALF(1e4).Proposed_Geoid     = [];

prevaccept=0;

if isempty(Continuing_Step) % beginning run
    %--------------------------------------------------------------------------
    % Acceptance and rejection count for the different Markov chain steps     %
    %--------------------------------------------------------------------------
    Markov_Chain.ABirth_count = 0;
    Markov_Chain.RBirth_count = 0;
    Markov_Chain.ADeath_count = 0;
    Markov_Chain.RDeath_count = 0;
    Markov_Chain.AValue_count = 0;
    Markov_Chain.RValue_count = 0;
    Markov_Chain.AMove_count  = 0;
    Markov_Chain.RMove_count  = 0;
    Markov_Chain.ANoise_count = 0;
    Markov_Chain.RNoise_count = 0;
    % -------------------------------------------------------------------------
    %   Initial State                                                         %
    % -------------------------------------------------------------------------
    Markov_Chain.InitialState.Zmax                = [Max_depth];        % Maximum allowed depth (surface)
    Markov_Chain.InitialState.Zmin                = [Min_depth];        % Minimum allowed depth (CMB)
    Markov_Chain.InitialState.SigmaInit           = [5];                % Initial Noise parameter
    Markov_Chain.InitialState.Iteration           = [1];                % start of sampling
    Markov_Chain.InitialState.Current_Layers      = [(CMB+1000)/R];     % Initial prescribed depth
    Markov_Chain.InitialState.Current_Viscosity   = log([1e1]);         % Initial prescribed viscosity contrast
    Markov_Chain.InitialState.Fixed_Visc_Contrast = fixed_visc;         % Fixed Viscosity contrast 
    Markov_Chain.InitialState.Fixed_Visc_depth    = fixed_depth;        % Fixed Viscosity contrast Layer
    Markov_Chain.InitialState.Min_Visc_value      = min_visc;           % Minimum acceptable viscosity value
    Markov_Chain.InitialState.Max_Visc_value      = max_visc;           % Maximum acceptable viscosity value

    % -------------------------------------------------------------------------
    %   Estimate likelihood of the initial state and store at position 1        
    % -------------------------------------------------------------------------
    
    RJMCMC_Chain_Array_Updates = emptystate_MCMC();
    RJMCMC_Chain_Array_Updates(Markov_Chain.InitialState.Iteration) = initializechain_MCMC(Markov_Chain);

    % -------------------------------------------------------------------------
    % Initialized allocation to store MCMC simulation for posterior analysis
    % -------------------------------------------------------------------------

    IT = 1;

else % load from saved file 
    Markov_Chain.Nsample = Step;
    RJMCMC_Chain_Array_Updates = Continuing_Step.RJMCMC_Chain_Array_Updates;		
    Markov_Chain.Iteration = RJMCMC_Chain_Array_Updates.Iteration;		
    IT = Markov_Chain.Iteration+1;
    
    for nn = 1:IT-1
        PROPOSALSGLOBAL(nn).LogLikelihood      = cont_Proposal.PROPOSALSGLOBAL(nn).LogLikelihood;
        PROPOSALSGLOBAL(nn).Misfit             = cont_Proposal.PROPOSALSGLOBAL(nn).Misfit;
        PROPOSALSGLOBAL(nn).Proposed_Layers    = cont_Proposal.PROPOSALSGLOBAL(nn).Proposed_Layers;
        PROPOSALSGLOBAL(nn).Proposed_Viscosity = cont_Proposal.PROPOSALSGLOBAL(nn).Proposed_Viscosity;
        PROPOSALSGLOBAL(nn).Accpetance         = cont_Proposal.PROPOSALSGLOBAL(nn).Acceptance;
        PROPOSALSGLOBAL(nn).RMS_Error          = cont_Proposal.PROPOSALSGLOBAL(nn).RMS_Error;
        PROPOSALSGLOBAL(nn).SigmaProposed      = cont_Proposal.PROPOSALSGLOBAL(nn).SigmaProposed;
        PROPOSALSGLOBAL(nn).NProposed_Layers   = cont_Proposal.PROPOSALSGLOBAL(nn).NProposed_Layers;
        
        prevaccept=prevaccept+cont_Proposal.PROPOSALSGLOBAL(nn).Acceptance;
    end
end

% Sample Prior 
Markov_Chain.Proposal_sample = true; %(false = Prior sample)
% start saving sampling
Start_save = 9e5;
Model_varify = true; %(if false only final solution will be output)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Do RJMCMC Transdimensional Bayesian Inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
Accept = 0;
arbC = 1/5; % equal probability of any Markov chain step
accepted = -1;
Iteration2 = 1;

% Print run parameters
fprintf('Start and end iterations: %d to %d\n', IT, Step)
fprintf('Tomography model: %s\n', Markov_Chain.Tomography)
fprintf('Region: %s\n', Markov_Chain.Region)
if ischar(Markov_Chain.Scaling) || isstring(Markov_Chain.Scaling)
    disp(['Scaling: ' Markov_Chain.Scaling])
else
    disp(['Scaling: ' num2str(Markov_Chain.Scaling)])
end
fprintf('Bandwidth: 2 to %d\n', Markov_Chain.Maximum_Degree)
fprintf('Lithosphere thickness: %d\n', Markov_Chain.Noden)
fprintf('Min and max depths: %f to %f\n', Min_depth, Max_depth)
fprintf('Max layers: %d\n', Max_num_layers)
fprintf('Viscosity range: %e to %e\n', exp(1)^min_visc, exp(1)^max_visc)
fprintf('Fixed viscosity Contrast: %e at %f\n', exp(1)^fixed_visc, fixed_depth)

for Iteration = IT:Markov_Chain.Nsample   
     while (accepted < 0)
        Markov_Chain.Iteration = Iteration;
        decision=rand(1); 
       if ((decision <= arbC))
           [RJ_MCMC_THBI_state,accepted]=birth_MCMC(RJMCMC_Chain_Array_Updates,Markov_Chain);
       elseif ((decision > arbC)) && ((decision <= 2*arbC))
           [RJ_MCMC_THBI_state,accepted]=death_MCMC(RJMCMC_Chain_Array_Updates,Markov_Chain);
        elseif ((decision > 2*arbC)) && ((decision <= 3*arbC))
            [RJ_MCMC_THBI_state,accepted]=value_MCMC(RJMCMC_Chain_Array_Updates,Markov_Chain);
        elseif ((decision > 3*arbC)) && ((decision <= 4*arbC))
            [RJ_MCMC_THBI_state,accepted]=move_MCMC(RJMCMC_Chain_Array_Updates,Markov_Chain);
       elseif ((decision > 4*arbC)) && ((decision <= 5*arbC))
           [RJ_MCMC_THBI_state,accepted]=noise_MCMC(RJMCMC_Chain_Array_Updates,Markov_Chain);
       end
     end
      RJMCMC_Chain_Array_Updates = RJ_MCMC_THBI_state;
% -------------------------------------------------------------------------
% Collect final results at specific intervals and save to final PROPOSALSGLOBAL       
% -------------------------------------------------------------------------
	if (Iteration > Start_save) && (mod(Iteration,10) == 0)
		PROPOSALSGLOBALF(Iteration2).LogLikelihood      = RJ_MCMC_THBI_state.LogLikelihood;
		PROPOSALSGLOBALF(Iteration2).Misfit             = RJ_MCMC_THBI_state.Misfit;
		PROPOSALSGLOBALF(Iteration2).Residual           = RJ_MCMC_THBI_state.Residual;
		PROPOSALSGLOBALF(Iteration2).Proposed_Layers    = RJ_MCMC_THBI_state.Proposed_Layers;
		PROPOSALSGLOBALF(Iteration2).Proposed_Viscosity = RJ_MCMC_THBI_state.Proposed_Viscosity;
		PROPOSALSGLOBALF(Iteration2).Accpetance         = RJ_MCMC_THBI_state.Accpetance;
		PROPOSALSGLOBALF(Iteration2).RMS_Error          = RJ_MCMC_THBI_state.RMS_Error;
		PROPOSALSGLOBALF(Iteration2).SigmaProposed      = RJ_MCMC_THBI_state.SigmaProposed;
		PROPOSALSGLOBALF(Iteration2).NProposed_Layers   = RJ_MCMC_THBI_state.NProposed_Layers;
        PROPOSALSGLOBALF(Iteration2).Proposed_Geoid     = RJ_MCMC_THBI_state.Proposed_Geoid;
		Iteration2 = Iteration2+1;
	end
% -------------------------------------------------------------------------
% Output simulation at specified intervals and analyse to see the
% evolution of the model
% -------------------------------------------------------------------------

    if (Model_varify) && (mod(Iteration,1) == 0)
        PROPOSALSGLOBAL(Iteration).LogLikelihood      = RJ_MCMC_THBI_state.LogLikelihood;
        PROPOSALSGLOBAL(Iteration).Misfit             = RJ_MCMC_THBI_state.Misfit;
        PROPOSALSGLOBAL(Iteration).Proposed_Layers    = RJ_MCMC_THBI_state.Proposed_Layers;
        PROPOSALSGLOBAL(Iteration).Proposed_Viscosity = RJ_MCMC_THBI_state.Proposed_Viscosity;
        PROPOSALSGLOBAL(Iteration).Accpetance         = RJ_MCMC_THBI_state.Accpetance;
        PROPOSALSGLOBAL(Iteration).RMS_Error          = RJ_MCMC_THBI_state.RMS_Error;
        PROPOSALSGLOBAL(Iteration).SigmaProposed      = RJ_MCMC_THBI_state.SigmaProposed;
        PROPOSALSGLOBAL(Iteration).NProposed_Layers   = RJ_MCMC_THBI_state.NProposed_Layers;
    end
        
% -------------------------------------------------------------------------
% Estimate the counts and rate of acceptance    
% -------------------------------------------------------------------------
    Accept = Accept+accepted+prevaccept; % Report every 1000th sampled 
    % save Markov chain state variables every 1000 iterations
    if mod(Iteration,1000) == 0
            disp(['Saving file for iteration ' num2str(Iteration)])
            save(['PROPOSALS_GLOBAL_' num2str(Iteration) '.mat'],'PROPOSALSGLOBAL')
            save(['RJMCMC_Chain_Array_Updates' num2str(Iteration) '.mat'],'RJMCMC_Chain_Array_Updates')
    end
    disp(['Iteration: ' num2str(Iteration) '; Acceptance rate: ' num2str(Accept/Iteration)])
    accepted = -1;
    RJ_MCMC_THBI_state =[];
 end
 toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the final inversion results to file for post processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	save(['PROPOSALSGLOBAL_' num2str(Iteration) '_final_.mat'],'PROPOSALSGLOBALF')
end
