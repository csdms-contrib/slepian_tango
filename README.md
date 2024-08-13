# slepian_tango

Version 0.1
The current commit is a tentative set of files that still need further commenting for readibility.

These functions allow the user to create a Markov chain Monte Carlo inversion for the viscosity structure of Earth's mantle using constraints from the observed geoid. These inversions can be either global or regional, with the regional techniques using Slepian functions created for repositories alpha through delta. 

The primary function in this repository is markovchain.m, which is where the user sets all of the initial parameters for the MCMC inversion and loads the file that contains the observed geoid. This function returns a set of ensemble variables that are used to create posterior statistics. This function can be started at the first iteration of the Markov chain or at an intermediate iteration if the inversion was interrupted. A sample script for running this function is provided below.

% set the Matlab path for the various local Slepian repositories at the top of the script

Step = 1000000; % set the maximum number of iterations for the inversion

% set the next iteration to be computed; jj=1 if starting new run
ii=[]; % continuing run  
jj=1;  % starting a new run

% load state variables if continuing previous run
if ~isempty(ii)
        RJMCMC_Model = ''; % load saved 'RJMCMC_Chain_Array_UpdatesXXXX.m' file
        cont_Proposal = ''; % load saved 'PROPOSALS_GLOBAL_XXXX.m' file
        Iteration = ii;
else
        Iteration = 0; % new run
end

while jj < 2
    if (jj == 1 && Iteration == 0) % new run
        Continuing_Step = [];
        cont_Proposal = [];
        PROPOSALSGLOBAL = markovchain(Step,Continuing_Step,cont_Proposal);
    elseif (jj>=1 && Iteration>0) % continuing from saved run
        Continuing_Step = RJMCMC_Model;
        PROPOSALSGLOBAL = markovchain(Step,Continuing_Step,cont_Proposal);
    end
    % save state variables every 1000 iterations in case inversion is interrupted
    save(['PROPOSALSGLOBALend' num2str(jj) '.mat'],'PROPOSALSGLOBAL')
    jj = jj+1;
    RJMCMC_Model = PROPOSALSGLOBAL;
end