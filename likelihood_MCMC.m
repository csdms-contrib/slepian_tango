function [LogLikelihood, Prefactor, Uncertainty, Misfit, Residual,RMS_Error] = likelihood_MCMC(RJ_MCMC_THBI_state)
%**************************************************************************
% p(d_obs|m)=1/sqrt(2pi)^N|Cd| exp{1/2(d_obs-d_mod )^t.Cd^-1.(d_obs-d_mod)}
% d_mod = g(m)
% Cd => data covariance matrix / uncertainty with its determinant as |Cd|
%**************************************************************************
    % Inputs:
            % rj_MCMC_THBI_state.Observed   :: observed data with length nlm
            % rj_MCMC_THBI_state.Predicted  :: predicted model data with length nlm
            % rj_MCMC_THBI_state.noiseType  :: uncertainty estimation applied


    % Output:
            % P_Log_Likelihood    :: Probability of likelihood
            % mistfit             :: Estimated mistfit from observed and Proposed data
if RJ_MCMC_THBI_state.Iteration == 1 && RJ_MCMC_THBI_state.NoiseType == 1
%**************************************************************************
%   A.    Estimate likelihood of initial current model with noiseType 1
%**************************************************************************
        Observed 	= RJ_MCMC_THBI_state.Observed_Geoid;
        Predicted	= RJ_MCMC_THBI_state.Proposed_Geoid;
        SigmaInit   = RJ_MCMC_THBI_state.SigmaInit;
        Nlm         = length(Observed);
        Residual 	= Observed - Predicted;
        Diff 		= norm(Residual);
        covariance_matrix   = SigmaInit*eye(length(Residual));
        Uncertainty         = chol(covariance_matrix,'lower');
%--------------------------------------------------------------------------
% Guassian likelihood %Neglect the constant term (1/(sqrt(2*pi)*sigma))^N as 
% it will pull down the likelihood value to zero for increasing value of N
%--------------------------------------------------------------------------
	G                       = Uncertainty\Residual;
	Mahalanobis_distance    = G'*G;
    Misfit                  = Diff;
	P_Likelihood            = -0.5*Mahalanobis_distance;
	Prefactor               = NaN;
    RMS_Error               = sqrt(sum((Residual(:).^2))/Nlm);

elseif RJ_MCMC_THBI_state.Iteration == 1 && RJ_MCMC_THBI_state.NoiseType == 2
%**************************************************************************
%   B.    Applying correlated data noise using normal distribution
%**************************************************************************

        Observed        = RJ_MCMC_THBI_state.Observed_Geoid;
        Predicted       = RJ_MCMC_THBI_state.Kernel_Geoid;
        SigmaInit       = RJ_MCMC_THBI_state.SigmaInit;
        Nlm             = length(Observed);
        
        Cov_Mat1      	= (eye(Nlm))*SigmaCurrent^2;
        Cov_Mat2      	= zeros(size(Cov_Mat1));
        Propo_CovMat  	= (Cov_Mat1 + Cov_Mat2);
        Propo_Cov_Det 	= det(Propo_CovMat);
        Uncertainty   	= Propo_CovMat;
        %Generate the Mahalanobis distance Cd
%       ====================================== 
        Residual        = Observed - Predicted;
%         MahalanobisDistance = -(0.5*(Residual'*(inv(Propo_CovMat))*Residual));
        Mahalanobis_distance = -(0.5*(Residual'/(Propo_CovMat)*Residual));
        logPhi = -log((sqrt(2*pi)^Nlm)*Propo_Cov_Det);
        P_Likelihood = logPhi + Mahalanobis_distance;

        %Estimate L2 norm misfit
%       =========================
        Misfit = 0;
        for ii = 1:length(Observed)
%         Misfit = Misfit + ((Observed(ii)-Predicted(ii))/diag(Propo_CovMat(ii)))^2;
        Misfit = Misfit + ((Observed(ii)-Predicted(ii))/SigmaCurrent)^2;
        end
        
        Misfit = Mahalanobis_distance;
        RMS_Error = sqrt(sum((Residual(:).^2))/Nlm);
        Residual = norm(Residual);
elseif RJ_MCMC_THBI_state.Iteration > 1 && RJ_MCMC_THBI_state.NoiseType == 1
%**************************************************************************
%   C.    Estimate likelihood of initial current model with noiseType 1
%**************************************************************************
        Observed        = RJ_MCMC_THBI_state.Observed_Geoid;
        Proposed        = RJ_MCMC_THBI_state.Proposed_Geoid;
        SigmaInit       = RJ_MCMC_THBI_state.SigmaInit;        
        SigmaProposed 	= RJ_MCMC_THBI_state.SigmaProposed;
        SigmaCurrent    = RJ_MCMC_THBI_state.SigmaCurrent;
        Nlm             = length(Observed);
        
        Residual            = Observed - Proposed;
        Diff                = norm(Residual);
        covariance_matrix   = SigmaInit*eye(length(Residual));
        Uncertainty         = (chol(covariance_matrix,'lower'))*sqrt(SigmaProposed);
%--------------------------------------------------------------------------
%      Guassian likelihood %Neglect the constant term 1/sqrt(2*pi)^N as 
%       it will pull %down the likelihood value to zero for increasing value of N
%--------------------------------------------------------------------------
        G                       = Uncertainty\Residual;
        Mahalanobis_distance    = G'*G;
        Misfit                  = Diff; 
        P_Likelihood            = -0.5*Mahalanobis_distance;
        Prefactor               = -0.5*Nlm*log(SigmaProposed/SigmaCurrent);
        RMS_Error               = sqrt(sum((Residual(:).^2))/Nlm);
elseif RJ_MCMC_THBI_state.Iteration > 1 && RJ_MCMC_THBI_state.NoiseType > 1
%**************************************************************************
%  D.     Applying correlated data noise using normal distribution
%**************************************************************************
        Observed = RJ_MCMC_THBI_state.Observed_Geoid;
        Proposed= RJ_MCMC_THBI_state.Proposed_Geoid;
        SigmaProp       = RJ_MCMC_THBI_state.SigmaProposed;
        SigmaCurrent    = RJ_MCMC_THBI_state.SigmaCurrent;
        Nlm             = length(Observed);
        
        Cov_Mat1        = (eye(Nlm))*SigmaProp^2;
        Cov_Mat2        = zeros(size(Cov_Mat1));
        Prop_CovMat     = (Cov_Mat1 + Cov_Mat2);
        Prop_CovDet     = det(Prop_CovMat);
        Uncertainty     = Prop_CovMat;
        %Generate the Mahalanobis distance Cd
%       ======================================        
        Residual             = Observed - Proposed;
        Mahalanobis_distance = (-0.5*(Residual'/(Prop_CovMat)*Residual));
        LogPhi = -log(sqrt(((2*pi)^Nlm)*Prop_CovDet)); 
        P_Likelihood = LogPhi + Mahalanobis_distance;
        
%         Prefactor = -0.5*Nlm*log(SigmaProp/SigmaCurrent)
%         L1 = chol(Prop_CovMat,'lower');
%         Z  = L1\Residual;
%         Mahalanobis_distance = -0.5*(Z'*Z);
        %Estimate L2 norm misfit and RMS Error
%       =========================         
        Misfit = Mahalanobis_distance;
        RMS_Error = sqrt(sum((Residual(:).^2))/Nlm);
        Residual = norm(Residual);
end
    LogLikelihood = P_Likelihood;
end
