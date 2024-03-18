function [Cfinal,Ctotal,Theta,Rxx] = MaxCap_RIS_passive(H,F,G,Theta,Pt,sigma2,UPA,varargin)

% Description: Finds the passive diagonal RIS, Theta, and the transmit covariance matrix, Rxx, that
% maximizes capacity in a MIMO link.
% It implements the algorithm in S. Zhang, R. Zhang, "Capacity
% characterization for IRS aided MIMO Communications", JSAC, vol. 38, no.
% 8, pp. 1823-1838, 2020.
% At each iteration the algorithm updates Rxx for fixed Theta, and then updates Theta
% (one element at a time) for fixed Rxx
%
% Input data:
% H,F,G: channels (direct, RIS->Tx, Tx->RIS, resp.)
% Theta: initial RIS matrix
% Pt, sigma2: Tx power, noise variance
% UPA: if UPA == 1 (uniform power allocation), otherwise optimal power allocation via waterfilling
% varargin: structure with optimization parameters
%
% Output data:
% Cfinal: Capacity final value (this is Ctotal(end))
% Ctotal: Convergence curve for capacity
% Theta: RIS
% Rxx: Transmit covariance matrix
%
% Ignacio Santamaria, UC, March. 2024


%% Initialization
[Nrx,~] = size(F);   % Matrix F is Nrx \times M
[Ntx,~] = size(G);   % Matrix G is Ntx \times M
%% Default values
opt_params = struct();
opt_params.niterRIS = 100;
opt_params.thresholdRIS = 1e-3;

if nargin < 7
    error(message('TooFewInputs'));
elseif nargin == 8
    params = varargin{1};
    for arg = fieldnames(params)'
        parameter = arg{1};
        param_value = params.(parameter);
        switch parameter
            case 'niterRIS'
                opt_params.niterRIS = param_value;
            case 'thresholdRIS'
                opt_params.thresholdRIS = param_value;
            case 'niterAO'
                opt_params.niterAO = param_value;
            case 'thresholdAO'
                opt_params.thresholdAO = param_value;
        end
    end
elseif nargin > 8
    error(message('TooManyInputs'));
end

%% Max-Cap optimization
disp('Start passive RIS optimization')
if UPA == 1       % Uniform power allocation
    Rxx = (Pt/Ntx)*eye(Ntx);
    [~, Ctotal, Theta] = OptimizeRIS(H,F,G,Theta,Rxx,sigma2,opt_params);
    iter = length(Ctotal);
else              % Optimal Power Allocation
    true = 1;
    niterAO = opt_params.niterAO;
    thresholdAO = opt_params.thresholdAO;
    iter = 1;
    Rxx = (Pt/Ntx)*eye(Ntx);  % Initial Rxx
    Ht = H + F*Theta*G';
    Ctotal(iter) = real(log(det(eye(Nrx) + Ht*Rxx*Ht'/sigma2)));
    while true == 1             % Alternating optimization loop
        iter = iter+1;
        %% Optimize Rxx: For a fixed Theta we optimize the Tx covariance matrix
        Ht = H + F*Theta*G';
        Si = sigma2*eye(Nrx);   % Noise covariance matrix
        [Rxx,~,~] = OptTransmitCovMatrix(Ht,Si,Pt);
        %%  Optimize Theta: For a fixed Rxx, we optimize Theta
        [Cfinalinner, Ctotalinner, Theta] = OptimizeRIS(H,F,G,Theta,Rxx,sigma2,opt_params); %#ok<*ASGLU>
        Ctotal(iter) = Cfinalinner;      % This is the final solution of the inner loop
        DeltaCap =  Ctotal(iter)- Ctotal(iter-1);
        %% Check convergence
        if (DeltaCap  < thresholdAO) || (iter==niterAO)
            true = 0;
        end
    end
end
disp('Passive RIS optimization finished')
Ctotal = Ctotal(1:iter)*log2(exp(1));      % Convergence of the outer loop, we've used natural logarithms
Cfinal = Ctotal(end);
