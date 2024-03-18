function [Cfinal, Ctotal, Theta] = OptimizeRIS(H,F,G,Theta,Rxx,sigma2,varargin)

% Description: This function optimizes a diagonal RIS to maximizae capacity
% in a MIMO link for a fixed transmit covariance.
%
% 
%   Detailed explanation goes here (TBC)

[Nrx,M] = size(F);   % Matrix F is Nrx \times M
%% Default values
opt_params = struct();
opt_params.niterRIS = 100;
opt_params.thresholdRIS = 1e-3;

if nargin < 7
    error(message('TooFewInputs'));
elseif nargin == 7
    params = varargin{1};
    for arg = fieldnames(params)'
        parameter = arg{1};
        param_value = params.(parameter);
        switch parameter
            case 'niterRIS'
                opt_params.niterRIS = param_value;
            case 'thresholdRIS'
                opt_params.thresholdRIS = param_value;
        end
    end
elseif nargin > 7
    error(message('TooManyInputs'));
end

%% Decompose Rxx 
[V,D,~] = svd(Rxx);
p = diag(D);
indp = find(p~=0);  % eliminate zeros. Indices with zero allocated power
p = p(indp);
V = V(:,indp);
Hq = H*V*diag(sqrt(p));    % transformed direct channel (with Rxx absorved)
Gq = diag(sqrt(p))*V'*G;   % transformed G (with Rxx absorved)

%% Algorithm begins
niter = opt_params.niterRIS;
threshold = opt_params.thresholdRIS;
Ctotal = zeros(1,niter);                    % Store capacity values
iter = 1;
Ht = H + F*Theta*G';
Ctotal(iter) = real(log(det(eye(Nrx) + Ht*Rxx*Ht'/sigma2)));
true = 1;
while true == 1             % Optimization loop
    iter = iter+1;
    %%  Optimize Theta: For a fixed Rxx, we optimize Theta
    % the RIS elements are updated one at a time
    theta = diag(Theta);
    for mm = 1:M  % loop to update the mth RIS element
        mindex = 1:M;
        thetam = theta;
        mindex(mm) = [];
        thetam(mm) = [];
        Thetam = diag(thetam);
        Fm = F(:,mindex);  % select the fixed columns
        Gm = Gq(:,mindex);
        fm = F(:,mm);      % select the column to update
        gm = Gq(:,mm);
        S = Hq + Fm*Thetam*Gm';    % fixed matrix (it does not depend on m)
        rm = S*gm;
        A = eye(Nrx)+ (1/sigma2)*(S*S'+ fm*fm'/(gm'*gm));
        angleopt = angle(fm'*(A\rm));  % optimal phase
        theta(mm) = exp(1i*angleopt);
    end
    Theta = diag(theta);           % new RIS
    Heq = Hq + F*Theta*Gq';        % equivalent channel
    Ctotal(iter) = log(real(det(eye(Nrx) + Heq*Heq'/sigma2)));  % This is the final solution of the inner loop
    DeltaCap =  Ctotal(iter)- Ctotal(iter-1);
    %% Check convergence
    if (DeltaCap  < threshold) || (iter==niter)
        true = 0;
    end
end
Ctotal = Ctotal(1:iter);
Cfinal = Ctotal(end);

