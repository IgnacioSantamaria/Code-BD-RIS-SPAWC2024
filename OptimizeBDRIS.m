function [Cfinal, Ctotal, Theta, Qt] = OptimizeBDRIS(H,F,G,Qt,Rxx,sigma2n,varargin)


% Description: This function optimizes a beyond-diagonal RIS (BDRIS) to maximize capacity
% in a MIMO link for a fixed transmit covariance. The BDRIS matrix is
% unitary+symmetric.
%
% It implements the algorithm in I. Santamaria et al, "MIMO capacity
% maximization with beyond-diagonal RIS" SPAWC 2024.
% The unitary +symmetric BD-RIS is factorized as Theta = Q.Q.', where Q is
% unitary
% 
%
% Input parameters:
% H,F,G: (direct, RIS->Tx, Tx->RIS, resp.)
% Qt : initial unitary matrix RIS is Qt*Qt.'
% Rxx : Tx covariance matrix
% sigma2 : noise variance
% varargin: structure with the algoritm parameters
%
% I. Santamaria, UC, Nov. 2023

[Nrx,~] = size(F);   % Matrix F is Nrx \times M

%% Default values
opt_params = struct();
opt_params.niterMM = 20;        % Maximum number of iterations for MM (inner loop)
opt_params.thresholdMM = 1e-3;  % To check convergence of the inner loop
opt_params.niterMO = 5000;      % maximum number of iterations for the manifold optimization (MO) algorithm
opt_params.thresholdMO = 1e-5;  % convergence threshold for the manifold optimization algorithm
opt_params.muMO = 1e-1;          % initial learning rate for the manifold optimization algorithm

if nargin < 7
    error(message('TooFewInputs'));
elseif nargin == 7
    params = varargin{1};
    for arg = fieldnames(params)'
        parameter = arg{1};
        param_value = params.(parameter);
        switch parameter
            case 'niterMM'
                opt_params.niterMM  = param_value;
            case 'thresholdMM'
                opt_params.thresholdMM  = param_value;
            case 'niterMO'
                opt_params.niterMO  = param_value;
            case 'thresholdMO'
                opt_params.thresholdMO  = param_value;
            case 'muMO'
                opt_params.muMO  = param_value;
        end
    end
elseif nargin > 7
    error(message('TooManyInputs'));
end

Theta = Qt*Qt.';
%% Decompose Rxx
[V,D,~] = svd(Rxx);
p = diag(D);
indp = find(p~=0);         % eliminate zeros. Indices with zero allocated power
p = p(indp);
V = V(:,indp);
Hq = H*V*diag(sqrt(p));    % transformed direct channel (with Rxx absorved)
Gq = diag(sqrt(p))*V'*G;   % transformed G (with Rxx absorved)
Ht = (Hq + F*Theta*Gq');   % initial equivalent channel
Sit = sigma2n*eye(Nrx);    % Noise Covariance matrix (fixed)
Ct = real(log(det(eye(Nrx) + Sit\(Ht*Ht'))));  % capacity

%% Algorithm begins
niterMM = opt_params.niterMM;
thresholdMM = opt_params.thresholdMM;
niterMO = opt_params.niterMO;
muMO = opt_params.muMO;
thresholdMO = opt_params.thresholdMO;

%% For a fixed transmit covariance matrix, we optimize the RIS elements
% minorizer-maximization (MM) inner loop
iterMM = 1;
trueMM = 1;
CapMM = zeros(1,niterMM);
CapMM(iterMM) = Ct;
while trueMM == 1  %MM (inner) loop
    iterMM = iterMM + 1;
    Rt = inv(Sit) - inv(Ht*Ht'+ Sit);
    [Uaux, Paux] = svd(Rt);
    Saux = sqrtm(Paux);
    sqrtRt = Uaux * Saux * Uaux';
    Ft = sqrtRt*F;
    Zt = Ht'/Sit - Hq'*Rt;
    A = Gq'*Zt*F;

    %% Solving a quadratic problem with a unitary constraint (this is the new part)
    Cost = zeros(1,niterMO);
    iterMO = 1;
    Cost(iterMO) = 2*real(trace(Zt*F*Theta*Gq')) - norm(Ft*Theta*Gq','fro')^2;  % Initial cost function
    trueMO = 1;
    while trueMO
        
        Grad = -((Ft'*Ft)*(Qt*Qt.')*(Gq'*Gq) - A')*conj(Qt);  % Unconstrained gradient (Changed Feb. 24)
        SkewHermitian = (Qt'*Grad-Grad'*Qt)/2;  % Projection onto the tangent plane
        Utnew = Qt*expm(muMO*SkewHermitian);    % Move along the geodesic
        Thetatnew = Utnew*Utnew.';
        Costnew  =  2*real(trace(Zt*F*Thetatnew*Gq')) - norm(Ft*Thetatnew*Gq','fro')^2;
        if Costnew>Cost(iterMO)
            iterMO = iterMO + 1;  %increase iteration
            Theta = Thetatnew;
            Qt = Utnew;
            Cost(iterMO) = Costnew;
            muMO = 1.1*muMO;     % increase step size
        else
            muMO = 0.9*muMO;     % decrease step size
        end

        %% Check convergence
        if (iterMO>2)&&((abs(Cost(iterMO)-Cost(iterMO-1))<thresholdMO)||(iterMO==niterMO))
            trueMO = 0;
        end
        
    end
    Ht = (Hq + F*Theta*Gq');         % New equivalent channel
    Ct = real(log(det(eye(Nrx) + Sit\(Ht*Ht'))));

    CapMM(iterMM) = Ct;
    DeltaCap = CapMM(iterMM)-CapMM(iterMM-1);
    %% Check convergence inner loop
    if (DeltaCap < thresholdMM) || (iterMM==niterMM)
        trueMM = 0;
    end

end
Ctotal = CapMM(1:iterMM);
Cfinal = Ctotal(end);  % This is the final solution of the inner loop
