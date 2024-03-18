
% This Script illustrates the achievable rate vs RIS position of the algorithm proposed in
% I. Santamaria, M. Soleymani, E. Jorswieck, J. Gutierrez, "MIMO capacity
% maximization with beyond-diagonal RIS", SPAWC 2024.This program maximize capacity in a RIS-assisted 
% Ntx \times Nrx MIMO link.
% For the capacity maximization algorithm we optimize alternatively the
% transmit covariance matrix and the RIS matrix.
%
% The different Passive RIS architectures correspond to different constraints on the RIS matrix:
% a) Passive diagonal RIS (with optimized and random phases)
% b) Passive BD-RIS 
% d) no RIS
%
% I. Santamaria, UC Jan. 2024
%
% Notes: You need cvx in the path. The inner loop uses a "parfor", change to "for" if you experiencie problems with the parallel computing toolbox.

format compact
clc; clear;
%addpath(genpath('D:\Usuarios\Ignacio\Matlab\cvx'));    % add all cvx subfolders to the path

%% Parameters: antennas, RIS elements, powers
Ntx = 4;                 % Number of transmit antennas
Nrx = 4;                 % Number of receive antennas
M = 100;                 % Number of RIS elements  (the size of the problem for BD-RIS architectures grows as M^2, so keep M<20)
Ptotal = 0.1;            % Total Power for Tx 
UPA = 0;                 % If 1 uniform power allocation Rx = (P/Ntx)*I, otherwise Rx is optimized

NsimMC = 1;              % Number of Monte Carlo simulations
show = 1;                % To plot intermediate convergence/rate results (only without the parfor)

%% Optimization parameters BDRIS
opt_paramsBDRIS = struct();
opt_paramsBDRIS.niterAO = 10;          % Maximum number of iterations for Alternating Optimization (outer loop)
opt_paramsBDRIS.niterMM = 100;         % Maximum number of iterations for MM (inner loop)
opt_paramsBDRIS.thresholdAO = 1e-3;    % To check convergence of the outer loop
opt_paramsBDRIS.thresholdMM = 1e-3;    % To check convergence of the inner loop
opt_paramsBDRIS.muMO = 1e-1;           % initial learning rate for the manifold optimization algorithm
opt_paramsBDRIS.niterMO = 1000;        % maximum number of iterations for the manifold optimization (MO) algorithm
opt_paramsBDRIS.thresholdMO = 1e-5;    % convergence threshold for the manifold optimization algorithm

%% Optimization parameters RIS
opt_paramsRIS = struct();
opt_paramsRIS.niterinner = 100;
opt_paramsRIS.thresholdinner = 1e-3;
opt_paramsRIS.niterAO = 10;
opt_paramsRIS.thresholdAO = 1e-3;

%% Parameters for figures
fs = 12;   % fontsize
lw = 1.5;  % linewidth
ms = 8;    % markersize
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Scenario (parameters)
ray_fading = 0;  % Set it "1" if all channels follow Rayleigh fading
blocked = 0;     % Set to 1 if all direct channels are blocked
RiceFactor = 3;  % Rician factor for links through RIS
% ===== Parameters for Large-scale Path Loss =============
pl_0      = -33.05 + 5; % Path loss at a reference distance (d_0)
alpha_RIS = 2;          % Path loss exponent for the RIS links
alpha_direct  = 3.75;   % Path loss exponent for the direct links
% ===== Position of the Tx/Rx/RIS (units in meters) ======
sqr_size = 50;                  % square of size sqr_size
PosTx_XYZ = [0 0 1.5];          % Position Tx
PosRx_XYZ = [sqr_size 0 1.5];   % Position Rx

posicion_RIS = 10:10:100; % We vary the x coordinate of RIS position

B = 20;   % Bandwidth MHz
NF = 0;   % Noise Factor in dBs
noiseVariancedBm = -174 + 10*log10(B*10^6) + NF;
sigma2n = 10^(noiseVariancedBm/10);       % additive noise variance

noiseVariancedBmRIS = -160 + 10*log10(B*10^6) + NF;
sigma2v = 10^(noiseVariancedBmRIS/10);    % variance of transmit noise (for active RISs)

%% Variables to store the rates (b/s/Hz)
CRIS_tot = zeros(size(posicion_RIS));               % Diagonal passive RIS (optimized)
CRISrnd_tot = zeros(size(posicion_RIS));            % Diagonal passive RIS (random phases)
CBDRIS_tot = zeros(size(posicion_RIS));             % Passive BD-RIS (unitary+symmetric)
CnoRIS_tot = zeros(size(posicion_RIS));             % Without RIS

for mm = 1:NsimMC

    disp(['Simulation:', num2str(mm)])

    for dd = 1:length(posicion_RIS)  % You can also remove the parfor and use for instead

        disp(['Position RIS:', num2str(posicion_RIS(dd))])

        %% Generate channels
        PosRIS_XYZ = [posicion_RIS(dd), 5, 5];     % (x,y,z) coordinates for the RIS position
        [H,G,F] = ChannelsMIMO(M,Nrx,Ntx,PosTx_XYZ,PosRx_XYZ,PosRIS_XYZ,ray_fading, RiceFactor, ...
            pl_0,alpha_RIS,alpha_direct,blocked);
        F = F.';  % F should be an Nrx \times M matrix
        G = G.';  % G should be an Ntx \times M matrix

        %% Solution for a diagonal passive RIS
        ThetaRISini = diag(exp(1i*rand(1,M)*2*pi));    % Initial diagonal RIS
        [CRIS, CtotalRIS,ThetaRIS,QRIS] = MaxCap_RIS_passive(H,F,G,ThetaRISini,Ptotal,sigma2n,UPA,opt_paramsRIS);
        [~,Paux,Vaux] = svd(QRIS);
        HeqRIS = (H + F*ThetaRIS*G')*Vaux*sqrt(Paux);  %Equivalent MIMO Channel with Theta and Q absorbed in the channel

        %% Solution for a diagonal (passive) RIS with random coefficients
        ThetaRISrnd = diag(exp(1i*2*pi*rand(M,1)));
        Hrnd = H + F*ThetaRISrnd*G';
        if UPA == 1
            QRISrnd = (Ptotal/Ntx)*eye(Ntx);
            Vaux = eye(Ntx);
            paux = (Ptotal/Ntx)*ones(1,Ntx);  % Q = V*diag(p)*V'
        else
            [QRISrnd,Vaux,paux] = OptTransmitCovMatrix(Hrnd,sigma2n*eye(Nrx),Ptotal);  % optimal Q
        end
        HeqRISrnd = (H + F*ThetaRISrnd*G')*Vaux*diag(sqrt(paux));

        %% Solution for a passive BD-RIS unitary+symmetric (init:Mao)
        % Initalization
        ThetaMao = (F'*H*G); % sqrt(M(dd))*(F'*H'*G)/norm(F'*H'*G, 'fro');
        ThetaMaosym = (ThetaMao + ThetaMao.')/2;
        [F_SVD,~,G_SVD] = svd(ThetaMaosym);
        dummy = F_SVD'*conj(G_SVD);
        % Eigenvectors post-correction
        phases = angle(diag(dummy))/2;
        F_SVD = F_SVD*diag(exp(1i*phases));
        Qtini = F_SVD;   %Initial Qt (Theta = Qt*Qt.')
        %Qtini = sqrtm(ThetaRISini);          % or random RIS initialization
        
        [CBDRIS, CtotalBDRIS, ThetaBDRIS, QBDRIS] = ...
            MaxCap_BDRIS_passive(H,F,G,Qtini,Ptotal,sigma2n,UPA,opt_paramsBDRIS);
        [~,Paux,Vaux] = svd(QBDRIS);
        HeqBDRIS = (H + F*ThetaBDRIS*G')*Vaux*sqrt(Paux);
        
        %% Solution without RIS
        if UPA == 1
            QnoRIS = (Pt/Ntx)*eye(Ntx);
            Vaux = eye(Ntx);
            paux = (Pt/Ntx)*ones(1,Ntx);  % Q = V*diag(p)*V'
        else
            [QnoRIS,Vaux,paux] = OptTransmitCovMatrix(H,sigma2n*eye(Nrx),Ptotal);   % optimal Q
        end
        HeqnoRIS = H*Vaux*diag(sqrt(paux));

        %% Store Rates
        CRIS_tot(dd) = CRIS_tot(dd) + real(log2(det(eye(Nrx) + HeqRIS*HeqRIS'./sigma2n)));
        CRISrnd_tot(dd) = CRISrnd_tot(dd) + real(log2(det(eye(Nrx) + HeqRISrnd*HeqRISrnd'./sigma2n))); 
        CBDRIS_tot(dd) = CBDRIS_tot(dd) +  real(log2(det(eye(Nrx) + (HeqBDRIS*HeqBDRIS')./sigma2n)));
        CnoRIS_tot(dd) = CnoRIS_tot(dd) +  real(log2(det(eye(Nrx) + (HeqnoRIS* HeqnoRIS')./sigma2n)));

    end
    if show ==1
        figure(10);clf; plot(posicion_RIS,CBDRIS_tot, 'r-o','MarkerSize',ms,'LineWidth',lw);
        hold on;
        plot(posicion_RIS,CRIS_tot,'b-d','MarkerSize',ms,'LineWidth',lw);
        plot(posicion_RIS,CRISrnd_tot,'c-^','MarkerSize',ms,'LineWidth',lw);
        plot(posicion_RIS,CnoRIS_tot,'k--s','MarkerSize',ms,'LineWidth',lw);
        legend('BDRIS (Proposed)', ...
            'RIS (opt.)', 'RIS (rand. phases)',...
            ' no RIS');
        xlabel('RIS position (m)');
        ylabel('Rate (b/s/Hz)')
        hold off
    end
end

%% Results
CRIS_tot = CRIS_tot/NsimMC;
CRISrnd_tot = CRISrnd_tot/NsimMC;
CBDRIS_tot = CBDRIS_tot/NsimMC;
CnoRIS_tot = CnoRIS_tot/NsimMC;

figure(10);clf; plot(posicion_RIS,CBDRIS_tot, 'r-o','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(posicion_RIS,CRIS_tot,'b-d','MarkerSize',ms,'LineWidth',lw);
plot(posicion_RIS,CRISrnd_tot,'c-^','MarkerSize',ms,'LineWidth',lw);
plot(posicion_RIS,CnoRIS_tot,'k--s','MarkerSize',ms,'LineWidth',lw);
legend('BDRIS (Proposed)', ...
            'RIS (opt.)', 'RIS (rand. phases)',...
            ' no RIS');
xlabel('RIS position (m)');
ylabel('Rate (b/s/Hz)')
hold off