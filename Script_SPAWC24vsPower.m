
% This program compares several passive RIS architectures that maximize capacity in 
% a RIS-assisted NtxNr MIMO link. We optimize the Tx covariance matrix and
% the BDRIS (or RIS) matrix.
% The different solutions correspond to different constraints on the RIS matrix:
% a) Conventional RIS (diagonal matrix, passive) and  b) BD-RIS passive
% This Script illustrates the achievable rate vs Tx power of the algorithm proposed in
% I. Santamaria, M- Soleymani, E. Jorswieck, J. Gutierrez, "MIMO capacity
% maximization with beyond-diagonal RIS", SPAWC 2024.
%
% I. Santamaria, UC Nov. 2023
%
% Notes: You need cvx in the path. The inner loop might use a "parfor", change to
% "for" if you experiencie problems with the parallel computing toolbox.

format compact
clc; clear;
%addpath(genpath('D:\Usuarios\Ignacio\Matlab\cvx'));    % add all cvx subfolders to the path

%% Parameters
Ntx = 4;                 % Number of transmit antennas
Nrx = 4;                 % Number of receive antennas
M = 100;                 % Number of RIS elements  (the size of the problem for BD-RIS architectures grows as M^2, so keep M<20)

PtdBm = 4:2:30;
Ptotal = 1e-3*10.^(PtdBm/10);   % Total Power for TX + RIS (if active)
UPA = 0;                        % If 1 uniform power allocation Rx = (P/Ntx)*I, otherwise Rx is optimized
NsimMC = 1;     % Number of Monte Carlo simulations
show = 1;       % To plot intermediate convergence/rate results (only without the parfor)

%% Optimization parameters BDRIS
opt_paramsBDRIS = struct();
opt_paramsBDRIS.niterAO = 50;          % Maximum number of iterations for Alternating Optimization (outer loop)
opt_paramsBDRIS.niterMM = 100;          % Maximum number of iterations for MM (inner loop)
opt_paramsBDRIS.thresholdAO = 1e-3;    % To check convergence of the outer loop
opt_paramsBDRIS.thresholdMM = 1e-3;    % To check convergence of the inner loop
opt_paramsBDRIS.muMO = 1e-1;            % initial learning rate for the manifold optimization algorithm
opt_paramsBDRIS.niterMO = 1000;        % maximum number of iterations for the manifold optimization (MO) algorithm
opt_paramsBDRIS.thresholdMO = 1e-5;    % convergence threshold for the manifold optimization algorithm

%% Optimization parameters RIS
opt_paramsRIS = struct();
opt_paramsRIS.niterinner = 100;
opt_paramsRIS.thresholdinner = 1e-3;
opt_paramsRIS.niterAO = 50;
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

posicion_RIS = 50; % We vary the x coordinate of RIS position

B = 20;   % Bandwidth MHz
NF = 0;   % Noise Factor in dBs
noiseVariancedBm = -174 + 10*log10(B*10^6) + NF;
sigma2n = 10^(noiseVariancedBm/10);       % additive noise variance

noiseVariancedBmRIS = -160 + 10*log10(B*10^6) + NF;
sigma2v = 10^(noiseVariancedBmRIS/10);    % variance of transmit noise (for active RISs)

%% Variables to store the rates (b/s/Hz)
CRIS_tot = zeros(size(Ptotal));               % Diagonal passive RIS (optimized)
CBDRIS_tot = zeros(size(Ptotal));             % Passive BD-RIS (unitary+symmetric, init: Mao)

%% Variables to store the number of streams
StreamsRIS_tot = zeros(size(Ptotal));          % Diagonal passive RIS (optimized)
StreamsBDRIS_tot = zeros(size(Ptotal));        % Passive BD-RIS (unitary+symmetric, init: Mao)

for mm = 1:NsimMC  

    disp(['Simulation:', num2str(mm)])
    for dd = 1:length(Ptotal)                  % Use for or parfor here
        disp(['Power:', num2str(Ptotal(dd))])

        %% Generate channels
        PosRIS_XYZ = [posicion_RIS, 5, 5];     % (x,y,z) coordinates for the RIS position
        [H,G,F] = ChannelsMIMO(M,Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,ray_fading, RiceFactor, ...
            pl_0,alpha_RIS,alpha_direct,blocked);
        F = F.';  % F should be an Nrx \times M matrix
        G = G.';  % G should be an Ntx \times M matrix

        %% Solution for a diagonal (passive) RIS
        ThetaRISini = diag(exp(1i*rand(1,M)*2*pi));    % Initial diagonal RIS
        [CRIS, CtotalRIS,ThetaRIS,QRIS] = MaxCap_RIS_passive(H,F,G,ThetaRISini,Ptotal(dd),sigma2n,UPA,opt_paramsRIS);
        [~,Paux,Vaux] = svd(QRIS);
        HeqRIS = (H + F*ThetaRIS*G')*Vaux*sqrt(Paux);  %Equivalent MIMO Channel with Theta and Q absorbed in the channel

        %% Solution for a passive BD-RIS unitary+symmetric (init:Mao)
        % Initialization
        ThetaBDRISini = F'*H*G;
        ThetaMaosym = (ThetaBDRISini + ThetaBDRISini.')/2;
        [Q,Sigma] = TakagiSVD(ThetaMaosym);
        Qtini = Q;   %Initial Q (Theta = Q*Q.')
        
        [CBDRIS, CtotalBDRIS, ThetaBDRIS, QBDRIS] =  ...
            MaxCap_BDRIS_passive(H,F,G,Qtini,Ptotal(dd),sigma2n,UPA,opt_paramsBDRIS);
        [~,Paux,Vaux] = svd(QBDRIS);
        HeqBDRIS = (H + F*ThetaBDRIS*G')*Vaux*sqrt(Paux);

        %% Store Rates
        CRIS_tot(dd) = CRIS_tot(dd) + real(log2(det(eye(Nrx) + HeqRIS*HeqRIS'./sigma2n)));
        CBDRIS_tot(dd) = CBDRIS_tot(dd) +  real(log2(det(eye(Nrx) + (HeqBDRIS*HeqBDRIS')./sigma2n)));

        StreamsRIS_tot(dd) = StreamsRIS_tot(dd) +rank(QRIS);               % Diagonal passive RIS (optimized)
        StreamsBDRIS_tot(dd) = StreamsBDRIS_tot(dd) + rank(QBDRIS);        % Passive BD-RIS (unitary+symmetric, init: Mao)
        
    end
    figure(20);clf; plot(10*log10(Ptotal*1000),CBDRIS_tot, 'r-o','MarkerSize',ms,'LineWidth',lw);
    hold on;
    plot(10*log10(Ptotal*1000),CRIS_tot,'b-d','MarkerSize',ms,'LineWidth',lw);
    legend('BDRIS','RIS');
    xlabel('Power (dBm)');
    ylabel('Rate (b/s/Hz)')
    hold off

    figure(30);clf; plot(10*log10(Ptotal*1000),StreamsBDRIS_tot, 'r-o','MarkerSize',ms,'LineWidth',lw);
    hold on;
    plot(10*log10(Ptotal*1000),StreamsRIS_tot,'b-d','MarkerSize',ms,'LineWidth',lw);
    legend('BDRIS','RIS');
    xlabel('Power (dBm)');
    ylabel('Number of streams')
    hold off

end

%% Plot results
CRIS_tot = CRIS_tot/NsimMC;
CBDRIS_tot = CBDRIS_tot/NsimMC;

StreamsRIS_tot = StreamsRIS_tot/NsimMC;
StreamsBDRIS_tot = StreamsBDRIS_tot/NsimMC;

figure(20);clf; plot(10*log10(Ptotal*1000),CBDRIS_tot, 'r-o','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(10*log10(Ptotal*1000),CRIS_tot,'b-d','MarkerSize',ms,'LineWidth',lw);
legend('BDRIS','RIS');
xlabel('Power (dBm)');
ylabel('Rate (b/s/Hz)')
hold off

figure(30);clf; plot(10*log10(Ptotal*1000),StreamsBDRIS_tot, 'r-o','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(10*log10(Ptotal*1000),StreamsRIS_tot,'b-d','MarkerSize',ms,'LineWidth',lw);
legend('BDRIS','RIS');
xlabel('Power (dBm)');
ylabel('Number of streams')
hold off