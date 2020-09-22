function gg_200921_jointfig_gen_1
clear;
close all;

rng(100);

%define 4 burst model kinetic parameters
gamma = 0.7; %degradation rate 
gammapar = gamma;
splic = 1.2; %splice rate
kinit = 0.1; %burst frequency
bs = 15; %burst size;


sde_params = zeros(3,3);

%%%%%%%%%%%% intrinsic
alpha = 0.01;
lambda = kinit;
kappa = lambda/alpha;
eta = 1/(bs*kappa);
sde_params(1,:) = [kappa, lambda, eta];
%%%%%%%%%%%% extrinsic
alpha = kinit/splic;
kappa = splic/10;
lambda = alpha*kappa;
eta = 1/bs;
sde_params(2,:) = [kappa, lambda, eta];
%%%%%%%%%%%% constitutive
MU = 5;
eta = 20;
alpha = MU*splic*eta;
lambda = kinit;
kappa = lambda/alpha;
sde_params(3,:) = [kappa, lambda, eta];


kpar = [kinit,splic,gamma];

%%%%%%%%%%%%%%
% set sim parameters
nCells = 10000;
nT = 300;
S = [1 0; -1 1; 0 -1; 0 0];

n_par_sets = 3;
% for k__ = 1:n_par_sets
Y = NaN(n_par_sets, nCells, nT,3);

T_ = [1,1,0.2];

for k__ = 1:3
    SDE = sde_params(k__,:);
    kappa = SDE(1);
    lambda = SDE(2);
    alpha = lambda/kappa;
    eta = SDE(3);
    
    Tmax=T_(k__)./min([kappa, gamma, lambda, eta]);
    tvec = linspace(0,Tmax,nT);
    t_matrix = repmat(tvec,nCells,1);

    [X,~,~]=...
        gg_200921_gillespie_gou_2step_8(kpar,t_matrix,S,nCells,SDE);
    Y(k__,:,:,:) = X(:,:,1:3);
    %%%%%%%%%%%%%%%%%%%%%%%
end
save('gg_200921_simresults_3.mat');
return
