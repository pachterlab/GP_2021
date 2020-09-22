function gg_200921_samplefig_viz_1
clear;
close all;

rng(100);

%%%%%%%%%%%%%% define meta parameters
nr = 2;
nc = 5;

nr_glob = 1;
nc_glob = 4;
% LINECOL = [0.85,0.33,0.10];
LINECOL_SIMPLE = [1 0 0];
LINECOL_SNB = [0 0 1];
HISTCOL = 0.6*[1 1 1];
WIDTH_MEAN = 2;
WIDTH_ONE = 0.5;
WIDTH_THEORY = 1.4;
LINECOL_MEAN = [0 0 0];
LINECOL_ONE = 0.5*[1 1 1];


load('gg_200921_simresults_samp_1.mat');

gamma=gammapar;
% kpar = [kinit,splic,gamma];


n_par_sets = 1;

f = figure(1);
f.Position =  [201.8000 253.8000 1132 394.4000];

offs=5;

IND_ = [51,6,4];


k__ = 1;
SDE = sde_params(k__,:);
kappa = SDE(1);
lambda = SDE(2);
alpha = lambda/kappa;
eta = SDE(3);
X = squeeze(Y(k__,:,:,:));

ind = IND_(k__);

figure(k__);

subplot(nr,nc,1);
plot(tvec,mean(X(:,:,3),1),'-','LineWidth',WIDTH_MEAN,'Color',LINECOL_MEAN);hold on;
plot(tvec,X(ind,:,3),'-','LineWidth',WIDTH_ONE,'Color',LINECOL_ONE);
plot([0,Tmax],[1 1]*alpha/eta,'r--','LineWidth',WIDTH_THEORY);
xlim([0,Tmax]);
xlabel('time');
ylabel('initiation rate');

subplot(nr,nc,nc+1);
H=histogram(X(:,end,3),300,'normalization','pdf',...
    'FaceColor',HISTCOL,'EdgeColor','none'); hold on;
EDGES_GAMMA = H.BinEdges;
x_gamma = EDGES_GAMMA(1:end-1) + diff(EDGES_GAMMA)/2;
x = x_gamma;
y = gampdf(x,lambda/kappa,1/eta);% y = y./sum(y);
plot(x,y,'-','LineWidth',WIDTH_THEORY,'Color','r');
xlabel('initiation rate');
ylabel('probability');
set(gca,'yscale','log')

subplot(nr,nc,2);
y = mean(X(:,:,1),1); hold on;
plot(tvec,y,'-','LineWidth',WIDTH_MEAN,'Color',LINECOL_MEAN);
plot(tvec,X(ind,:,1),'-','LineWidth',WIDTH_ONE,'Color',LINECOL_ONE);
plot([0, Tmax],[1 1]*alpha/eta/splic,'r--','LineWidth',WIDTH_THEORY);
xlim([0,Tmax]);
xlabel('time');
ylabel('# nascent');

subplot(nr,nc,nc+2);
H = histogram(X(:,nT,1),'binmethod','integer','normalization','probability',...
    'FaceColor',HISTCOL,'EdgeColor','none'); hold on;
EDGES_X = H.BinEdges;
x_nas = EDGES_X(1:end-1) + diff(EDGES_X)/2;
fit = nbinfit(X(:,nT,1));
y = nbinpdf(x_nas,fit(1),fit(2));
plot(x_nas,y,'-','LineWidth',WIDTH_THEORY,'Color','b');
set(gca,'yscale','log')
xlabel('# nascent');
ylabel('probability');

subplot(nr,nc,3);
y = mean(X(:,:,2),1); hold on;
plot(tvec,y,'-','LineWidth',WIDTH_MEAN,'Color',LINECOL_MEAN);
plot(tvec,X(ind,:,2),'-','LineWidth',WIDTH_ONE,'Color',LINECOL_ONE);
plot([0, Tmax],[1 1]*alpha/eta/gammapar,'r--','LineWidth',WIDTH_THEORY);
xlim([0,Tmax]);
xlabel('time');
ylabel('# mature');

subplot(nr,nc,nc+3);
H = histogram(X(:,nT,2),'binmethod','integer','normalization','probability',...
    'FaceColor',HISTCOL,'EdgeColor','none'); hold on;
EDGES_X = H.BinEdges;
x_mat = EDGES_X(1:end-1) + diff(EDGES_X)/2;
fit = nbinfit(X(:,nT,2));
y = nbinpdf(x_mat,fit(1),fit(2));
plot(x_mat,y,'-','LineWidth',WIDTH_THEORY,'Color','b');
set(gca,'yscale','log')
xlabel('# mature');
ylabel('probability');

disp('');

subplot(nr,nc,[4,5,9,10]); hold on;
noise = 1+randn(nCells,2)/30;
axis([x_nas(1),x_nas(end),x_mat(1),x_mat(end)]+1);
nas_ = (X(:,nT,1)+1).*noise(:,1);
filt = X(:,nT,1)==0 & nas_<1;
nas_(filt) = 2-nas_(filt);
mat_ = (X(:,nT,2)+1).*noise(:,2);
filt = X(:,nT,2)==0 & mat_<1;
mat_(filt) = 2-mat_(filt);

scatter(nas_,mat_,1,'k','filled',...
    'MarkerFaceAlpha',0.7);
xlabel('# nascent + 1');
ylabel('# mature + 1');
set(gca,'yscale','log','xscale','log');
    %%%%%%%%%%%%%%%%%%%%%%%
pos = get(gca,'position');
pos(1) = pos(1)+0.03;
set(gca,'position',pos);


exportgraphics(figure(1),'fig_200921/samp_fig_1.png','Resolution',600);
return