function gg_200921_jointfig_viz_1
clear;
close all;

rng(100);

%%%%%%%%%%%%%% define meta parameters
nr = 3;
nc = 2;

nr_glob = 3;
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


load('gg_200921_simresults_3.mat');

gamma=gammapar;
% kpar = [kinit,splic,gamma];


n_par_sets = 3;

for k__ = 1:3
    f = figure(k__);
    f.Position = [381.8000 121.8000 573.6000 588];
end

f=figure(n_par_sets+1);
f.Position = [339.4000 122.6000 809.6000 612];

% T_ = [1,1,0.2];
offs=5;

IND_ = [51,6,4];

titles_joint = {'K(t=\infty) distribution','Nascent distribution','Mature distribution','Joint distribution'};
systems_joint = {'Intrinsic','Extrinsic','Constitutive'};

for k__ = 1:3
    figure(k__)
    subplot(nr,nc,1);hold on;
    title(sprintf('%s time-series',systems_joint{k__}),'FontWeight','normal');
    subplot(nr,nc,2);hold on;
    title(sprintf('%s distribution',systems_joint{k__}),'FontWeight','normal');
end

figure(n_par_sets+1);
for i = 1:4
    subplot(nr_glob,nc_glob,i); hold on;
    title(titles_joint{i},'fontweight','normal');
end
for k__ = 1:3
    subplot(nr_glob,nc_glob,(k__-1)*nc_glob+1); hold on;
    ylabel(systems_joint{k__});
end

for k__ = 1:3
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

    subplot(nr,nc,2);
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

    subplot(nr,nc,3);
    y = mean(X(:,:,1),1); hold on;
    plot(tvec,y,'-','LineWidth',WIDTH_MEAN,'Color',LINECOL_MEAN);
    plot(tvec,X(ind,:,1),'-','LineWidth',WIDTH_ONE,'Color',LINECOL_ONE);
    plot([0, Tmax],[1 1]*alpha/eta/splic,'r--','LineWidth',WIDTH_THEORY);
    xlim([0,Tmax]);
    xlabel('time');
    ylabel('# nascent');

    subplot(nr,nc,4);
    H = histogram(X(:,nT,1),'binmethod','integer','normalization','probability',...
        'FaceColor',HISTCOL,'EdgeColor','none'); hold on;
    EDGES_X = H.BinEdges;
    x_nas = EDGES_X(1:end-1) + diff(EDGES_X)/2;
    y = compute_nas(k__,x_nas,lambda,eta,splic,gammapar,kappa);
    plot(x_nas,y,'-','LineWidth',WIDTH_THEORY,'Color',LINECOL_SIMPLE);
    set(gca,'yscale','log')
    xlabel('# nascent');
    ylabel('probability');

    subplot(nr,nc,5);
    y = mean(X(:,:,2),1); hold on;
    plot(tvec,y,'-','LineWidth',WIDTH_MEAN,'Color',LINECOL_MEAN);
    plot(tvec,X(ind,:,2),'-','LineWidth',WIDTH_ONE,'Color',LINECOL_ONE);
    plot([0, Tmax],[1 1]*alpha/eta/gammapar,'r--','LineWidth',WIDTH_THEORY);
    xlim([0,Tmax]);
    xlabel('time');
    ylabel('# mature');

    subplot(nr,nc,6);
    H = histogram(X(:,nT,2),'binmethod','integer','normalization','probability',...
        'FaceColor',HISTCOL,'EdgeColor','none'); hold on;
    % x = H.; 
    EDGES_X = H.BinEdges;
    x_mat = EDGES_X(1:end-1) + diff(EDGES_X)/2;
    y = compute_mat(k__,x_mat,lambda,eta,splic,gammapar,kappa);
    plot(x_mat,y,'-','LineWidth',WIDTH_THEORY,'Color',LINECOL_SIMPLE);
    set(gca,'yscale','log')
    xlabel('# mature');
    ylabel('probability');

    disp('');
    
    %%%%%%%%%%%%%%
    figure(n_par_sets+1);
    subplot(nr_glob,nc_glob,(k__-1)*nc_glob+1);
    H=histogram(X(:,end,3),300,'normalization','pdf',...
        'FaceColor',HISTCOL,'EdgeColor','none'); hold on;
    EDGES_GAMMA = H.BinEdges;
    x_gamma = EDGES_GAMMA(1:end-1) + diff(EDGES_GAMMA)/2;
    y = gampdf(x_gamma,lambda/kappa,1/eta);% y = y./sum(y);
    plot(x_gamma,y,'-','LineWidth',WIDTH_THEORY,'Color','r');
    xlim([x_gamma(1), x_gamma(end)]);
    xlabel('initiation rate');
    set(gca,'yscale','log')
    
    subplot(nr_glob,nc_glob,(k__-1)*nc_glob+2);
    H = histogram(X(:,nT,1),'binmethod','integer','normalization','probability',...
        'FaceColor',HISTCOL,'EdgeColor','none'); hold on;
    EDGES_X = H.BinEdges;
    x_nas = EDGES_X(1:end-1) + diff(EDGES_X)/2;
    x_nas = 0:(x_nas(end)+offs);
    y = compute_nas(k__,x_nas,lambda,eta,splic,gamma,kappa);
    plot(x_nas,y,'-','LineWidth',WIDTH_THEORY,'Color',LINECOL_SIMPLE);
    xlim([0, x_nas(end)]);
    xlabel('# nascent');
    set(gca,'yscale','log')
%     ylabel('probability');
    
    subplot(nr_glob,nc_glob,(k__-1)*nc_glob+3);
    H = histogram(X(:,nT,2),'binmethod','integer','normalization','probability',...
        'FaceColor',HISTCOL,'EdgeColor','none'); hold on;
    EDGES_X = H.BinEdges;
    x_mat = EDGES_X(1:end-1) + diff(EDGES_X)/2;    
    x_mat = 0:(x_mat(end)+offs);
    y = compute_mat(k__,x_mat,lambda,eta,splic,gamma,kappa);
    plot(x_mat,y,'-','LineWidth',WIDTH_THEORY,'Color',LINECOL_SIMPLE);
    xlim([0,x_mat(end)]);
    xlabel('# mature');
    set(gca,'yscale','log')
%     ylabel('probability');
    
    subplot(nr_glob,nc_glob,(k__-1)*nc_glob+4); hold on;
%     noise = zeros(nCells,2);
    noise = 1+randn(nCells,2)/20;
    z = compute_joint(k__,x_nas,x_mat,lambda,eta,splic,gamma,kappa);
    s=pcolor(x_nas+.5,x_mat+.5,(z')); s.EdgeColor='none';hold on;
    colormap summer
    axis([x_nas(1),x_nas(end),x_mat(1),x_mat(end)]+1);
    nas_ = (X(:,nT,1)+1).*noise(:,1);
    filt = X(:,nT,1)==0 & nas_<1;
    nas_(filt) = 2-nas_(filt);
    mat_ = (X(:,nT,2)+1).*noise(:,2);
    filt = X(:,nT,2)==0 & mat_<1;
    mat_(filt) = 2-mat_(filt);
    
    scatter(nas_,mat_,0.3,'k','filled',...
        'MarkerFaceAlpha',0.7);
    xlabel('# nascent + 1');
    ylabel('# mature + 1');
    set(gca,'yscale','log','xscale','log');
    %%%%%%%%%%%%%%%%%%%%%%%
end
for J = 1:4
    exportgraphics(figure(J),sprintf('fig_200921/fig%i.png',J),'Resolution',600);
end
return

function y = compute_nas(k__,x,lambda,eta,splic,gamma,kappa)

switch k__
    case 1 %intrinsic
        kinit = lambda;
        bs = 1/(kappa*eta);
        y = gg_200921_numint_burst_1(kinit,bs,gamma,splic,x(end)+1,1,inf,'nascent');
    case 2
        alpha = lambda/kappa;
        y = nbinpdf(x,alpha, splic*eta/(splic*eta+1));
    case 3
        mu = lambda/kappa/eta;
        y = poisspdf(x,mu/splic);
end
return

function y = compute_mat(k__,x,lambda,eta,splic,gamma,kappa)

switch k__
    case 1 %intrinsic
        kinit = lambda;
        bs = 1/(kappa*eta);
        y = gg_200921_numint_burst_1(kinit,bs,gamma,splic,1,x(end)+1,inf,'mature');
    case 2
        alpha = lambda/kappa;
        y = nbinpdf(x,alpha,gamma*eta/(gamma*eta+1));
    case 3
        mu = lambda/kappa/eta;
        y = poisspdf(x,mu/gamma);
end
return


function y = compute_joint(k__,x1,x2,lambda,eta,splic,gamma,kappa)

switch k__
    case 1 %intrinsic
        kinit = lambda;
        bs = 1/(kappa*eta);
        y = log(gg_200921_numint_burst_1(kinit,bs,gamma,splic,x1(end)+1,x2(end)+1,inf,'none'));
    case 2
        a = lambda/kappa;
        C = 1+(1/splic+1/gamma)/eta;
        x1 = x1';
        y = gammaln(a+x1+x2)-gammaln(a)-gammaln(x1+1)-gammaln(x2+1)...
            +a.*log(C)-x1.*log(C.*eta.*splic)-x2.*log(C.*eta.*gamma);
    case 3
        mu = lambda/kappa/eta;
        y = log(poisspdf(x1,mu/splic)'.*poisspdf(x2,mu/gamma));
end
return

