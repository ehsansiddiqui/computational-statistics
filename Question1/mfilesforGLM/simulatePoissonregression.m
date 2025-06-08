
%if exist('samplesize') ~= 1, samplesize = 1000; end
%if exist('fullmodelsize') ~= 1, fullmodelsize = 15; end
if exist('realrandom') ~= 1, realrandom = false; end
if exist('studnumber') ~= 1, studnumber = 12345; end
if exist('betamax') ~= 1, betamax = 10; end

if exist('maxit') ~= 1, maxit = 100; end

if exist('figsize') ~= 1, figsize = [1280 420]; end
if exist('papersize') ~= 1, papersize = figsize/96; end
if exist('figpaperpos') ~= 1, figpaperpos = [0 0 papersize]; end

forestgreen = [34 139 34]/255; fg = forestgreen;
linedotopt = {'-o','linewidth',3,'markersize',8,'color'};
linedotoptg = {linedotopt{:},fg,'markerfacecolor','w'};
linedotoptr = {linedotopt{:},'r','markerfacecolor','w'};
linedotoptb = {linedotopt{:},'b','markerfacecolor','w'};

bulletcolor = 'b';
bullet = {'o','markersize',6,'linewidth',2,'color',bulletcolor,...
          'markerfacecolor',bulletcolor};

seed = studnumber;
if realrandom==false, randn('state',seed); rand('state',seed); end

n = 10000; m = 50;
x = sort(rand(n,1));
mu = (cos(5*x.^2)+2)*5; Y = poissonnoise(mu); % DGP
X = ones(n,m); for k=1:m-1, X(1:n,k+1) = x.^k; end
% alternative: DGP within specified model
% beta = round(rand(m,1)*2*betamax-betamax);
% theta = X*beta; mu = exp(theta); Y = poissonnoise(mu); % DGP
theta = log(mu);

% Model selection
submodels = cell(1,m); submodels{1} = 1;
for k=2:m, submodels{k} = (1:k); end
pp = (1:m);

% CpGLM = icGLM(X,Y,'Poisson','Cp',submodels,maxit)*stdev^2/n;
AICGLM = icGLM(X,Y,'Poisson','AIC',submodels,maxit);
BICGLM = icGLM(X,Y,'Poisson','BIC',submodels,maxit);

plot(pp,AICGLM,linedotoptb{:})


% added values

figure(1)
plot(pp,AICGLM,linedotoptb{:})
title('AIC'); 
xlabel('modelsize');              
ylabel('ic'); 

figure(2)
plot(pp,AICGLM,linedotoptb{:})
title('AIC'); 
xlabel('modelsize');              
ylabel('ic'); 
xlim([1 15]);
ylim([-7.0 -5.0]);

figure(3)
plot(pp,BICGLM,linedotoptb{:})
title('BIC'); 
xlabel('modelsize');              
ylabel('ic');  

figure(4)
plot(pp,BICGLM,linedotoptb{:})
title('BIC'); 
xlabel('modelsize');              
ylabel('ic');  
xlim([1 12]);
ylim([-7.0 -5.0]);
% end values


figure(5)
plot(pp,AICGLM,linedotoptr{:}, 'DisplayName', 'AIC');
hold on;
plot(pp,BICGLM,linedotoptb{:}, 'DisplayName', 'BIC');
title('n = 100 : AIC VS BIC'); 
xlabel('modelsize');              
ylabel('ic'); 
legend('show');
figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos);





% [minCp pCp] = min(CpGLM);
[maxAIC pAIC] = max(AICGLM);
[maxBIC pBIC] = max(BICGLM);

% pCp
pAIC
pBIC

% IRLS with BIC
S = submodels{pBIC};
Xp = X(1:n,S);


[betahat muhat thetahat niter] = IRLS(Xp,Y,'Poisson',maxit);
% thetahat = X*betahat; muhat = exp(thetahat);
% plot(x,mu,'r','linewidth',3)
% hold on
% plot(x,muhat,'m','linewidth',3)
% plot(x,Y,bullet{:})
% hold off
% figpos = get(gcf,'position'); figpos(3:4) = figsize;
% set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos)



% added values
figure(6)
plot(x,mu,'r','linewidth',3)
hold on
plot(x,muhat,'m','linewidth',3)
plot(x,Y,bullet{:})
hold off 
xlabel('x');              
ylabel('mu');  
legend('true model', 'IRLS estimation', 'data')
figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos)

%IRLS with AIC
S = submodels{pAIC};
Xp = X(1:n,S);
[betahat, muhat, thetahat, niter] = IRLS(Xp,Y,'Poisson',maxit);
% thetahat = X*betahat; phat = 1./(1+exp(-thetahat));
figure(7)
plot(x,mu,'r','linewidth',3)
hold on
plot(x,muhat,'m','linewidth',3)
plot(x,Y,bullet{:})
hold off
xlabel('x');              
ylabel('mu');  
legend('true model', 'IRLS estimation', 'data')
figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos)
% end values