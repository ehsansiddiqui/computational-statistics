screensize = get(0,'screensize');
figwidth0 = screensize(3); figsize = [figwidth0 420];
papersize = figsize/96; figpaperpos = [0 0 papersize];
% if exist('samplesize') ~= 1, samplesize = 20; end
if exist('realrandom') ~= 1, realrandom = false; end
if realrandom==false, randn('state',0); rand('state',0); end
n = 20; % samplesize
ns = 1000; % number of samples

M = 3;
p = rand(1,M);
r = 2.3;
X = zeros(n,ns);
for m=1:M, X = X+randnegbin(size(X),p(m),r); end
d = sum((1./p).*(1./p-1))/sum(1./p-1);
mu = r*sum(1./p-1);

B = 10000;
tstaralfa = zeros(ns,2);
ustaralfa = zeros(ns,2);
alfa = 0.05;
qq = round([alfa/2 1-alfa/2]*B);
muhat = mean(X); T = muhat';
S = sqrt(muhat'*d);
for s=1:ns,
   Xs = X(1:n,s);
   r = ceil(rand(n,B)*n);
   Xstar = Xs(r);
   mustar = T(s);
   muhatstar = sort(mean(Xstar));
   tstaralfa(s,1:2) = muhatstar(qq);
   Sstar = sqrt(sum((1-p)./(p.^2))/n); % Standard error of bootstrap mean
   Ustar = (muhatstar - mustar) ./ Sstar; % Studentized statistic
   ustaralfa(s,1:2) = sort(Ustar)(qq); % Quantiles of studentized statistic
end

basicCI = 2*[T T]-tstaralfa(1:ns,[2 1]);
bootstraptCI = [T T] + S .* ustaralfa(1:ns,[1 2]); % Bootstrap-t CI

% Debugging: Check inputs to plotasblocks
disp('Size of bootstraptCI(:,1):'); disp(size(bootstraptCI(1:ns,1)));
disp('Size of bootstraptCI(:,2):'); disp(size(bootstraptCI(1:ns,2)));

figure(1)
try
    [i1 b1] = plotasblocks(basicCI(1:ns,1));
    [i2 b2] = plotasblocks(basicCI(1:ns,2));
    nn = length(i2);
    ii = [i1; i2(nn:-1:1)];
    bb = [b1; b2(nn:-1:1)];
    fill(ii,bb,'y')
catch e
    disp('Error in plotasblocks for Figure 1:');
    disp(e.message);
    % Fallback: Simple plot
    plot(1:ns, basicCI(1:ns,1), 'b-', 1:ns, basicCI(1:ns,2), 'b-');
    
end
xlabel('sample')
ylabel('mean');
hold on
plotasblocks(T,'k')
plot([1 ns],[mu mu],'r')
hold off
figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos)

figure(2)
try
    [i1 b1] = plotasblocks(bootstraptCI(1:ns,1));
    [i2 b2] = plotasblocks(bootstraptCI(1:ns,2));
    nn = length(i2);
    ii = [i1; i2(nn:-1:1)];
    bb = [b1; b2(nn:-1:1)];
    fill(ii,bb,[1 0.7 1])
catch e
    disp('Error in plotasblocks for Figure 2:');
    disp(e.message);
    % Fallback: Simple plot
    plot(1:ns, bootstraptCI(1:ns,1), 'b-', 1:ns, bootstraptCI(1:ns,2), 'b-');
    
end
xlabel('sample')
ylabel('mean');
hold on
plotasblocks(T,'k')
plot([1 ns],[mu mu],'r')
hold off
figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'papersize',papersize,'paperposition',figpaperpos)

coverbasicCI = sum((basicCI(1:ns,1)<mu)&(mu<basicCI(1:ns,2)))/ns
coverbootstraptCI = sum((bootstraptCI(1:ns,1)<mu)&(mu<bootstraptCI(1:ns,2)))/ns
widthbasicCI = basicCI(1:ns,2)-basicCI(1:ns,1);
widthbootstraptCI = bootstraptCI(1:ns,2)-bootstraptCI(1:ns,1);
meanwidthbasicCI = mean(widthbasicCI)
meanwidthbootstraptCI = mean(widthbootstraptCI)