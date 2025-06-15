n = 1000;
x = linspace(-1, 1, n);
f = exp(-1./sqrt(1 - x.^2));
K = 1 / (trapz(x, f)); % Numerical integration
fX = K * f;
g = 1 - x.^2;
L = 3/4;
gX = L * g;
M = (4 * K / 3) * exp(-1);
MgX = M * gX;

figure(1);
plot(x, fX, 'b-', 'LineWidth', 2, 'DisplayName', 'f_X(x)');
hold on;
plot(x, gX, 'r--', 'LineWidth', 2, 'DisplayName', 'g_X(x)');
plot(x, MgX, 'k-', 'LineWidth', 2, 'DisplayName', 'M \cdot g_X(x)');
legend('show');
xlabel('x');
ylabel('Density');
title('Rejection Sampling: f_X(x) \leq M \cdot g_X(x)');
hold off;



n = 10000;
X = randexpinvsqrt(n,1);
h = 0.1; kerneltype = 'cos';
u = (0.5+(0:n-1))/n;
x = 2*u-1;
fXhat = kernelestimation(x,X,h,kerneltype);
figure(2)
f = exp(-1./(1-x.^2));
dx = initmom(x,0);
K = 1/(f*dx);
fX = K*f;
plot(x,fXhat,'b','linewidth',3);
hold on
plot(x,fX,'r','linewidth',3);
xlabel('x');
ylabel('Density');
title('Kernel Density Estimation vs. True Density');
hold off
axy = axis; axy(2) = 1.1; axy(4) = axy(4)*1.1; axis(axy);
fticks = bestticks(0,axy(4)*1.1,10);
axis off
drawaxes
drawgrids(xticks,fticks)

varX = (fX.*x.^2)*dx
varXhat = mean(X.^2)

% The variance estimator may not be that accurate, but the procedure allows us
% to generate pseudo-random numbers.




varXhat = mean(X.^2);
fprintf('Estimated variance: %f\n', varXhat);

% For comparison, compute true variance numerically
x = linspace(-1, 1, 1000);
f = exp(-1./sqrt(1 - x.^2));
K = 1 / (trapz(x, f));
fX = K * f;
dx = x(2) - x(1);
varX = sum(x.^2 .* fX) * dx;
fprintf('True variance: %f\n', varX);
