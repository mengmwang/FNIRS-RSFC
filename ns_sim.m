%% RSFC Non-stationary correction
% work from my phd thesis
% effect of non-stationary signals to the statistical tests for correlation
% results presented in sfnirs 2021

%% white signal * non-stationary weights

% simulate correlated random noise
T=3000;
N=2000;
rho=0.2;

% simulate random non-stationary signals from a inverse-gamma distribution
nbin=200;
gpar = [2 1/2];  % gamma parameters
ig = 1./gamrnd(gpar(1),gpar(2),[T,1]);
fia = fitvariance(ig,'figure','off');
binnsa = [0:max(ig)/nbin:max(ig)];
inda = find(strcmp({fia.name}, 'inversegamma')==1);
co_iga = fia(inda).par;
pdfiga = (co_iga(2)^co_iga(1)/gamma(co_iga(1)).*binnsa.^(-co_iga(1)-1).*exp(-co_iga(2)./binnsa));

% simulate awgn first
xg1 = normrnd(0,1,[T,N]);
xg3 = normrnd(0,1,[T,N]);
si = 10*rand(2,1); % sigma_x and sigma_y
xg2 = rho.*xg1+sqrt(1-rho^2)*xg3;

xg1 = xg1*sqrt(si(1));
xg2 = xg2*sqrt(si(2));

ig2 = ig; % weight from the same ig distribution
xgn1 = xg1.*sqrt(ig);
xgn2 = xg2.*sqrt(ig2);

xg = [xg1 xg2];
xgn = [xgn1 xgn2];

corrg = diag(corr(xg1,xg2))';
corrgn = diag(corr(xgn1,xgn2))';

% correction
vg1 = var(xgn1,[],2);
vg2 = var(xgn2,[],2);
vg = var(xgn,[],2);

xg1c = xgn1./sqrt(vg1);
xg2c = xgn2./sqrt(vg2);

corrgc = diag(corr(xg1c,xg2c))';

% plot results
figure(1)
nbinhist = 100;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2 0 0.54 1]);
bin = [min([corrg corrgn corrgc]):(max([corrg corrgn corrgc])-min([corrg corrgn corrgc]))/nbinhist:max([corrg corrgn corrgc])];
ax1 = subplot(221);
histogram(corrg,bin,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
xlabel('Original Sample Correlation')
ylabel('Probability Density')
ax2 = subplot(222);
histogram(corrgn,bin,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
xlabel('Uncorrected Sample Correlation')
ylabel('Probability Density')
ax3 = subplot(223);
histogram(corrgc,bin,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
xlabel('Corrected Sample Correlation')
ylabel('Probability Density')
linkaxes([ax1 ax2 ax3],'xy');
subplot(224);
scatter(corrg,corrgn,'k+')
hold on
scatter(corrg,corrgc,'r.')
axis equal
xlim([min([corrg corrgn corrgc])-0.01 max([corrg corrgn corrgc])]+0.01)
ylim([min([corrg corrgn corrgc])-0.01 max([corrg corrgn corrgc])]+0.01)
legend('Non-stationary timeseries','Corrected timeseries','Location','northwest')
xlabel('Stationary Timeseries Correlation')
ylabel('Correlation')
suptitle('Correction for White Gaussian Signals with Non-Stationary weights')


%%  ar * non-stationary weights

% ar(p) p=1
phi1 = 0.8;
phi2 = 0.8;
rho = 0.2;
cov_ep = rho*(1-phi1*phi2)/sqrt((1-phi1^2)*(1-phi2^2));

mu = [0 0];
sigma = [1 cov_ep; cov_ep 1];

for i = 1:1:N
    eps = mvnrnd(mu,sigma,T);
    
    x = zeros(T,1);
    y = zeros(T,1);
    
    x(1) = eps(1,1);
    y(1) = eps(1,2);
    
    for t=2:T
        x(t) = phi1*x(t-1)+eps(t,1);
        y(t) = phi2*y(t-1)+eps(t,2);
    end
    
    xs1(:,i) = x;
    xs2(:,i) = y;
    
end

% weight
xs1 = xs1*sqrt(si(1));
xs2 = xs2*sqrt(si(2));

xn1 = xs1.*sqrt(ig);
xn2 = xs2.*sqrt(ig2);

xs = [xs1 xs2];
xn = [xn1 xn2];

corrs = diag(corr(xs1,xs2))';
corrn = diag(corr(xn1,xn2))';

% correction
v1 = var(xn1,[],2);
v2 = var(xn2,[],2);
v = var(xn,[],2);

x1c = xn1./sqrt(v1);
x2c = xn2./sqrt(v2);

corrc = diag(corr(x1c,x2c))';

figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2 0 0.54 1]);
bin = [min([corrs corrn corrc]):(max([corrs corrn corrc])-min([corrs corrn corrc]))/nbinhist:max([corrs corrn corrc])];
ax1 = subplot(221);
histogram(corrs,bin,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
xlabel('AR Sample Correlation')
ylabel('Probability Density')
ax2 = subplot(222);
histogram(corrn,bin,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
xlabel('Uncorrected Sample Correlation')
ylabel('Probability Density')
ax3 = subplot(223);
histogram(corrc,bin,'Normalization','probability','FaceColor','k','EdgeColor','k','FaceAlpha',1)
xlabel('Corrected Sample Correlation')
ylabel('Probability Density')
linkaxes([ax1 ax2 ax3],'xy');
subplot(224);
scatter(corrs,corrn,'k+')
hold on
scatter(corrs,corrc,'r.')
axis equal
xlim([min([corrs corrn corrc])-0.01 max([corrs corrn corrc])+0.01])
ylim([min([corrs corrn corrc])-0.01 max([corrs corrn corrc])+0.01])
legend('AR non-stationary timeseries','AR corrected timeseries','Location','northwest')
xlabel('AR Stationary Timeseries Correlation')
ylabel('Correlation')
suptitle('Correction for AR Signals with Non-Stationary weights')

vg = [var(corrg) var(corrgn) var(corrgc)];
v = [var(corrs) var(corrn) var(corrc)];
%
