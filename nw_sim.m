%% RSFC Non-white correction
% work from my phd thesis
% effect of non-white signals to the statistical tests for correlation
% results presented in OHBM 2021

%% define parameters

% define signal parameters
Fs = 50; % sampling frequency
T = 2000; % T: length of the signal
N=1000; % N: number of iterations

% define the values for corr and var
mu = [0,0];
var_def = 1;
corr_def = 0;
cor = corr_def;
sigma = var_def;
cov_mat = [sigma cor; cor sigma];

% define the filter
order = 3;   %order of filter
fcutlow = 0.01;   %low cut frequency in Hz
fcuthigh = 0.1;   %high cut frequency in Hz
[bf,af] = butter(order,[fcutlow/(Fs/2) fcuthigh/(Fs/2)],'bandpass');

% zero matrices
corr_0 = zeros(1,N);
corr_b = zeros(1,N);
corr_f = zeros(1,N);
s1 = zeros(T,N);
s2 = zeros(T,N);
s3 = zeros(T,N);
s1_f = zeros(T,N);
s2_f = zeros(T,N);
b1 = zeros(T,N);
b2 = zeros(T,N);
b3 = zeros(T,N);

% coloured spectrum parameters
c_def = 1;
a_def = 0.8;

%% generate simulated signals

% simulate non-white power spectrum
k = 1:T/2;
psh = zeros(T/2,1);
psh(k) = (c_def./k.^(a_def));

% generate simulated signals
for i = 1:1:N
    p = zeros(T/2,1);
    p1 = zeros(T/2,1);
    
    p(k) = normrnd(0,sqrt(psh(k)));
    p1(k) = normrnd(0,sqrt(psh(k)));
    
    p_all = [0;p;flipud(p(1:end-1))];
    p = (p_all).^2; % psd is positive
    
    p1_all = [0;p1;flipud(p1(1:end-1))];
    p1 = (p1_all).^2; % make sure psd is positive
    
    % bivariate normal
    t = T;
    r = mvnrnd(mu,cov_mat,t);
    sig1 = r(:,1);
    sig2 = r(:,2);
    sig1 = sig1-mean(sig1);
    sig2 = sig2-mean(sig2);
    b1(:,i) = sig1;
    b2(:,i) = sig2;
    
    % fnirs signals with decay psd
    [t1,~] = signalpsd(p,Fs,T/Fs);
    [temp,~] = signalpsd(p1,Fs,T/Fs);
    t1 = t1-mean(t1);
    temp = temp-mean(temp);
    t2 = corr_def * t1 + sqrt(1-corr_def^2)*temp;
    t2 = t2-mean(t2);
    s1(:,i) = t1;
    s2(:,i) = t2;
    s3(:,i) = temp;
    
    % fnirs signals after bpf
    t1_f = filter(bf,af,t1);
    t2_f = filter(bf,af,t2);
    s1_f(:,i) = t1_f;
    s2_f(:,i) = t2_f;
    
    corr_0(i) = corr(sig1,sig2);
    corr_b(i) = corr(t1,t2);
    corr_f(i) = corr(t1_f,t2_f);
    
end

T_all = [round((1-corr_def^2)^2/var(corr_0)) round((1-corr_def^2)^2/var(corr_b)) round((1-corr_def^2)^2/var(corr_f))];

% fft
fb1 = fft(b1);
fb2 = fft(b2);
fs1 = fft(s1);
fs2 = fft(s2);
fs1_f = fft(s1_f);
fs2_f = fft(s2_f);

%% correction

fs12 = [fs1,fs2];
va_re_emp = var(real(fs12),[],2);
va_im_emp = var(imag(fs12),[],2);

va2_emp = va_re_emp + va_im_emp;
va2_emp = va2_emp./max(va2_emp);
va2 = va2_emp;

% bpf fk2=dk2*bk2
[h,w] = freqz(bf,af,T/2,Fs);
h(1)=0;
hhalf = h.*conj(h);
fva2_emp = va2_emp.*([hhalf;flipud(hhalf)]);

% DoF for z
K_b = (sum(va2_emp).^2)/sum(va2_emp.^2);
K_f = (sum(fva2_emp).^2)/sum(fva2_emp.^2);
% DoF for t-distribution with circular convolution
K_b_t = abs((sum(cconv(va2_emp,va2_emp,T)))^2/(2*sum(cconv(va2_emp,va2_emp,T).^2)));
K_f_t = abs((sum(cconv(fva2_emp,fva2_emp,T)))^2/(2*sum(cconv(fva2_emp,fva2_emp,T).^2)));


% correct z-transform and t-test
% null hypothesis H0: rho=0
% define critical p-value
p_ci_z = 0.05;
p_ci_t = 0.05;
thre_z = norminv(1-p_ci_z/2);
thre_t_b_0 = tinv(1-p_ci_t/2,T-1);
thre_t_b = tinv(1-p_ci_t/2,K_b_t);
thre_t_f = tinv(1-p_ci_t/2,K_f_t);

% z-transform
z_0 = atanh(corr_0);
z_b = atanh(corr_b);
z_f = atanh(corr_f);

% standard normal
z_0_n = atanh(corr_0)*sqrt(T-1);
z_b_n = atanh(corr_b)*sqrt(T-1);
z_bc_n = atanh(corr_b)*sqrt(K_b);
z_f_n = atanh(corr_f)*sqrt(T-1);
z_fc_n = atanh(corr_f)*sqrt(K_f);

% test against standard normal
pz_0 = (1-normcdf(abs(z_0_n),atanh(corr_def)*sqrt(T-1),1))*2;
pz = (1-normcdf(abs(z_b_n),atanh(corr_def)*sqrt(T-1),1))*2;
pz_f = (1-normcdf(abs(z_f_n),atanh(corr_def)*sqrt(T-1),1))*2;
pz_fc = (1-normcdf(abs(z_fc_n),atanh(corr_def)*sqrt(K_f),1))*2;

% t-score
t_0_n = corr_0./sqrt(1-corr_0.^2)*sqrt(T-1);
t_b_n = corr_b./sqrt(1-corr_b.^2)*sqrt(T-1);
t_bc_n = corr_b./sqrt(1-corr_b.^2)*sqrt(K_b);
t_f_n = corr_f./sqrt(1-corr_f.^2)*sqrt(T-1);
t_fc_n = corr_f./sqrt(1-corr_f.^2)*sqrt(K_f);

corr_b_0 = corr_b - corr_def;
corr_f_0 = corr_f - corr_def;

t_b_n_0 = corr_b_0./sqrt(1-corr_b_0.^2)*sqrt(K_b);
t_f_n_0 = corr_f_0./sqrt(1-corr_f_0.^2)*sqrt(K_b);
t_fc_n_0 = corr_f_0./sqrt(1-corr_f_0.^2)*sqrt(K_f);

% test against t-distribution with correct DoF
pt = (1-tcdf(abs(t_b_n_0),K_b_t))*2;
pt_f = (1-tcdf(abs(t_f_n_0),K_b_t))*2;
pt_fc = (1-tcdf(abs(t_fc_n_0),K_f_t))*2;

% count the number of correlation that is statistically significant
ss_z = [sum(pz_0(:)<=p_ci_z) sum(pz(:)<=p_ci_z) sum(pz_f(:)<=p_ci_z) sum(pz_fc(:)<=p_ci_z)];
ss_t = [sum(pt(:)<=p_ci_t) sum(pt_f(:)<=p_ci_t) sum(pt_fc(:)<=p_ci_t)];


%% generate figures

% simulated coloured spectrum before and after filtering
figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2 0.2 0.6 0.6]);
fi = [0:Fs/(T-1):Fs];
ax1 = subplot(221);
plot(fi,va2_emp,'k')
title('d_k^2')
xlim([0 50])
xlabel('Hz')
ax3 = subplot(223);
plot(fi,fva2_emp,'k')
title('f_k^2 = d_k^2 b_k^2')
xlim([0 50])
xlabel('Hz')
linkaxes([ax1 ax3],'xy');
ax2 = subplot(222);
plot(fi(1:50),va2_emp(1:50),'k')
title('d_k^2 (Zoomed)')
xline(0.01,'-.r');
xline(0.1,'-.r');
xlabel('Hz')
ax4 = subplot(224);
plot(fi(1:50),fva2_emp(1:50),'k')
title('f_k^2 = d_k^2 b_k^2 (Zoomed)')
xline(0.01,'-.r');
xline(0.1,'-.r');
xlabel('Hz')
linkaxes([ax2 ax4],'xy');
suptitle('Illustrations of Non-white Power Spectrum')

% change in dof
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2 0.2 0.6 0.6]);
bar(T_all)
set(gca,'XTickLabel',{'White','Unfiltered Coloured','Filtered Coloured'});
text(1:3,T_all,num2str(T_all'),'vert','bottom','horiz','center');
title('Changes in Degrees of Freedom for White and Coloured Noise')

% statistical correction
figure(3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.95 1]);
bin = [-50:0.5:50];
ax1 = subplot(241);
histogram(z_0_n,bin,'Normalization','probability','FaceColor','none','EdgeColor','k','FaceAlpha',1);
xline(atanh(corr_def)*sqrt(T-1)-thre_z,'--r');
xline(atanh(corr_def)*sqrt(T-1)+thre_z,'--r');
xlim([-10 10])
xlabel('z-score for white gaussian noise')
ylabel('Probability Density')
ax2 = subplot(242);
histogram(z_b_n,bin,'Normalization','probability','FaceColor','none','EdgeColor','k','FaceAlpha',1);
xline(atanh(corr_def)*sqrt(T-1)-thre_z,'--r');
xline(atanh(corr_def)*sqrt(T-1)+thre_z,'--r');
xlim([-10 10])
xlabel('z-score for unfiltered coloured noise')
ylabel('Probability Density')
ax3 = subplot(243);
histogram(z_f_n,bin,'Normalization','probability','FaceColor','none','EdgeColor','k','FaceAlpha',1);
xline(atanh(corr_def)*sqrt(T-1)-thre_z,'--r');
xline(atanh(corr_def)*sqrt(T-1)+thre_z,'--r');
xlim([-40 40])
xlabel('z-score for filtered coloured noise')
ylabel('Probability Density')
ax4 = subplot(244);
histogram(z_fc_n,bin,'Normalization','probability','FaceColor','none','EdgeColor','k','FaceAlpha',1);
xline(atanh(corr_def)*sqrt(K_f)-thre_z,'--r');
xline(atanh(corr_def)*sqrt(K_f)+thre_z,'--r');
xlim([-10 10])
xlabel('z-score for filtered coloured noise (corrected)')
ylabel('Probability Density')

ax5 = subplot(245);
histogram(t_0_n,bin,'Normalization','probability','FaceColor','none','EdgeColor','k','FaceAlpha',1);
xline(corr_def./sqrt(1-corr_def.^2)*sqrt(T-1)-thre_t_b_0,'--r');
xline(corr_def./sqrt(1-corr_def.^2)*sqrt(T-1)+thre_t_b_0,'--r');
xlim([-10 10])
xlabel('t-score for white gaussian noise')
ylabel('Probability Density')
ax6 = subplot(246);
histogram(t_b_n,bin,'Normalization','probability','FaceColor','none','EdgeColor','k','FaceAlpha',1);
xline(corr_def./sqrt(1-corr_def.^2)*sqrt(T-1)-thre_t_b_0,'--r');
xline(corr_def./sqrt(1-corr_def.^2)*sqrt(T-1)+thre_t_b_0,'--r');
xlim([-10 10])
xlabel('t-score for unfiltered coloured noise')
ylabel('Probability Density')
ax7 = subplot(247);
histogram(t_f_n,bin,'Normalization','probability','FaceColor','none','EdgeColor','k','FaceAlpha',1);
xline(corr_def./sqrt(1-corr_def.^2)*sqrt(T-1)-thre_t_b_0,'--r');
xline(corr_def./sqrt(1-corr_def.^2)*sqrt(T-1)+thre_t_b_0,'--r');
xlim([-40 40])
xlabel('t-score for filtered coloured noise')
ylabel('Probability Density')
ax8 = subplot(248);
histogram(t_fc_n,bin,'Normalization','probability','FaceColor','none','EdgeColor','k','FaceAlpha',1);
xline(corr_def./sqrt(1-corr_def.^2)*sqrt(K_f)-thre_t_f,'--r');
xline(corr_def./sqrt(1-corr_def.^2)*sqrt(K_f)+thre_t_f,'--r');
xlim([-10 10])
xlabel('t-score for filtered coloured noise (corrected)')
ylabel('Probability Density')
linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7 ax8],'y');
suptitle('Statistical Correction for Simulated Signals')
%
