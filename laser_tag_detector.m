%% Analysis for the initial detector
%% create the <H> matrix first
freqs = [1471, 1724, 2000, 2273, 2632, 2941, 3333, 3571, 3846, 4167];
% declare <H>
H = zeros(2000, 20);
for i = 1:1:10
for j = 0:1:1999 %vector length 2000 starting at t=0
%cosine vector
H(j+1,i*2-1) = cos(2*pi*freqs(i)*(j*.0001));
%sine vector
H(j+1,i*2) = sin(2*pi*freqs(i)*(j*.0001));
end
end
%% beta vs SNR
%declare variable
snr = linspace(0,10,10000);
alpha = [.1 .01 .001 .0001];
%compute thresholds
p = 20; %h is an Nxp matrix
N = 2000; % number of samples we care about.
figure(1)
for iter = 1:1:4
fo = icdf('F',1-alpha(iter),p,N-p);
beta = 1-cdf('Noncentral F',fo,p,N-p,snr.^2);
semilogy(snr,beta);
hold on; grid on;
end
hold off;
legend('\alpha = .1','\alpha = .01','\alpha = .001','\alpha = .0001','Location','SouthEast');
ylabel('\beta');
xlabel('SNR');
title('Probability of detection for specific values of \alpha');
%% beta vs alpha
alpha = linspace(0,1,10000);
snr = [.8 1.6 2.4 3.2 4.8 6.4]; %length 6
fo = icdf('F',1-alpha,p,N-p);
figure(2)
for iter = 1:1:6
beta = 1-cdf('Noncentral F',fo,p,N-p,snr(iter)^2);
plot(alpha,beta);
hold on; grid on;
end
hold off;
legend('SNR = .8','SNR = 1.6','SNR = 2.4','SNR = 3.2','SNR = 4.8','SNR = 6.4','Location','SouthEast');
ylabel('\beta');
xlabel('\alpha');
title('\beta vs \alpha');
%% s_m'*s_m
%generate random phis on a scale from zero to 2*pi
phi = rand(10,1)*2*pi;
freqs = [1471, 1724, 2000, 2273, 2632, 2941, 3333, 3571, 3846, 4167];
t = linspace(.0001,.2,2000);
results = zeros(10);
for iter = 1:1:10
for k = 1:1:10
sm = cos(2*pi*freqs(iter)*t+phi(k))';
results(iter,k) = sm'*sm;
end
end
%% Analysis for the Second Detector
% While it has been said that the use of the chisquare distribution is
% possible, I will use the F distribution because overall compuation
% ability is not my concern, and it is easier to follow the guide given in
% the handouts, since technically we should be using the F distribuiton for
% this project
% dist = zeros(10000,1);
% for iter = 1:1:10000
%% set up <s> = <H><theta>
freqs = [ 1471,1724, 2000, 2273, 2632, 2941, 3333, 3571, 3846, 5167];
% declare <H>
H = zeros(2000, 20);
for i = 1:1:10
for j = 0:1:1999 %vector length 2000 starting at t=0
%cosine vector
H(j+1,i*2-1) = cos(2*pi*freqs(i)*(j*.0001));
%sine vector
H(j+1,i*2) = sin(2*pi*freqs(i)*(j*.0001));
end
end
%% find threshold for given alpha
alpha = .001; % 0.1 percent probability of error
p = 20; %h is an Nxp matrix
N = 2000; % number of samples we care about.
%determine threshold fo
fo = icdf('F',1-alpha,p,N-p);
%% create a random vector <x>
x = zeros(2000,1);
%create signal
a = 1; %zero for no signal 1 for signal
frq = 3; %frequency number
for j = 0:1:1999
x(j+1) = a*cos(2*pi*freqs(frq)*(j*.0001));
end
%add noise
sigma2= 1;
x = x+sigma2*randn(2000,1); %you MUST have noise in the received signal for this detector to work
%% determine if received signal is a hit or miss
Ph = H*(H'*H)^-1*H';
numerator = x'*Ph*x/(p);
I = eye(N);
denominator = x'*(I-Ph)*x/(N-p);
f = numerator/denominator;% see handout page 12, eq 59
%now print the comparison
if f>fo
output = 'Hit';
else
output = 'Miss';
end
if a == 0 && strcmp(output,'Hit')
result = 'incorrect detect';
elseif a == 1 && strcmp(output,'Miss')
result = 'incorrect miss';
else
fprintf('correct decision: a = %d\n',a);
end
%where noise is not super gaussian like since we only have 2000samples)
% u2 = 0.01 this is mu square; sigma2 = 1 this is simga square. this way
% to make signal to noise ratio smaller like range from 1:10
M = cell(10,1);% modes
Pm= cell(10,1);%projection of each modes
E = zeros(10,1); %Engergy of all 10
SNR = zeros(10,1);%SNR of all 10
xx = 0:0.01:100;
figure;
for ii = 1:10
M{ii} = H(:,2*ii-1:2*ii);
Pm{ii} = M{ii}*(M{ii}'*M{ii})^-1*M{ii}';
E(ii)=x'*Pm{ii}*x;
SNR(ii) = 0.01*E(ii)/sigma2;
plot(xx,pdf('Noncentral F',xx,2,20,SNR(ii)));
hold on;
end
xlim([0 200]);
ylim([0 1]);
legend('f1','f2','f3','f4','f5','f6','f7','f8','f9','f10');
title('when frequency'+ string(frq)+' is in force, SNR ='+ string(SNR(frq)));
xlabel('x');
ylabel('F_{2,1980}(x)');
hold off;
%fd = cdf('Noncentral F',113.3,2,20,SNR(2)); % probability of false detection(2.4003e-9) super small
% fot the correct detection probability we can say it's almost 1
%% for different u2/sigma2, we make u2/sigma2 = const
% M = cell(10,1);% modes
% Pm= cell(10,1);%projection of each modes
% E = zeros(10,1); %Engergy of all 10
noise = randn(2000,1);
figure;
for kk = 1:10
threshold = zeros(10,1);
fd = zeros(10,1);
SNR_k = zeros(10,1);
SNR_m= zeros(10,1);
for jj = 1:10
x = zeros(2000,1);
for j = 0:1:1999
x(j+1) = cos(2*pi*freqs(kk)*(j*.0001));
end
%add noise
sigma2= 1;
u2 =jj*0.003*sigma2;
x = x+sigma2*noise; %you MUST have noise in the received signal for this detector to work
%figure;
for ii = 1:10
M{ii} = H(:,2*ii-1:2*ii);
Pm{ii} = M{ii}*(M{ii}'*M{ii})^-1*M{ii}';
E(ii)=x'*Pm{ii}*x;
SNR_m(ii) = u2*E(ii)/sigma2;
%f = pdf('Noncentral F',xx,2,20,SNR_m(ii));
% plot(xx,f);
end
for ii = 1:10
F_func = @(t) pdf('Noncentral F',t,2,20,SNR_m(kk))- pdf('Noncentral F',t,2,20,SNR_m(ii));
threshold(ii) = fzero( F_func, 4);
end
SNR_k(jj) = SNR_m(kk) ;
threshold(kk)= 0;
fd(jj) = cdf('Noncentral F',max(threshold),2,20,SNR_m(kk));
end
plot(SNR_k,(1-fd)*100);
hold on;
end
xlabel('SNR');
ylabel('Probability of Correct Detection (%)');
title('The correct detection probability given SNR')
legend('f1','f2','f3','f4','f5','f6','f7','f8','f9','f10');