
figure; 
subplot(2,2,1); 
tmax = 400;
tnew = linspace(0,tmax,40000);

ss =1; %prs.ss; 
sr = 1; %prs.sr;
k = ss;
% beta = 1/sr;
% 
% theta = 1/beta;

invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
plot(tnew, invgam(tnew,k,sr));
 
subplot(2,2,2); 
tmax = 1;
tnew = linspace(0,tmax,1000);

ss =2; %prs.ss; 
sr = 0.005; %prs.sr;
k = ss;
% beta = 1/sr;
% 
% theta = 1/beta;

invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
plot(tnew, invgam(tnew,k,sr));


subplot(2,2,3); 
tmax = 4;
tnew = linspace(0,tmax,1000);

ss =100; %prs.ss; 
sr = 1; %prs.sr;
k = ss;
% beta = 1/sr;
% 
% theta = 1/beta;

invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
plot(tnew, invgam(tnew,k,sr));



subplot(2,2,4); 
tmax = 1;
tnew = linspace(0,tmax,1000);

ss = 70; %prs.ss; 
sr = 0.1; %prs.sr;
k = ss; 
% beta = 1/sr;
% 
% theta = 1/beta;

invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
plot(tnew, invgam(tnew,k,sr));
