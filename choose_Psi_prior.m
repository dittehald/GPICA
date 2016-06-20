
Eval =  mean(val); % (e(i,:)*e(i,:)' )
Psi_pri = 0.6; %200;
Psi_pos = 0.8; %325;
% N


s = (Psi_pos*(0.5*N+1) -0.5*Eval - Psi_pri ) / (Psi_pri - Psi_pos);

r = Psi_pri*(s+1);


disp(s)
disp(r)

SigmaMedian = median(Sigmaone,2);
%%
figure; subplot(4,1,1); 

if prs.std
tmax = 1;    
else
tmax = 500; %max(SigmaMedian);
end

histogram(SigmaMedian,100,'BinLimits',[0 tmax]);
title(['ss: ' num2str(prs.ss) ', sr: ' num2str(prs.sr) ', std: ' num2str(prs.std)])
axis tight 


% tnew = linspace(0,tmax,1000);
% xax = gca;
% ss =prs.ss; 
% sr = prs.sr;
% k = ss;
% 
% invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
% plot(tnew, invgam(tnew,k,sr));
% xlim(xax.XLim)
% axis tight 



e = X_0mean-WZ;

% histogram(abs(e),'Normalization','pdf','BinLimits',[0 tmax])
% axis tight 
% 
% %%
% figure
% subplot(2,1,1)
tnew = linspace(0,tmax,1000);
sstest = s;
srtest = r;

subplot 413
histogram(1./gamrnd( sstest, 1/srtest, [ d, 1 ] ),'Normalization','pdf','BinLimits',[0 tmax])
axis tight 

val = zeros(1,d);
subplot 414
for i = 1:d
%     subplot(1,2,1)
% N = 1;
 SigmaNew(i) = 1/gamrnd( sstest + 0.5*N, ...
                     1/( srtest+ 0.5*( e(i,:)*e(i,:)' ) ) );
                 
                 
%                  var1 = prs.ss + 0.5*N;
%                  var2 = ( prs.sr + 0.5*( (e(i,:)*e(i,:)' ) ));
% invpdftest = invgam(tnew, var1 , var2);
%    plot(tnew, invpdftest);
%        hold all
       val(i) = (e(i,:)*e(i,:)' );
end
histogram(SigmaNew,100,'Normalization','pdf','BinLimits',[0 tmax]);
axis tight 

% figure;
subplot(4,1,2)
histogram(val,'BinLimits',[0 max(val)])
% rmsval = rms(rms(e))

