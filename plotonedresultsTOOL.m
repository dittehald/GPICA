
set(0,'defaultfigurecolor',[1 1 1])
% % set(0, 'DefaultFigurePosition', [100 100 1000 700 ]);
% set(0, 'DefaultFigurePosition', [50 50 1700 900 ]);
set(0,'DefaultAxesFontSize',14)
set(0,'defaultlinelinewidth',1)

set(0,'defaulttextinterpreter','latex')
set(0,'defaulttextfontname','cmu serif')
% set(0,'defaultaxesfontname','cmu serif')
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultfigurecolor',[1 1 1])
set(0,'DefaultFigureUnits','normalized')
set(0,'DefaultFigurePosition',[0 0 0.95 0.9])

if isfield(prs,'filename') == 1
    filename = prs.filename;
end

% load Dip3711;
% X = X3711';
% l1 = 82; l2 = 68;

if prs.std == 1
X_0mean = zscore( X, 0, 2 );
else
X_0mean = bsxfun( @minus, X , mean( X , 2 ) );
end



mask = zeros(l1,l2);
[row,col] = find(Mask3711);
mask(row) = 1;

% TR = 333e-3;
TR = prs.TR;
dt = TR;
[ d ,N ] = size( X );
t = (0:dt:N*dt-dt)';


%% Stimuli
% Ref. function
pa=zeros(N,1);
for j=1:10,
    pa((30+(j-1)*121):(60+(j-1)*121))=1;
end


%%

bigZ = reshape( res.hZ, prs.m*N, prs.nsamples, prs.chains);
bigW = reshape( res.hW, prs.m*d, prs.nsamples, prs.chains);
bigSigma = reshape( res.hSigma, d, prs.nsamples+prs.skip-1,prs.chains);
bighell = reshape( res.hell, prs.m, prs.nsamples, prs.chains);
bighellC = reshape( res.hellC, prs.m, prs.nsamples+prs.skip-1, prs.chains);
bighaccr = reshape( res.haccr, prs.m, prs.nsamples+prs.skip-1, prs.chains);
endval = squeeze(bighaccr(:,end,:));
bigalpha = reshape( res.halpha, prs.m*d, prs.nsamples+prs.skip-1, prs.chains);


if  isfield(res,'hloglikelihood')
    bighloglikelihood = reshape( res.hloglikelihood, prs.nsamples+prs.skip-1, prs.chains);
    bighloglike = bighloglikelihood;
    
else
    
    bighloglike = reshape( res.hloglike, prs.m, prs.nsamples+prs.skip-1, prs.chains);

end


for cc = 1%:prs.chains
    
    Zone = bigZ(:,:,cc);
    Wone = bigW(:,:,cc);
    Sigmaone = bigSigma(:,:,cc);
    hellone = bighell(:,:,cc);
    hellCone = bighellC(:,:,cc);
    haccrone = bighaccr(:,:,cc);
    if cc>1
        haccrone = bsxfun(@minus,haccrone,endval(:,cc-1));
    end
    hloglikeone = bighloglike(:,cc);
end
%%



Wall = reshape(Wone, d ,  prs.m , prs.nsamples);
Zall = reshape(Zone, prs.m,  N  , prs.nsamples);

%% energy
Energy_over_all_samples;

% Energy_over_mean_error

%%






%%


%%

who2turn = double(max(Z') > abs(min(Z')));
who2turn(who2turn == 0) = -1;

Z = bsxfun(@times,Z,who2turn');
Zu = bsxfun(@times,Zu,who2turn');
Zl = bsxfun(@times,Zl,who2turn');

W = bsxfun(@times,W,who2turn);

%%


hell = hellone.*N.*TR;


f4 = figure('name','hist over ell');
for i=1:prs.m
    subplot( prs.m, 1, i )
    h(i)=histogram(hell(i,:)','Normalization','pdf');
    hold off
    if i==1
        BinLimits = h(1).BinLimits;
    else
        BinLimits(1) = min(BinLimits(1),h(i).BinLimits(1));
        BinLimits(2) = max(BinLimits(2),h(i).BinLimits(2));
    end
    
%     if i == 2
%         keyboard
%     end
    
end
for i=1:prs.m
    subplot( prs.m, 1, i )
    h(i)=histogram(hell(i,:)','Normalization','pdf','BinLimits',BinLimits);
    hold on
    if i == prs.m
        xlabel('Time scale $[\# \ samples]$, $\ell*N$')
    end
    
    if i == ceil(prs.m/2)+1
        ylabel('$ \ \ \ \ \ \ \ \ \ \ p(\ell)$')
    end
    
end
close(f4);

%%

convergence_diag

%%


minW = min(min(min(W)));
maxW = min(max(max(W)));

spat_map = mean(X,2);


Cmap = colormap(gray(128));
Cmap(1,:) = [1 0 0];
Cmap(end,:) = [0 0 1];

plotall = 0;


for mm = 1:prs.m
    
    sig = 2.5; % 2.5
    SMmax = W(:,mm) < quantile( W(:,mm), sig/100   , 1 );
    SMmin = W(:,mm) > quantile( W(:,mm), 1-sig/100 , 1 );
    
    SM = spat_map;
    SM(SMmax) = max(spat_map)+100;
    SM(SMmin) = min(spat_map)-100;
    
    
    mM=min(spat_map)*ones(l1*l2,1);
    mM(logical(mask))=SM;
    mM=reshape(mM,[l1 l2]);
    
    MAPS(:,:,mm) = mM;
    
end




%%

%% All together
mainfig = figure('Name',[filename '_chain' num2str(cc)]);
% figure
Cmap = colormap(gray(128));
Cmap(1,:) = [1 0 0];
Cmap(end,:) = [0 0 1];

ss = 1;
for mm = 1:prs.m
    
    subplot(prs.m,10,ss)
    imagesc(MAPS(:,:,mm)); %
    axis image;
    colormap(Cmap)
    caxis([min(spat_map)-1 max(spat_map)+1])
    
    
    %         if length(es) == 1
    %             ylabel(['MCMC ' num2str(mm)])
    %         else
    %             ylabel(['$\ell$ '  num2str(es(mm)/er(mm)*N)])
    %         end
    
    % ylabel([num2str(PartVar(mm),3) ' \%'])
    
    if mm ==1
        axMAP = gca;
        text(axMAP.XLim(1)-1.5*axMAP.XLim(2),axMAP.YLim(2)/2,...
            {['\bf{E ' num2str(ExpEng(mm),3) ' \%}']}); %,...
    else
        axMAP = gca;
        text(axMAP.XLim(1)-1.5*axMAP.XLim(2),axMAP.YLim(2)/2,...
            {['\bf{E ' num2str(ExpEng(mm),3) ' \%}'],['C ' num2str(ex(mm)*100,3) ...
            ' \%'] }); %,...
        %      'HorizontalAlignment','right')
    end
    
    
    if mm == 1;
        title('(a)')
    end
    
    % axis off;
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    
    subplot(prs.m,10,ss+1:ss+5)
    plot(t, Z(mm,:), 'b-' )
    hold on
    plot(t, Zu(mm,:), 'r-.' )
    plot(t, Zl(mm,:), 'r-.' )
    plot(t,pa,'k')
    
    axis tight
    if mm == prs.m
        xlabel('Time (s)')
    end
    
    if mm == 1;
        title('(b)')
    end
    
    
    
    subplot(prs.m,10,ss+6:ss+7)
    x =  Z(mm,:);
    freqspec
    axis tight
    
    if mm == prs.m
        xlabel('Frequency (Hz)')
    end
    
    if mm == 1;
        title('(c)')
    end
    
    set(gca, 'YTick', []);
    
    % for i=1:prs.m
    subplot(prs.m,10,ss+8:ss+9)
    
    if mm == 1
        
        for mmm = 1:prs.m;
            h(mmm)=histogram(hell(mmm,:)','Normalization','pdf');
            hold off
            
            if mmm==1
                BinLimits = h(1).BinLimits;
            else
                BinLimits(1) = min(BinLimits(1),h(mmm).BinLimits(1));
                BinLimits(2) = max(BinLimits(2),h(mmm).BinLimits(2));
            end
            
        end
    end
    
    % for i=1:prs.m
    subplot(prs.m,10,ss+8:ss+9)
    
    if hell(mm,1) == hell(mm,end)
        h(mm) = 1;
    else
        h(mm)=histogram(hell(mm,:)','Normalization','pdf','BinLimits',BinLimits);
        
    end
    
    
    hold on
    plot([1/N 1/N]*(N*TR),[0 max(h(mm).Values)],'r')
    
    if mm == 1
        TE = text(TR,max(h(mm).Values)+max(h(mm).Values)*0.15,'TR');
    end
    TE.HorizontalAlignment = 'center';
    
    axis tight
    ax = gca;
    ax.XMinorTick = 'on';
    %         ax.XTick = 1:BinLimits(2);
    %         plot([1/N 1/N]*(N*TR),[0 0.6])
    % ax.XTickLabel = {'1','','','',5,'','','','',10,'',''};
    
    %         ylim([0 0.6])
    if mm == prs.m
        xlabel('Time scale [s] , $\ell*N*TR $')
    end
    
    %         set(gca, 'YTick', []);
    set(gca, 'YAxisLocation', 'right');
    
    if mm == 1;
        title('(d)')
    end
    
    % end
    
    
    
    ss = ss+10;
    
    
end
%     end


% prompt = 'Save images (and overwrite old)? - then press 1 : ';
% savefigure = input(prompt);
% % y = x*10
% 
savefigure  = 0;
if savefigure == 1
    export_fig(mainfig, 'mainfig','-pdf')
    export_fig(logfig, ['logfig'],'-pdf')
export_fig(alphafig, 'alpha','-pdf')

end


%%
% figure;
% %    SigmaOne
% subplot(121)
% 
% plot_quantiles = [5 25 50 75 95];
% for pq = 1:length(plot_quantiles)
%     Sigma_quan=  quantile( Sigmaone, plot_quantiles(pq)/100, 1 );
%     
%     plot(Sigma_quan)
%     hold all
% end
% axinfo = gca;
% axinfo.YLim = [0,axinfo.YLim(2)];
% 
% subplot(122)
% tnew = linspace(0,axinfo.YLim(2),1000);
%%
% ss = prs.ss; 
% sr =prs.sr;
% k = ss;
% sub1 = 1; sub2 = 3;
% invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
% figure;subplot(sub1,sub2,1);
% plot(1./invgam(tnew,k,sr),tnew);
% subplot(sub1,sub2,2);
% plot(invgam(tnew,k,sr),tnew);
% % plot(invgam(1./tnew,k,sr),tnew);
% subplot(sub1,sub2,3);
% ynew = pdf('Gamma',1./tnew,prs.ss,1/prs.sr);
% plot(ynew,tnew)

%%

SigmaMedian = median(Sigmaone,2);
%%
figure; subplot(4,1,1); 
% for i = 1:size(Sigmaone,2)
% SigmaHistOne = histogram(Sigmaone(:,i),100,'BinLimits',[0 2]);
% SigmaHist(:,i) = SigmaHistOne.Values';
% end
% surf(SigmaHist)
if prs.std
tmax = 1;    
else
tmax = 500; %max(SigmaMedian);
end

histogram(SigmaMedian,100,'BinLimits',[0 tmax]);
title(['ss: ' num2str(prs.ss) ', sr: ' num2str(prs.sr) ', std: ' num2str(prs.std)])
axis tight 

subplot(4,1,2)
tnew = linspace(0,tmax,1000);
xax = gca;
ss =prs.ss; 
sr = prs.sr;
k = ss;

invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
plot(tnew, invgam(tnew,k,sr));
xlim(xax.XLim)
axis tight 



e = X_0mean-W*Z;

% histogram(abs(e),'Normalization','pdf','BinLimits',[0 tmax])
% axis tight 
% 
% %%
% figure
% subplot(2,1,1)
tnew = linspace(0,tmax,1000);
sstest =prs.ss;
srtest = prs.sr;

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
% figure;
%  histogram(val)
% rmsval = rms(rms(e))

histogram(SigmaNew,100,'Normalization','pdf','BinLimits',[0 tmax]);
axis tight 

%%
% hold all;
% pdf('Gamma',tnew,prs.ss(1),1/prs.sr(1))
% plot(1./invgam(tnew,k,sr),tnew);
%xlim(xax.XLim)

% ylim([axinfo.YLim])

%%
% % % figure; subplot(2,1,1); 
% % % SigmaMedian = median(Sigmaone,2);
% % % 
% % % for i = 1:size(Sigmaone,2)
% % % SigmaHistOne = histogram(Sigmaone(:,i),100,'BinLimits',[0 2]);
% % % SigmaHist(:,i) = SigmaHistOne.Values';
% % % end
% % % surf(SigmaHist)
% % % % histogram(SigmaMedian,100);
% % % 
% % % subplot(2,1,2)
% % % tmax = max(SigmaMedian);
% % % tnew = linspace(0,2,1000);
% % % xax = gca;
% % % title('')
% % % 
% % % ss =prs.ss; 
% % % sr = prs.sr;
% % % k = ss;
% % % 
% % % invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
% % % plot(tnew, 1./invgam(tnew,k,sr));
% % % xlim(xax.XLim)
%%
% %%
% % figure; 
% subplot(2,2,3)
% histogram(1./gamrnd( prs.ss, 1/prs.sr, [ d, 1 ] ),'Normalization','pdf')
% xlim([0 1])
% 
% i = 1000;
% e = X_0mean-WZ;
%  Sigma = 1/gamrnd( prs.ss + 0.5*N, ...
%                     1/( prs.sr + 0.5*( (e(i,:)*e(i,:)' ) )));
% subplot(2,2,4)
% 
% % sstest = N/2;
% % srtest = 200;
% % 
% % figure
% % tnew = linspace(0,1,10000);
% invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
% 
% figure
% subplot(2,1,1)
% tnew = linspace(0,70,1000);
% for i = 1:d
% %     subplot(1,2,1)
% N = 1;
%  SigmaNew(i) = 1/gamrnd( prs.ss + 0.5*N, ...
%                      1/( prs.sr+ 0.5*( e(i,:)*e(i,:)' ) ) );
%                  
%                  
%                  var1 = prs.ss + 0.5*N;
%                  var2 = ( prs.sr + 0.5*( (e(i,:)*e(i,:)' ) ));
% invpdftest = invgam(tnew, var1 , var2);
%    plot(tnew, invpdftest);
%        hold all
%        val(i) = (e(i,:)*e(i,:)' );
% end
% 
% subplot(2,1,2)
% histogram(SigmaNew,100);


% %%
% figure;
% tnew = linspace(0,500,10000);
% plot(tnew, invgam(tnew, prs.ss + 0.5*N,1/( prs.sr + 0.5*( (e(i,:)*e(i,:)' ) ))));

% subplot(1,3,3); 
% SigmaMedian = median(Sigmaone,2);
% histogram(SigmaMedian,100);



%%
% figure
% tnew = linspace(0,2,1000);
% ss = .1*N/2;
% sr = rms(rms(X_0mean))/rms(rms(X_0mean-WZ));
% ynew2 = (pdf('Gamma',tnew,ss,1/sr));
% plot(tnew,ynew2)
% %%
% figure
% tnew = linspace(0,10,1000);
% ss = 1;
% sr = 1;
% ynew2 = 1./(pdf('Gamma',tnew,ss,1/sr));
% % ynew = (pdf('Gamma',325,ss,1/sr))
% plot(tnew,ynew2)
% hold all
% 
% invgamval = 1./gamrnd(ss, sr, [ 1000, 1 ] );
%  histogram(invgamval,'Normalization','pdf')
% 
% hold all
% invgamval = 1./gamrnd(ss, 1/sr, [ 1000, 1 ] );
% figure; histogram(invgamval,'Normalization','pdf')
% %%
% tmax =100;
% tnew = linspace(0,tmax,10000);
% 
% ss =100; sr = 1;
% 
% sz =[10000 1];
% 
% alpha = ss;
% beta = 1/sr;
% 
% x = gamrnd(alpha,beta,sz);
% k= alpha;
% theta = 1/beta;
% r = 1./gamrnd(k,beta,sz);
% 
% figure; 
% subplot(2,2,1);
% histogram(x,'Normalization','pdf')
% xlim([0 tmax])
% ax1 = gca;
% subplot(2,2,2)
% histogram(r,100000,'Normalization','pdf')
% xlim([0 tmax])
% ax2 = gca;
% 
% subplot(2,2,3)
% ynew2 = pdf('Gamma',tnew,alpha,beta);
% plot(tnew,ynew2)
% xlim(ax1.XLim);
% 
% subplot(2,2,4)
% % ynew2 = 1./(pdf('Gamma',tnew,k,1/theta));
% % plot(tnew,ynew2)
% invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
% plot(tnew, invgam(tnew,k,theta));
% xlim(ax2.XLim);


%%

% tnew = linspace(300,400,1000);
% ss = 1;
% sr = 1;
% ynew2 = (pdf('InvGamma',tnew,ss,1/sr));
% % ynew = (pdf('Gamma',325,ss,1/sr))
% plot(tnew,ynew2)

%ylim([0 1])
%%

% ynew2 = (pdf('Gamma',tnew,prs.ss(1),1/1));

% 
% Gmean = 275;
% GVar = 200;
% 
% ynew2 = (pdf('Gamma',tnew,Gmean^2/GVar,1/(Gmean/GVar )  ));

% 
%%

alphaall = reshape(bigalpha(:,:,1), d ,  prs.m , prs.nsamples+prs.skip-1);
 alphaall = alphaall(:,indx(K:-1:1),:);

alphafig = figure;
plot_quantiles = [50 5 95];
mcor = 'rgbcmk';
linsty{1}= '-'; linsty{2}='--'; linsty{3}=':';
for mm = 1: prs.m
for pq = 1:length(plot_quantiles)
    alpha_quan=  quantile( squeeze(alphaall(:,mm,:)), plot_quantiles(pq)/100, 1 );
    
    subplot(2,6,1:3)
    plothandle(mm,pq) = plot(squeeze(alpha_quan),'color',mcor(mm), 'LineStyle' , (linsty{pq}));
    hold all
    
   
   
    
   
end

    
    
hh = subplot(2,6,6+mm);
pp = get(hh,'pos');
pp(2) = pp(2) + 0.04;
set(hh, 'pos', pp);


alphapic = median(alphaall(:,mm,prs.skip:end),3);
% alphapic = W(:,mm);
% alphapic = bsxfun(@times,alphapic,who2turn(mm));
 mM=min(alphapic)*ones(l1*l2,1);
    mM(logical(mask))=alphapic;
    mM=reshape(mM,[l1 l2]);

 imagesc(mM); axis image; axis off;
 titlehandle = title(['GP ' num2str(mm)]);
 titlehandle.FontSize = 12;
% imghandel = gca;
% if mm == 1
% text(imghandel.XLim(1),imghandel.YLim(2), '(c) $\alpha$ median')
% end

end

subplot(2,6,1:3)
leghandle = legend(plothandle(:,1),'GP 1','GP 2','GP 3','GP 4','GP 5','GP 6');


subplot(2,6,1:3)
axinfo = gca;
axinfo.XLim = [-50,axinfo.XLim(2)];
axinfo.YLim = [0,axinfo.YLim(2)];
xlabel('Iterations (incl warm-up)') 
ylabel('$\alpha$')
title('(a)')

subplot(2,6,4:6)
tnew = linspace(0,axinfo.YLim(2),1000);
ynew2 = (pdf('Gamma',tnew,prs.ts(1),1/prs.tr(1)));
 plot(ynew2,tnew)
ylim([axinfo.YLim])
xlabel('Prior $\alpha$ probability') 
% ylabel('$\alpha$')
title('(b)')

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

text(0.5, 0.5,'(c)','HorizontalAlignment' ...
    ,'center','VerticalAlignment', 'top','FontSize',16)
return

% % %%
% % figure;
% % Z1 = Z(1,:);
% % triallength = 121;
% % trialnr = 10;
% % totallength = triallength*trialnr;
% % 
% % pa1 = [zeros(1,30) ones(1,31) zeros(1,60)];
% % pa2 = repmat(pa1,1,10);
% % 
% % subplot(2,1,1)
% % % Z1 = Z1(1:totallength);
% % Boldimage = reshape(Z1,triallength,trialnr);
% % Boldimage = Boldimage';
% % imagesc(Boldimage./max(max(Boldimage)))
% % colorbar
% % gca;
% % subplot(2,1,2)
% % 
% % Boldmean = mean(Boldimage);
% % 
% % Baselinemean = mean(Boldmean(1:30));
% % 
% % plot((Boldmean-Baselinemean)./max((Boldmean-Baselinemean)))
% % 
% % % zeros(length(hrf),1);
% % acti1c = conv(  pa2(1:triallength),hrf);
% % acti1c = acti1c(1:triallength);
% % 
% % % acti1c2 = conv( hrf, [zeros(length(hrf),1); pa(1:triallength)] );
% % % acti1c2 = acti1c2((length(hrf)):(length(hrf)-2)+length(pa(1:triallength)));
% % 
% % hold all;
% % plot(acti1c./max(acti1c))
% % axis tight
% % 
% % % plot(acti1c2./max(acti1c2))
% % 
% % plot(pa2(1:triallength))
% % %% =========================================
% % 
% % return
% % showvideo =0;
% % 
% % if showvideo == 1
% %     figure
% %     for iter = 1:size(Sigmaone,2)
% %         
% %         mM=ones(l1*l2,1);
% %         mM(logical(mask))=Sigmaone(:,iter);
% %         %         mM=reshape(mM,[l1 l2]);
% %         
% %         imagesc(reshape(mM,l1,l2))
% %         axis image
% %         axis off
% %         caxis([min(min(Sigmaone))-100 max(max(Sigmaone))])
% %         drawnow
% %     end
% % end
% % 
% % %
% % %    figure
% % %    mM=min(spat_map)*ones(l1*l2,1);
% % %         mM(logical(mask))=spat_map;
% % %         mM=reshape(mM,[l1 l2]);
% % %    imagesc(mM); axis image; axis off; colormap gray
% % 
% % %%
% % shape = prs.es(1);
% % scale = 1/prs.er(1);
% % % scale = 1/800
% % tnew = (0:1/(N*1000):1)'; %N*dt-dt)';
% % 
% % ell = bighellC(:,1);  %gamrnd( prs.es, 1./prs.er)';
% % 
% % PPplot = figure;
% % tsek = tnew.*(TR*N);
% % alpha = TR*N;
% % ynew = (alpha).*pdf('Gamma',tsek,shape,scale*alpha);
% % plot(tsek,ynew,'linewidth',2.5)
% % hold all
% % 
% % for mm= 1: length(ell)
% %     ellPsek(:,mm)  = alpha.*pdfrG(tsek,ell(mm)*alpha,prs.ess*alpha);
% %     
% %     plot(tsek,ellPsek(:,mm),'color',[1 0 0 ],'linewidth',2)
% % end
% % xlim([0 0.0035*alpha])
% % % ylim([0 1000])
% % 
% % 
% % xlabel('$\ell$ [s]','FontSize',28)
% % ylabel('Probability', 'FontSize',28)
% % title('Prior for $\ell$ and the initial proposal distributions', 'FontSize',28)
% % lege = legend(['$p( \ell_p ) \sim  \Gamma ( $' num2str(prs.es(1)) ', 1/' num2str(prs.er(1)) '$ )$'],...
% %     ['$q(\ell_p) \sim  \mathcal{N}_{Ref} (\ell,$' num2str(prs.ess) ' )']);
% % lege.FontSize = 30;
% % 
% % %%
% % 
% % NN = 1;
% % mumean = shape*scale;
% % % samples per second
% % %    dtt = N*TR/NN;                     % seconds per sample
% % %    tpts = (0:dtt:(N*TR-TR)/10);
% % %  dss = sq_dist( tpts, tpts );
% % % K = exp( -dss/(2*mumean^2) );
% % %
% % %
% % % f2 = figure('Name','Kernel for GP');
% % % surf(tpts,tpts,K)
% % % title('Kernel for GP, $\ell$ = $s_e / r_e$')
% % % xlabel('$t_i$')
% % % ylabel('$t_j$')
% % %%
% % % close all
