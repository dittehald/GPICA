
set(0,'defaultfigurecolor',[1 1 1])
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
% TR = 333e-3;

if prs.std == 1
    X_0mean = zscore( X, 0, 2 );
else
    X_0mean = bsxfun( @minus, X , mean( X, 2 ) );
end

mask = zeros(l1,l2);
[row,col] = find(Mask3711);
mask(row) = 1;

dt = prs.TR;
[ d ,N ] = size( X );
t = (0:dt:N*dt-dt)';

% Ref. function
pa=zeros(N,1);
for j=1:10,
    pa((30+(j-1)*121):(60+(j-1)*121))=1;
end


bigZ = reshape( res.hZ, prs.m*N, prs.nsamples, prs.chains);
bigW = reshape( res.hW, prs.m*d*( prs.lags+1), prs.nsamples, prs.chains);

bighloglikelihood = reshape( res.hloglikelihood, prs.nsamples+prs.skip-1, prs.chains);
bighloglike = bighloglikelihood;
bigSigma = reshape( res.hSigma, d, prs.nsamples+prs.skip-1, prs.chains);
bighell = reshape( res.hell, prs.m, prs.nsamples, prs.chains);
bighellC = reshape( res.hellC, prs.m, prs.nsamples+prs.skip-1, prs.chains);
bighaccr = reshape( res.haccr, prs.m, prs.nsamples+prs.skip-1, prs.chains);
endval = squeeze(bighaccr(:,end,:));

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



Energy_conv_allsamples


%% No turning??

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

%%

%
% logfig = figure('Position',[0 0.2 1 0.6]);
%
% hh2=subplot(1,4,2);
% pp = get(hh2,'pos')
% pp(3) = pp(3)*0.5;
% pp(1) = pp(1)*1.2;
% set(hh2, 'pos', pp);
%
% tnew = linspace(0,0.1,1000);
% ynew2 = pdf('Gamma',tnew,prs.es(1),1/prs.er(1));
%  plot(ynew2,tnew)
% % ylim(yl)
% set(gca,'XTickLabel','')
%  set(gca,'YTickLabel','')
%  title('Prior')
% xlabel('Prop.')
%
% hh1=subplot(1,4,1);
% pp1 = get(hh1,'pos')
% pp1(3) = pp1(3)*1.5;
% set(hh1, 'pos', pp1);
%
%
% plot(hellCone')
% yl=get(gca,'ylim');
% title('Developement of $\ell$')
%
% xlabel(['Iterations (incl burn-in of ' num2str(prs.skip) ')'])
% ylabel(['$\ell$'])
%
%
%  set(hh2, 'ylim', yl);
%
% hh =subplot(1,4,3);
% pp = get(hh,'pos');
% % pp(3) = pp(3) + 0.05;
% % pp(1) = pp(1) - 0.05;
% set(hh, 'pos', pp);
%
% accrates = round(haccrone./(repmat(1:(prs.skip+prs.nsamples-1),prs.m,1)).*100);
%
% plot(accrates')
% title('Acceptance rate developement')
%
% accrate = round(haccrone(:,end)./(prs.skip+prs.nsamples-1).*100);
%
% xlabel(['Iterations (incl burn-in of ' num2str(prs.skip) ')'])
% ylabel(['Acceptance rate in \%'])
%
%
%             subplot(1,4,4)
%             plot(hloglikeone); axis tight;
%             title('Loglikelihood')
%             xlabel(['Iterations (incl burn-in of ' num2str(prs.skip) ')'])
%             ylabel('$p(z_{[p]}|\ell_p)$')
%             xlim([-10 prs.nsamples+prs.skip-1])
%             ylim([min(min(hloglikeone(200:end)))-100 max(max(hloglikeone(prs.skip:end)))+100])
% %          return
%         if  isfield(res,'hloglikelihood')
%             ylabel('$p(X|WZ,\Psi)$')
%         end


convergence_diag

%%

Zalllags = lags_shift(Z,prs.lags);

Wmean = zeros(prs.m,prs.lags+1);
for mm = 1%:prs.m
    f5(mm) = figure('Name','MCMCica maps');
    
    for la = 1:prs.lags+1
        
        
        Wextract = W(:,mm,la);
        
        minW = min(min(min(Wextract)));
        maxW = max(max(max(Wextract)));
        
        
        m1=minW*ones(l1*l2,1);
        m1(logical(mask))=Wextract;
        m11=reshape(m1,[l1 l2]);
        figure(f5(mm))
        subplot(2,prs.lags+1,la)
        imagesc(m11); title(['Lag ' num2str(la-1)]); axis off; axis image;
        caxis([minW maxW])
        
        Wmean(mm,la) = norm(Wextract*Zalllags(mm,:,la))^2;
        
    end
    figure(f5(mm))
    meanplot =subplot(2,prs.lags+1,la+1:la+1+prs.lags);
     plot(0:prs.lags,(Wmean(mm,:))./Wmean(mm,1),'-*');
    
     ax= gca;
     
    hold all
    
    xlim([-0.5 prs.lags+0.5])
    meanplot.XTick = 0:prs.lags;
    meanplot.XTickLabel = {'0','1','2','3','4'};
    
    xlabel('Lags ($\tau $)')
ylabel('$E [ X_\tau ] / E [ X_0 ]$')
    
%     ax.XLabel = 'Lags ($\tau $)';
%    ax.YLabel = '$X_\tau $';
    % histogram(Wextract(:))
    % ylim([0 100])
    
    
end

%%


minW = min(min(min(W)));
maxW = min(max(max(W)));

spat_map = mean(X,2);


Cmap = colormap(gray(128));
Cmap(1,:) = [1 0 0];
Cmap(end,:) = [0 0 1];

plotall = 0;


for mm = 1:prs.m
    
    for la = 1:prs.lags+1
        sig = 2.5; % 2.5
        SMmax = W(:,mm,la) < quantile( W(:,mm,la), sig/100   , 1 );
        SMmin = W(:,mm,la) > quantile( W(:,mm,la), 1-sig/100 , 1 );
        
        SM = spat_map;
        SM(SMmax) = max(spat_map)+100;
        SM(SMmin) = min(spat_map)-100;
        
        
        mM=min(spat_map)*ones(l1*l2,1);
        mM(logical(mask))=SM;
        mM=reshape(mM,[l1 l2]);
        
        MAPS(:,:,mm,la) = mM;
        
    end
end


%%
% red = colormap(autumn(64));
% blue = flipud(colormap(winter(64)));
% Twocolors =[red; blue];
% Twocolors(1,:) = [0.8 0.8 0.8];
% Twocolors(end,:) = [0 0 0 ];
% for mm = 1:prs.m
%
%     for la = 1 :prs.lags+1
% figure
%
% minW = min(min(min(W)));
% maxW = max(max(max(W)));
%
% minW = min(W(:,mm,la));
% maxW = max(W(:,mm,la));
%
% m1=minW*ones(l1*l2,1);
% m1(logical(mask))=W(:,mm,la);
% m11=reshape(m1,[l1 l2]);
% figure(f5(mm))
% subplot(1,prs.lags+1,la)
% imagesc(m11); title(['MCMC ' num2str(mm)]); axis off; axis image;
% caxis([minW maxW])
% figure
% imagesc(MAPS(:,:,mm,la))
% axis off; axis image;
%
%         SMmax = W(:,mm,la) < quantile( W(:,mm,la), sig/100   , 1 );
%         SMmin = W(:,mm,la) > quantile( W(:,mm,la), 1-sig/100 , 1 );
%
%         Wsingle = W(:,mm,la);
%         Part1 = Wsingle(SMmax);
%         Part1 = Part1./(max(abs(Part1)));
%         Part2 = Wsingle(SMmin);
%         Part2 = Part2./(max(abs(Part2)));
%
%         if abs(sum(sign(Part1))) ~= length(Part1) || ...
%                 abs(sum(sign(Part2))) ~= length(Part2)
%             keyboard
%         end
%
%         SM = ones(d,1).*min(Part1).*1.10;
%         SM(SMmax) =Part1;
%         SM(SMmin) =Part2;
%
%         mM=ones(l1*l2,1).*max(Part2).*1.10;
%         mM(logical(mask))=SM;
%         mM=reshape(mM,[l1 l2]);
%
%          MAPS(:,:,mm,la) = mM;
%         figure; imagesc(mM)
%         axis off; axis image;
%         colormap(Twocolors)
%     end
% end





%% All together
mainfig = figure('Name',[filename '_chain' num2str(cc)]);
% figure
Cmap = colormap(gray(128));
Cmap(1,:) = [1 0 0];
Cmap(end,:) = [0 0 1];

ss = 1;
numbersubplot = prs.lags+8;
for mm = 1:prs.m
    for la = 1:prs.lags+1
        subplot(prs.m,numbersubplot,ss)
        imagesc(MAPS(:,:,mm,la)); %
        axis image;
        %         colormap(Twocolors)
        colormap(Cmap)
        caxis([min(spat_map)-1 max(spat_map)+1])
        
        
        %         if length(es) == 1
        %             ylabel(['MCMC ' num2str(mm)])
        %         else
        %             ylabel(['$\ell$ '  num2str(es(mm)/er(mm)*N)])
        %         end
        if la == 1
            % ylabel({[num2str(ExpEng(mm),3) ' \%'],[num2str(ex(mm)*100,3) ...
            %      ' \%'] , ['(+' num2str(dex(mm)*100,3) ' \% )' ]})
            % %                 axMAP = gca;
            % % axMAP.YTickLabelRotation=90;
            % % keyboard
            % set(get(gca,'YLabel'),'Rotation',0);
            
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
            
        end
        
        
        if mm == prs.m;
            xlabel(['lag ' num2str(la-1)])
        end
        
        if mm == 1 && la == 3;
            title('(a)')
        end
        
        % axis off;
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        
        ss = ss + 1;
    end % end lags
    ss = ss-1;
    %        keyboard
    stime = subplot(prs.m,numbersubplot,ss+1:ss+5);
    
    pp = get(stime,'pos');
pp(1) = pp(1) + 0.01;
pp(3) = pp(3) - 0.01;
set(stime, 'pos', pp);
    
    plot(t, Z(mm,:), 'b-' )
    hold on
    plot(t, Zu(mm,:), 'r-.' )
    plot(t, Zl(mm,:), 'r-.' )
    %         plot(t,pa,'k')
    
    axis tight
    if mm == prs.m
        xlabel('Time (s)')
    end
    
    if mm == 1;
        title('(b)')
    end
    
    
    subplot(prs.m,numbersubplot,ss+6:ss+7)
    
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
    subplot(prs.m,numbersubplot,ss+6:ss+7)
    
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
        title('(c)')
    end
    
    % end
    
    
    %         keyboard
    ss = ss+8;
    
    
end
%%

prompt = 'Save images (and overwrite old)? - then press 1 : ';
savefigure = input(prompt);
% y = x*10

if savefigure == 1
    export_fig(mainfig, 'mainfig_conv_std','-pdf')
    export_fig(logfig, 'logfig_conv_std','-pdf')
end

%%
singleplot = figure;
ss = 1;
numbersubplot = 5;
for mm = 1
    for la = 1:prs.lags+1
       
        
 subplot(2,numbersubplot,ss)
        imagesc(MAPS(:,:,mm,la)); %
        axis image;
        %         colormap(Twocolors)
        colormap(Cmap)
        caxis([min(spat_map)-1 max(spat_map)+1])

      
%         if la == 1
%             
%             if mm ==1
%                 axMAP = gca;
%                 text(axMAP.XLim(1)-1.5*axMAP.XLim(2),axMAP.YLim(2)/2,...
%                     {['\bf{E ' num2str(ExpEng(mm),3) ' \%}']}); %,...
%             else
%                 axMAP = gca;
%                 text(axMAP.XLim(1)-1.5*axMAP.XLim(2),axMAP.YLim(2)/2,...
%                     {['\bf{E ' num2str(ExpEng(mm),3) ' \%}'],['C ' num2str(ex(mm)*100,3) ...
%                     ' \%'] });
%             end
%             
%         end
        
        
        if mm == 1;
            xlabel(['lag ' num2str(la-1)])
        end
        
        if mm == 1 && la == 3;
            title('(a)')
        end
        
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        
        ss = ss + 1;
    end % end lags
    ss = ss-1;
    
    stime = subplot(2,numbersubplot,ss+1:ss+5);
    
%     pp = get(stime,'pos');
% pp(1) = pp(1) + 0.01;
% pp(3) = pp(3) - 0.01;
% set(stime, 'pos', pp);
    
    plot(t, Z(mm,:), 'b-' )
    hold on
    plot(t, Zu(mm,:), 'r-.' )
    plot(t, Zl(mm,:), 'r-.' )
             plot(t,pa.*0.2,'k')
    
    axis tight
    if mm == 1
        xlabel('Time (s)')
    end
    
    if mm == 1;
        title('(b)')
    end
end


%%
figure;
%    SigmaOne
subplot(121)
plot_quantiles = [5 25 50 75 95];
for pq = 1:length(plot_quantiles)
    Sigma_quan=  quantile( Sigmaone, plot_quantiles(pq)/100, 1 );
    
    plot(Sigma_quan)
    hold all
end

axinfo = gca;
axinfo.YLim = [0,axinfo.YLim(2)];

subplot(122)
tnew = linspace(0,axinfo.YLim(2),1000);
ynew2 = (pdf('Gamma',tnew,prs.ss(1),1/prs.sr(1)));
 plot(ynew2,tnew)
ylim([axinfo.YLim])

figure; 
plot(tnew,ynew2)

%%
showvideo =0;

if showvideo == 1
    figure
    for iter = 1:size(Sigmaone,2)
        
        mM=ones(l1*l2,1);
        mM(logical(mask))=Sigmaone(:,iter);
        %         mM=reshape(mM,[l1 l2]);
        
        imagesc(reshape(mM,l1,l2))
        axis image
        axis off
        caxis([min(min(Sigmaone))-100 max(max(Sigmaone))])
        drawnow
    end
end

%
%    figure
%    mM=min(spat_map)*ones(l1*l2,1);
%         mM(logical(mask))=spat_map;
%         mM=reshape(mM,[l1 l2]);
%    imagesc(mM); axis image; axis off; colormap gray

%%
shape = prs.es(1);
scale = 1/prs.er(1);

tnew = (0:1/(N*1000):1)'; %N*dt-dt)';

ell = bighellC(:,1);  %gamrnd( prs.es, 1./prs.er)';

PPplot = figure;
tsek = tnew.*(TR*N);
alpha = TR*N;
ynew = (alpha).*pdf('Gamma',tsek,shape,scale*alpha);
plot(tsek,ynew)
hold all

for mm= 1: length(ell)
    ellPsek(:,mm)  = alpha.*pdfrG(tsek,ell(mm)*alpha,prs.ess*alpha);
    
    plot(tsek,ellPsek(:,mm),'color',[1 0 0 ])
end
xlim([0 0.003*alpha])
% ylim([0 1000])


xlabel('$\ell$ [s]')
ylabel('Probability')
title('Prior for $\ell$ and the initial proposal distributions')
legend(['$p( \ell_p ) \sim  \Gamma ( $' num2str(prs.es(1)) ', 1/' num2str(prs.er(1)) '$ )$'],...
    ['$q(\ell_p) \sim  \mathcal{N}_{Ref} (\ell,$' num2str(prs.ess) ' )'])
%%

NN = 1;
mumean = shape*scale;
% samples per second
%    dtt = N*TR/NN;                     % seconds per sample
%    tpts = (0:dtt:(N*TR-TR)/10);
%  dss = sq_dist( tpts, tpts );
% K = exp( -dss/(2*mumean^2) );
%
%
% f2 = figure('Name','Kernel for GP');
% surf(tpts,tpts,K)
% title('Kernel for GP, $\ell$ = $s_e / r_e$')
% xlabel('$t_i$')
% ylabel('$t_j$')
%%
% close all
return

% Twocolors
red = colormap(autumn(64));
blue = flipud(colormap(winter(64)));
Twocolors =[red; blue];
Twocolors(1,:) = [0.8 0.8 0.8];
Twocolors(end,:) = [0 0 0 ];
for mm = 1
    
    for la = 1 :prs.lags+1
        % figure
        
        % minW = min(min(min(W)));
        % maxW = max(max(max(W)));
        
        minW = min(W(:,mm,la));
        maxW = max(W(:,mm,la));
        
        m1=minW*ones(l1*l2,1);
        m1(logical(mask))=W(:,mm,la);
        m11=reshape(m1,[l1 l2]);
        % figure(f5(mm))
        % subplot(1,prs.lags+1,la)
        % imagesc(m11); title(['MCMC ' num2str(mm)]); axis off; axis image;
        % caxis([minW maxW])
        % figure
        % imagesc(MAPS(:,:,mm,la))
        % axis off; axis image;
        
        SMmax = W(:,mm,la) < quantile( W(:,mm,la), sig/100   , 1 );
        SMmin = W(:,mm,la) > quantile( W(:,mm,la), 1-sig/100 , 1 );
        
        Wsingle = W(:,mm,la);
        Part1 = Wsingle(SMmax);
        Part1 = Part1./(max(abs(Part1)));
        Part2 = Wsingle(SMmin);
        Part2 = Part2./(max(abs(Part2)));
        
        if abs(sum(sign(Part1))) ~= length(Part1) || ...
                abs(sum(sign(Part2))) ~= length(Part2)
            keyboard
        end
        
        SM = ones(d,1).*min(Part1).*1.10;
        SM(SMmax) =Part1;
        SM(SMmin) =Part2;
        
        mM=ones(l1*l2,1).*max(Part2).*1.10;
        mM(logical(mask))=SM;
        mM=reshape(mM,[l1 l2]);
        
        figure; imagesc(mM)
        axis off; axis image;
        colormap(Twocolors)
    end
end


