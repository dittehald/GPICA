
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

plotall = 0;

% addpath(genpath([startpath '/home/ditha/network_perm_test/spm8']))

% X_0mean = bsxfun( @minus, X , mean( X, 2 ) );
% X_0mean = zscore( X, 0, 2 );
if prs.std == 1
X_0mean = zscore( X, 0, 2 );
else
X_0mean = bsxfun( @minus, X , mean( X, 2 ) );
end

[ d ,N ] = size( X );

dt = TR;
t = (0:dt:N*dt-dt)';

mask1 = zeros(l1*l2*l3,1);
mask1(Indx) = 1;
mas = mask1;
mask1(mask1 == 0) = nan;


bigZ = reshape( res.hZ, prs.m*N, prs.nsamples, prs.chains);
bigW = reshape( res.hW, prs.m*d, prs.nsamples, prs.chains);

bighloglikelihood = reshape( res.hloglikelihood, prs.nsamples+prs.skip-1, prs.chains);
bighloglike = bighloglikelihood;
if isfield(res, 'hSigma')
    if size(res.hSigma,2)>prs.nsamples
        bigSigma = reshape( res.hSigma, d, prs.nsamples + prs.skip-1, prs.chains);
    else
        bigSigma = reshape( res.hSigma, d, prs.nsamples, prs.chains);
        
    end
else
    bigSigma = 0;
end
bighell = reshape( res.hell, prs.m, prs.nsamples, prs.chains);
bighellC = reshape( res.hellC, prs.m, prs.nsamples+prs.skip-1, prs.chains);
bighaccr = reshape( res.haccr, prs.m, prs.nsamples+prs.skip-1, prs.chains);
endval = squeeze(bighaccr(:,end,:));
if isfield(res, 'halpha')
    if d > 10000
        bigalpha = reshape( res.halpha, prs.m*d, prs.chains);
        %ones( prs.m*d, prs.nsamples+prs.skip-1);
    else
        bigalpha = reshape( res.halpha, prs.m*d, prs.nsamples+prs.skip-1, prs.chains);
    end
end

for cc = 1%:prs.chains
    
    Zone = bigZ(:,:,cc);
    Wone = bigW; %(:,:,cc)
    Sigmaone = bigSigma(:,:,cc);
    hellone = bighell(:,:,cc);
    hellCone = bighellC(:,:,cc);
    haccrone = bighaccr(:,:,cc);
    if cc>1
        haccrone = bsxfun(@minus,haccrone,endval(:,cc-1));
    end
    hloglikeone = bighloglike(:,cc);
    
end

if prs.nsamples > 800
    clear bigZ bigW bigSigma bighell bighellC bighaccr bighloglike
    %res= rmfield(res,{'hZ','hW','hSigma'});
end

%%


Wall = reshape(Wone, d ,  prs.m , prs.nsamples);
Zall = reshape(Zone, prs.m,  N  , prs.nsamples);

%% energy
% Energy_over_all_samples;

Energy_over_mean_error



%%

% Zone = reshape(Zone, N*prs.m, prs.nsamples );
% Wone = reshape(Wone, d*prs.m, prs.nsamples );
% 
% Z = reshape( median( Zone, 2 ), prs.m, N );
% 
% % Errorbars
% Zu = reshape( quantile( Zone, 0.95, 2 ), prs.m, N );
% Zl = reshape( quantile( Zone, 0.05, 2 ), prs.m, N );
% 
% 
% clear Zone
% 
% if  prs.nsamples>500
%     Wr1 =  median( Wone( 1:d*prs.m/4,:), 2 );
%     Wr2 = median( Wone( d*prs.m/4+1:2*d*prs.m/4,:), 2 );
%     Wr3 = median( Wone( 2*d*prs.m/4+1:3*d*prs.m/4,:), 2 );
%     Wr4 = median( Wone( 3*d*prs.m/4+1:end,:), 2 );
%     W = [Wr1; Wr2; Wr3; Wr4];
%     W = reshape( W, d, prs.m);
% else
%     W = reshape( median( Wone, 2 ), d, prs.m);
% end
% 
% 
% %% energy
% [d , N] = size(X_0mean);
% Xback = zeros(d,N,prs.m);
% EnergyX = norm(X_0mean,'fro')^2;
% 
% for m=1:prs.m
%     Xback(:,:,m) = (W(:,m)*Z(m,:));
%     Residual(m) = norm(X_0mean-Xback(:,:,m),'fro')^2; % all samples
% end
% 
% [vals,ordr] = sort(Residual);
% 
% ExpEng = 100-100*Residual./EnergyX;
% [vals5,ordr5] = sort(ExpEng);
% 
% 
% K = prs.m;
% 
% W = W(:,ordr5(K:-1:1)); % same as ordr
% Z = Z(ordr5(K:-1:1),:);
% ExpEng = vals5(K:-1:1);
% Zl = Zl(ordr5(K:-1:1),:);
% Zu = Zu(ordr5(K:-1:1),:);
% 
% hellone = hellone(ordr5(K:-1:1),:);
% hellCone = hellCone(ordr5(K:-1:1),:);
% haccrone = haccrone(ordr5(K:-1:1),:);
% 
% if length(prs.es) ~= 1
%     es = prs.es(ordr5(K:-1:1));
%     er = prs.er(ordr5(K:-1:1));
% else
%     es = prs.es;
%     er = prs.er;
% end
% 
% dex = zeros(prs.m,1);
% ex = zeros(prs.m,1);
% 
% WZ = zeros(d,N);
% 
% for m=1:prs.m
%     WZ = WZ + (W(:,m)*Z(m,:));
%     
%     ex(m) = 1 -   ( norm(X_0mean-WZ,'fro')^2 / EnergyX ); % note it is the variance explained by sources 1:m
%     if m==1
%         dex(1) = ex(1);
%     else
%         dex(m) = ex(m)-ex(m-1);
%     end
% end
% 


%%

who2turn = double(max(Z') > abs(min(Z')));
who2turn(who2turn == 0) = -1;

Z = bsxfun(@times,Z,who2turn');
Zu = bsxfun(@times,Zu,who2turn');
Zl = bsxfun(@times,Zl,who2turn');

W = bsxfun(@times,W,who2turn);
%     Wu = bsxfun(@times,Wu,who2turn);
%     Wl = bsxfun(@times,Wl,who2turn);



% return

% % Sort according to variance
% S = Z;
% K=prs.m;
% Wr = reshape( median( res.hW, 2 ), d, prs.m);
% A = Wr;
%
% Avar=diag(A'*A)/d;
% Svar=diag(S*S')/N;
% vS=var(S');
% sig=Avar.*Svar;
% [a,indx]=sort(sig);
%
% Z=S(indx(K:-1:1),:);
% W=A(:,indx(K:-1:1));
% Zu=Zu(indx(K:-1:1),:);
% Zl=Zl(indx(K:-1:1),:);
% res.hell = res.hell(indx(K:-1:1),:);
% res.accr = res.accr(indx(K:-1:1),:);
%



%%


hell = hellone.*N.*TR;
hell = sort(hell);
f4 = figure('name','hist over ell');
for i=1:prs.m
    subplot( prs.m, 1, i )
    h(i)=histogram(hell(i,:)','Normalization','probability');
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
    h(i)=histogram(hell(i,:)','Normalization','probability','BinLimits',BinLimits);
    hold on
    if i == prs.m
        xlabel('Time scale $[\# \ samples]$, $\ell*N$')
    end
    
    if i == ceil(prs.m/2)+1
        ylabel('$ \ \ \ \ \ \ \ \ \ \ p(\ell)$')
    end
    xlim([1 1.8])
end


%%

convergence_diag

%
% logfig = figure('Position',[0 0.2 1 0.6]);
%
%
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
% % hold all
%
%
%  set(hh2, 'ylim', yl);
%
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
% %      val_ell = es./er;
% %     legtex = [repmat('GP ',prs.m,1) num2str([1:prs.m]')  ...
% %         repmat(' - $\ell$ ',prs.m,1) num2str(val_ell')];
% %     legend(cellstr(legtex))
% %
%
%
% % legtex = [repmat('GP ',prs.m,1) num2str([1:prs.m]')  ...
% %     repmat(' - acc rate ',prs.m,1) num2str(accrate) repmat('\%',prs.m,1)];
% % legend1 = legend(cellstr(legtex),'location','Northwest');
% % set(legend1,...
% %     'Position',[0.19 0.19 0.24 0.27]);
%
%
% xlabel(['Iterations (incl burn-in of ' num2str(prs.skip) ')'])
% ylabel(['Acceptance rate in \%'])
%
%
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
%
%


%%


for mm = 1; % choose number
    
    sig = 2.5; % 2.5
    SMmax = W(:,mm) < quantile( W(:,mm), sig/100   , 1 );
    SMmin = W(:,mm) > quantile( W(:,mm), 1-sig/100 , 1 );
    
    SM1 = zeros(size(W(:,mm)));%spat_map;
    SM2 = SM1;
    
    SM1(SMmax) = -1;
    SM2(SMmin) = -2;
    
    m1=zeros(l1*l2*l3,1);
    m1(logical(mas))=SM1;
    m11=reshape(m1,[l1 l2 l3]);
    
    m1=zeros(l1*l2*l3,1);
    m1(logical(mas))=SM2;
    m12=reshape(m1,[l1 l2 l3]);
    
    
    %% MIP
    
    % im = m11;
    %
    % im_mip = (im);
    % im_mip(isnan(im_mip)) = 0;
    
    
    mas_mip = (reshape(mas,l1,l2,l3));
    
    if plotall == 1
        f1 = figure;
        f2 = figure;
        for i = 1:3
            
            mip1 = squeeze(sum(m11,i));
            mip1 = double(mip1<0);
            
            mip2 = squeeze(sum(m12,i));
            mip2 = double(mip2<0);
            
            mip_M = squeeze(sum(mas_mip,i));
            mip_M = double(mip_M>0);
            
            figure(f1)
            subplot(2,2,i)
            MIP = mip_M + mip1;
            
            MIP = imrotate(MIP, 90);
            imagesc(MIP); axis image
            
            figure(f2)
            subplot(2,2,i)
            MIP = mip_M + mip2;
            
            MIP = imrotate(MIP, 90);
            imagesc(MIP); axis image
            
        end
    end
    
    % return
    
    
    %% 3D image
    
    
    im = abs(SM1 + SM2);
    m1=zeros(l1*l2*l3,1);
    m1(logical(mas))=im;
    im=reshape(m1,[l1 l2 l3]);
    
    
    
    
    nrCnum = unique(im);
    nrC = length(nrCnum)-1;
    
    
    im = im.*reshape(mask1,l1,l2,l3);
    
    
    
    
    
    % if plotall == 1
    % f1 = figure;
    % clear legendInfo
    % c = distinguishable_colors(nrC);
    %
    % for i = 1:nrC;
    %
    %     ii = nrCnum(i+1);
    %
    %     clear clu
    %     clu = im==ii;
    %     clu = double(clu);
    %
    %     p2 = patch(isosurface(clu,.5),...
    %         'FaceColor',c(i,:),'EdgeColor','none');
    %     isonormals(clu,p2);
    %
    %     hold on
    %
    %     legendInfo{i} = ['Cluster ' num2str(i)]; % or whatever is appropriate
    %
    % end
    %
    % p2 = patch(isosurface(reshape(mas,l1,l2,l3),.5),...
    %     'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
    % isonormals(reshape(mas,l1,l2,l3),p2)
    % view(-31,4); daspect([1 1 1]); axis tight
    % grid on
    %
    % set(p2,'FaceColor',[0.5 0.5 0.5])
    % alpha(0.3)
    %
    % camlight;  camlight(-80,-10); lighting phong;
    % % title([ForG ' clusters - ' nameeffekt ' hand activation'])
    % axis([0 l2 0 l1 0 l3])
    % hold on
    %
    % legendInfo{i +1} = 'Mask'; % or whatever is appropriate
    %
    % legend(legendInfo)
    %
    
    % end
    
end

if plotall == 1
    slices = 1:2:l3;
    for i = 1:prs.m
        figure('name',['component ' num2str(i)]);
        WWmax = max(W(:,i));
        WWmin = min(W(:,i));
        
        co = 1;
        for l = slices
            subplot(4,ceil(length(slices)/4),co)
            
            m1=WWmin.*ones(l1*l2*l3,1);
            m1(logical(mas))=W(:,i);
            WW=reshape(m1,[l1 l2 l3]);
            
            WW = imrotate(WW, 90);
            imagesc(squeeze(WW(:,:,l))); axis image; axis off;
            caxis([WWmin WWmax])
            co = co+1;
        end
    end
end

%%

spat_map = mean(X,2);

Cmap = colormap(gray(128));
Cmap(1,:) = [1 0 0];
Cmap(end,:) = [0 0 1];

plotall = 0;

for mm = 19
    
    sig = 2.5; % 2.5
    SMmax = W(:,mm) < quantile( W(:,mm), sig/100   , 1 );
    SMmin = W(:,mm) > quantile( W(:,mm), 1-sig/100 , 1 );
    
    SM = spat_map;
    SM(SMmax) = max(spat_map)+2000;
    SM(SMmin) = min(spat_map)-2000;
    
    
    mM=min(spat_map)*ones(l1*l2*l3,1);
    mM(logical(mas))=SM;
    mM=reshape(mM,[l1 l2 l3]);
    
    MAPS = mM;
    
end

%%

if plotall == 1
    %%
    slices =  round(linspace(1,l3-5,20));
    for i = 19 %1:prs.m
        onecomp = figure('name',['Collected component ' num2str(i)]);
        
        
        
        
        
        hh = subplot(3,10,6:10);
        h(i)=histogram(hell(i,:)','Normalization','probability','BinLimits',[0 5]); %BinLimits);
        hold on
%         if i == prs.m
            xlabel('Time scale $[s]$, $\ell*N*TR$')
%         end
        
%         if i == ceil(prs.m/2)+1
%             ylabel('$ \ \ \ \ \ \ \ \ \ \ p(\ell)$')
%         end
title('(b)')
         pp = get(hh,'pos');
    pp(1) = pp(1) + 0.1;
%     pp(2) = pp(2) - 0.02;
%     pp(4) = pp(4) + 0.04;
    pp(3) = pp(3) - 0.1;
    set(hh, 'pos', pp);
        
        
        
        
        
        hh = subplot(3,10,1:5);
        
         pp = get(hh,'pos');
%     pp(1) = pp(1) + 0.06;
%     pp(2) = pp(2) - 0.02;
%     pp(4) = pp(4) + 0.04;
    pp(3) = pp(3) + 0.08;
    set(hh, 'pos', pp);
        
        
        plot(t, Z(i,:), 'b-' )
        hold on
        plot(t, Zu(i,:), 'r:' )
        plot(t, Zl(i,:), 'r:' )
        %     plot(t,pa)
        
       stiml = plot(t,stim_left.*max(Z(i,:)),'k:');
%     hold on
    stimr = plot(t,stim_right.*max(Z(i,:)),'--','color',[0.8 0.8 0.8]);
    
    leg = legend([stiml, stimr],'Left','Right');
   leg.Location = 'eastoutside';
   leg.Box = 'off';
leg.FontName =  'cmu serif';
    
%         plot(t,stim_left,'--')
%         % hold on
%         plot(t,stim_right,'--r')
        hold off
        axis tight
%         if i == prs.m
            xlabel('Time (s)')
%         end
        title('(a)'); %['Component ' num2str(i)])
        
        
%         x =  Z(i,:);
%         Fs = TR;
%         subplot(4,10,11:15)
%         freqspec
        
        
        
        WWmax = max(W(:,i));
        WWmin = min(W(:,i));
        
        co = 11;
        for l = slices
            
            if co == 16
                title('$\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \  $ (c)')
            end 
           hh = subplot(3,10,co);
            if co > 20
                 pp = get(hh,'pos');
%     pp(1) = pp(1) + 0.06;
    pp(2) = pp(2) + 0.15;
%     pp(4) = pp(4) - 0.15;
%     pp(3) = pp(3) - 0.15;
    set(hh, 'pos', pp);
            end
                
                plotmaps = imrotate(fliplr(MAPS(:,:,l)),270);
            imagesc(plotmaps); %
    axis image; axis off;
    colormap(Cmap)
    caxis([min(spat_map)-100 max(spat_map)+100])
    
    
    
%             m1=WWmin.*ones(l1*l2*l3,1);
%             m1(logical(mas))=W(:,i);
%             WW=reshape(m1,[l1 l2 l3]);
%             
% %             WW = imrotate(WW, 270);
%             imagesc(squeeze(WW(:,:,l))); axis image; axis off;
%             caxis([WWmin WWmax])
%             colormap(jet)
            co = co+1;
        end
    end
    %%
end


%%

figure('Name', filename)

ss = 1;

% l1 = 53; l2 = 63; l3 = 46;
% indx = Xall.indx;
%
% mask1 = zeros(l1*l2*l3,1);
% mask1(indx) = 1;
% mas = mask1;
% mask1(mask1 == 0) = nan;

for mm = 1:prs.m
    
    sig = 2.5; % 2.5
    SMmax = W(:,mm) < quantile( W(:,mm), sig/100   , 1 );
    SMmin = W(:,mm) > quantile( W(:,mm), 1-sig/100 , 1 );
    
    SM1 = zeros(size(W(:,mm)));%spat_map;
    SM2 = SM1;
    
    SM1(SMmax) = -1;
    SM2(SMmin) = -2;
    
    % m1=zeros(l1*l2*l3,1);
    % m1(logical(mas))=SM1;
    % m11=reshape(m1,[l1 l2 l3]);
    %
    % m1=zeros(l1*l2*l3,1);
    % m1(logical(mas))=SM2;
    % m12=reshape(m1,[l1 l2 l3]);
    
    
    %
    %     x =  Z(mm,:);
    %     Fs = 2.5;
    %     subplot(prs.m,7,ss+5:ss+6)
    %     freqspec
    %
    %
    %     axis tight
    %     if mm == prs.m
    %         xlabel('Frequency (Hz)')
    %     end
    %
    
    %     subplot(prs.m,7,ss+1:ss+4)
    hh = subplot(prs.m/2+1,2,mm);
    
    pp = get(hh,'pos');
    pp(4) = pp(4) + 0.04;
    pp(3) = pp(3) + 0.05;
    set(hh, 'pos', pp);
    plot(t, Z(mm,:), 'b-' )
    hold on
    plot(t, Zu(mm,:), 'r:' )
    plot(t, Zl(mm,:), 'r:' )
    %     plot(t,pa)
    hold off
    axis tight
    %     if mm == prs.m
    %         xlabel('Time (s)')
    %     end
    
    %     axis off
    %     subplot(prs.m,7,ss)
    %     plot3Dbrain
    %     zoom(1.5)
    %     ss = ss+7;
    
    %     set(gca,'xtick',[],'ytick',[],'box','off')
    %     ylabel(num2str(round(PartVar(mm))))
    axis off
    text(-0.15,0.5,{['Comp. ' num2str(mm) ] ;...
        [' ExpEnergy: ' num2str(ExpEng(mm),2)] ; ...
        [' CumEnergy: ' num2str(ex(mm),2)]},'units','normalized')
    
end


%
%     stimuli_l = zeros(1,N);
%      stimuli_r = zeros(1,N);
%
%     s_left =  [12.0000
%   36.0000
%   60.0000
%   84.0000
%  108.0000
%  132.0000
%  156.0000
%  180.0000
%  204.0000
%  228.0000];
%
% s_right = [ 24.0000
%   48.0000
%   72.0000
%   96.0000
%  120.0000
%  144.0000
%  168.0000
%  192.0000
%  216.0000
%  240.0000];
%
% for i = 1:length(s_left)
%  stimuli_l(s_left(i):s_left(i)+6) = 1;
%  stimuli_r(s_right(i):s_right(i)+6) = 1;
% end
%
%
% hh =subplot(prs.m/2+1,2,mm+1);
% pp = get(hh,'pos');
% % pp(4) = pp(4) + 0.04;
% pp(3) = pp(3) + 0.05;
% set(hh, 'pos', pp);
%
% plot(t,stim_left,'--')
%     hold on
%     plot(t,stim_right,'--r')
%     axis off
%
%
%
% % figure; plot(t,stimuli_l)
% % hold all
% % plot(t,stim_left,'r')
%
%   hh = subplot(prs.m/2+1,2,mm+2);
%     pp = get(hh,'pos');
% % pp(4) = pp(4) + 0.04;
% pp(3) = pp(3) + 0.05;
% set(hh, 'pos', pp);
% plot(t,stim_left,'--')
%     hold on
%     plot(t,stim_right,'--r')
%
axis off


%%

load('KasperStim')

mainfig = figure('Name',filename);
ss = 1;

for mm = 1:prs.m
    
    sig = 2.5; % 2.5
    SMmax = W(:,mm) < quantile( W(:,mm), sig/100   , 1 );
    SMmin = W(:,mm) > quantile( W(:,mm), 1-sig/100 , 1 );
    
    SM1 = zeros(size(W(:,mm)));%spat_map;
    SM2 = SM1;
    
    SM1(SMmax) = -1;
    SM2(SMmin) = -2;
    
    
    hh1 = subplot(prs.m/2+1,4,ss);
    
    plot3Dbrain
    zoom(1.5)
    
    if mm ==1
        axMAP = gca;
        text(axMAP.XLim(1)-axMAP.XLim(2)*1.25, axMAP.YLim(2)*1.25, axMAP.ZLim(2)/2, ...
            {['\bf{1: E ' num2str(ExpEng(mm),3) ' \%}']}); %,...
        %         keyboard
    else
        axMAP = gca;
        text(axMAP.XLim(1)-axMAP.XLim(2)*1.25, axMAP.YLim(2)*1.25, axMAP.ZLim(2)/2,...
            {['\bf{' num2str(mm) ': E ' num2str(ExpEng(mm),3) ' \%}'],['C ' num2str(ex(mm)*100,3) ...
            ' \%'] }); %,...
        %      'HorizontalAlignment','right')
    end
    
    
    hh = subplot(prs.m/2+1,4,ss+1);
    
    pp = get(hh,'pos');
    pp(1) = pp(1) - 0.06;
    pp(2) = pp(2) - 0.02;
    pp(4) = pp(4) + 0.04;
    pp(3) = pp(3) + 0.07;
    set(hh, 'pos', pp);
    plot(t, Z(mm,:), 'b-' )
    hold on
    plot(t, Zu(mm,:), 'r:' )
    plot(t, Zl(mm,:), 'r:' )
    %     plot(t,pa)
    hold off
    axis tight
    
    ss = ss+2;
    
    
    axis off
    
end


hh =subplot(prs.m/2+1,4,ss+1);
pp = get(hh,'pos');
pp(1) = pp(1) - 0.06;
pp(2) = pp(2) - 0.06;
pp(4) = pp(4) + 0.04;
pp(3) = pp(3) + 0.07;
set(hh, 'pos', pp);

plot(t,stim_left,'k:')
    hold on
    plot(t,stim_right,'--','color',[0.8 0.8 0.8])
    axis off
    ylim([-0.5 1.5])
    
    
   leg = legend('Left','Right');
   leg.Position = [0.18 0.0708 0.0637 0.0572];
   leg.Box = 'off';
   leg.FontName =  'cmu serif';
%     ,'Location','westoutside','box','off')
    
 hh =subplot(prs.m/2+1,4,ss+3);
pp = get(hh,'pos');
pp(1) = pp(1) - 0.06;
pp(2) = pp(2) - 0.06;
pp(4) = pp(4) + 0.04;
pp(3) = pp(3) + 0.07;
set(hh, 'pos', pp);

plot(t,stim_left,'k:')
    hold on
    plot(t,stim_right,'--','color',[0.8 0.8 0.8])
    axis off
       ylim([-0.5 1.5])
       
   leg = legend('Left','Right');
   leg.Position = [0.6 0.0708 0.0637 0.0572];
   leg.Box = 'off';
leg.FontName =  'cmu serif';


%
% prompt = 'Save images (and overwrite old)? - then press 1 : ';
% savefigure = input(prompt);
% % y = x*10
%
savefigure =0;
if savefigure == 1
    export_fig(mainfig, 'mainfig_3D','-pdf')
    export_fig(logfig, 'logfig_3D','-pdf')
     export_fig(onecomp, 'comp_3D','-pdf')
end



%%
shape = prs.es(1);
scale = 1/prs.er(1);
%scale = 200

tnew = (0:1/(N*1000):1)'; %N*dt-dt)';

ell = bighellC(:,1);  %gamrnd( prs.es, 1./prs.er)';

PPplot = figure;
tsek = tnew.*(TR*N);
alphaa = TR*N;
ynew = (alphaa).*pdf('Gamma',tsek,shape,scale*alphaa);
plot(tsek,ynew)
hold all

for mm= 1: length(ell)
    ellPsek(:,mm)  = alphaa.*pdfrG(tsek,ell(mm)*alphaa,prs.ess*alphaa);
    
    plot(tsek,ellPsek(:,mm),'color',[1 0 0 ])
end
xlim([0 0.05*alphaa])
% ylim([0 1000])


xlabel('$\ell$ [s]')
ylabel('Probability')
title('Prior for $\ell$ and the initial proposal distributions')
% legend(['$p( \ell_p ) \sim  \Gamma ( $' num2str(prs.es(1)) ', 1/' num2str(prs.er(1)) '$ )$'],...
%     ['$q(\ell_p) \sim  \mathcal{N}_{Ref} (\ell,$' num2str(prs.ess) ' )'])
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

ss =prs.ss; 
sr = prs.sr;
k = ss;

invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
plot(invgam(tnew,k,sr),tnew);
% plot(invgam(1/tnew,k,sr),tnew);



%%
tnew = linspace(0,1.2,1000);
ss =1 %prs.ss; 
sr = 0.5 %prs.sr;
k = ss;
sub1 = 1; sub2 = 3;
invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
figure;subplot(sub1,sub2,1);
plot(1./invgam(tnew,k,sr),tnew);
subplot(sub1,sub2,2);
plot(invgam(tnew,k,sr),tnew);
% plot(invgam(1./tnew,k,sr),tnew);
subplot(sub1,sub2,3);
ynew = pdf('Gamma',1./tnew,prs.ss,1/prs.sr);
plot(ynew,tnew)

%% Sigmaone(:,end); %

SigmaMedian = median(Sigmaone,2);
%%
figure; subplot(4,1,1); 
% for i = 1:size(Sigmaone,2)
% SigmaHistOne = histogram(Sigmaone(:,i),100,'BinLimits',[0 2]);
% SigmaHist(:,i) = SigmaHistOne.Values';
% end
% surf(SigmaHist)
tmax =1; % max(SigmaMedian);
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



e = X_0mean-WZ;

% histogram(abs(e),'Normalization','pdf','BinLimits',[0 tmax])
% axis tight 
% 
% %%
% figure
% subplot(2,1,1)
tnew = linspace(0,tmax,1000);
sstest =  8.5714%90 %100;
srtest =  9.5714%60 %30;

subplot 413
histogram(1./gamrnd( sstest, 1/srtest, [ d, 1 ] ),'Normalization','pdf','BinLimits',[0 tmax])
axis tight 


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
%  histogram(val(i),100,'Normalization','pdf','BinLimits',[0 tmax])
rmsval = rms(rms(e))

histogram(SigmaNew,100,'Normalization','pdf','BinLimits',[0 tmax]);
axis tight 

%%
% % figure; subplot(2,1,1); 
% % % SigmaMedian = median(Sigmaone,2);
% % histogram(SigmaMedian,100);
% % 
% % subplot(2,1,2)
% % % tmax = max(SigmaMedian);
% % tnew = linspace(0,tmax,1000);
% % 
% % ss =20; %prs.ss; 
% % sr = 10; %prs.sr;
% % k = ss;
% % % beta = 1/sr;
% % % 
% % % theta = 1/beta;
% % 
% % invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
% % plot(tnew, invgam(tnew,k,sr));
% % 

%%
% % figure;
% % tnew = linspace(0,tmax,1000);
% % 
% % ss =20; %1/N; %prs.ss; 
% % sr = 10; %1*0.8; %prs.sr;
% % 
% % invgam = @(x,a,b) b^a/gamma(a).*(1./x).^(a+1).*exp(-b./x);
% % 
% % % for ss =0.0000001
% % plot(tnew, invgam(tnew,ss,sr));
% % hold all
% % % end

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
% ynew2 = (pdf('Gamma',tnew,prs.ss(1),1/prs.sr(1)));
% 
% plot(ynew2,tnew)
% ylim([axinfo.YLim])
% %
%%


    alphaall = reshape(bigalpha, d ,  prs.m);
     alphaall = alphaall(:,ordr5(K:-1:1));


figure;
plot_quantiles = [50 5 95];
mcor = distinguishable_colors(prs.m); %'rgbcmk';
linsty{1}= '-'; linsty{2}='--'; linsty{3}=':';
for mm = 1: prs.m
    for pq = 1:length(plot_quantiles)
        alpha_quan=  quantile( squeeze(alphaall(:,mm,:)), plot_quantiles(pq)/100, 1 );
        
        subplot(1,prs.m,1:prs.m/2)
        plothandle(mm,pq) = plot([alpha_quan alpha_quan],'color',mcor(mm,:), 'LineStyle' , (linsty{pq}));
        hold all
    end
end

subplot(1,prs.m,1:prs.m/2)
axinfo = gca;
axinfo.XLim = [1,axinfo.XLim(2)];
axinfo.YLim = [0,axinfo.YLim(2)];

xlabel('End sample')
axinfo.XTick = [];
ylabel('$\alpha$')
title('(a)')

subplot(1,prs.m,prs.m/2+1:prs.m)
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

%%
figure('position',[0.1 0.1 0.9 0.5])
for mm = 1:prs.m
hh = subplot(2,prs.m/2,mm);
    pp = get(hh,'pos');
    pp(2) = pp(2) + 0.04;
    set(hh, 'pos', pp);
    
        alphapic = alphaall(:,mm);
    
    sig = 5; % 2.5
    SMmax = alphapic < quantile( alphapic, sig/100   , 1 );
    %SMmin = alphapic > quantile( alphapic, 1-sig/100 , 1 );
    
    SM1 = zeros(size(alphapic));%spat_map;
    SM2 = SM1;
    
    SM1(SMmax) = -1;
    SM2(1) = 0;
    
    plot3Dbrain
   zoom(1.5)
    %
    %  imagesc(mM); axis image; axis off;
    titlehandle = title(['GP ' num2str(mm)]);
    titlehandle.FontSize = 12;
    % imghandel = gca;
    % if mm == 1
    % text(imghandel.XLim(1),imghandel.YLim(2), '(c) $\alpha$ median')
    % end
end


%%
