

logfig =figure('Position',[0 0 0.95 0.5]');

hh = subplot(1,3,1);
ellcurves = plot(hellCone');
yl=get(gca,'ylim');
title({'Developement of $\ell$','    '})

xlabel(['Iterations (incl warm-up of ' num2str(prs.skip) ')'])
ylabel(['$\ell$'])

hold all
ax1 = gca; % current axes
areacolor = [0.8 0.8 0.8];
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'YAxisLocation','right',...
    'XAxisLocation','top',...
    'Color','none');




% ax2.XTick = [];
% ax2.XTickLabel = 'hej';
% ax2.TickDir = 'both'
% ax2.YTick = [];

ax2.XColor = areacolor;
ax2.YColor = areacolor;
ax2.XDir = 'reverse';
tnew = (linspace(0,yl(2),1000));
ynew2 = (pdf('Gamma',tnew,prs.es(1),1/prs.er(1)));
%  plot(ynew2,tnew)
hold all
% figure
a = line(ynew2,tnew,'Parent',ax2,'color',areacolor,'linewidth',4);
te = text(ax2.XLim(2)/2,ax2.YLim(2),'Prior $\ell$ propability');
te.HorizontalAlignment = 'center';
te.VerticalAlignment = 'top';
te.FontSize = 18;
te.Color =  [0.6 0.6 0.6];
 drawnow
a.Parent.Position = ax1.Position;

% return

% set (gca, 'Yscale', 'log');

% a.FaceColor = areacolor;
% a.EdgeColor = areacolor;

% 
% ax1_pos = ax1.Position; % position of first axes
% ax3 = axes('Position',ax1_pos,...
%     'Color','none');
% plot(hellCone','Parent',ax3)



%  uistack(ellcurves, 'top')
% set(ax1,'children',flipud(get(ax1,'children')))
% plot(hellCone','Parent',ax1)

% a.FaceAlpha = 0.05;
% alpha(a,.5)
%  set(hh2, 'ylim', yl);

% figure
% aaa = area(1:10);
%  aaa.FaceAlpha = 0.05;

%%
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
 
hh =subplot(1,3,2);
% pp = get(hh,'pos');
% pp(3) = pp(3) + 0.05;
% pp(1) = pp(1) - 0.05;
% set(hh, 'pos', pp);

accrates = round(haccrone./(repmat(1:(prs.skip+prs.nsamples-1),prs.m,1)).*100);

plot(accrates')
title({'Acceptance rate developement','  '})

accrate = round(haccrone(:,end)./(prs.skip+prs.nsamples-1).*100);

xlabel(['Iterations (incl warm-up of ' num2str(prs.skip) ')'])
ylabel(['Acceptance rate in \%'])


            subplot(1,3,3)
            plot(hloglikeone); axis tight;
            title({'Loglikelihood','  '})
            xlabel(['Iterations (incl warm-up of ' num2str(prs.skip) ')'])
            ylabel('$p(z_{[p]}|\ell_p)$')
            xlim([-10    prs.nsamples+prs.skip-1])
            
            if size(hloglikeone,1) < 200
                ylim([min(min(hloglikeone(:)))-100 max(max(hloglikeone(:)))+100])
            
            else
            ylim([min(min(hloglikeone(200:end)))-100 max(max(hloglikeone(:)))+100])
            end
            %          return        
        if  isfield(res,'hloglikelihood')
            ylabel('log ( $p(X|WZ,\Psi)$ )')
        end