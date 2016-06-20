

figu =  findall(0,'type','figure');
fignames = get(figu,'name');

if sum(strcmp(fignames,'ell_running'))
    figure(ellfig)
    clf
else
    ellfig = figure('Name','ell_running');    
end



            
            ELLPLOT = reshape( hellC, prs.m, count, prs.chains);
            
            ellcurves = plot(ELLPLOT');
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
line(ynew2,tnew,'Parent',ax2,'color',areacolor,'linewidth',4);
te = text(ax2.XLim(2)/2,ax2.YLim(2),'Prior $\ell$ propability');
te.HorizontalAlignment = 'center';
te.VerticalAlignment = 'top';
te.FontSize = 18;
te.Color =  [0.6 0.6 0.6];
            