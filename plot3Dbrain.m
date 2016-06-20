

% SM1, SM2, mas,mask1

im = abs(SM1 + SM2);
m1=zeros(l1*l2*l3,1);
m1(logical(mas))=im;
im=reshape(m1,[l1 l2 l3]);




nrCnum = unique(im);
nrC = length(nrCnum)-1;


im = im.*reshape(mask1,l1,l2,l3);


% f1 = figure;
clear legendInfo 




c = distinguishable_colors(nrC);

for i = 1:nrC;
    
    ii = nrCnum(i+1);
    
    clear clu
    clu = im==ii;
    clu = double(clu);
    
    p2 = patch(isosurface(clu,.5),...
        'FaceColor',c(i,:),'EdgeColor','none');
    isonormals(clu,p2);
    
    hold on
    
    legendInfo{i} = ['Cluster ' num2str(i)]; % or whatever is appropriate
    
end

p2 = patch(isosurface(reshape(mas,l1,l2,l3),.5),...
    'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
isonormals(reshape(mas,l1,l2,l3),p2)
view(-31,4); daspect([1 1 1]); axis tight
grid on

set(p2,'FaceColor',[0.5 0.5 0.5])
alpha(0.3)

camlight;  camlight(-80,-10); lighting phong;
% title([ForG ' clusters - ' nameeffekt ' hand activation'])
axis([0 l2 0 l1 0 l3])
hold on

legendInfo{i +1} = 'Mask'; % or whatever is appropriate

% legend(legendInfo)
axis off



