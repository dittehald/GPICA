
global h_main ica

 %% For GUI   
if ~isfield(ica,'FileNames') || isempty(ica.FileNames);
    errordlg('No files have been selected')
    return 
end 

      %% data    
       status = 'Loading data...';
    set(h_main.settings.status.text, 'String', status), drawnow
    
        rdata = load(ica.FileNames);
        
        
        
        
        
   status = 'Done loading data';
    set(h_main.settings.status.text, 'String', status), drawnow
    

	

% TR = 333e-3;

TR = eval(ica.TR{1});
L = eval(ica.Lsize{1});

dataformat = length(L);

if dataformat == 2  % 2D data
    l1 = L(1);
    l2 = L(2);
    % l1 = 82; l2 = 68;
    if isfield(rdata, 'mask1')
        Mask3711 = rdata.mask1;
        X = rdata.X1;
   
    else
         Mask3711 = rdata.Mask3711;
         X = rdata.X3711';
%          X = X(:,1:363);
    end
else                % 3D data
    l1 = L(1);
    l2 = L(2);
    l3 = L(3);
    %l1 = 53; l2 = 63; l3 = 46; TR= 2.49; [53, 63, 46]
     
    X = rdata.X';
    Indx = rdata.indx;

% mask1 = zeros(l1*l2*l3,1);
% mask1(indx) = 1;
% mas = mask1;
% mask1(mask1 == 0) = nan;
%    Mask3711 = mask1; 
end

 [ d ,N ] = size( X );

Fs = 1/TR;                      % samples per second
   dt = 1/Fs;                     % seconds per sample
   t = (0:dt:N*dt-dt)';
   
     %% parameters    
    
       status = 'Performing initialization...';
    set(h_main.settings.status.text, 'String', status), drawnow
    

    prs.chains = 1;
    prs.stride = 1;
%     prs.std = 0;        % standardize with z-score!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
prs.std = str2double(get(h_main.settings.ica.algorithm.std,'String'));


    prs.eu = 1;         % 1 if you want to optimize the parameters of K
prs.jit = 1e-4;     % add to the diagonal of K for numerical stability
prs.tpts = (1:N)/N; %t; % time points, it is normalized time right now
prs.drawfollow = 0;
prs.TR = TR;
    
    prs.skip = str2double(get(h_main.settings.ica.algorithm.skip,'String'));
    prs.nsamples = str2double(get(h_main.settings.ica.algorithm.nsamples,'String'));
    prs.m = str2double(get(h_main.settings.ica.algorithm.m,'String'));
    
    prs.ts = str2double(get(h_main.settings.ica.algorithm.ts,'String'));
    prs.tr = str2double(get(h_main.settings.ica.algorithm.tr,'String'));
   
    
    prs.ss = str2double(get(h_main.settings.ica.algorithm.ss,'String'));
    %prs.ss = prs.ss/N;
    
    prs.sr = str2double(get(h_main.settings.ica.algorithm.sr,'String'));
    
%     prs.es = str2double(get(h_main.settings.ica.algorithm.es,'String'));
%     prs.er = str2double(get(h_main.settings.ica.algorithm.er,'String'));
prs.emean = str2double(get(h_main.settings.ica.algorithm.emean,'String'));

prs.es = 1; 
prs.emean = str2double(get(h_main.settings.ica.algorithm.emean,'String')); %TR*1.2;

prs.er = prs.es/ (prs.emean/(N*TR));
% prs.er = 1000

prs.ess = str2double(get(h_main.settings.ica.algorithm.ess,'String'));
   
    if length(prs.er) ~= prs.m 

    prs.es = [prs.es repmat(prs.es(end),1,prs.m-length(prs.es))];
    prs.er = [prs.er repmat(prs.er(end),1,prs.m-length(prs.er))];
%     prs.ess = [prs.ess repmat(prs.ess(end),1,prs.m-length(prs.ess))]';
    end
 
filename = ['GUIRUN_m' num2str(prs.m) '_' num2str(N) '_' num2str(prs.nsamples) 'burn' num2str(prs.skip) '.mat'];

if exist(filename,'file') == 2
   filename = [filename(1:end-4) datestr(datetime,'ddmmyyHHMM') '.mat'];
end

prs.filename = filename;
    
  
    
    %%
    if h_main.settings.ica.checkconv.Value == 1% CONV
        prs.lags = str2double(get(h_main.settings.ica.algorithm.convlags,'String'));
        status = 'Running model...';
        set(h_main.settings.status.text, 'String', status), drawnow
        
        
        % run model
        res = GPICA4fmriCONV( X, prs );
    else
        status = 'Running model...';
        set(h_main.settings.status.text, 'String', status), drawnow
        
        
        % run model
        res = GPICA4fmri( X, prs ); %remember
    end
    %% 
    
%     prs = res.prs;
    
    status = 'Done';
    set(h_main.settings.status.text, 'String', status), drawnow
    
    ica.res = res;
    ica.prs = prs;
    
    save(filename,'ica','prs','-v7.3');
    

    if h_main.settings.ica.checkconv.Value == 1% CONV
        if dataformat == 2
            plotconvresultsTOOL % 
        else
            plot3dresultsTOOL % not working?
        end
    else
         if dataformat == 2
            plotonedresultsTOOL
        else
            plot3dresultsTOOL
        end
    end
    


%                 
%                 if  k == 1 && mod( 10*l, prs.skip ) == 0
%                     status = ['Running model... (burn in: ' num2str(procc*10) '%)'];
%                     set(h_main.settings.status.text, 'String', status), drawnow
%                     fprintf( ':' )
%                     procc = procc + 1;
%                     
%                 elseif l == 1 && mod( 10*k, prs.nsamples ) == 0
%                     procc1 = procc1 + 1;
%                     status = ['Running model... ' num2str(procc1*10) '%'];
%                     set(h_main.settings.status.text, 'String', status), drawnow
%                     fprintf( '-' )
%                 end
         
