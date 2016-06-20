
Wall_withlags = reshape(Wone, d ,  prs.m ,(prs.lags +1), prs.nsamples);
Zall = reshape(Zone, prs.m,  N  , prs.nsamples);


%%
Xbackpro = zeros(d,N);
for la = 1:prs.lags
    
    Zlag =  lags_shift_with_samples(Zall,la-1);
    Xbackpro =Xbackpro+ Wall_withlags(:,:,la,end)*Zlag(:,:,end);
    
    
end

% Zall_withlags = lags_shift(Zone,prs.lags);  % Update all lags

% Wall = squeeze(sum(Wall_withlags,3));
% Zall = reshape(Zone, prs.m,  N  , prs.nsamples);

%%
% Xlag_proj = zeros(d,N,prs.m);

[d N] = size(X_0mean);
Xback = zeros(d,N,prs.m);
EnergyX = norm(X_0mean,'fro')^2;
Residual = zeros(prs.m,1);
pvaf = zeros(prs.m,1);

%% Slow all-sample calculations
Residual = zeros(prs.m,1);

tic
for m=1:prs.m
    tic
    for r=1:prs.nsamples
        e = X_0mean;
    for la=1:prs.lags+1
        
        Zlag =  lags_shift_with_samples(Zall(:,:,r),la-1);
        Wlag =  squeeze(Wall_withlags(:,:,la,:));
    
    
         e = e - squeeze(Wlag(:,m,r))*squeeze(Zlag(m,:));
    end
    Residual(m) = Residual(m) + norm(e,'fro')^2;
%          counter(m) = counter(m)+1;
   
    
    end
    
    chain_time(m) = toc;
   
    
    disp(chain_time(m)) % = 1s pr prs.nsamples
end
time_residual = toc; % One source = 16 min (1s pr sample iteration)

Residual2 = Residual / (prs.nsamples);
ExpEng = 100-100*Residual2./EnergyX

%% Cummulative

[vals5,ordr5] = sort(ExpEng);

K = prs.m;

Word = Wall_withlags(:,ordr5(K:-1:1),:,:); 
Zord = Zall(ordr5(K:-1:1),:,:);
ExpEng = vals5(K:-1:1);

% Calc. cummulated energy
dex = zeros(prs.m,1);
ex = zeros(prs.m,1);
[d, n] = size(X_0mean);
WZ = zeros(d,n);

Residual = zeros(prs.m,1);
for m=1:prs.m
    tic
    
    for r=1:prs.nsamples
    e = X_0mean;
    for la=1:prs.lags+1
     
    Zlag =  lags_shift_with_samples(Zord(:,:,r),la-1);
        Wlag =  squeeze(Word(:,:,la,:));
    
     
          e = e - squeeze(Wlag(:,1:m,r))*squeeze(Zlag(1:m,:));
          
     end
     Residual(m) = Residual(m) + norm(e,'fro')^2;
    end
    
 chain_time(m) = toc;
 disp('time for one source')
 disp(chain_time(m))
     ex(m) = 1 -  Residual(m) / (prs.nsamples) / EnergyX; 
     if m==1
          dex(1) = ex(1);
     else
          dex(m) = ex(m)-ex(m-1);
    end
end

% EnergyX = norm(X_0mean,'fro')^2;
% 
% for m = 1:prs.m
%     for la=1:prs.lags+1
%         
%         Zlag =  lags_shift_with_samples(Zall,la-1);
%         Wlag =  squeeze(Wall_withlags(:,:,la,:));
%         
%         
%         Xlag_proj(:,:,m) = Xlag_proj(:,:,m)+ squeeze(Wlag(:,m,:))*squeeze(Zlag(m,:,:))'./prs.nsamples;
%         %
%     end
%     Residual3(m) = norm(X_0mean-Xlag_proj(:,:,m),'fro')^2; % all samples
%     
% end
% 
% [vals,ordr] = sort(Residual3);
% 
% ExpEng3 = 100-100*Residual3./EnergyX;
% [vals5,ordr5] = sort(ExpEng3);
% % 
% %%
% K = prs.m;
% 
% Word = Wall_withlags(:,ordr5(K:-1:1),:,:); % same as ordr
% Zord = Zall(ordr5(K:-1:1),:,:);
% ExpEng = vals5(K:-1:1);
% 
% dex = zeros(prs.m,1);
% ex = zeros(prs.m,1);
% [d, n] = size(X_0mean);
% WZ = zeros(d,n);
% Xlag_proj = zeros(d,N,prs.m);
% 
% for m=1:prs.m
%     
%     for la=1:prs.lags+1
%         
%         Zlag =  lags_shift_with_samples(Zord,la-1);
%         Wlag =  squeeze(Word(:,:,la,:));
%         
%         
%         Xlag_proj(:,:,m) = Xlag_proj(:,:,m)+ squeeze(Wlag(:,m,:))*squeeze(Zlag(m,:,:))'./prs.nsamples;
%         %
%     end
%     
%     WZ = WZ + Xlag_proj(:,:,m);
%     
%     ex(m) = 1 -   ( norm(X_0mean-WZ,'fro')^2 / EnergyX ); % note it is the variance explained by sources 1:m
%     if m==1
%         dex(1) = ex(1);
%     else
%         dex(m) = ex(m)-ex(m-1);
%     end
% end

%%


Zone = reshape(Zord, N*prs.m, prs.nsamples );

Word = Wall_withlags(:, ordr5(K:-1:1),:,:);
Wone = reshape(Word, d*prs.m*(prs.lags +1), prs.nsamples );


Z = reshape( median( Zone, 2 ), prs.m, N );


% Errorbars
Zu = reshape( quantile( Zone, 0.95, 2 ), prs.m, N );
Zl = reshape( quantile( Zone, 0.05, 2 ), prs.m, N );

clear Zone

partofW =  floor(d* prs.m* ( prs.lags +1)/4);
if  prs.nsamples>2000
    Wr1 =  median( Wone( 1:partofW,:), 2 );
    Wr2 = median( Wone( partofW+1:2*partofW,:), 2 );
    Wr3 = median( Wone( 2*partofW+1:3*partofW,:), 2 );
    Wr4 = median( Wone( 3*partofW+1:end,:), 2 );
    W = [Wr1; Wr2; Wr3; Wr4];
    W = reshape( W, d, prs.m , prs.lags +1);
else
    W = reshape( median( Wone, 2 ), d, prs.m, prs.lags +1);
end


indx = ordr5;
hellone = hellone(indx(K:-1:1),:);
hellCone = hellCone(indx(K:-1:1),:);
haccrone = haccrone(indx(K:-1:1),:);

if length(prs.es) ~= 1
    es = prs.es(indx(K:-1:1));
    er = prs.er(indx(K:-1:1));
else
    es = prs.es;
    er = prs.er;
end