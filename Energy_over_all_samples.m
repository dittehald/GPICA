%% Energy over all samples

[d N] = size(X_0mean);
% Xback = zeros(d,N,prs.m);
EnergyX = norm(X_0mean,'fro')^2;



Residual3 = zeros(prs.m,1);
pvaf = zeros(prs.m,1);

tic
for m=1:prs.m
    Wsq = squeeze(Wall(:,m,:));
    Zsq = squeeze(Zall(m,:,:));
    WZ = (Wsq * Zsq');
    Residual3(m) = EnergyX -  2*sum(sum(X_0mean .* WZ)) ./ prs.nsamples + ...
        (sum( (Wsq).^2 )) * (sum( (Zsq).^2 ))'./ prs.nsamples;
    pvaf(m) = 100-100*mean(var(X_0mean - WZ./prs.nsamples))/mean(var(X_0mean));
end
time_residual = toc;

ExpEng = 100-100*Residual3./EnergyX;

         
%% Cummulative


[vals5,ordr5] = sort(ExpEng);

K = prs.m;

Word = Wall(:,ordr5(K:-1:1),:); 
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
          e = X_0mean - squeeze(Word(:,1:m,r))*squeeze(Zord(1:m,:,r));
          Residual(m) = Residual(m) + norm(e,'fro')^2;
     end
 chain_time(m) = toc;
 disp('time for one source')
 disp(chain_time(m))
     ex(m) = 1 -  Residual(m) / prs.nsamples / EnergyX; 
     if m==1
          dex(1) = ex(1);
     else
          dex(m) = ex(m)-ex(m-1);
    end
end
% ex = ex.*100;



Zone = reshape(Zord, N*prs.m, prs.nsamples );
Wone = reshape(Word, d*prs.m, prs.nsamples );

Z = reshape( median( Zone, 2 ), prs.m, N );

% Errorbars
Zu = reshape( quantile( Zone, 0.95, 2 ), prs.m, N );
Zl = reshape( quantile( Zone, 0.05, 2 ), prs.m, N );



if  d>4000
    Wr1 =  median( Wone( 1:d*prs.m/4,:), 2 );
    Wr2 = median( Wone( d*prs.m/4+1:2*d*prs.m/4,:), 2 );
    Wr3 = median( Wone( 2*d*prs.m/4+1:3*d*prs.m/4,:), 2 );
    Wr4 = median( Wone( 3*d*prs.m/4+1:end,:), 2 );
    W = [Wr1; Wr2; Wr3; Wr4];
    W = reshape( W, d, prs.m);
else
    W = reshape( median( Wone, 2 ), d, prs.m);
end


% W = reshape( median( Wone, 2 ), d, prs.m);


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
