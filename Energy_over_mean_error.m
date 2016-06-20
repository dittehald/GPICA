%% Energy over mean error

Zone = reshape(Zall, N*prs.m, prs.nsamples );
Wone = reshape(Wall, d*prs.m, prs.nsamples );

Z = reshape( median( Zone, 2 ), prs.m, N );

% Errorbars
Zu = reshape( quantile( Zone, 0.95, 2 ), prs.m, N );
Zl = reshape( quantile( Zone, 0.05, 2 ), prs.m, N );

W = reshape( median( Wone, 2 ), d, prs.m);


%% energy
[d , N] = size(X_0mean);
Xback = zeros(d,N,prs.m);
EnergyX = norm(X_0mean,'fro')^2;

for m=1:prs.m
    Xback(:,:,m) = (W(:,m)*Z(m,:));
    Residual(m) = norm(X_0mean-Xback(:,:,m),'fro')^2; % all samples
end

[vals,ordr] = sort(Residual);

ExpEng = 100-100*Residual./EnergyX;
[vals5,ordr5] = sort(ExpEng);


K = prs.m;

W = W(:,ordr5(K:-1:1)); % same as ordr
Z = Z(ordr5(K:-1:1),:);
ExpEng = vals5(K:-1:1);
Zl = Zl(ordr5(K:-1:1),:);
Zu = Zu(ordr5(K:-1:1),:);

dex = zeros(prs.m,1);
ex = zeros(prs.m,1);

WZ = zeros(d,N);

for m=1:prs.m
    WZ = WZ + (W(:,m)*Z(m,:));
    
    ex(m) = 1 -   ( norm(X_0mean-WZ,'fro')^2 / EnergyX ); % note it is the variance explained by sources 1:m
    if m==1
        dex(1) = ex(1);
    else
        dex(m) = ex(m)-ex(m-1);
    end
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