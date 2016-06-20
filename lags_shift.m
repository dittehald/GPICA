function Zny = lags_shift(Z,lags)

[m,N] = size(Z);
Znul = [Z 0.5.*randn(m, lags)];

for i = 1:lags+1;
    Zs(:,:,i) = circshift(Znul,[0 i-1]);
    
end

Zny(:,:,:) = Zs(1:m,1:N,:);
end