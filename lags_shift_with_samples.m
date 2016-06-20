 function Zny = lags_shift_with_samples(Z,lags)

if lags == 0
    Zny = Z;
else
[m,N,s] = size(Z);

    perZ = permute(Z,[1 3 2] );
    newsamples = 0.5.*randn(m, s,lags);
    Zekstra = cat(3,perZ, newsamples);
    Zekstra = permute(Zekstra,[1 3 2] );
    Zs = circshift(Zekstra,[ 0 lags]);
    Zny(:,:,:) = Zs(1:m,1:N,:);
end
end
