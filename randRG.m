function ell = randRG(mu,sigma,sizeD ) 

if nargin == 2
    sizeD = [1 1];
end



ell =  sigma.*randn(sizeD)+mu;
ell = ell.*sign(ell);

while any(ell(:) == 0)
    ell(ell == 0) = sigma.*randn([sum(sum(ell == 0)) 1])+mu;
    ell = ell.*sign(ell);
end

end
