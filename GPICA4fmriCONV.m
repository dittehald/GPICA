function [ result ] = GPICA4fmriCONV( X, prs )
%% MCMC for fMRI
% MCMC algorithms work by taking random samples from a given
% distribution. A Bayesian framework is used to determine from
% which distribution to sample from.
%
% Assumes the following model for the data: X = W*Z + epsilon (noise)
% Where W is the mixing matrix (the spatial dispersion) and Z is the
% corresponding timeseries. IID is assumed on W, and ID is assumed on Z.
% The noise is assumed to be normal distributed with a covariance decribed
% by Sigma.
%
% Input:        X       - data set
%               prs     - structure containing all inputs - see
%                         runMCMC4fmri
% Output        results - strucure containing the mixing matrix W, the time
%                         series Z, Sigma, upsilon, the 'accuracy' accr and
%                         the elapsed time etime
%
% -------------------------------------------------------------------------
% DATE       :  20.05.2010
% AUTHORS    :  Ricardo Henoa
%               Simon Grützmeier, Ditte Hald, Stine Harder and Christian
%               Blücher
% COURSE     :  02459 - Machine Learning for Signal Processing, DTU
% SUPERVISOR :  Ricardo Henoa and Ole Winther IMM DTU
%
% =========================================================================


% sizes
m = prs.m;
[ d ,N ] = size( X );

% standarize
if prs.std
    [X,mu_X,sig_X] = zscore( X, 0, 2 );
    
    prs.tr = prs.tr/mean(sig_X)^2;
    prs.sr = prs.sr/mean(sig_X)^2;
    
else
    X = bsxfun( @minus, X , mean( X, 2 ) );
end



itime = cputime;

% array allocation
hlen = prs.nsamples*prs.chains;
% hlenall = (prs.nsamples+prs.skip-1)*prs.chains;
hSigma = zeros( d, hlen );
hW = zeros( d*m*(prs.lags+1), hlen);
% hupsilon = zeros( m, hlen );
hell = zeros( m, hlen );
% haccr = zeros(m, hlenall);
hZ = zeros( m*N, hlen );


K = zeros( N, N, m );
R = ones( d, m );
accr = zeros( m, 1 );
count = 0;

% Calculates the squared distance
ds = sq_dist( prs.tpts, prs.tpts );

% Preallocation
AA = zeros(d,prs.lags+1);
W = zeros(d,m,prs.lags +1);
Y = zeros(d,N,prs.lags+1);
a = zeros(1,prs.lags +1);

for r=1:prs.chains                                      % Chain loop
    
    % set initial values
    weight = exp(linspace(0,-3,prs.lags+1));     
    for s = 1:prs.lags+1
        W(:,:,s) = weight(s)*randn( d, m);              % Normal
    end
    Sigma = 1./gamrnd( prs.ss, 1/prs.sr, [ d, 1 ] );    % Gamma
    Z = 0.5*randn( m, N );                              % Normal
    %upsilon = prs.upsilon;
%     upsilon = gamrnd( prs.us, 1./prs.ur, [ m 1 ]);      % Gamma 
    
if length(prs.es) > 1
    ell = gamrnd( prs.es, 1./prs.er)'; 
else
     ell = gamrnd( prs.es, 1./prs.er, [ m 1 ]); 
end
    
    D = ones([d m prs.lags+1]);
    
    Z = lags_shift(Z,prs.lags);                         % Lags shift Z
    
    % Creating the K matrix - use jitter for numerical stability
    for j=1:m
        K(:,:,j) = exp( -ds/(2*ell(j)^2)) + prs.jit*eye( N );
    end
    
    for k=1:prs.nsamples                                % Sample loop
        
        disp(k)
        
        % burn-in loop
        for l=1:prs.skip*( k == 1 ) + prs.stride*( k > 1 )
            
            disp(l)
            
            % ========================= Gibbs for Z ===================================
            for j = randperm( m )                       % Random order
                
                Z(j,:,:) = 0;
                
                for s = 1:prs.lags+1
                    Y(:,:,s) =  W(:,:,s)*Z(:,:,s);
                end
                e = X - sum(Y,3);
                    
                % Mean and covariance calculated from the Bayesian
                % framework
% prior middel for sum w^2 skal være en konstant
                AA(:,:) = W(:,j,:);
                xs = AA ./ Sigma(:,ones(prs.lags+1,1)); % d x prs.lags+1
                
                % construct A to replace cc     
                for s = 1:prs.lags+1                     
                    a(s) = sum(sum( xs(:,1:prs.lags+1-(s-1)) .* AA(:,s:end) ));
                end

                % Zeros since 0>T (ej 2T)
                A = toeplitz([a zeros(1,N-(prs.lags+1))]);

                % Cholesky  decomposition
                [L, p] = chol( K(:,:,j) + A\eye(N), 'lower' ); 
                if p > 0
                    sendmail('ditte.hald@gmail.com','MATLAB keyboard','MATLAB keyboard')
                    keyboard
                end
                
                V = L\K(:,:,j);
                Sj = K(:,:,j) - V'*V;
                
                
                
                [L,p] = chol( Sj, 'lower' );            % standard deviation
                if p > 0
                    sendmail('ditte.hald@gmail.com','MATLAB keyboard','MATLAB keyboard')
                    keyboard
                end
                
                
                mj = Sj*sum(e'*xs,2);               % mean, note we do the convolution here!
                
                Z(j,:,1) = mj + L*randn( N, 1);     % Update Z_j
                
                Z = lags_shift(Z(:,:,1),prs.lags);  % Update all lags
                
            end
            
            % ========================== Gibbs for W ==================================
            for s = 1:prs.lags+1
                Y(:,:,s) =  W(:,:,s)*Z(:,:,s);
            end
            e = X - sum(Y,3);
            
            for j = randperm( m )
               
               for s = 1:prs.lags+1
                   e = e + W(:,j,s)*Z(j,:,s);
                   
                   B = 1./( Z(j,1:N-(s-1),s)*Z(j,1:N-(s-1),s)' + D(:,j,s).*Sigma );
                   
                   mj = B.*( Z(j,:,s)*e' )';
                   
                   Sj = B.*Sigma;
                   
                   W(:,j,s) = normrnd( mj, sqrt( Sj ) ).*R(:,j);
                   
                   e = e - W(:,j,s)*Z(j,:,s);
               end
                
            end
            
            
            
            % ========================= Gibbs for Sigma ===============================
            % sample from conditional for alpha
            for s = 1:prs.lags+1
                D(:,:,s) = gamrnd( ( prs.ts + 0.5 )*ones( d, m ), 1./( prs.tr + 0.5*W(:,:,s).^2 ) );
% D(:,:,s) = ones( d, m );
           
            end
            
            % sample from conditional for Sigma
            for i=1:d
                Sigma(i) = 1/gamrnd( prs.ss + 0.5*N, ...
                    1/( prs.sr + 0.5*( e(i,:)*e(i,:)' ) ) );
            end
            
     % ============================ Update K ===================================
            if prs.eu
                
                % Gamma distribution for upsilon
%                 upsilonp = gamrnd( 1/prs.uss^2, prs.uss^2*upsilon ); % proposal has mean upsilon and standard deviation prs.uss*upsilon
%                 ellp = gamrnd( 1/prs.ess^2, prs.ess^2*ell );
                ellp = randRG(ell,prs.ess,[prs.m 1]);
                
                for j=1:m
                    % New K
%                     Kp = exp( -upsilonp(j)*ds ) + prs.jit*eye( N );
                    Kp = exp( -ds/(2*ellp(j)^2) ) + prs.jit*eye( N );
                    
                    num = logmvnpdf( Z(j,:,1), zeros( 1, N ), Kp ); %New
                    den = logmvnpdf( Z(j,:,1), zeros( 1, N ), K(:,:,j) );  %Old
                    
%                     ratio = exp( num - den ) * ...
%                         gampdf( upsilonp(j), prs.us, 1/prs.ur) / ...
%                         gampdf( upsilon(j), prs.us, 1/prs.ur) * ...
%                         gampdf( upsilon(j), 1/prs.uss^2, prs.uss^2*upsilonp(j) ) / ...
%                         gampdf( upsilonp(j), 1/prs.uss^2, prs.uss^2*upsilon(j) ) ;
%                         gampdf( upsilon(j), upsilonp(j).^2/prs.uss^2, prs.uss^2/upsilonp(j) ) / ...
%                         gampdf( upsilonp(j), upsilon(j).^2/prs.uss^2, prs.uss^2/upsilon(j) ) ;
%                     


                       ratio = exp( num - den ) * ...
                        gampdf( ellp(j), prs.es(j), 1/prs.er(j)) / ...
                        gampdf( ell(j), prs.es(j), 1/prs.er(j));  % * ...
%                         gampdf( ell(j), 1/prs.ess^2, prs.ess^2*ellp(j) ) / ...
%                         gampdf( ellp(j), 1/prs.ess^2, prs.ess^2*ell(j) ) ;
                    
                    

                    if rand() < min( ratio, 1 )     % Metropolis
                        accr(j) = accr(j) + 1;      % accuracy 
%                         upsilon(j) = upsilonp(j);   % Update upsilon
                        ell(j) = ellp(j);
                        K(:,:,j) = Kp;              % Update K
                    end 
                end 
                
%                 count = count +1;
                for s = 1:prs.lags+1
                Y(:,:,s) =  W(:,:,s)*Z(:,:,s);
                end
            xc = X - sum(Y,3);
            
                 const = -0.5 * d * N * log(2*pi);
%                     xc = X- W*Z;
%                     SIGMA = diag(Sigma);
%                     keyboard
                    term1 = -0.5*sum( sum(xc.^2,2)' ./ Sigma');
%                     term2 = const - 0.5 * logdet(SIGMA)   % scalar
                    term2 = const - 0.5 * N * sum(log(Sigma)) ;   % scalar
                    
                    loglikelihood = term1' + term2;
                    
                
                
            end
            count = count +1;
            disp(accr./count)
        hellC(:,count) = ell;
        haccr(:,count) = accr;
%         hloglike(:,count) = loglike;
%         hlikeratio(:,count) = likeratio;
        hloglikelihood(:,count) = loglikelihood;
%             accr/count
             hSigma(:,count) = Sigma;
             halpha(:,count) = D(:);
             
            if ~rem(l,floor(prs.skip/4))
                fprintf('%g%% of warm-up completed - accr %2.1f%%+/-%2.1f%%\n',...
                    l/prs.skip*100, mean(100*accr/l), std(100*accr/l));
            
%             HELLplot
%             drawnow
%             pause(0.001)
            end
            
            
            
            
            
        end % end burn-in loop
       
        if ~rem(k,floor(prs.nsamples/4))
            fprintf('%g%% of samples completed - accr %2.1f%%+/-%2.1f%%\n',...
                k/prs.nsamples*100, mean(100*accr/(k+prs.skip)), std(100*accr/(k+prs.skip)));
%         HELLplot
%         drawnow
%             pause(0.001)
        end
        
        
        idx = k + ( r - 1 )*prs.nsamples;
        
        % saves
%         hSigma(:,idx) = Sigma;
%         hW(:,idx) = W(:);
%         Z1 = Z(:,:,1);
%         
%         hZ(:,idx) = Z1(:);
%         Accr(:,idx) = accr(:);  
%         hell(:,idx) = ell;
        
           
        hW(:,idx) = W(:);
        Z1 = Z(:,:,1);
        hZ(:,idx) = Z1(:);
%         hupsilon(:,idx) = upsilon;
        hell(:,idx) = ell;
        
        
%         hellC(:,idx) = 
%         hupsilon(:,idx) = upsilon;
           
    end % end sample loop
end % chain loop

etime = cputime - itime;

fprintf( 'Done.\n' )

% return results
result.prs = prs;
result.hSigma = hSigma;
result.halpha = halpha;
result.hW = hW;
result.hZ = hZ;
result.etime = etime;
% result.hupsilon = hupsilon;
result.hell = hell;
result.accr = accr;
result.hellC = hellC;
result.haccr = haccr;
% result.hloglike = hloglike;
% result.hlikeratio = hlikeratio;
result.hloglikelihood = hloglikelihood;
       
end

function C = sq_dist(a, b, Q)
%% Computes the squared distance between a and b

if nargin < 1 || nargin > 3 || nargout > 1
    error('Wrong number of arguments.');
end

if nargin == 1 || isempty(b)         % input arguments are taken to be
    b = a;                           % identical if b is missing or empty
end

[D, n] = size(a);
[d, m] = size(b);
if d ~= D
    error('Error: column lengths must agree.');
end

if nargin < 3
    C = zeros(n,m);
    for d = 1:D
        C = C + (repmat(b(d,:), n, 1) - repmat(a(d,:)', 1, m)).^2;
    end
    % C = repmat(sum(a.*a)',1,m)+repmat(sum(b.*b),n,1)-2*a'*b could be used to
    % replace the 3 lines above; it would be faster, but numerically less stable.
else
    if [n m] == size(Q) %#ok<BDSCA>
        C = zeros(D,1);
        for d = 1:D
            C(d) = sum(sum((repmat(b(d,:), n, 1) - repmat(a(d,:)', 1, m)).^2.*Q));
        end
    else
        error('Third argument has wrong size.');
    end
end

end


function [logp] = logmvnpdf(x,mu,Sigma)
% outputs log likelihood array for observations x  where x_n ~ N(mu,Sigma)
% x is NxD, mu is 1xD, Sigma is DxD

[~,D] = size(x);
const = -0.5 * D * log(2*pi);

xc = bsxfun(@minus,x,mu);

term1 = -0.5 * sum((xc / Sigma) .* xc, 2); % N x 1
term2 = const - 0.5 * logdet(Sigma);    % scalar
logp = term1' + term2;

end

function y = logdet(A)

U = chol(A);
y = 2*sum(log(diag(U)));

end
