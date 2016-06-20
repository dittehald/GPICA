function [ result ] = GPICA4fmri( X, prs )
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
% =========================================================================


% sizes
m = prs.m;
[ d N ] = size( X );

% standarize
if prs.std == 1
     [X,mu_X,sig_X] = zscore( X, 0, 2 );
    
    prs.tr = prs.tr/mean(sig_X)^2;
    prs.sr = prs.sr/mean(sig_X)^2;
else
    X = bsxfun( @minus, X , mean( X, 2 ) );
end

itime = cputime;

% array allocation
hlen = prs.nsamples*prs.chains;
hSigma = zeros( d, hlen );
hW = zeros( d*m, hlen );
hell = zeros( m, hlen );
hZ = zeros( m*N, hlen );
K = zeros( N, N, m );
R = ones( d, m );
accr = zeros( m, 1 );
count = 0;
D = ones([d m]);

% Calculates the squared distance
ds = sq_dist( prs.tpts, prs.tpts );

for r=1:prs.chains                                      % Chain loop
    
    % set initial values
    W = 0.5*randn( d, m );                              % Normal
    Sigma = 1./gamrnd( prs.ss, 1/prs.sr, [ d, 1 ] );    % Gamma
    Z = 0.5*randn( m, N );                              % Normal
    
    
    if length(prs.es ) > 1
        ell = gamrnd( prs.es, 1./prs.er)';
    else
        ell = gamrnd( prs.es, 1./prs.er, [ prs.m 1 ]);
    end
    
    
    
    % Creating the K matrix - use jitter for numerical stability
    for j=1:m
        K(:,:,j) = exp( -ds/(2*(ell(j))^2)) + prs.jit*eye( N );
    end
    
    
    for k=1:prs.nsamples                                % Sample loop
        
        k
        
        % burn-in loop
        for l=1:prs.skip*( k == 1 ) + prs.stride*( k > 1 )
            
            l
            % ========================= Gibbs for Z ===================================
            for j = randperm( m )                       % Random order
                
                Z(j,:) = 0;
                
                e = X - W*Z;
                
                % Mean and covariance calculated from the Bayesian
                % framework
                
                xs = W(:,j)./Sigma;
                cc = xs'*W(:,j);
                
                % Cholesky  decomposition
                L = chol( K(:,:,j) + eye( N )/cc, 'lower' );
                V = L\K(:,:,j);
                Sj = K(:,:,j) - V'*V;
                
                
                L = chol( Sj, 'lower' );      %chol to draw samples with that covar
                
                mj = Sj*(e'*xs);
                
                Z(j,:) = mj + L*randn( N, 1 );
            end
            
            % ========================== Gibbs for W ==================================
            e = X - W*Z;
            
            for j = randperm( m )
                
                e = e + W(:,j)*Z(j,:);
                
                B = 1./( Z(j,:)*Z(j,:)' + D(:,j).*Sigma );
                
                mj = B.*( Z(j,:)*e' )';
                
                Sj = B.*Sigma;
                
                W(:,j) = normrnd( mj, sqrt( Sj ) ).*R(:,j);
                
                e = e - W(:,j)*Z(j,:);
            end
            
            % ====================== Update Alpha (W prior) ===========================
            
            D = gamrnd( ( prs.ts + 0.5 )*ones( d, m ), 1./( prs.tr + 0.5*W(:,:).^2 ) );
            
            
            % ========================= Gibbs for Sigma ===============================
            
            % sample from conditional for Sigma
            for i=1:d
                Sigma(i) = 1/gamrnd( prs.ss + 0.5*N, ...
                    1/( prs.sr + 0.5*( e(i,:)*e(i,:)' ) ) );
            end
            
            % ============================ Update K ===================================
            if prs.eu
                
                ellp = randRG(ell,prs.ess,[prs.m 1]);
                
                for j=1:m
                    
                    % New K
                    Kp = exp( -ds/(2*ellp(j)^2) ) + prs.jit*eye( N );
                    
                    num = logmvnpdf( Z(j,:), zeros( 1, N ), Kp ); %New
                    den = logmvnpdf( Z(j,:), zeros( 1, N ), K(:,:,j) );  %Old
                    
                    ratio = exp( num - den ) * ...
                        gampdf( ellp(j), prs.es(j), 1/prs.er(j)) / ...
                        gampdf( ell(j), prs.es(j), 1/prs.er(j));
                    
                    
                    if rand() < min( ratio, 1 )     % Metropolis
                        accr(j) = accr(j) + 1;      % accuracy %Hvor mange gange der bliver accepteret
                        %                         upsilon(j) = upsilonp(j);   % Update upsilon
                        ell(j) = ellp(j);
                        K(:,:,j) = Kp;              % Update K
                        
                    end
                    
                    
                    
                end % end source loop
                
                % Calc loglikehihood
                const = -0.5 * d * N * log(2*pi);
                xc = X- W*Z;
                term1 = -0.5*sum( sum(xc.^2,2)' ./ Sigma');
                term2 = const - 0.5 * N * sum(log(Sigma)) ;   % scalar
                
                loglikelihood = term1 + term2;
                
                
            end % end if K update loop
            
            count = count +1;
            disp(accr./count)
            hellC(:,count) = ell;
            haccr(:,count) = accr;
            hloglikelihood(:,count) = loglikelihood;
            
             hSigma(:,count) = Sigma;
             
             if d < 10000
              halpha(:,count) = D(:);
             end
             
             
             
                if ~rem(l,floor(prs.skip/4))
%             fprintf('%g%% of samples completed - accr %2.1f%%+/-%2.1f%%\n',...
%                 k/prs.nsamples*100, mean(100*accr/(k+prs.skip)), std(100*accr/(k+prs.skip)));
        HELLplot
        drawnow
            pause(0.001)
                end
             
        end % end burn-in loop
        
         if ~rem(k,floor(prs.nsamples/4)) %&& prs.nsamples>500
%             fprintf('%g%% of samples completed - accr %2.1f%%+/-%2.1f%%\n',...
%                 k/prs.nsamples*100, mean(100*accr/(k+prs.skip)), std(100*accr/(k+prs.skip)));
        HELLplot
        drawnow
            pause(0.001)
        end
        
        idx = k + ( r - 1 )*prs.nsamples;
        
        % saves
         hSigma(:,idx) = Sigma;
        hW(:,idx) = W(:);
        hZ(:,idx) = Z(:);
        hell(:,idx) = ell;
        
    end % end sample loop
    
end % end chain loop

if d > 10000
    halpha = D(:);
end

etime = cputime - itime;

fprintf( 'Done.\n' )

% return results

result.prs = prs;
result.halpha = halpha;
result.hSigma = hSigma;
result.hW = hW;
result.hZ = hZ;
result.etime = etime;
result.hell = hell;
result.accr = accr;
result.hellC = hellC;
result.haccr = haccr;
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

[N,D] = size(x);
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