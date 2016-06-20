
function f = pdfrG(x,mu,sigma)
%%
% x = t;

y1 = exp(-0.5 .* ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);

y2 = exp(-0.5 .* ((-x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);

y3 = y2+y1;

y3(x<=0) = 0;

% if length(x)>100
% figure;
% hold all
% plot(-x,y1)
% plot(x,y1)
% plot(x,y3)
% end

f = y3;
%%

% end
% %%
% close all
% mu = 4;
% sigma = 3;
% prs.m = 10000;
% 
% x1 = sigma.*randn(prs.m,1)+mu;
%                    x = x1.*sign(x1);
%                    
% 
% H = histfit(x1,100);
% 
% figure;
% histogram(x,100)
% hold all
% 
% xx = H(2).XData;
% yy = H(2).YData;
% plot(xx,yy)
% 
% %%
% figure;
% 
% x = t;
% y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
% 
% if any(x<=0)
%     y2 = exp(-0.5 * ((x(x<=0) - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
%     plot(-x(x<=0),y2)
% end
% 
% x1 = [-x(x<=0) x(x>0)];
% x2 = unique(sort(x1));
% 
% y3 = y;
% y3(x<=0) = [];
% 
% y4 = y3 + y2;
% 
% hold all
% plot(x,y4)

% y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);

%%

% for t<0
%     z = y(-x)+y(-x)
% end

% 
% t = -5:0.1:20;
% y2 = pdf('Normal',t,mu,sigma);
% plot(t,y2)
%%

end



% Phi = cdf('normal',(-mu/sigma),0,1);
% DI = dirac(x);
% stepfunc = double(x>0);
% 
% f = Phi.*DI + (1/sqrt(2*pi*sigma^2)).*exp(-(x-mu).^2/(2*sigma^2)).*stepfunc;
% 
% end
% 
% x = -10:0.1:10; f = pdfrG(x, 5, 2); figure; plot(x,f)
% hold all; pf = pdf('Normal',x,5,2); plot(x,pf)