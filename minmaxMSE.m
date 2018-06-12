% min-max MSE
clear
clc
%% basic environment
PL = [40 60];
alpha = db2lin(-PL);
sigma = sqrt(db2lin(-174)*180000);
c = 1;
Ptot = 250;
tau_d = 6;
Nr = 8;

% alpha = [1 1];
% sigma = 1;
% c = 1;
% Ptot = 1;
% tau_d = 1;

%%
% y = minmaxylb(tau_d,sigma,Ptot,alpha,c)
% y=10^(-10)
% f = minmaxfun(sigma,Ptot,tau_d,c,alpha,y)
%  yl = minmaxylb(tau_d,sigma,Ptot,alpha,c)

% y_list = linspace(3+2*sqrt(2)+0.01,10,100);
% for i = 1:length(y_list)
%     y = y_list(i);
%     df(i) = minmaxfun(sigma,Ptot,tau_d,c,alpha,y);
%     f(i) = minmaxobj(y,sigma,Ptot,tau_d,c,alpha);
% end
% plot(y_list,df);
% figure,plot(y_list,f)

% syms y
% calphap2 = c.*alpha.^2;
% term3 = sqrt(tau_d^2+c^2.*alpha.^4.*Ptot^2*y^2-2*tau_d*y.*(c.*alpha.^2*Ptot+2*sigma^2));
% term3 = sum(term3);
% f = y*(sigma^2+Ptot/(2*tau_d)*sum(calphap2))-1/(2*tau_d)*term3
%% lower bound of y
yl = minmaxylb(tau_d,sigma,Ptot,alpha,c);
for i = 1:2
    coe = [c^2*alpha(i)^4*Ptot^2,-2*tau_d*(c*alpha(i)^2*Ptot+2*sigma^2),tau_d^2];
    r(i,:) = sort(roots(coe),'ascend');
end
ylb_complex = max(r(end,:));%the bounds that avoid complex solution
yup_complex = min(r(1,:));
yl = max(ylb_complex,yl)+1;% +1 ensures the derivative does not go to infinity
yl_order = floor(log(abs(yl))./log(10));
%%
% y_list = linspace(yl+1,4e8,10000);
% for i = 1:length(y_list)
%     y = y_list(i);
%     df(i) = minmaxfun(sigma,Ptot,tau_d,c,alpha,y);
%     f(i) = minmaxobj(y,sigma,Ptot,tau_d,c,alpha);
% end
% plot(y_list,df);
% figure,plot(y_list,f);
% [val,pos] = min(f);
% y_opt = y_list(pos);

%% upper bound
i = 0;
step = 10^(yl_order/2);% step length of search
while true
    i = i+1;
    yu = yl+step;
    if minmaxfun(sigma,Ptot,tau_d,c,alpha,yl)*minmaxfun(sigma,Ptot,tau_d,c,alpha,yu) < 0
        break
    else
        yl = yu;
    end
end
%% bisection search
threshold = 0.1;
j = 1;
while abs(yl-yu)>threshold
    j = j+1;
    y = (yl+yu)/2;
    if minmaxfun(sigma,Ptot,tau_d,c,alpha,y)>0
        yu = y;
    else
        yl = y;
    end
end
opty = (yl+yu)/2;% optimal y
optP = minmaxBP(opty,tau_d,c,Ptot,alpha,sigma);% optimal data power of all users
%% Optimal MSE
Pgamma = c^2.*optP.*alpha.^4.*(Ptot-tau_d.*optP)./(sigma^2+c.*alpha.^2.*(Ptot-tau_d.*optP));
mu = -1+(sum(alpha.^2.*optP.*c)+sigma^2)./Pgamma;
for i = 1:2
    MSE(i) = MSEmu(Nr,mu(i));
end
% load('minmaxMSE');
% MSE = [MSE;tempMSE];
% save('minmaxMSE','MSE')
