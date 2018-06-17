function Pd = minmaxMSEalgorithm(PL)
%min-max MSE
%% basic environment
alpha = db2lin(-PL);
sigma = sqrt(db2lin(-174)*180000);
c = 1;
Ptot = 250;
tau_d = 6;
Nr = 8;
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
Pd = minmaxBP(opty,tau_d,c,Ptot,alpha,sigma);% optimal data power of all users

end

