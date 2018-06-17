function f = minmaxobj(y,sigma,Ptot,tau_d,c,alpha)
%The objective function

calphap2 = c.*alpha.^2;
term3 = sqrt(tau_d^2+c^2.*alpha.^4.*Ptot^2*y^2-2*tau_d*y.*(c.*alpha.^2*Ptot+2*sigma^2));
term3 = sum(term3);
f = y*(sigma^2+Ptot/(2*tau_d)*sum(calphap2))-1/(2*tau_d)*term3;

end

