function yl = minmaxylb(tau_d,sigma,Ptot,alpha,c)
%Compute the lower bound on the variable y

ylb_nu = tau_d.*(2*sigma.*(sigma^2+Ptot*c.*alpha.^2).^(3/2)+Ptot^2.*alpha.^4.*c^2+2*sigma^4+...
    3*Ptot*c*sigma^2.*alpha.^2);
ylb_de = Ptot^2*c^2.*alpha.^4.*(sigma^2+Ptot*c.*alpha.^2);
y = ylb_nu./ylb_de;
yl = max(ylb_nu./ylb_de);

end

