function n_p = Gnoise( Ns,N,sigma )
% Generate Gaussian noise
%   Ns:number of samples
%   N:number of anttenas
%   sigma:covariance

cov_noise=sigma^2.*eye(N);
mvncov_noise=[0.5*real(cov_noise) -0.5*imag(cov_noise);0.5*imag(cov_noise) 0.5*real(cov_noise)];
mu_noise=zeros(1,2*N);
r_noise=mvnrnd(mu_noise,mvncov_noise,Ns);
x_noise=r_noise(:,1:N);
y_noise=r_noise(:,N+1:2*N);
n_p=complex(x_noise,y_noise);
% n_p=transpose(n_p);%%%
end

