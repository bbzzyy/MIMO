function n_p = Gmultinoise( N,tau_p,sigma )
%spatially and temporally additive white Gaussian noise generator
% N:numbers of attenna
% taup:length of pilot sequence
% sigma:variance

mvncov_noise=repmat(0.5*sigma^2,1,2*N);
mu_noise=zeros(1,2*N);
r_noise=mvnrnd(mu_noise,mvncov_noise,tau_p);
x_noise=r_noise(:,1:N);
y_noise=r_noise(:,N+1:2*N);
n_p=transpose(complex(x_noise,y_noise));

end

