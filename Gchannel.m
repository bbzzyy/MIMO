function h = Gchannel( Ns,N,cov )
% Generate Gaussian channel
%   Ns:number of samples
%   N:number of anttenas
%   cov:covariance matrix

if size(cov)~=N
    error('Incorrect dimension of covariance matrix');
end

mu_channel = zeros(1,2*N);%mean value
mvncov_channel = [0.5*real(cov) -0.5*imag(cov);0.5*imag(cov) 0.5*real(cov)];
r_channel = mvnrnd(mu_channel,mvncov_channel,Ns);
x_channel = r_channel(:,1:N);
y_channel = r_channel(:,N+1:2*N);
h = complex(x_channel,y_channel);%channel h

% h = transpose(h);%%%

end

