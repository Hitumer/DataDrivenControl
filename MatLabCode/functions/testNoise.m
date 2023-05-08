mu = [2,3];
sigma = [1 1.5; 1.5 3];
[V,D] = eig(sigma);
W = repmat(mu',1,100)+(V*sqrt(D))*randn(2,100);
