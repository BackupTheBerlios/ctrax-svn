function v = systematic_ml_distr(x)

mu = mean(x,2);
sig = std(x,1,2);
sig(all(isinf(x),2)) = 0;

v = [mu,sig]';
v = v(:);