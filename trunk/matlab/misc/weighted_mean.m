function mu = weighted_mean(x,w)

[n,d] = size(x);

z = sum(w);
mu = sum(x.*repmat(w,1,d),1)/z;