function [dnose2ell,t0s,t1s] = distfromnose2otherfly(trx,fly1,t0,t1)

nflies = length(trx);
i0 = trx(fly1).f2i(t0);
i1 = trx(fly1).f2i(t1);
xnose = trx(fly1).x(i0:i1) + 2*trx(fly1).a(i0:i1).*cos(trx(fly1).theta(i0:i1));
ynose = trx(fly1).y(i0:i1) + 2*trx(fly1).a(i0:i1).*sin(trx(fly1).theta(i0:i1));

t0s = ones(1,nflies);
t1s = zeros(1,nflies);
dnose2ell = cell(1,nflies);

for fly2 = 1:nflies,
  if fly2 == fly1, continue; end
  fprintf('Computing distances from noses to fly %d for all frames, other flies.\n',fly2);
  % loop through each frame
  t0c = max(t0,trx(fly2).firstframe);
  t1c = min(t1,trx(fly2).endframe);
  for t = t0c:t1c,
    
    i2 = trx(fly2).f2i(t);
    
    % compute distance to fly2's ellipse from each nose
    docompute = t >= t0s & t <= t1s;
    
    d(:) = nan;
    i1 = f2i(t);
    d(docompute) = ellipsedist_hack(trx(fly2).x(i2),trx(fly2).y(i2),trx(fly2).a(i2),...
      trx(fly2).b(i2),trx(fly2).theta(i2),xnose(docompute,i1),ynose(docompute,i1));
  
    % store
    for fly1 = flies1(docompute),
      if fly1 == fly2, continue; end
      trx(fly1).dnose2ell(fly2,trx(fly1).f2i(t)) = d(fly1);
    end

  end % end loop over frame number

end % end loop over fly2