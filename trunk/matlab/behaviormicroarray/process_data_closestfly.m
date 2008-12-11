function trx = process_data_closestfly(trx,fov)

nflies = length(trx);

for fly1 = 1:nflies,
  
  fprintf('Computing closest fly statistics for fly %d\n',fly1);

  % position of nose
  xnose = trx(fly1).x + 2*trx(fly1).a.*cos(trx(fly1).theta);
  ynose = trx(fly1).y + 2*trx(fly1).a.*sin(trx(fly1).theta);
 
  % initialize
  dcenter = nan(nflies,trx(fly1).nframes);
  dnose2ell = nan(nflies,trx(fly1).nframes);
  dell2nose = nan(nflies,trx(fly1).nframes);
  magveldiff = nan(nflies,trx(fly1).nframes);
  veltoward = nan(nflies,trx(fly1).nframes-1);
  absthetadiff = nan(nflies,trx(fly1).nframes);
  absphidiff = nan(nflies,trx(fly1).nframes-1);
  absanglefrom2to1 = nan(nflies,trx(fly1).nframes);
  anglesub = nan(nflies,trx(fly1).nframes);
  
  % loop over other flies
  for fly2 = 1:nflies,
    if fly2 == fly1, continue; end
    
    fprintf('Other fly = %d\n',fly2);
    
    % get start and end frames of overlap
    t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
    t1 = min(trx(fly1).endframe,trx(fly2).endframe);
    
    % no overlap
    if t1 < t0, continue; end
    
    % indices for these frames
    offi = trx(fly1).firstframe-1;
    offj = trx(fly2).firstframe-1;
    i0 = t0 - offi;
    i1 = t1 - offi;
    j0 = t0 - offj;
    j1 = t1 - offj;

    % centroid distance
    dx = trx(fly2).x(j0:j1)-trx(fly1).x(i0:i1);
    dy = trx(fly2).y(j0:j1)-trx(fly1).y(i0:i1);
    z = sqrt(dx.^2 + dy.^2);
    dcenter(fly2,i0:i1) = z;
    % direction to other fly
    dx = dx ./ z;
    dy = dy ./ z;
    
    % distance from fly1's nose to fly2
    for t = t0:t1,
      i = t - offi;
      j = t - offj;
      dnose2ell(fly2,i) = ellipsedist_hack(trx(fly2).x(j),trx(fly2).y(j),...
        trx(fly2).a(j),trx(fly2).b(j),trx(fly2).theta(j),...
        xnose(i),ynose(i));
    end
    
    % distance from fly2's nose to fly1
    for t = t0:t1,
      i = t - offi;
      j = t - offj;    
      xnose2 = trx(fly2).x(j) + 2*trx(fly2).a(j).*cos(trx(fly2).theta(j));
      ynose2 = trx(fly2).y(j) + 2*trx(fly2).a(j).*sin(trx(fly2).theta(j));
      dell2nose(fly2,i) = ellipsedist_hack(trx(fly1).x(i),trx(fly1).y(i),...
        trx(fly1).a(i),trx(fly1).b(i),trx(fly1).theta(i),...
        xnose2,ynose2);
    end
    
    % angle of fly1's vision subtended by fly2
    for t = t0:t1,
      i = t - offi;
      j = t - offj;
      anglesub(fly2,i) = anglesubtended(...
        trx(fly1).x(i),trx(fly1).y(i),trx(fly1).a(i),trx(fly1).b(i),trx(fly1).theta(i),...
        trx(fly2).x(j),trx(fly2).y(j),trx(fly2).a(j),trx(fly2).b(j),trx(fly2).theta(j),...
        fov);
    end
    
    % velocity difference
    magveldiff(fly2,i0:i1-1) = ...
      sqrt( (diff(trx(fly1).x(i0:i1))-diff(trx(fly2).x(j0:j1))).^2 + ...
      (diff(trx(fly1).y(i0:i1))-diff(trx(fly2).y(j0:j1))).^2 );
    
    % velocity in direction of other fly
    veltoward(fly2,i0:i1-1) = dx(1:end-1).*diff(trx(fly1).x(i0:i1)) + ...
      dy(1:end-1).*diff(trx(fly2).y(j0:j1));
    
    % orientation of fly2 relative to orientation of fly1
    absthetadiff(fly2,i0:i1) = abs(modrange(trx(fly2).theta(j0:j1) - trx(fly1).theta(i0:i1),-pi,pi));

    % velocity direction of fly2 relative to fly1's velocity direction
    absphidiff(fly2,i0:i1-1) = abs(modrange(trx(fly2).phi(j0:j1-1)-trx(fly1).phi(i0:i1-1),-pi,pi));
  
    % direction to fly2 from fly1
    absanglefrom2to1(fly2,i0:i1) = abs(modrange(atan2(dy,dx)-trx(fly1).theta(i0:i1),-pi,pi));

  end % end loop over fly2

  % closest fly according to centroid distance
  [trx(fly1).dcenter,trx(fly1).closestfly_center] = min(dcenter,[],1);
  
  % closest fly according to dnose2ell
  [trx(fly1).dnose2ell,trx(fly1).closestfly_nose2ell] = min(dnose2ell,[],1);
  
  % closest fly according to dell2nose
  [trx(fly1).dell2nose,trx(fly1).closestfly_ell2nose] = min(dell2nose,[],1);
  
  % closest fly according to angle subtended
  [trx(fly1).anglesub,trx(fly1).closestfly_anglesub] = max(anglesub,[],1);
  
  % magveldiff for each of these measured closest flies
  idx = 1:trx(fly1).nframes-1;
  trx(fly1).magveldiff_center = magveldiff(sub2ind(size(magveldiff),trx(fly1).closestfly_center(idx),idx));
  trx(fly1).magveldiff_nose2ell = magveldiff(sub2ind(size(magveldiff),trx(fly1).closestfly_nose2ell(idx),idx));
  trx(fly1).magveldiff_ell2nose = magveldiff(sub2ind(size(magveldiff),trx(fly1).closestfly_ell2nose(idx),idx));
  trx(fly1).magveldiff_anglesub = magveldiff(sub2ind(size(magveldiff),trx(fly1).closestfly_anglesub(idx),idx));
  
  % veltoward for each of these measured closest flies
  idx = 1:trx(fly1).nframes-1;
  trx(fly1).veltoward_center = veltoward(sub2ind(size(veltoward),trx(fly1).closestfly_center(idx),idx));
  trx(fly1).veltoward_nose2ell = veltoward(sub2ind(size(veltoward),trx(fly1).closestfly_nose2ell(idx),idx));
  trx(fly1).veltoward_ell2nose = veltoward(sub2ind(size(veltoward),trx(fly1).closestfly_ell2nose(idx),idx));
  trx(fly1).veltoward_anglesub = veltoward(sub2ind(size(veltoward),trx(fly1).closestfly_anglesub(idx),idx));

  % absthetadiff for each of these measured closest flies
  idx = 1:trx(fly1).nframes;
  trx(fly1).absthetadiff_center = absthetadiff(sub2ind(size(absthetadiff),trx(fly1).closestfly_center(idx),idx));
  trx(fly1).absthetadiff_nose2ell = absthetadiff(sub2ind(size(absthetadiff),trx(fly1).closestfly_nose2ell(idx),idx));
  trx(fly1).absthetadiff_ell2nose = absthetadiff(sub2ind(size(absthetadiff),trx(fly1).closestfly_ell2nose(idx),idx));
  trx(fly1).absthetadiff_anglesub = absthetadiff(sub2ind(size(absthetadiff),trx(fly1).closestfly_anglesub(idx),idx));

  % absphidiff for each of these measured closest flies
  idx = 1:trx(fly1).nframes-1;
  trx(fly1).absphidiff_center = absphidiff(sub2ind(size(absphidiff),trx(fly1).closestfly_center(idx),idx));
  trx(fly1).absphidiff_nose2ell = absphidiff(sub2ind(size(absphidiff),trx(fly1).closestfly_nose2ell(idx),idx));
  trx(fly1).absphidiff_ell2nose = absphidiff(sub2ind(size(absphidiff),trx(fly1).closestfly_ell2nose(idx),idx));
  trx(fly1).absphidiff_anglesub = absphidiff(sub2ind(size(absphidiff),trx(fly1).closestfly_anglesub(idx),idx));

  % anglefrom2to1 for each of these measured closest flies
  idx = 1:trx(fly1).nframes;
  trx(fly1).absanglefrom2to1_center = absanglefrom2to1(sub2ind(size(absanglefrom2to1),trx(fly1).closestfly_center(idx),idx));
  trx(fly1).absanglefrom2to1_nose2ell = absanglefrom2to1(sub2ind(size(absanglefrom2to1),trx(fly1).closestfly_nose2ell(idx),idx));
  trx(fly1).absanglefrom2to1_ell2nose = absanglefrom2to1(sub2ind(size(absanglefrom2to1),trx(fly1).closestfly_ell2nose(idx),idx));
  trx(fly1).absanglefrom2to1_anglesub = absanglefrom2to1(sub2ind(size(absanglefrom2to1),trx(fly1).closestfly_anglesub(idx),idx));

  % change in various parameters
  trx(fly1).ddcenter = diff(trx(fly1).dcenter);
  trx(fly1).ddnose2ell = diff(trx(fly1).dnose2ell);
  trx(fly1).ddell2nose = diff(trx(fly1).dell2nose);
  trx(fly1).danglesub = diff(trx(fly1).anglesub);
  
end