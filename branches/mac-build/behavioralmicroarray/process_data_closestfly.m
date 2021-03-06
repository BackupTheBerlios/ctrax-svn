function trx = process_data_closestfly(trx,fov)

nflies = length(trx);

for fly1 = 1:nflies,
  
  fprintf('Computing closest fly statistics for fly %d\n',fly1);

  % position of nose
  xnose = trx(fly1).x_mm + 2*trx(fly1).a_mm.*cos(trx(fly1).theta);
  ynose = trx(fly1).y_mm + 2*trx(fly1).a_mm.*sin(trx(fly1).theta);
  xtail1 = trx(fly1).x_mm - 2*trx(fly1).a_mm.*cos(trx(fly1).theta);
  ytail1 = trx(fly1).y_mm - 2*trx(fly1).a_mm.*sin(trx(fly1).theta);
  
 
  % initialize
  dcenter = nan(nflies,trx(fly1).nframes);
  dnose2ell = nan(nflies,trx(fly1).nframes);
  dell2nose = nan(nflies,trx(fly1).nframes);
  dnose2tail = nan(nflies,trx(fly1).nframes);
  dtail2nose = nan(nflies,trx(fly1).nframes);
  magveldiff = nan(nflies,trx(fly1).nframes-1);
  veltoward = nan(nflies,trx(fly1).nframes-1);
  absthetadiff = nan(nflies,trx(fly1).nframes);
  absphidiff = nan(nflies,trx(fly1).nframes-1);
  anglefrom1to2 = nan(nflies,trx(fly1).nframes);
  absanglefrom1to2 = nan(nflies,trx(fly1).nframes);
  anglefrom2to1 = nan(nflies,trx(fly1).nframes);
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
    dx = trx(fly2).x_mm(j0:j1)-trx(fly1).x_mm(i0:i1);
    dy = trx(fly2).y_mm(j0:j1)-trx(fly1).y_mm(i0:i1);
    z = sqrt(dx.^2 + dy.^2);
    dcenter(fly2,i0:i1) = z;
    % direction to other fly
    dx = dx ./ z;
    dy = dy ./ z;
    
    % distance from fly1's nose to fly2
    for t = t0:t1,
      i = t - offi;
      j = t - offj;
      dnose2ell(fly2,i) = ellipsedist_hack(trx(fly2).x_mm(j),trx(fly2).y_mm(j),...
        trx(fly2).a_mm(j),trx(fly2).b_mm(j),trx(fly2).theta(j),...
        xnose(i),ynose(i));
    end
    
    % distance from fly2's nose to fly1
    xnose2 = trx(fly2).x_mm(j0:j1) + 2*trx(fly2).a_mm(j0:j1).*cos(trx(fly2).theta(j0:j1));
    ynose2 = trx(fly2).y_mm(j0:j1) + 2*trx(fly2).a_mm(j0:j1).*sin(trx(fly2).theta(j0:j1));

    for t = t0:t1,
      i = t - offi;
      j = t - t0 + 1;
      dell2nose(fly2,i) = ellipsedist_hack(trx(fly1).x_mm(i),trx(fly1).y_mm(i),...
        trx(fly1).a_mm(i),trx(fly1).b_mm(i),trx(fly1).theta(i),...
        xnose2(j),ynose2(j));
    end
    
    % distance from fly1's nose to fly2's tail
    xtail2 = trx(fly2).x_mm(j0:j1) - 2*trx(fly2).a_mm(j0:j1).*cos(trx(fly2).theta(j0:j1));
    ytail2 = trx(fly2).y_mm(j0:j1) - 2*trx(fly2).a_mm(j0:j1).*sin(trx(fly2).theta(j0:j1));
    dnose2tail(fly2,i0:i1) = sqrt((xnose(i0:i1)-xtail2).^2 + (ynose(i0:i1)-ytail2).^2);

    % distance from fly1's tail to fly2's nose
    dtail2nose(fly2,i0:i1) = sqrt((xnose2-xtail1(i0:i1)).^2 + (ynose2-ytail1(i0:i1)).^2);
    
    % angle of fly1's vision subtended by fly2
    for t = t0:t1,
      i = t - offi;
      j = t - offj;
      anglesub(fly2,i) = anglesubtended(...
        trx(fly1).x_mm(i),trx(fly1).y_mm(i),trx(fly1).a_mm(i),trx(fly1).b_mm(i),trx(fly1).theta(i),...
        trx(fly2).x_mm(j),trx(fly2).y_mm(j),trx(fly2).a_mm(j),trx(fly2).b_mm(j),trx(fly2).theta(j),...
        fov);
    end
    
    % velocity difference
    magveldiff(fly2,i0:i1-1) = ...
      sqrt( (diff(trx(fly1).x_mm(i0:i1))-diff(trx(fly2).x_mm(j0:j1))).^2 + ...
      (diff(trx(fly1).y_mm(i0:i1))-diff(trx(fly2).y_mm(j0:j1))).^2 ).*trx(fly1).fps;
    if i1 < trx(fly1).nframes && i1 -1 >= i0,
      magveldiff(fly2,i1) = magveldiff(fly2,i1-1);
    end
    
    % velocity in direction of other fly
    if i1 < trx(fly1).nframes,
      if i1 - i0 + 1 > 0,
        veltoward(fly2,i0:i1) = (dx(1:end).*diff(trx(fly1).x_mm(i0:i1+1)) + ...
          dy(1:end).*diff(trx(fly1).y_mm(i0:i1+1))).*trx(fly1).fps;
      end
    else
      if i1 - i0 > 0,
        veltoward(fly2,i0:i1-1) = (dx(1:end-1).*diff(trx(fly1).x_mm(i0:i1)) + ...
          dy(1:end-1).*diff(trx(fly1).y_mm(i0:i1))).*trx(fly1).fps;
      end
    end
    
    % orientation of fly2 relative to orientation of fly1
    absthetadiff(fly2,i0:i1) = abs(modrange(trx(fly2).theta(j0:j1) - trx(fly1).theta(i0:i1),-pi,pi));

    % velocity direction of fly2 relative to fly1's velocity direction
    absphidiff(fly2,i0:i1-1) = abs(modrange(trx(fly2).phi(j0:j1-1)-trx(fly1).phi(i0:i1-1),-pi,pi));
    if i1 < trx(fly1).nframes && i1 -1 >= i0,
      absphidiff(fly2,i1) = absphidiff(fly2,i1-1);
    end
  
    % direction to fly2 from fly1
    anglefrom1to2(fly2,i0:i1) = modrange(atan2(dy,dx)-trx(fly1).theta(i0:i1),-pi,pi);
    absanglefrom1to2(fly2,i0:i1) = abs(anglefrom1to2(fly2,i0:i1));

    % direction to fly1 from fly2
    anglefrom2to1(fly2,i0:i1) = modrange(atan2(-dy,-dx)-trx(fly2).theta(i0:i1),-pi,pi);
    absanglefrom2to1(fly2,i0:i1) = abs(anglefrom2to1(fly2,i0:i1));

    
  end % end loop over fly2

  % closest fly according to centroid distance
  [trx(fly1).dcenter,trx(fly1).closestfly_center] = min(dcenter,[],1);
  trx(fly1).units.dcenter = parseunits('mm');
  trx(fly1).units.closestfly_center = parseunits('unit');
  
  % closest fly according to dnose2ell
  [trx(fly1).dnose2ell,trx(fly1).closestfly_nose2ell] = min(dnose2ell,[],1);
  trx(fly1).units.dnose2ell = parseunits('mm');
  trx(fly1).units.closestfly_nose2ell = parseunits('unit');
  
  % closest fly according to dell2nose
  [trx(fly1).dell2nose,trx(fly1).closestfly_ell2nose] = min(dell2nose,[],1);
  trx(fly1).units.dell2nose = parseunits('mm');
  trx(fly1).units.closestfly_ell2nose = parseunits('unit');

  % closest fly according to dnose2tail
  [trx(fly1).dnose2tail,trx(fly1).closestfly_nose2tail] = min(dnose2tail,[],1);
  trx(fly1).units.dnose2tail = parseunits('mm');
  trx(fly1).units.closestfly_nose2tail = parseunits('unit');

  % closest fly according to dtail2nose
  [trx(fly1).dtail2nose,trx(fly1).closestfly_tail2nose] = min(dtail2nose,[],1);
  trx(fly1).units.dtail2nose = parseunits('mm');
  trx(fly1).units.closestfly_tail2nose = parseunits('unit');
  
  % closest fly according to angle subtended
  [trx(fly1).anglesub,trx(fly1).closestfly_anglesub] = max(anglesub,[],1);
  trx(fly1).units.anglesub = parseunits('rad');
  trx(fly1).units.closestfly_anglesub = parseunits('unit');
  
  % magveldiff for each of these measured closest flies
  idx = 1:trx(fly1).nframes-1;
  trx(fly1).magveldiff_center = magveldiff(sub2ind(size(magveldiff),trx(fly1).closestfly_center(idx),idx));
  trx(fly1).units.magveldiff_center = parseunits('mm/s');
  trx(fly1).magveldiff_nose2ell = magveldiff(sub2ind(size(magveldiff),trx(fly1).closestfly_nose2ell(idx),idx));
  trx(fly1).units.magveldiff_nose2ell = parseunits('mm/s');
  trx(fly1).magveldiff_ell2nose = magveldiff(sub2ind(size(magveldiff),trx(fly1).closestfly_ell2nose(idx),idx));
  trx(fly1).units.magveldiff_ell2nose = parseunits('mm/s');
  trx(fly1).magveldiff_nose2tail = magveldiff(sub2ind(size(magveldiff),trx(fly1).closestfly_nose2tail(idx),idx));
  trx(fly1).units.magveldiff_nose2tail = parseunits('mm/s');
  trx(fly1).magveldiff_tail2nose = magveldiff(sub2ind(size(magveldiff),trx(fly1).closestfly_tail2nose(idx),idx));
  trx(fly1).units.magveldiff_tail2nose = parseunits('mm/s');
  trx(fly1).magveldiff_anglesub = magveldiff(sub2ind(size(magveldiff),trx(fly1).closestfly_anglesub(idx),idx));
  trx(fly1).units.magveldiff_anglesub = parseunits('mm/s');
  
  % veltoward for each of these measured closest flies
  idx = 1:trx(fly1).nframes-1;
  trx(fly1).veltoward_center = veltoward(sub2ind(size(veltoward),trx(fly1).closestfly_center(idx),idx));
  trx(fly1).units.veltoward_center = parseunits('mm/s');
  trx(fly1).veltoward_nose2ell = veltoward(sub2ind(size(veltoward),trx(fly1).closestfly_nose2ell(idx),idx));
  trx(fly1).units.veltoward_nose2ell = parseunits('mm/s');
  trx(fly1).veltoward_ell2nose = veltoward(sub2ind(size(veltoward),trx(fly1).closestfly_ell2nose(idx),idx));
  trx(fly1).units.veltoward_ell2nose = parseunits('mm/s');
  trx(fly1).veltoward_nose2tail = veltoward(sub2ind(size(veltoward),trx(fly1).closestfly_nose2tail(idx),idx));
  trx(fly1).units.veltoward_nose2tail = parseunits('mm/s');
  trx(fly1).veltoward_tail2nose = veltoward(sub2ind(size(veltoward),trx(fly1).closestfly_tail2nose(idx),idx));
  trx(fly1).units.veltoward_tail2nose = parseunits('mm/s');
  trx(fly1).veltoward_anglesub = veltoward(sub2ind(size(veltoward),trx(fly1).closestfly_anglesub(idx),idx));
  trx(fly1).units.veltoward_anglesub = parseunits('mm/s');

  % absthetadiff for each of these measured closest flies
  idx = 1:trx(fly1).nframes;
  trx(fly1).absthetadiff_center = absthetadiff(sub2ind(size(absthetadiff),trx(fly1).closestfly_center(idx),idx));
  trx(fly1).units.absthetadiff_center = parseunits('rad');
  trx(fly1).absthetadiff_nose2ell = absthetadiff(sub2ind(size(absthetadiff),trx(fly1).closestfly_nose2ell(idx),idx));
  trx(fly1).units.absthetadiff_nose2ell = parseunits('rad');
  trx(fly1).absthetadiff_ell2nose = absthetadiff(sub2ind(size(absthetadiff),trx(fly1).closestfly_ell2nose(idx),idx));
  trx(fly1).units.absthetadiff_ell2nose = parseunits('rad');
  trx(fly1).absthetadiff_nose2tail = absthetadiff(sub2ind(size(absthetadiff),trx(fly1).closestfly_nose2tail(idx),idx));
  trx(fly1).units.absthetadiff_nose2tail = parseunits('rad');
  trx(fly1).absthetadiff_tail2nose = absthetadiff(sub2ind(size(absthetadiff),trx(fly1).closestfly_tail2nose(idx),idx));
  trx(fly1).units.absthetadiff_tail2nose = parseunits('rad');
  trx(fly1).absthetadiff_anglesub = absthetadiff(sub2ind(size(absthetadiff),trx(fly1).closestfly_anglesub(idx),idx));
  trx(fly1).units.absthetadiff_anglesub = parseunits('rad');

  % absphidiff for each of these measured closest flies
  idx = 1:trx(fly1).nframes-1;
  trx(fly1).absphidiff_center = absphidiff(sub2ind(size(absphidiff),trx(fly1).closestfly_center(idx),idx));
  trx(fly1).units.absphidiff_center = parseunits('rad');
  trx(fly1).absphidiff_nose2ell = absphidiff(sub2ind(size(absphidiff),trx(fly1).closestfly_nose2ell(idx),idx));
  trx(fly1).units.absphidiff_nose2ell = parseunits('rad');
  trx(fly1).absphidiff_ell2nose = absphidiff(sub2ind(size(absphidiff),trx(fly1).closestfly_ell2nose(idx),idx));
  trx(fly1).units.absphidiff_ell2nose = parseunits('rad');
  trx(fly1).absphidiff_nose2tail = absphidiff(sub2ind(size(absphidiff),trx(fly1).closestfly_nose2tail(idx),idx));
  trx(fly1).units.absphidiff_nose2tail = parseunits('rad');
  trx(fly1).absphidiff_tail2nose = absphidiff(sub2ind(size(absphidiff),trx(fly1).closestfly_tail2nose(idx),idx));
  trx(fly1).units.absphidiff_tail2nose = parseunits('rad');
  trx(fly1).absphidiff_anglesub = absphidiff(sub2ind(size(absphidiff),trx(fly1).closestfly_anglesub(idx),idx));
  trx(fly1).units.absphidiff_anglesub = parseunits('rad');

  % anglefrom1to2 for each of these measured closest flies
  idx = 1:trx(fly1).nframes;
  trx(fly1).absanglefrom1to2_center = absanglefrom1to2(sub2ind(size(absanglefrom1to2),trx(fly1).closestfly_center(idx),idx));
  trx(fly1).anglefrom1to2_center = anglefrom1to2(sub2ind(size(anglefrom1to2),trx(fly1).closestfly_center(idx),idx));
  trx(fly1).units.absanglefrom1to2_center = parseunits('rad');
  trx(fly1).units.anglefrom1to2_center = parseunits('rad');
  trx(fly1).absanglefrom1to2_nose2ell = absanglefrom1to2(sub2ind(size(absanglefrom1to2),trx(fly1).closestfly_nose2ell(idx),idx));
  trx(fly1).anglefrom1to2_nose2ell = anglefrom1to2(sub2ind(size(anglefrom1to2),trx(fly1).closestfly_nose2ell(idx),idx));
  trx(fly1).units.absanglefrom1to2_nose2ell = parseunits('rad');
  trx(fly1).units.anglefrom1to2_nose2ell = parseunits('rad');
  trx(fly1).absanglefrom1to2_ell2nose = absanglefrom1to2(sub2ind(size(absanglefrom1to2),trx(fly1).closestfly_ell2nose(idx),idx));
  trx(fly1).anglefrom1to2_ell2nose = anglefrom1to2(sub2ind(size(anglefrom1to2),trx(fly1).closestfly_ell2nose(idx),idx));
  trx(fly1).units.absanglefrom1to2_ell2nose = parseunits('rad');
  trx(fly1).units.anglefrom1to2_ell2nose = parseunits('rad');
  trx(fly1).absanglefrom1to2_nose2tail = absanglefrom1to2(sub2ind(size(absanglefrom1to2),trx(fly1).closestfly_nose2tail(idx),idx));
  trx(fly1).anglefrom1to2_nose2tail = anglefrom1to2(sub2ind(size(anglefrom1to2),trx(fly1).closestfly_nose2tail(idx),idx));
  trx(fly1).units.absanglefrom1to2_nose2tail = parseunits('rad');
  trx(fly1).units.anglefrom1to2_nose2tail = parseunits('rad');
  trx(fly1).absanglefrom1to2_tail2nose = absanglefrom1to2(sub2ind(size(absanglefrom1to2),trx(fly1).closestfly_tail2nose(idx),idx));
  trx(fly1).anglefrom1to2_tail2nose = anglefrom1to2(sub2ind(size(anglefrom1to2),trx(fly1).closestfly_tail2nose(idx),idx));
  trx(fly1).units.absanglefrom1to2_tail2nose = parseunits('rad');
  trx(fly1).units.anglefrom1to2_tail2nose = parseunits('rad');
  trx(fly1).absanglefrom1to2_anglesub = absanglefrom1to2(sub2ind(size(absanglefrom1to2),trx(fly1).closestfly_anglesub(idx),idx));
  trx(fly1).anglefrom1to2_anglesub = anglefrom1to2(sub2ind(size(anglefrom1to2),trx(fly1).closestfly_anglesub(idx),idx));
  trx(fly1).units.absanglefrom1to2_anglesub = parseunits('rad');
  trx(fly1).units.anglefrom1to2_anglesub = parseunits('rad');
  
  % anglefrom2to1 for each of these measured closest flies
  idx = 1:trx(fly1).nframes;
  trx(fly1).absanglefrom2to1_center = absanglefrom2to1(sub2ind(size(absanglefrom2to1),trx(fly1).closestfly_center(idx),idx));
  trx(fly1).anglefrom2to1_center = anglefrom2to1(sub2ind(size(anglefrom2to1),trx(fly1).closestfly_center(idx),idx));
  trx(fly1).units.absanglefrom2to1_center = parseunits('rad');
  trx(fly1).units.anglefrom2to1_center = parseunits('rad');
  trx(fly1).absanglefrom2to1_nose2ell = absanglefrom2to1(sub2ind(size(absanglefrom2to1),trx(fly1).closestfly_nose2ell(idx),idx));
  trx(fly1).anglefrom2to1_nose2ell = anglefrom2to1(sub2ind(size(anglefrom2to1),trx(fly1).closestfly_nose2ell(idx),idx));
  trx(fly1).units.absanglefrom2to1_nose2ell = parseunits('rad');
  trx(fly1).units.anglefrom2to1_nose2ell = parseunits('rad');
  trx(fly1).absanglefrom2to1_ell2nose = absanglefrom2to1(sub2ind(size(absanglefrom2to1),trx(fly1).closestfly_ell2nose(idx),idx));
  trx(fly1).anglefrom2to1_ell2nose = anglefrom2to1(sub2ind(size(anglefrom2to1),trx(fly1).closestfly_ell2nose(idx),idx));
  trx(fly1).units.absanglefrom2to1_ell2nose = parseunits('rad');
  trx(fly1).units.anglefrom2to1_ell2nose = parseunits('rad');
  trx(fly1).absanglefrom2to1_nose2tail = absanglefrom2to1(sub2ind(size(absanglefrom2to1),trx(fly1).closestfly_nose2tail(idx),idx));
  trx(fly1).anglefrom2to1_nose2tail = anglefrom2to1(sub2ind(size(anglefrom2to1),trx(fly1).closestfly_nose2tail(idx),idx));
  trx(fly1).units.absanglefrom2to1_nose2tail = parseunits('rad');
  trx(fly1).units.anglefrom2to1_nose2tail = parseunits('rad');
  trx(fly1).absanglefrom2to1_tail2nose = absanglefrom2to1(sub2ind(size(absanglefrom2to1),trx(fly1).closestfly_tail2nose(idx),idx));
  trx(fly1).anglefrom2to1_tail2nose = anglefrom2to1(sub2ind(size(anglefrom2to1),trx(fly1).closestfly_tail2nose(idx),idx));
  trx(fly1).units.absanglefrom2to1_tail2nose = parseunits('rad');
  trx(fly1).units.anglefrom2to1_tail2nose = parseunits('rad');
  trx(fly1).absanglefrom2to1_anglesub = absanglefrom2to1(sub2ind(size(absanglefrom2to1),trx(fly1).closestfly_anglesub(idx),idx));
  trx(fly1).anglefrom2to1_anglesub = anglefrom2to1(sub2ind(size(anglefrom2to1),trx(fly1).closestfly_anglesub(idx),idx));
  trx(fly1).units.absanglefrom2to1_anglesub = parseunits('rad');
  trx(fly1).units.anglefrom2to1_anglesub = parseunits('rad');

  % change in various parameters
  trx(fly1).units.ddcenter = parseunits('mm/s');
  trx(fly1).ddnose2ell = diff(trx(fly1).dnose2ell).*trx(fly1).fps;
  trx(fly1).units.ddnose2ell = parseunits('mm/s');
  trx(fly1).ddell2nose = diff(trx(fly1).dell2nose).*trx(fly1).fps;
  trx(fly1).units.ddell2nose = parseunits('mm/s');
  trx(fly1).ddnose2tail = diff(trx(fly1).dnose2tail).*trx(fly1).fps;
  trx(fly1).units.ddnose2tail = parseunits('mm/s');
  trx(fly1).ddtail2nose = diff(trx(fly1).dtail2nose).*trx(fly1).fps;
  trx(fly1).units.ddtail2nose = parseunits('mm/s');
  trx(fly1).danglesub = diff(trx(fly1).anglesub).*trx(fly1).fps;
  trx(fly1).units.danglesub = parseunits('rad/s');
  
end