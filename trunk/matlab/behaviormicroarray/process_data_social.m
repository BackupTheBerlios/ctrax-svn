function trk = process_data_social(trx,fly1,t0,t1)

fov = pi;

if ~exist('t0','var'),
  t0 = trx(fly1).firstframe;
end
if ~exist('t1','var'),
  t1 = trx(fly1).endframe;
end

nflies = length(trx);
t0 = max(t0,trx(fly1).firstframe);
t1 = min(t1,trx(fly1).endframe);

trk = GetPartOfTrack(trx(fly1),t0,t1);

% don't get last frame
t1 = t1 - 1;
nframes = max(0,t1 - t0 + 1);

% allocate
trk.dcenter = nan(nflies,nframes);
trk.units.dcenter = parseunits('mm');
trk.dnose2ell = nan(nflies,nframes);
trk.units.dnose2ell = parseunits('mm');
trk.dell2nose = nan(nflies,nframes);
trk.units.dell2nose = parseunits('mm');
trk.magveldiff = nan(nflies,nframes);
trk.units.magveldiff = parseunits('mm/s');
trk.veltoward = nan(nflies,nframes);
trk.units.veltoward = parseunits('mm/s');
trk.thetadiff = nan(nflies,nframes);
trk.units.thetadiff = parseunits('rad');
trk.phidiff = nan(nflies,nframes);
trk.units.phidiff = parseunits('rad');
trk.minvelmag = nan(nflies,nframes);
trk.units.minvelmag = parseunits('mm/s');
trk.maxvelmag = nan(nflies,nframes);
trk.units.maxvelmag = parseunits('mm/s');
trk.anglefrom2to1 = nan(nflies,nframes);
trk.units.anglefrom2to1 = parseunits('rad');
trk.anglesub = nan(nflies,nframes);
trk.units.anglesub = parseunits('rad');

% position of nose
xnose = trx(fly1).x + 2*trx(fly1).a.*cos(trx(fly1).theta);
ynose = trx(fly1).y + 2*trx(fly1).a.*sin(trx(fly1).theta);

for fly2 = 1:nflies,
  
  if fly1 == fly2,
    continue;
  end

  % get start and end frames of overlap
  t0curr = max(t0,trx(fly2).firstframe);
  t1curr = min(t1,trx(fly2).endframe);
  if t1curr < t0curr, continue; end
  offi = trx(fly1).firstframe-1;
  offj = trx(fly2).firstframe-1;
  offk = trk.firstframe-1;
  i0 = t0curr - offi;
  i1 = t1curr - offi;
  j0 = t0curr - offj;
  j1 = t1curr - offj;
  k0 = t0curr - offk;
  k1 = t1curr - offk;
    
  % distance from fly1's nose to fly2
  for t = t0curr:t1curr,
    i = t - offi;
    j = t - offj;
    k = t - offk;
    trk.dnose2ell(fly2,k) = ...
      ellipsedist_hack(trx(fly2).x(j),trx(fly2).y(j),...
      trx(fly2).a(j),trx(fly2).b(j),trx(fly2).theta(j),...
      xnose(i),ynose(i));
  end % end loop over frame number
    
  % distance from fly2's nose to fly1
  for t = t0curr:t1curr,
    i = t - offi;
    j = t - offj;
    k = t - offk;
    
    xnose2 = trx(fly2).x(j) + 2*trx(fly2).a(j).*cos(trx(fly2).theta(j));
    ynose2 = trx(fly2).y(j) + 2*trx(fly2).a(j).*sin(trx(fly2).theta(j));
    trk.dell2nose(fly2,k) = ...
      ellipsedist_hack(trx(fly1).x(i),trx(fly1).y(i),...
      trx(fly1).a(i),trx(fly1).b(i),trx(fly1).theta(i),...
      xnose2,ynose2);
  end
  
  % magnitude of difference in velocity vectors
  trk.magveldiff(fly2,k0:k1-1) = ...
    sqrt( (diff(trx(fly1).x(i0:i1))-diff(trx(fly2).x(j0:j1))).^2 + ...
    (diff(trx(fly1).y(i0:i1))-diff(trx(fly2).y(j0:j1))).^2 );
    
  % direction to other fly
  dx = trx(fly2).x(j0:j1)-trx(fly1).x(i0:i1);
  dy = trx(fly2).y(j0:j1)-trx(fly1).y(i0:i1);
  z = sqrt(dx.^2 + dy.^2);
  dx = dx ./ z;
  dy = dy ./ z;
  
  % distance between centers
  trk.dcenter(fly2,k0:k1) = z;  

  % velocity in direction of other fly  
  trk.veltoward(fly2,k0:k1-1) = dx(1:end-1).*diff(trx(fly1).x(i0:i1)) + ...
    dy(1:end-1).*diff(trx(fly2).y(j0:j1));
  
  % orientation of fly2 relative to orientation of fly1
  trk.thetadiff(fly2,k0:k1) = modrange(trx(fly2).theta(j0:j1) - trx(fly1).theta(i0:i1),-pi,pi);

  % velocity direction of fly2 relative to fly1's velocity direction
  trk.phidiff(fly2,k0:k1-1) = modrange(trx(fly2).phi(j0:j1-1)-trx(fly1).phi(i0:i1-1),-pi,pi);
  
  % minimum velocity magnitude of the two
  trk.minvelmag(fly2,k0:k1-1) = min(trx(fly1).velmag(i0:i1-1),trx(fly2).velmag(j0:j1-1));
  trk.maxvelmag(fly2,k0:k1-1) = max(trx(fly1).velmag(i0:i1-1),trx(fly2).velmag(j0:j1-1));

  % direction to fly2 from fly1
  trk.anglefrom2to1(fly2,k0:k1) = modrange(atan2(dy,dx)-trx(fly1).theta(i0:i1),-pi,pi);

  % angle of fly1's vision subtended by fly2
  for t = t0:t1,
    i = t - offi;
    j = t - offj;
    k = t - offk;
    trk.anglesub(fly2,k) = anglesubtended(...
      trx(fly1).x(i),trx(fly1).y(i),trx(fly1).a(i),trx(fly1).b(i),trx(fly1).theta(i),...
      trx(fly2).x(j),trx(fly2).y(j),trx(fly2).a(j),trx(fly2).b(j),trx(fly2).theta(j),...
      fov);
  end
  
end

trk.absthetadiff = abs(trk.thetadiff);
trk.units.absthetadiff = parseunits('rad');
trk.absphidiff = abs(trk.phidiff);
trk.units.absphidiff = parseunits('rad');
trk.absanglefrom2to1 = abs(trk.anglefrom2to1);
trk.units.absanglefrom2to1 = parseunits('rad');

trk.dcenter = trk.dcenter / trk.pxpermm;
trk.dnose2ell = trk.dnose2ell / trk.pxpermm;
trk.dell2nose = trk.dell2nose / trk.pxpermm;
trk.magveldiff = trk.magveldiff / trk.pxpermm * trk.fps;
trk.veltoward = trk.veltoward / trk.pxpermm * trk.fps;
trk.minvelmag = trk.minvelmag / trk.pxpermm * trk.fps;
trk.maxvelmag = trk.maxvelmag / trk.pxpermm * trk.fps;