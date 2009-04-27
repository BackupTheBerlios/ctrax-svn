function pairtrx = process_data_pairs(trx,fly1,t0,t1)

nflies = length(trx);

% make a trajectory for each pair
pairtrx = repmat(trx(fly1),[nflies,1]);

% overlapping frames
fns = fieldnames(trx);
for fly2 = 1:nflies,
  pairtrx(fly2).fly1 = fly1;
  t0curr = max(t0,trx(fly2).firstframe);
  t1curr = min(t1,trx(fly2).endframe);

  for i = 1:length(fns),
    fn = fns{i};
    i0 = trx(fly1).f2i(t0curr);
    i1 = trx(fly1).f2i(t1curr);
    if length(pairtrx(fly2).(fn)) == trx(fly1).nframes-1,
      pairtrx(fly2).(fn) = trx(fly1).(fn)(i0:i1-1);
    elseif length(pairtrx(fly2).(fn)) == trx(fly1).nframes,
      pairtrx(fly2).(fn) = trx(fly1).(fn)(i0:i1);
    end
  end
  
  pairtrx(fly2).firstframe = t0curr;
  pairtrx(fly2).endframe = t1curr;
  pairtrx(fly2).f2i = @(f) f - pairtrx(fly2).firstframe + 1;
  pairtrx(fly2).nframes = t1curr - t0curr + 1;
  
end

% get rid of fly1
fns = fieldnames(pairtrx);
for i = 1:length(fns),
  if length(pairtrx(fly1).(fns{i})) >= pairtrx(fly1).nframes - 2,
    pairtrx(fly1).(fns{i}) = [];
  end
end
pairtrx(fly1).firstframe = 1;
pairtrx(fly1).endframe = 0;
pairtrx(fly1).nframes = 0;

% location of the fly 1's nose
i0 = trx(fly1).f2i(t0);
i1 = trx(fly1).f2i(t1);
xnose = trx(fly1).x_mm(i0:i1) + 2*trx(fly1).a_mm(i0:i1).*cos(trx(fly1).theta(i0:i1));
ynose = trx(fly1).y_mm(i0:i1) + 2*trx(fly1).a_mm(i0:i1).*sin(trx(fly1).theta(i0:i1));

% units
units = trx(1).units;
units.distnose2ell = struct('num',{{'mm'}},'den',{{}});
units.dcenter = struct('num',{{'mm'}},'den',{{}});
units.magveldiff = struct('num',{{'mm'}},'den',{{'s'}});
units.veltoward = struct('num',{{'mm'}},'den',{{'s'}});
units.thetadiff = struct('num',{{'rad'}},'den',{{}});
units.absthetadiff = struct('num',{{'rad'}},'den',{{}});
units.phidiff = struct('num',{{'rad'}},'den',{{}});
units.absphidiff = struct('num',{{'rad'}},'den',{{}});
units.minvelmag = struct('num',{{'mm'}},'den',{{'s'}});
units.maxvelmag = struct('num',{{'mm'}},'den',{{'s'}});
units.anglefrom1to2 = struct('num',{{'rad'}},'den',{{}});
units.absanglefrom1to2 = struct('num',{{'rad'}},'den',{{}});
units.anglesub = struct('num',{{'rad'}},'den',{{}});
units.minabsanglefrom1to2 = struct('num',{{'rad'}},'den',{{}});
units.danglefrom1to2 = struct('num',{{'rad'}},'den',{{'s'}});
units.danglesub = struct('num',{{'rad'}},'den',{{'s'}});

for fly2 = 1:nflies,
  
  pairtrx(fly2).units = units;
  
  if fly1 == fly2,
    continue;
  end
  
  % statistics of pairs of flies
  i0 = trx(fly1).f2i(pairtrx(fly2).firstframe);
  i1 = trx(fly1).f2i(pairtrx(fly2).endframe);
  j0 = trx(fly2).f2i(pairtrx(fly2).firstframe);
  j1 = trx(fly2).f2i(pairtrx(fly2).endframe);
  
  % distance from fly1's nose to fly2
  pairtrx(fly2).distnose2ell = zeros(1,pairtrx(fly2).nframes);
  for t = pairtrx(fly2).firstframe:pairtrx(fly2).endframe,
    i_pair = pairtrx(fly2).f2i(t);
    i_nose1 = t - t0 + 1;
    i_fly2 = trx(fly2).f2i(t);
    pairtrx(fly2).distnose2ell(i_pair) = ...
      ellipsedist_hack(trx(fly2).x_mm(i_fly2),trx(fly2).y_mm(i_fly2),trx(fly2).a_mm(i_fly2),...
      trx(fly2).b_mm(i_fly2),trx(fly2).theta(i_fly2),xnose(i_nose1),ynose(i_nose1));
  end
  
  % magnitude of difference in velocity vectors
  pairtrx(fly2).magveldiff = sqrt( (diff(trx(fly1).x_mm(i0:i1))-diff(trx(fly2).x_mm(j0:j1))).^2 + ...
    (diff(trx(fly1).y_mm(i0:i1))-diff(trx(fly2).y_mm(j0:j1))).^2 ) *trx(fly1).fps;
  
  dx = trx(fly2).x_mm(j0:j1)-trx(fly1).x_mm(i0:i1);
  dy = trx(fly2).y_mm(j0:j1)-trx(fly1).y_mm(i0:i1);
  z = sqrt(dx.^2 + dy.^2);
  
  % distance between centers
  pairtrx(fly2).dcenter = z;
  
  % velocity in direction of other fly 
  dx = dx ./ z;
  dy = dy ./ z;
  pairtrx(fly2).veltoward = (dx(1:end-1).*diff(trx(fly1).x_mm(i0:i1)) + ...
    dy(1:end-1).*diff(trx(fly2).y_mm(j0:j1))).*trx(fly1).fps;
  
  % orientation of fly2 relative to orientation of fly1
  pairtrx(fly2).thetadiff = modrange(trx(fly2).theta(j0:j1) - trx(fly1).theta(i0:i1),-pi,pi);
  pairtrx(fly2).absthetadiff = abs(pairtrx(fly2).thetadiff);
  
  % velocity direction of fly2 relative to fly1's velocity direction
  pairtrx(fly2).phidiff = modrange(trx(fly2).phi(j0:j1-1)-trx(fly1).phi(i0:i1-1),-pi,pi);
  pairtrx(fly2).absphidiff = abs(pairtrx(fly2).phidiff);
  
  % minimum velocity magnitude of the two
  pairtrx(fly2).minvelmag = min(trx(fly1).velmag(i0:i1-1),trx(fly2).velmag(j0:j1-1));
  pairtrx(fly2).maxvelmag = max(trx(fly1).velmag(i0:i1-1),trx(fly2).velmag(j0:j1-1));
  
  % direction to fly2 from fly1
  pairtrx(fly2).anglefrom1to2 = modrange(atan2(dy,dx)-trx(fly1).theta(i0:i1),-pi,pi);
  pairtrx(fly2).absanglefrom1to2 = abs(pairtrx(fly2).anglefrom1to2);
  
  % angle subtended and minimum angle from fly1 to some point on fly2
  pairtrx(fly2).minabsanglefrom1to2 = zeros(1,pairtrx(fly2).nframes);
  pairtrx(fly2).anglesub = zeros(1,pairtrx(fly2).nframes);
  for t = pairtrx(fly2).firstframe:pairtrx(fly2).endframe,
    i_pair = pairtrx(fly2).f2i(t);
    i1 = trx(fly1).f2i(t);
    i2 = trx(fly2).f2i(t);
    % compute angles to fly2
    [pairtrx(fly2).anglesub(i_pair),psi1,psi2] = tangentangles(trx(fly1).x(i1),trx(fly1).y(i1),...
      trx(fly1).a(i1)*2,trx(fly1).b(i1)*2,trx(fly1).theta(i1),...
      trx(fly2).x(i2),trx(fly2).y(i2),trx(fly2).a(i2)*2,...
      trx(fly2).b(i2)*2,trx(fly2).theta(i2));
    if sign(psi1) ~= sign(psi2),
      pairtrx(fly2).minabsanglefrom1to2(i_pair) = 0;
    else
      pairtrx(fly2).minabsanglefrom1to2(i_pair) = min(abs(psi1),abs(psi2));
    end
  end
  
  pairtrx(fly2).danglesub = diff(pairtrx(fly2).anglesub) * trx(fly1).fps;
  pairtrx(fly2).danglefrom1to2 = modrange(diff(pairtrx(fly2).anglefrom1to2),-pi,pi) * trx(fly1).fps;
  
end