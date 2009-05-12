function [trk,f] = FixTrackFly(fly,f0,f1,handles)

trk = handles.trx(fly);
boxrad = handles.maxjump;
for f = f0+1:f1
  
  i = trk.f2i(f);
  [isfore,xpred,ypred,thetapred,r0,r1,c0,c1,im] = FixBgSub(fly,f,handles);

  [cc,ncc] = bwlabel(isfore);
  isdeleted = [];
  for fly2 = 1:handles.nflies,
    if fly2 == fly, continue; end
    i2 = handles.trx(fly2).f2i(f);
    bw = ellipsepixels([handles.trx(fly2).x(i2),handles.trx(fly2).y(i2),...
      handles.trx(fly2).a(i2)*4,handles.trx(fly2).b(i2)*4,handles.trx(fly2).theta(i2)],...
      [r0,r1,c0,c1]);
    for j = 1:ncc,
      if ismember(j,isdeleted), continue; end
      fracoverlap = nnz((cc==j) & bw) / nnz(cc==j);
      if fracoverlap > .5
        isfore(cc==j) = false;
        isdeleted(end+1) = j;
      end
    end
  end
  % choose the closest connected component
  if ~any(isfore(:)),
    msgbox('Frame %d: Could not find the selected fly. Quitting',f);
    return;
  end
  [cc,ncc] = bwlabel(isfore);
  xfit = zeros(1,ncc);
  yfit = zeros(1,ncc);
  thetafit = zeros(1,ncc);
  afit = zeros(1,ncc);
  bfit = zeros(1,ncc);
  for j = 1:ncc,
    [y,x] = find(cc==j);
    w = dfore(cc==j);
    [mu,S] = weighted_mean_cov([x;y]',w(:));
    xfit(j) = mu(1);
    yfit(j) = mu(2);
    [afit(j),bfit(j),thetafit(j)] = cov2ell(S);
  end
  afit = afit / 2;
  bfit = bfit / 2;
  
  xfit = xfit - c0 + 1;
  yfit = yfit - r0 + 1;
  
  if ncc == 1,
    trk.x(i) = xfit;
    trk.y(i) = yfit;
    trk.theta(i) = thetafit;
    trk.a(i) = afit;
    trk.b(i) = bfit;
  else
    err = (xpred - xfit).^2 + (ypred - yfit).^2 + handles.ang_dist_wt*(modrange(theta - thetafit,-pi/2,pi/2)).^2;
    j = argmin(err);
    trk.x(i) = xfit(j);
    trk.y(i) = yfit(j);
    trk.theta(i) = thetafit(j);
    trk.a(i) = afit(j);
    trk.b(i) = bfit(j);
  end
  dtheta = modrange(trk.theta(i)-trk.theta(i-1),-pi/2,pi/2);
  trk.theta(i) = trk.theta(i-1)+dtheta;
  
end
