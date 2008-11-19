function [isfore,dfore,xpred,ypred,thetapred,r0,r1,c0,c1,im] = FixBgSub(fly,f,handles)

trk = handles.trx(fly);
boxrad = handles.maxjump;

i = trk.f2i(f);
x1 = trk.x(i-1);
y1 = trk.y(i-1);
theta1 = trk.theta(i-1);
if i == 2,
  xpred = x1;
  ypred = y1;
  thetapred = theta1;
else
  x2 = trk.x(i-2);
  y2 = trk.y(i-2);
  theta2 = trk.theta(i-2);
  [xpred,ypred,thetapred] = cvpred(x2,y2,theta2,x1,y1,theta1);
end

r0 = max(floor(ypred-boxrad),1); r1 = min(ceil(ypred+boxrad),handles.nr);
c0 = max(floor(xpred-boxrad),1); c1 = min(ceil(xpred+boxrad),handles.nc);
im = handles.readframe(f);
im = im(r0:r1,c0:c1);

bg = handles.bgcurr(r0:r1,c0:c1);
dfore = im - bg;
if handles.lighterthanbg == 1
  dfore = max(dfore,0);
elseif handles.lighterthanbg == -1
  dfore = max(-dfore,0);
else
  dfore = abs(dfore);
end
isfore = dfore >= handles.bgthresh;
se = strel('disk',1);
isfore = imclose(isfore,se);
isfore = imopen(isfore,se);