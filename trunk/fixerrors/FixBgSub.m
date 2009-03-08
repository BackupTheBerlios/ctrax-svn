function [isfore,dfore,xpred,ypred,thetapred,r0,r1,c0,c1,im] = FixBgSub(flies,f,handles)

trx = handles.trx(flies);
nflies = length(flies);
boxrad = handles.maxjump;

xpred = zeros(1,nflies);
ypred = zeros(1,nflies);
thetapred = zeros(1,nflies);
for j = 1:nflies,
  i = trx(j).f2i(f);
  x1 = trx(j).x(i-1);
  y1 = trx(j).y(i-1);
  theta1 = trx(j).theta(i-1);
  if i == 2,
    xpred(j) = x1;
    ypred(j) = y1;
    thetapred(j) = theta1;
  else
    x2 = trx(j).x(i-2);
    y2 = trx(j).y(i-2);
    theta2 = trx(j).theta(i-2);
    [xpred(j),ypred(j),thetapred(j)] = cvpred(x2,y2,theta2,x1,y1,theta1);
  end
end

r0 = max(floor(min(ypred)-boxrad),1); r1 = min(ceil(max(ypred)+boxrad),handles.nr);
c0 = max(floor(min(xpred)-boxrad),1); c1 = min(ceil(max(xpred)+boxrad),handles.nc);
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