function [x1,y1,x2,y2] = rfrac2center(trk,rfrac)

x1 = trk.x(1:end-1) + rfrac(1,:).*trk.a(1:end-1).*2.*cos(trk.theta(1:end-1)) - rfrac(2,:).*trk.b(1:end-1).*2.*sin(trk.theta(1:end-1));
y1 = trk.y(1:end-1) + rfrac(1,:).*trk.a(1:end-1).*2.*sin(trk.theta(1:end-1)) + rfrac(2,:).*trk.b(1:end-1).*2.*cos(trk.theta(1:end-1));
x2 = trk.x(2:end) + rfrac(1,:).*trk.a(2:end).*2.*cos(trk.theta(2:end)) - rfrac(2,:).*trk.b(2:end).*2.*sin(trk.theta(2:end));
y2 = trk.y(2:end) + rfrac(1,:).*trk.a(2:end).*2.*sin(trk.theta(2:end)) + rfrac(2,:).*trk.b(2:end).*2.*cos(trk.theta(2:end));
