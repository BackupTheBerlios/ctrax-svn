function newdata = process_data_crabwalks(data)

nflies = length(data);
for fly = 1:nflies,
  newdata(fly) = process_data_crabwalks_main(data(fly));
end

function data = process_data_crabwalks_main(data)

% location of tail
tailx = data.x_mm + 2*cos(-data.theta).*data.a_mm;
taily = data.y_mm + 2*sin(-data.theta).*data.a_mm;
% location of nose
nosex = data.x_mm + 2*cos(data.theta).*data.a_mm;
nosey = data.y_mm + 2*sin(data.theta).*data.a_mm;

dx = diff(tailx);
dy = diff(taily);

% project onto body coords
data.du_tail = dx.*cos(data.theta(1:end-1)) + dy.*sin(data.theta(1:end-1))*data.fps;
data.units.du_tail = parseunits('mm/s');
data.dv_tail = dx.*cos(data.theta(1:end-1)+pi/2) + dy.*sin(data.theta(1:end-1)+pi/2)*data.fps;
data.units.dv_tail = parseunits('mm/s');
data.absdu_tail = abs(data.du_tail);
data.units.absdu_tail = parseunits('mm/s');
data.absdv_tail = abs(data.dv_tail);
data.units.absdv_tail = parseunits('mm/s');

% compute the rotation of nose around mean tail location
meantailx = (tailx(1:end-1)+tailx(2:end))/2;
meantaily = (taily(1:end-1)+taily(2:end))/2;
anglenose1 = atan2(nosey(1:end-1)-meantaily,nosex(1:end-1)-meantailx);
anglenose2 = atan2(nosey(2:end)-meantaily,nosex(2:end)-meantailx);
data.dtheta_tail = modrange(anglenose2-anglenose1,-pi,pi)*data.fps;
data.units.dtheta_tail = parseunits('rad/s');
data.absdtheta_tail = abs(data.dtheta_tail);
data.units.absdtheta_tail = parseunits('rad/s');

% how sideways is the velocity direction?
phi = atan2(dy,dx);
data.phisideways = modrange(phi-data.theta(1:end-1),-pi/2,pi/2);
data.units.phisideways = parseunits('rad');
data.absphisideways = abs(data.phisideways);
data.units.absphisideways = parseunits('rad');
