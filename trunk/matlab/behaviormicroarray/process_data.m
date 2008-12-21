function data = process_data(data,matname,moviename,annname)

nflies = length(data);

if ~exist('matname','var')
  matname = '';
end
if ~exist('moviename','var')
  moviename = '';
end
if ~exist('annname','var')
  annname = '';
end

% set moviename
for fly = 1:nflies,
  data(fly).moviename = moviename;
  data(fly).matname = matname;
  data(fly).annname = annname;
end

if isfield(data,'xpred'),
  data = rmfield(data,'xpred');
end
if isfield(data,'ypred'),
  data = rmfield(data,'ypred');
end
if isfield(data,'thetapred'),
  data = rmfield(data,'thetapred');
end

% read arena dimensions
if ~isempty(annname) && exist(annname,'file'),
  [arena.x,arena.y,arena.r] = arena_params(annname,moviename);
else
  arena.x = nan; arena.y = nan; arena.r = nan;
end
for fly = 1:nflies,
  data(fly).arena = arena;
end

thetafil = [1     4     6     4     1]/16;

% compute velocities in the canonical coordinates of the fly
for fly = 1:nflies,
  
  % change in body orientation
  data(fly).dtheta = modrange(diff(data(fly).theta),-pi,pi)*data(fly).fps;
  
  % change in center position
  dx = diff(data(fly).x_mm);
  dy = diff(data(fly).y_mm);
  
  % forward motion of body center
  data(fly).du_ctr = dx.*cos(data(fly).theta(1:end-1)) + dy.*sin(data(fly).theta(1:end-1))*data(fly).fps;
  % sideways motion of body center
  data(fly).dv_ctr = dx.*cos(data(fly).theta(1:end-1)+pi/2) + dy.*sin(data(fly).theta(1:end-1)+pi/2)*data(fly).fps;
  
  % find the center of rotation
  [data(fly).corfrac,data(fly).corisonfly] = center_of_rotation2(data(fly),false);
  [x_cor_curr,y_cor_curr,x_cor_next,y_cor_next] = rfrac2center(data(fly),data(fly).corfrac);

  % change in center of rotation
  dx_cor = x_cor_next - x_cor_curr;
  dy_cor = y_cor_next - y_cor_curr;
  
  % forward motion of center of rotation
  data(fly).du_cor = dx_cor.*cos(data(fly).theta(1:end-1)) + dy_cor.*sin(data(fly).theta(1:end-1))*data(fly).fps;
  % sideways motion of body center
  data(fly).dv_cor = dx_cor.*cos(data(fly).theta(1:end-1)+pi/2) + dy_cor.*sin(data(fly).theta(1:end-1)+pi/2)*data(fly).fps;
  
  % magnitude of velocity
  data(fly).velmag_ctr = sqrt(dx.^2 + dy.^2)*data(fly).fps;
  data(fly).velmag = sqrt(dx_cor.^2 + dy_cor.^2)*data(fly).fps;
  
  % acceleration magnitude
  tmp = sqrt(diff(dx).^2 + diff(dy).^2)*data(fly).fps^2;
  data(fly).accmag = [0,tmp];
  
  % flipped sign dv, dtheta
  data(fly).signdtheta = sign(data(fly).dtheta);
  data(fly).absdv_cor = abs(data(fly).dv_cor);
  data(fly).flipdv_cor = data(fly).dv_cor.*data(fly).signdtheta;
  %data(fly).realabsdv_cor = abs(data(fly).dv_cor);
  data(fly).absdtheta = abs(data(fly).dtheta);
  data(fly).d2theta = [0,modrange(diff(data(fly).dtheta),-pi,pi)]*data(fly).fps;
  data(fly).absd2theta = abs(data(fly).d2theta);
  
  % smoothed orientation
  data(fly).smooththeta = imfilter(unwrap(data(fly).theta),thetafil);
  data(fly).smoothdtheta = diff(data(fly).smooththeta)*data(fly).fps;
  data(fly).smooththeta = modrange(data(fly).smooththeta,-pi,pi);
  data(fly).abssmoothdtheta = abs(data(fly).smoothdtheta);
  data(fly).smoothd2theta = [0,modrange(diff(data(fly).smoothdtheta),-pi,pi)]*data(fly).fps;
  data(fly).abssmoothd2theta = abs(data(fly).smoothd2theta);

  data(fly).f2i = @(f) f - data(fly).firstframe + 1;
  
  % velocity direction
  dy1 = [data(fly).y(2)-data(fly).y(1),(data(fly).y(3:end)-data(fly).y(1:end-2))/2,data(fly).y(end)-data(fly).y(end-1)];
  dx1 = [data(fly).x(2)-data(fly).x(1),(data(fly).x(3:end)-data(fly).x(1:end-2))/2,data(fly).x(end)-data(fly).x(end-1)];
  data(fly).phi = atan2(dy1,dx1);
  
  % difference between velocity direction and orientation
  data(fly).yaw = modrange(data(fly).phi - data(fly).theta,-pi,pi);
  data(fly).absyaw = abs(data(fly).yaw);
end