function data = process_data_arena(data)

nflies = length(data);
for i = 1:nflies,
  data(i).dist2wall = max(0,data(i).arena.r-sqrt((data(i).arena.x-data(i).x).^2 + (data(i).arena.y-data(i).y).^2));
  angle2wall = atan2(data(i).y-data(i).arena.y,data(i).x-data(i).arena.x);
  % orientation relative to angle to wall
  data(i).theta2wall = modrange(data(i).theta-angle2wall,-pi,pi);
  data(i).ddist2wall = diff(data(i).dist2wall);
  data(i).absdangle2wall = abs(modrange(diff(data(i).theta2wall),-pi,pi));
end
