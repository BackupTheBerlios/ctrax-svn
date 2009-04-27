%% simple_diagnostics

% this script inputs a mat file containing trajectories, and histograms
% various simple properties. 

%% set up path, if necessary

setuppath;

%% set all defaults

matname = '';
matpath = '';

%% load settings

pathtocomputeperframestats = which('simple_diagnostics');
savedsettingsfile = strrep(pathtocomputeperframestats,'simple_diagnostics.m','.simplediagnosticssrc.mat');
if exist(savedsettingsfile,'file')
  load(savedsettingsfile);
end

%% choose a mat file to analyze
fprintf('Choose mat file to analyze.\n');
matname = [matpath,matname];
[matname,matpath] = uigetfile('*.mat','Choose mat file to analyze',matname);
if isnumeric(matname) && matname == 0,
  return;
end
fprintf('Matfile: %s%s\n\n',matpath,matname);

if exist(savedsettingsfile,'file'),
  save('-append',savedsettingsfile,'matname','matpath');
else
  save(savedsettingsfile,'matname','matpath');
end
matnameonly = matname;

%% get conversion to mm, seconds

[convertunits_succeeded,matname] = convert_units_f('isautomatic',true,'matname',matnameonly,'matpath',matpath);

%% load in the data

[trx,matname,loadsucceeded] = load_tracks(matname);
if ~loadsucceeded,
  msgbox('Could not load trx from file %s\n',matname);
  return;
end

%% plot all the trajectories

figure(1);

nflies = length(trx);

% put all the flies on one figure
n1 = round(sqrt(nflies));
n2 = ceil(nflies/n1);
hax = createsubplots(n1,n2,[[.05,.01];[.05,.05]]);
hax = reshape(hax,[n1,n2])'; % column major
if n1*n2 > nflies,
  delete(hax(nflies+1:end));
  hax = hax(1:nflies);
end

% get bounds for all flies
minx = min([trx.x_mm]);
maxx = max([trx.x_mm]);
miny = min([trx.y_mm]);
maxy = max([trx.y_mm]);
dx = maxx - minx;
dy = maxy - miny;

for fly = 1:nflies,
  
  % plot just the trajectory
  axes(hax(fly));
  plot(trx(fly).x_mm,trx(fly).y_mm,'k.-','markersize',3,'linewidth',.5);
  title(sprintf('Fly %d, frames %d to %d',fly,trx(fly).firstframe,trx(fly).endframe));
  axis([minx-.025*dx,maxx+.025*dx,miny-.025*dy,maxy+.025*dy]);
  axis equal;
  
  % onlt put an x-axis on the lowest plots
  [c,r] = ind2sub([n2,n1],fly);
  if r ~= n1 && r*n2+c <= nflies,
    set(hax(fly),'xticklabel',{});
  end
  if c ~= 1,
    set(hax(fly),'yticklabel',{});
  end
end

linkaxes(hax);

%% histogram position in the arena for all flies

figure(2); 
clf;

nbinsx = 100;
nbinsy = max(1,round(nbinsx*dy/dx));
edgesx = linspace(minx,maxx,nbinsx+1);
edgesy = linspace(miny,maxy,nbinsy+1);
centersx = (edgesx(1:end-1)+edgesx(2:end))/2;
centersy = (edgesy(1:end-1)+edgesy(2:end))/2;
counts = hist3([[trx.x_mm];[trx.y_mm]]',{centersx,centersy});
freq_position = counts / sum(counts(:));

hax = createsubplots(1,2,.05);
axes(hax(1));
imagesc([centersx(1),centersx(end)],[centersy(1),centersy(end)],freq_position');
axis image;
title('Position heat map, frequency');
colorbar;
axes(hax(2));
imagesc([centersx(1),centersx(end)],[centersy(1),centersy(end)],log(freq_position'));
title('Position heat map, log frequency');
axis image;
colorbar;

%% histogram speed

figure(3);
clf;
hax = createsubplots(1,3,[.05,.1]);

% compute speed for all flies
speed = [];
for fly = 1:nflies,
  speedcurr = sqrt(diff(trx(fly).x_mm).^2 + diff(trx(fly).y_mm).^2)*trx(fly).fps;
  speed = [speed,speedcurr];
end

% choose bins
nbins = min(100,length(speed)/20);
prctlastbin = 1;
edges = [linspace(0,prctile(speed,100-prctlastbin),nbins),max(speed)];
centers = (edges(1:end-1)+edges(2:end))/2;

% histogram for all flies
counts = histc(speed,edges);
counts(end-1) = counts(end-1) + counts(end);
counts = counts(1:end-1);
freq_speed_allflies = counts / sum(counts);

% do each fly individually
freq_speed_perfly = zeros(nflies,nbins);
for fly = 1:nflies,
  speedcurr = sqrt(diff(trx(fly).x_mm).^2 + diff(trx(fly).y_mm).^2)*trx(fly).fps;
  counts = histc(speedcurr,edges);
  counts(end-1) = counts(end-1) + counts(end);
  counts = counts(1:end-1);
  freq_speed_perfly(fly,:) = counts / sum(counts);
end

axes(hax(1));
xplot = centers;
nskip = nbins/10;
xplot(end) = (1+nskip)*centers(end-1)-nskip*centers(end-2);
plot(xplot,freq_speed_perfly,'-','linewidth',.25);
hold on;
plot(xplot,freq_speed_allflies,'k.-','linewidth',5);
dx = xplot(end)-xplot(1);
dy = max(freq_speed_allflies);
ax = [-dx/20,xplot(end)+dx/20,0,dy+dy/20];
axis(ax);

% reset xticks
xtick = get(gca,'xtick');
xticklabel = cellstr(get(gca,'xticklabel'));
xticklabel(xtick > xplot(end-1)) = [];
xtick(xtick > xplot(end-1)) = [];
xtick(end+1) = xplot(end);
xticklabel{end+1} = sprintf('> %.1f',edges(end-1));
set(gca,'xtick',xtick,'xticklabel',xticklabel);
xlabel('speed (mm/s)');
ylabel('frequency');
title('histogram of speed for all flies');

axes(hax(2));
semilogy(xplot,freq_speed_perfly,'-','linewidth',.25);
hold on;
h2 = semilogy(xplot,freq_speed_allflies,'k.-','linewidth',5);
% reset xticks
xtick = get(gca,'xtick');
xticklabel = cellstr(get(gca,'xticklabel'));
xticklabel(xtick > xplot(end-1)) = [];
xtick(xtick > xplot(end-1)) = [];
xtick(end+1) = xplot(end);
xticklabel{end+1} = sprintf('> %.1f',edges(end-1));
set(gca,'xtick',xtick,'xticklabel',xticklabel);
ax = [-dx/20,xplot(end)+dx/20,min(freq_speed_perfly(:)),dy+dy/20];
axis(ax);
title('log frequency plot');
xlabel('speed (mm/s)');
legend(h2,'All flies');

axes(hax(3));
loglog(xplot,freq_speed_perfly,'-','linewidth',.25);
hold on;
loglog(xplot,freq_speed_allflies,'k.-','linewidth',5);
% reset xticks
xtick = get(gca,'xtick');
xticklabel = cellstr(get(gca,'xticklabel'));
xticklabel(xtick > xplot(end-1)) = [];
xtick(xtick > xplot(end-1)) = [];
xtick(end+1) = xplot(end);
xticklabel{end+1} = sprintf('> %.1f',edges(end-1));
set(gca,'xtick',xtick,'xticklabel',xticklabel);
ax = [centers(1),xplot(end)+dx/20,min(freq_speed_perfly(:)),dy+dy/20];
axis(ax);
title('log-log frequency plot');
xlabel('speed (mm/s)');

%% histogram change in orientation

figure(4);
clf;
hax = createsubplots(1,3,[.05,.1]);

% compute speed for all flies
speed = [];
for fly = 1:nflies,
  speedcurr = sqrt(modrange(diff(trx(fly).theta),-pi,pi).^2 + ...
    modrange(diff(trx(fly).y_mm),-pi,pi).^2)*trx(fly).fps*180/pi;
  speed = [speed,speedcurr];
end

% choose bins
nbins = min(100,length(speed)/20);
prctlastbin = 1;
edges = [linspace(0,prctile(speed,100-prctlastbin),nbins),max(speed)];
centers = (edges(1:end-1)+edges(2:end))/2;

% histogram for all flies
counts = histc(speed,edges);
counts(end-1) = counts(end-1) + counts(end);
counts = counts(1:end-1);
freq_angularspeed_allflies = counts / sum(counts);

% do each fly individually
freq_angularspeed_perfly = zeros(nflies,nbins);
for fly = 1:nflies,
  speedcurr = sqrt(modrange(diff(trx(fly).theta),-pi,pi).^2 + ...
    modrange(diff(trx(fly).y_mm),-pi,pi).^2)*trx(fly).fps*180/pi;
  counts = histc(speedcurr,edges);
  counts(end-1) = counts(end-1) + counts(end);
  counts = counts(1:end-1);
  freq_angularspeed_perfly(fly,:) = counts / sum(counts);
end

axes(hax(1));
xplot = centers;
nskip = nbins/10;
xplot(end) = (1+nskip)*centers(end-1)-nskip*centers(end-2);
plot(xplot,freq_angularspeed_perfly,'-','linewidth',.25);
hold on;
plot(xplot,freq_angularspeed_allflies,'k.-','linewidth',5);
dx = xplot(end)-xplot(1);
dy = max(freq_speed_allflies);
ax = [-dx/20,xplot(end)+dx/20,0,dy+dy/20];
axis(ax);

% reset xticks
xtick = get(gca,'xtick');
xticklabel = cellstr(get(gca,'xticklabel'));
xticklabel(xtick > xplot(end-1)) = [];
xtick(xtick > xplot(end-1)) = [];
xtick(end+1) = xplot(end);
xticklabel{end+1} = sprintf('> %.1f',edges(end-1));
set(gca,'xtick',xtick,'xticklabel',xticklabel);
xlabel('angular speed (deg/s)');
ylabel('frequency');
title('histogram of angular speed for all flies');

axes(hax(2));
semilogy(xplot,freq_angularspeed_perfly,'-','linewidth',.25);
hold on;
h2 = semilogy(xplot,freq_angularspeed_allflies,'k.-','linewidth',5);
% reset xticks
xtick = get(gca,'xtick');
xticklabel = cellstr(get(gca,'xticklabel'));
xticklabel(xtick > xplot(end-1)) = [];
xtick(xtick > xplot(end-1)) = [];
xtick(end+1) = xplot(end);
xticklabel{end+1} = sprintf('> %.1f',edges(end-1));
set(gca,'xtick',xtick,'xticklabel',xticklabel);
ax = [-dx/20,xplot(end)+dx/20,min(freq_angularspeed_perfly(:)),dy+dy/20];
axis(ax);
title('log frequency plot');
xlabel('angular speed (deg/s)');
legend(h2,'All flies');

axes(hax(3));
loglog(xplot,freq_angularspeed_perfly,'-','linewidth',.25);
hold on;
loglog(xplot,freq_angularspeed_allflies,'k.-','linewidth',5);
% reset xticks
xtick = get(gca,'xtick');
xticklabel = cellstr(get(gca,'xticklabel'));
xticklabel(xtick > xplot(end-1)) = [];
xtick(xtick > xplot(end-1)) = [];
xtick(end+1) = xplot(end);
xticklabel{end+1} = sprintf('> %.1f',edges(end-1));
set(gca,'xtick',xtick,'xticklabel',xticklabel);
ax = [centers(1),xplot(end)+dx/20,min(freq_angularspeed_perfly(:)),dy+dy/20];
axis(ax);
title('log-log frequency plot');
xlabel('angular speed (deg/s)');