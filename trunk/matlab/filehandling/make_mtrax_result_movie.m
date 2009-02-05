function make_mtrax_result_movie(moviename,trx,aviname,varargin)

[colors,zoomflies,nzoomr,nzoomc,boxradius,taillength,bitrate,NFRAMESBUFFER,fps,maxnframes] = ...
  myparse(varargin,'colors',[],'zoomflies',[],'nzoomr',5,'nzoomc',3,...
  'boxradius',20,'taillength',100,'bitrate',20000,'nframesbuffer',2000,'fps',20,...
  'maxnframes',inf);

addpath /home/kristin/FLIES/code/exploremtraxresults;

[readframe,nframes,fid] = get_readframe_fcn(moviename);
nframes = min(nframes,maxnframes);
im = readframe(1);
[nr,nc] = size(im);
nids = length(trx);

if isempty(colors),
  colors = jet(nids);
elseif size(colors,1) ~= nids,
  colors = colors(modrange(0:nids-1,size(colors,1))+1,:);
end

% choose some random flies to zoom in on
nzoom = nzoomr*nzoomc;
if isempty(zoomflies),
  nframesperfly = getstructarrayfield(trx,'nframes');
  fliesmaybeplot = find(nframesperfly > 1000);
  if length(fliesmaybeplot) <= nzoom,
    if nids >= nzoom,
      [tmp,zoomflies] = sort(-nframesperfly);
      zoomflies = zoomflies(1:nzoom);
    else
      zoomflies = [1:nids,nan(1,nzoom-length(fliesmaybeplot))];
      warning('Not enough flies to plot');
    end
  else
    zoomflies = sort(randsample(fliesmaybeplot,nzoom));
  end
end
zoomflies = reshape(zoomflies,[nzoomr,nzoomc]);
rowszoom = floor(nr/nzoomr);

figure(1);
clf;
hold on;
hax = gca;
set(hax,'position',[0,0,1,1]);
axis off;

% corners of zoom boxes in plotted image coords
x0 = nc+(0:nzoomc-1)*rowszoom+1;
y0 = (0:nzoomr-1)*rowszoom+1;
x1 = x0 + rowszoom - 1;
y1 = y0 + rowszoom - 1;

% relative frame offset
nframesoff = getstructarrayfield(trx,'firstframe') - 1;

% pre-allocate
himzoom = zeros(nzoomr,nzoomc);
htail = zeros(1,nids);
htri = zeros(1,nids);
scalefactor = rowszoom / (2*boxradius+1);
hzoom = zeros(nzoomr,nzoomc);

% for compression
isfirst = true;

for frame = 1:nframes,
  fprintf('frame %d\n',frame);
  
  % relative frame
  idx = frame - nframesoff;

  isalive = frame >= getstructarrayfield(trx,'firstframe') & ...
    frame <= getstructarrayfield(trx,'endframe');
  
  % draw the unzoomed image
  im = uint8(readframe(frame));
  if frame == 1,
    him = image([1,nc],[1,nr],repmat(im,[1,1,3]));
    axis image;
    axis([.5,x1(end)+.5,.5,y1(end)+.5]);
    axis off;
  else
    set(him,'cdata',repmat(im,[1,1,3]));
  end
  
  % draw the zoomed image
  for i = 1:nzoomr,
    for j = 1:nzoomc,
      fly = zoomflies(i,j);
      
      % fly not visible?
      if isnan(fly) || ~isalive(fly),
        if frame == 1,
          himzoom(i,j) = imagesc([x0(j),x1(j)],[y0(i),y1(i)],123);
        else
          set(himzoom(i,j),'cdata',repmat(uint8(123),[1,1,3]));
        end
        continue;
      end
      
      % grab a box around (x,y)
      x = round(trx(fly).x(idx(fly)));
      y = round(trx(fly).y(idx(fly)));
      boxradx1 = min(boxradius,x-1);
      boxradx2 = min(boxradius,size(im,2)-x);
      boxrady1 = min(boxradius,y-1);
      boxrady2 = min(boxradius,size(im,1)-y);
      box = uint8(zeros(2*boxradius+1));
      box(boxradius+1-boxrady1:boxradius+1+boxrady2,...
        boxradius+1-boxradx1:boxradius+1+boxradx2) = ...
        im(y-boxrady1:y+boxrady2,x-boxradx1:x+boxradx2);
      if frame == 1,
        himzoom(i,j) = image([x0(j),x1(j)],[y0(i),y1(i)],repmat(box,[1,1,3]));
      else
        set(himzoom(i,j),'cdata',repmat(box,[1,1,3]));
      end
      
    end
  end;

  % plot the zoomed out position
  if frame == 1,
    for fly = 1:nids,
      if isalive(fly),
        i0 = max(1,idx(fly)-taillength);
        htail(fly) = plot(trx(fly).x(i0:idx(fly)),trx(fly).y(i0:idx(fly)),...
          '-','color',colors(fly,:));
        htri(fly) = drawflyo(trx(fly),idx(fly));
        set(htri(fly),'color',colors(fly,:));
      else
        htail(fly) = plot([],[],'-','color',colors(fly,:));
        htri(fly) = plot([],[],'-','color',colors(fly,:));
      end
    end
  else
    for fly = 1:nids,
      if isalive(fly),
        i0 = max(1,idx(fly)-taillength);
        set(htail(fly),'xdata',trx(fly).x(i0:idx(fly)),...
          'ydata',trx(fly).y(i0:idx(fly)));
        updatefly(htri(fly),trx(fly),idx(fly));
      else
        set(htail(fly),'xdata',[],'ydata',[]);
        set(htri(fly),'xdata',[],'ydata',[]);
      end
    end
  end
  
  % plot the zoomed views
  for i = 1:nzoomr,
    for j = 1:nzoomc,
      fly = zoomflies(i,j);
      if ~isnan(fly) && isalive(fly),
        x = trx(fly).x(idx(fly));
        y = trx(fly).y(idx(fly));
        x = boxradius + (x - round(x))+.5;
        y = boxradius + (y - round(y))+.5;
        x = x * scalefactor;
        y = y * scalefactor;
        x = x + x0(j) - 1;
        y = y + y0(i) - 1;
        a = trx(fly).a(idx(fly))*scalefactor;
        b = trx(fly).b(idx(fly))*scalefactor;
        theta = trx(fly).theta(idx(fly));
        if frame == 1,
          hzoom(i,j) = drawflyo(x,y,theta,a,b);
          set(hzoom(i,j),'color',colors(fly,:));
        else
          updatefly(hzoom(i,j),x,y,theta,a,b);
        end
      else
        if frame == 1,
          hzoom(i,j) = plot([],[],'-','color',colors(fly,:));
        else
          set(hzoom(i,j),'xdata',[],'ydata',[]);
        end
      end
    end
  end
  
  if frame == 1,
    input('Resize figure 1 to the desired size, hit enter when done.');
    set(1,'visible','off');
    aviobj = avifile(aviname,'fps',fps);
  end
  
  if frame == 1,
    fr = getframe_invisible(hax);
    [height,width,tmp] = size(fr);
  else
    fr = getframe_invisible(hax,[height,width]);
  end
  aviobj = addframe(aviobj,fr);
    
end

aviobj = close(aviobj);
fclose(fid);