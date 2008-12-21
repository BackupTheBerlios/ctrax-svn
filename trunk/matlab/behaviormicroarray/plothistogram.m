function plotstuff = plothistogram(histstuff,varargin)
%fighandles has member variables {'fig','haxes'}
%plotmode in {'Mean Fraction per Fly','Mean Count per Fly','Total Count'}
%plotindivs is true/false
%plotstd in {'off','standard deviation','standard error'}
%flycolors

[plotmode,plotstuff,plotindivs,plotstd] = ...
  myparse(varargin,'plotmode','Mean Fraction per Fly','fighandles',struct,...
  'plotindivs',false,'plotstd','standard deviation');

% total count has no single-fly stuff
if strcmpi(plotmode,'total count'),
  plotindivs = false;
  plotstd = 'off';
end

nprops = length(histstuff.centers);
if nprops == 1,
  nflies = size(histstuff.countsperfly,1);
else
  nflies = size(histstuff.countsperfly,3);
end
if ~isfield(plotstuff,'fig'),
  plotstuff.fig = figure;
else
  figure(plotstuff.fig);
  clf;
end

switch lower(plotmode),
  case 'mean count per fly',
    perflydata = histstuff.countsperfly;
    meandata = histstuff.countsmean;
    stddata = histstuff.countsstd;
    stderrdata = histstuff.countsstderr;
  case 'mean fraction per fly',
    perflydata = histstuff.fracperfly;
    meandata = histstuff.fracmean;
    stddata = histstuff.fracstd;
    stderrdata = histstuff.fracstderr;
  case 'total count',
    perflydata = nan;
    meandata = histstuff.countstotal;
    stddata = nan;
    stderrdata = nan;
end

% plus and minus one standard deviation or standard error
if strcmpi(plotstd,'standard error')
  above = meandata + stderrdata;
  below = meandata - stderrdata;
elseif strcmpi(plotstd,'standard deviation')
  above = meandata + stddata;
  below = meandata - stddata;
end

% plot one property histogram
if nprops == 1,
  
  plotstuff.haxes = gca;

  hold on;
  % plot standard deviation/standard error
  if ~strcmpi(plotstd,'off'),
    px = [histstuff.centers{1},fliplr(histstuff.centers{1})];
    py = [above,fliplr(below)];
    plotstuff.hsig = patch(px,py,[.8,.8,.8],'linestyle','none');
  end
  
  % plot individual flies
  if plotindivs,
    if ~isfield(plotstuff,'flycolors') || isempty(plotstuff.flycolors),
      % find projection on principal component
      coeff = princomp(perflydata');
      coeff = coeff(:,1);
      [tmp,order] = sort(coeff);
      [tmp,order] = sort(order);
      plotstuff.flycolors = (3+jet(nflies))/4*.7;
      plotstuff.flycolors = plotstuff.flycolors(order,:);
    end
    plotstuff.hindiv = zeros(1,nflies);
    for fly = 1:nflies,
      plotstuff.hindiv(fly) = plot(histstuff.centers{1},perflydata(fly,:),'-','color',plotstuff.flycolors(fly,:));
    end
  end

  % plot mean
  plotstuff.hmu = plot(histstuff.centers{1},meandata,'ko-','linewidth',2,'markerfacecolor','k','markersize',6);
  axisalmosttight;
   
  xlabel(histstuff.propnames{1});
  ylabel(plotmode);

else
  % two properties
  
  % create axes

  % if we're plotting individuals, compute a good number of rows, columns
  plotstuff.haxes = struct;
  if plotindivs,
    pos = get(plotstuff.fig,'position');
    w = pos(3); h = pos(4);
    nindivcol = max(1,round(sqrt(2*nflies*w/h)));
    nindivrow = ceil(nflies/nindivcol);
  end
  if strcmpi(plotstd,'off') && ~plotindivs,
    % mean only  
    plotstuff.haxes.mean = axes;
  elseif ~strcmpi(plotstd,'off') && ~plotindivs,
    % mean, plus and minus std
    plotstuff.haxes.minus = subplot(1,3,1);
    plotstuff.haxes.mean = subplot(1,3,2);
    plotstuff.haxes.plus = subplot(1,3,3);
  elseif strcmpi(plotstd,'off') && plotindivs,
    % mean, indivs
    plotstuff.haxes.mean = subplot(2,1,1);
    plotstuff.haxes.indiv = zeros(1,nflies);
    for i = 1:nflies,
      plotstuff.haxes.indiv(i) = subplot(nindivrow*2,nindivcol,nindivrow*nindivcol+i);
    end
  elseif ~strcmpi(plotstd,'off') && plotindivs,
    % mean, plus and minus std, indivs
    plotstuff.haxes.minus = subplot(2,3,1);
    plotstuff.haxes.mean = subplot(2,3,2);
    plotstuff.haxes.plus = subplot(2,3,3);
    plotstuff.haxes.indiv = zeros(1,nflies);
    plotstuff.haxes.indiv = zeros(1,nflies);
    for i = 1:nflies,
      plotstuff.haxes.indiv(i) = subplot(nindivrow*2,nindivcol,nindivrow*nindivcol+i);
    end
  end
  
  % plot the mean
  axes(plotstuff.haxes.mean);
  x = [histstuff.centers{2}(1),histstuff.centers{2}(end)];
  y = [histstuff.centers{1}(1),histstuff.centers{1}(end)];
  maxv = max(meandata(:));
  minv = min(meandata(:));
  plotstuff.him.mean = imagesc(x,y,meandata,[minv,maxv]);
  xlabel(histstuff.propnames{2});
  colormap jet;
  title(plotmode);
  if ~strcmpi(plotstd,'off'),
    minv = min(minv,min(below(:)));
    maxv = max(maxv,max(above(:)));
    axes(plotstuff.haxes.minus);
    plotstuff.him.minus = imagesc(x,y,below,[minv,maxv]);
    plotstuff.hcolorbar = colorbar('east');
    ylabel(histstuff.propnames{1});
    title(['-1 ',plotstd]);
    axes(plotstuff.haxes.plus);
    plotstuff.him.plus = imagesc(x,y,above,[minv,maxv]);
    set(plotstuff.haxes.mean,'clim',[minv,maxv]);
    title(['+1 ',plotstd]);
  else
    ylabel(histstuff.propnames{1});
    plotstuff.hcolorbar = colorbar;
  end
  set(plotstuff.hcolorbar,'edgecolor','w');
  if plotindivs,
    for i = 1:nflies,
      axes(plotstuff.haxes.indiv(i));
      imagesc(x,y,perflydata(:,:,i),[minv,maxv]);
      title(sprintf('Fly %d',i));
      axis off;
    end
  end
  % link the color maps of all the axes
  allaxes = struct2cell(plotstuff.haxes);
  allaxes = [allaxes{:}];
  plotstuff.hlink = linkprop(allaxes,'clim');
  
  % add a button for changing colormap
  plotstuff.hcolormapbutton = uicontrol('style','pushbutton','string',...
    'Edit Colormap...','position',[785     7   107    26],...
    'callback',@(hObject, eventdata, handles) colormapeditor);
  
end