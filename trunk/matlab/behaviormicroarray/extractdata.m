function [data,histstuff] = extractdata(trx,props,nbins,lb,ub,transforms,isseg,averaging,doinvert)

nflies = length(trx);
nprops = length(props);
if isseg && ~exist('averaging','var'),  
  averaging = 'None (per-frame)';
end
if ~exist('doinvert','var'),
  doinvert = false;
end

% invert seg, if necessary
if isseg && doinvert,
  for i = 1:length(trx),
    trx(i).seg = invertseg(trx(i).seg);
    trx(i).duration = (trx(i).seg.t2 - trx(i).seg.t1)/trx(i).fps;
  end
end

for j = 1:nprops,
  for i = 1:nflies,
    if strcmpi(transforms{j},'absolute value'),
      trx(i).(props{j}) = abs(trx(i).(props{j}));
    elseif strcmpi(transforms{j},'log absolute value'),
      trx(i).(props{j}) = log(abs(trx(i).(props{j})));
    end
  end
end
  

% data
data = cell(1,nflies);
for i = 1:nflies,  
  
  n = trx(i).nframes - 1;
  
  % extract out segments, if necessary
  if isseg,
    if strcmpi(averaging,'none (per-frame)'),
      bw = set_interval_ends(trx(i).seg.t1,trx(i).seg.t2-1,n);
      data{i} = zeros(nnz(bw),nprops);
      for j = 1:nprops,
        if strcmpi(props{j},'duration'),
          off = 0;
          for k = 1:length(trx(i).seg.t1),
            ncurr = trx(i).seg.t2(k) - trx(i).seg.t1(k);
            data{i}(off+(1:ncurr),j) = trx(i).duration(k);
            off = off + ncurr;
          end
        else
          data{i}(:,j) = trx(i).(props{j})(bw);
        end
      end
    else
      data{i} = zeros(length(trx(i).seg.t1),nprops);
      for j = 1:nprops,
        if strcmpi(props{j},'duration'),
          data{i}(:,j) = trx(i).duration;
        else
          switch lower(averaging),
            case 'interval mean',
              if isanglename(props{j}),
                % convert to per-frame
                if isvelocityname(props{j}),
                  datacurr = trx(i).(props{j})/trx(i).fps;
                elseif isaccelname(props{j}),
                  datacurr = trx(i).(props{j})/trx(i).fps^2;
                else
                  datacurr = trx(i).(props{j});
                end
                data{i}(:,j) = interval_mean_angle(datacurr',trx(i).seg.t1,trx(i).seg.t2-1);
                % convert back to per-second
                if isvelocityname(props{j}),
                  data{i}(:,j) = data{i}(:,j)*trx(i).fps;
                elseif isaccelname(props{j}),
                  data{i}(:,j) = data{i}(:,j)*trx(i).fps^2;
                else
                  datacurr = trx(i).(props{j});
                end
              else
                data{i}(:,j) = interval_mean(trx(i).(props{j})',trx(i).seg.t1,trx(i).seg.t2-1);
              end
            case 'interval median',
              if isanglename(props{j}),
                % convert to per-frame
                if isvelocityname(props{j}),
                  datacurr = trx(i).(props{j})/trx(i).fps;
                elseif isaccelname(props{j}),
                  datacurr = trx(i).(props{j})/trx(i).fps^2;
                end
                data{i}(:,j) = interval_median_angle(datacurr',trx(i).seg.t1,trx(i).seg.t2-1);
                % convert back to per-second
                if isvelocityname(props{j}),
                  data{i}(:,j) = data{i}(:,j)*trx(i).fps;
                elseif isaccelname(props{j}),
                  data{i}(:,j) = data{i}(:,j)*trx(i).fps^2;
                end
              else
                data{i}(:,j) = interval_median(trx(i).(props{j})',trx(i).seg.t1,trx(i).seg.t2-1);
              end
            case 'interval start',
              data{i}(:,j) = trx(i).(props{j})(trx(i).seg.t1);
            case 'interval end',
              data{i}(:,j) = trx(i).(props{j})(trx(i).seg.t2);
          end
        end
      end
    end
  else
    data{i} = zeros(n,nprops);
    for j = 1:nprops,
      data{i}(:,j) = trx(i).(props{j})(1:n);
    end
  end

end

% bin edges
edges = cell(1,nprops);
for i = 1:nprops,
  edges{i} = linspace(lb(i),ub(i),nbins(i)+1);
end
% bin centers
centers = cell(1,nprops);
for i = 1:length(centers),
  centers{i} = (edges{i}(1:end-1)+edges{i}(2:end))/2;
end

% change to degrees if necessary
for i = 1:length(centers),
  if isanglename(props{i}),
    fprintf('Converting %s from radians to degrees\n',props{i});
    centers{i} = centers{i}*180/pi;
  end
end

% allocate
histstuff.propnames = cell(1,nprops);

% histogram
if nprops == 1,
  % allocate per-fly stuff
  histstuff.countsperfly = zeros(nflies,nbins);
  histstuff.fracperfly = zeros(nflies,nbins);

  for i = 1:nflies,
    n = size(data{i},1);
    countscurr = histc(data{i},edges{1});
    countscurr(end) = [];
    histstuff.countsperfly(i,:) = countscurr;
    histstuff.fracperfly(i,:) = countscurr / n;
  end
elseif nprops == 2,
  % allocate per-fly stuff
  histstuff.countsperfly = zeros([nbins(1),nbins(2),nflies]);
  histstuff.fracperfly = zeros([nbins(1),nbins(2),nflies]);

  for i = 1:nflies,
    n = size(data{i},1);
    countscurr = hist3(data{i},'edges',edges);
    countscurr(:,end) = [];
    countscurr(end,:) = [];
    histstuff.countsperfly(:,:,i) = countscurr;
    histstuff.fracperfly(:,:,i) = countscurr / n;
  end
else
  fprintf('Cannot histogram with more than 2 properties\n');
  return;
end

for i = 1:nprops,
  if strcmpi(transforms{i},'absolute value'),
    props{i} = ['absolute value of ',lower(props{i})];
  elseif strcmpi(transforms{i},'log absolute value'),
    props{i} = ['log absolute value of ',lower(props{i})];
  end

  if isseg,
    if strcmpi(props{i},'duration'),
      histstuff.propnames{i} = 'Duration';
    elseif strcmpi(averaging,'none (per-frame)'),
      if doinvert,
        histstuff.propnames{i} = sprintf('%s not during behavior',props{i});
      else
        histstuff.propnames{i} = sprintf('%s during behavior',props{i});
      end
    else
      if doinvert,
        histstuff.propnames{i} = sprintf('%s not behavior %s',props{i});
      else
        histstuff.propnames{i} = sprintf('%s behavior %s',props{i});
      end
    end
  else
    histstuff.propnames{i} = sprintf('%s during all frames',props{i});
  end
end

if nprops == 1,
  dim = 1;
else
  dim = 3;
end
histstuff.countstotal = sum(histstuff.countsperfly,dim);
histstuff.countsmean = mean(histstuff.countsperfly,dim);
histstuff.countsstd = std(histstuff.countsperfly,1,dim);
histstuff.countsstderr = histstuff.countsstd / sqrt(nflies);
histstuff.fractotal = histstuff.countstotal / sum(histstuff.countstotal(:));
histstuff.fracmean = mean(histstuff.fracperfly,dim);
histstuff.fracstd = std(histstuff.fracperfly,1,dim);
histstuff.fracstderr = histstuff.fracstd / sqrt(nflies);
histstuff.centers = centers;