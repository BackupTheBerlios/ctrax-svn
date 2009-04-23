function [data,histstuff] = extractdata(trx,props,nbins,lb,ub,transforms,isseg,averaging,doinvert,flytype,dataname)

nflies = length(trx);
nprops = length(props);
ndata = length(isseg);
if ~exist('averaging','var'),
  averaging = repmat({'None (per-frame)'},[1,ndata]);
end
if ~exist('doinvert','var'),
  doinvert = false(1,ndata);
end


% invert seg, if necessary
for j = 1:ndata,
  if isseg(j) && doinvert(j),
    for i = 1:length(trx),
      trx(i).seg{j} = invertseg(trx(i).seg{j});
      trx(i).duration{j} = (trx(i).seg{j}.t2 - trx(i).seg{j}.t1)/trx(i).fps;
    end
  end
end

istype = false(ndata,nflies);
if ~isfield(trx,'type'),
  istype(:) = true;
else
  for i = 1:ndata,
    for j = 1:nflies,
      istype(i,j) = ismember(trx(j).type,flytype{i});
    end
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
data = cell(ndata,nflies);
for datai = 1:ndata,

  for i = 1:nflies,  
  
    if ~istype(datai,i), 
      data{datai,i} = []; 
      continue;
    end
    
    n = trx(i).nframes - 1;
    
  
    % extract out segments, if necessary
    if isseg(datai),
      if strcmpi(averaging{datai},'none (per-frame)'),
        bw = set_interval_ends(trx(i).seg{datai}.t1,trx(i).seg{datai}.t2-1,n);
        data{datai,i} = zeros(nnz(bw),nprops);
        for j = 1:nprops,
          if strcmpi(props{j},'duration'),
            off = 0;
            for k = 1:length(trx(i).seg{datai}.t1),
              ncurr = trx(i).seg{datai}.t2(k) - trx(i).seg{datai}.t1(k);
              data{datai,i}(off+(1:ncurr),j) = trx(i).duration{datai}(k);
              off = off + ncurr;
            end
          else
            data{datai,i}(:,j) = trx(i).(props{j})(bw);
          end
        end
      else
        data{datai,i} = zeros(length(trx(i).seg{datai}.t1),nprops);
        for j = 1:nprops,
          if strcmpi(props{j},'duration'),
            data{datai,i}(:,j) = trx(i).duration{datai};
          else
            switch lower(averaging{datai}),
             case 'interval mean',
              if any(strcmpi(trx(1).units.(props{j}).num,'rad')),
                % convert to per-frame
                n = nnz(strcmpi(trx(1).units.(props{j}).den,'s'));
                datacurr = trx(i).(props{j}) / trx(i).fps^n;
                data{datai,i}(:,j) = interval_mean_angle(datacurr',trx(i).seg{datai}.t1,...
                                                         trx(i).seg{datai}.t2-1);
                % convert back to per-second
                data{datai,i}(:,j) = data{datai,i}(:,j)*trx(i).fps^n;
              else
                data{datai,i}(:,j) = interval_mean(trx(i).(props{j})',trx(i).seg{datai}.t1,...
                                                   trx(i).seg{datai}.t2-1);
              end
             case 'interval median',
              if any(strcmpi(trx(1).units.(props{j}).num,'rad')),
                % convert to per-frame
                n = nnz(strcmpi(trx(1).units.(props{j}).den,'s'));
                datacurr = trx(i).(props{j}) / trx(i).fps^n;
                data{datai,i}(:,j) = interval_median_angle(datacurr',trx(i).seg{datai}.t1,...
                                                           trx(i).seg{datai}.t2-1);
                % convert back to per-second
                data{datai,i}(:,j) = data{datai,i}(:,j)*trx(i).fps^n;
              else
                data{datai,i}(:,j) = interval_median(trx(i).(props{j})',trx(i).seg{datai}.t1,...
                                                     trx(i).seg{datai}.t2-1);
              end
             case 'interval start',
              data{datai,i}(:,j) = trx(i).(props{j})(trx(i).seg{datai}.t1);
             case 'interval end',
              data{datai,i}(:,j) = trx(i).(props{j})(trx(i).seg{datai}.t2-1);
            end
          end
        end
      end
    else
      data{datai,i} = zeros(n,nprops);
      for j = 1:nprops,
        data{datai,i}(:,j) = trx(i).(props{j})(1:n);
      end
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
  if any(strcmpi(trx(1).units.(props{i}).num,'rad')),
    fprintf('Converting %s from radians to degrees\n',props{i});
    centers{i} = centers{i}*180/pi;
  end
end

% allocate
if nprops == 1,  
  % allocate per-fly stuff
  histstuff.countsperfly = zeros([nflies,nbins,ndata]);
  histstuff.fracperfly = zeros([nflies,nbins,ndata]);
else
  histstuff.countsperfly = zeros([nbins(1),nbins(2),nflies,ndata]);
  histstuff.fracperfly = zeros([nbins(1),nbins(2),nflies,ndata]);
end

% histogram
for datai = 1:ndata,

  if nprops == 1,
  
    for i = 1:nflies,
      if ~istype(datai,i),
        continue;
      end
      countscurr = histc(data{datai,i},edges{1});
      countscurr(end) = [];
      histstuff.countsperfly(i,:,datai) = countscurr;
      histstuff.fracperfly(i,:,datai) = countscurr / sum(countscurr);
    end
  elseif nprops == 2,
    
    for i = 1:nflies,
      if ~istype(datai,i),
        continue;
      end
      countscurr = hist3(data{datai,i},'edges',edges);
      countscurr(:,end) = [];
      countscurr(end,:) = [];
      histstuff.countsperfly(:,:,i,datai) = countscurr;
      histstuff.fracperfly(:,:,i,datai) = countscurr / sum(countscurr(:));
    end
  else
    fprintf('Cannot histogram with more than 2 properties\n');
    return;
  end
  
end

% labels of plot axes
for i = 1:nprops,
  if strcmpi(transforms{i},'absolute value'),
    props{i} = ['absolute value of ',lower(props{i})];
  elseif strcmpi(transforms{i},'log absolute value'),
    props{i} = ['log absolute value of ',lower(props{i})];
  end

%  for datai = 1:ndata,
%    
%    if isseg(datai),
%      if strcmpi(averaging{datai},'none (per-frame)'),
%        if doinvert(datai),
%          histstuff.propnames{datai,i} = sprintf('%s not during behavior',props{i});
%        else
%          histstuff.propnames{datai,i} = sprintf('%s during behavior',props{i});
%        end
%      else
%        if doinvert(datai),
%          histstuff.propnames{datai,i} = sprintf('%s not behavior %s',props{i},lower(averaging{datai}));
%        else
%          histstuff.propnames{datai,i} = sprintf('%s behavior %s',props{i},lower(averaging{datai}));
%        end
%      end
%    else
%      histstuff.propnames{datai,i} = sprintf('%s during all frames',props{i});
%    end
%    
%    %if ndata > 1,
%    %  histstuff.propnames{datai,i} = [histstuff.propnames{datai,i},', ',dataname{datai}];
%    %end
%  
%  end

end
histstuff.propnames = props;

if nprops == 1,

  % allocate
  % countsperfly = nflies x nbins x ndata
  histstuff.countstotal = zeros([1,nbins,ndata]);
  histstuff.countsmean = zeros([1,nbins,ndata]);
  histstuff.countsstd = zeros([1,nbins,ndata]);
  histstuff.countsstderr = zeros([1,nbins,ndata]);
  histstuff.fractotal = zeros([1,nbins,ndata]);
  histstuff.fracmean = zeros([1,nbins,ndata]);
  histstuff.fracstd = zeros([1,nbins,ndata]);
  histstuff.fracstderr = zeros([1,nbins,ndata]);

  for datai = 1:ndata,
    nfliescurr = nnz(istype(datai,:));
    histstuff.countstotal(:,:,datai) = sum(histstuff.countsperfly(istype(datai,:),:,datai),1);
    histstuff.countsmean(:,:,datai) = histstuff.countstotal(:,:,datai) / nfliescurr;
    histstuff.countsstd(:,:,datai) = std(histstuff.countsperfly(istype(datai,:),:,datai),1,1);
    histstuff.countsstderr(:,:,datai) = histstuff.countsstd(:,:,datai) / sqrt(nfliescurr);
    histstuff.fractotal(:,:,datai) = histstuff.countstotal(:,:,datai) / ...
        sum(histstuff.countstotal(:,:,datai));
    histstuff.fracmean(:,:,datai) = mean(histstuff.fracperfly(istype(datai,:),:,datai),1);
    histstuff.fracstd(:,:,datai) = std(histstuff.fracperfly(istype(datai,:),:,datai),1,1);
    histstuff.fracstderr(:,:,datai) = histstuff.fracstd(:,:,datai) / sqrt(nfliescurr);
  end
  
else
  
  % allocate
  % countsperfly = nbins1 x nbins2 x nflies x ndata
  histstuff.countstotal = zeros([nbins(1),nbins(2),1,ndata]);
  histstuff.countsmean = zeros([nbins(1),nbins(2),1,ndata]);
  histstuff.countsstd = zeros([nbins(1),nbins(2),1,ndata]);
  histstuff.countsstderr = zeros([nbins(1),nbins(2),1,ndata]);
  histstuff.fractotal = zeros([nbins(1),nbins(2),1,ndata]);
  histstuff.fracmean = zeros([nbins(1),nbins(2),1,ndata]);
  histstuff.fracstd = zeros([nbins(1),nbins(2),1,ndata]);
  histstuff.fracstderr = zeros([nbins(1),nbins(2),1,ndata]);

  for datai = 1:ndata,
    nfliescurr = nnz(istype(datai,:));
    histstuff.countstotal(:,:,:,datai) = sum(histstuff.countsperfly(:,:,istype(datai,:),datai),3);
    histstuff.countsmean(:,:,:,datai) = histstuff.countstotal(:,:,:,datai) / nfliescurr;
    histstuff.countsstd(:,:,:,datai) = std(histstuff.countsperfly(:,:,istype(datai,:),datai),1,3);
    histstuff.countsstderr(:,:,:,datai) = histstuff.countsstd(:,:,:,datai) / sqrt(nfliescurr);
    histstuff.fractotal(:,:,:,datai) = histstuff.countstotal(:,:,:,datai) / ...
        sum(sum(histstuff.countstotal(:,:,:,datai)));
    histstuff.fracmean(:,:,:,datai) = mean(histstuff.fracperfly(:,:,istype(datai,:),datai),3);
    histstuff.fracstd(:,:,:,datai) = std(histstuff.fracperfly(:,:,istype(datai,:),datai),1,3);
    histstuff.fracstderr(:,:,:,datai) = histstuff.fracstd(:,:,:,datai) / sqrt(nfliescurr);
  end  
end
histstuff.centers = centers;
histstuff.dataname = dataname;
histstuff.istype = istype;