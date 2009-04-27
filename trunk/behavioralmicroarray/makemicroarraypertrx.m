function [microarray,succeeded] = makemicroarraypertrx(trxname,segnames,basename,microarray,behredo)

isinput = true;
if ~exist('microarray','var'),
  microarray = struct;
  isinput = false;
end
if ~exist('behredo','var')
  behredo = [];
end
succeeded = false;

statmsg = sprintf('Computing microarray for trxfile %s\n',basename);
fprintf([statmsg,'\n']);

if isstruct(trxname),
  % actually input trx instead of trxname
  trx = trxname;
else

  % load in current trx
  try
    load(trxname,'trx');
  catch
    fprintf('Could not load trx from %s\n',trxname);
    return;
  end
end
if ~isfield(trx,'units'),
  fprintf('No variable trx.units in %s\n',trxname);
  return;
end

nflies = length(trx);
nbehaviors = length(segnames);

% nothing to do
if isinput && nbehaviors == length(microarray(1).duration) && isempty(behredo),
  return;
end

% per-frame properties
props = getperframepropnames(trx);
nprops = length(props);

hstatus = waitbar(0,statmsg);

% allocate
if ~isinput,
  ma = struct;
  nbehaviors0 = 0;
  ma.duration = zeros(1,nbehaviors);
  ma.fractime = zeros(1,nbehaviors);
  ma.frequency = zeros(1,nbehaviors);
  durings = {'during','invert'};
  averagings = {'none','mean','median','start','end'};
  for propi = 1:nprops,
    prop = props{propi};
    for duringi = 1:length(durings),
      during = durings{duringi};
      for avei = 1:length(averagings),
        averaging = averagings{avei};
        s = sprintf('%s_%s_%s',prop,during,averaging);
        ma.(s) = zeros(1,nbehaviors);
      end
    end
    s = sprintf('%s_allframes',prop);
    ma.(s) = 0;
  end
  microarray = repmat(ma,[1,nflies]);
  clear ma;
else
  nbehaviors0 = length(microarray(1).duration);
  % nbehaviors0 < nbehaviors
  for fly = 1:nflies,
    microarray(fly).duration = [microarray(fly).duration,zeros(1,nbehaviors-nbehaviors0)];
    microarray(fly).fractime = [microarray(fly).fractime,zeros(1,nbehaviors-nbehaviors0)];
    microarray(fly).frequency = [microarray(fly).frequency,zeros(1,nbehaviors-nbehaviors0)];
    durings = {'during','invert'};
    averagings = {'none','mean','median','start','end'};
    for propi = 1:nprops,
      prop = props{propi};
      for duringi = 1:length(durings),
        during = durings{duringi};
        for avei = 1:length(averagings),
          averaging = averagings{avei};
          s = sprintf('%s_%s_%s',prop,during,averaging);
          microarray(fly).(s) = [microarray(fly).(s),zeros(1,nbehaviors-nbehaviors0)];
        end
      end
    end
  end
end

behdo = union([nbehaviors0+1:nbehaviors],behredo);

statuscount = 1;
nstatus = 1 + nflies*(length(behdo)+double(~isinput));
waitbar(statuscount/nstatus,hstatus);

% loop over each behavior
for behi = behdo,
  
  try
    load(segnames{behi},'seg');
  catch
    fprintf('Could not load seg from %s\n',segnames{behi});
    delete(hstatus);
    return;
  end
  if length(seg) ~= nflies,
    fprintf('Length of seg = %d for %s\ndoes not match length of trx = %d for %s\n',...
            length(seg),segnames{behi},nflies,trxname);
    delete(hstatus);
    return;
  end
  
  % for now, we will only compute mean per fly, not standard deviation over intervals
  for fly = 1:nflies,
  
    nbeh = length(seg(fly).t1);
    nframesbeh = sum(seg(fly).t2 - seg(fly).t1);
    % mean duration in seconds
    microarray(fly).duration(behi) = (nframesbeh / nbeh) / trx(fly).fps;
    % fractime is unitless
    microarray(fly).fractime(behi) = nframesbeh / trx(fly).nframes;
    % frequency in nbeh / second
    microarray(fly).frequency(behi) = nbeh / trx(fly).nframes * trx(fly).fps;
    
    % compute per-frame property stuff
    bw_during = set_interval_ends(seg(fly).t1,seg(fly).t2-1,trx(fly).nframes-1);
    bw_invert = ~bw_during;
    seg_invert = invertseg(seg(fly),trx(fly).nframes-1);
    for propi = 1:nprops,
      prop = props{propi};
      isangle = any(strcmpi(trx(fly).units.(prop).num,'rad'));
      if isangle,
        ns = nnz(strcmpi(trx(fly).units.(prop).den,'s'));
        perframedata = trx(fly).(prop) / trx(fly).fps^ns;
      else
        perframedata = trx(fly).(prop);
      end

      % remove nans if necessary
      perframedata = removenansfromperframedata(perframedata,prop);
      
      % during
      
      % averaging = none
      s = sprintf('%s_during_none',prop);
      
      % averaging angles is annoying
      if isangle,
        % perframedata is converted to per-frame
        mu = set_mean_angle(perframedata(bw_during));
        % convert back to per-second
        microarray(fly).(s)(behi) = mu * trx(fly).fps^ns;
      else
        microarray(fly).(s)(behi) = mean(perframedata(bw_during));
      end
      
      % averaging = mean
      s = sprintf('%s_during_mean',prop);
      
      if isangle,
        % compute mean of each interval
        mus = interval_mean_angle(perframedata',seg(fly).t1,seg(fly).t2-1);
        % compute mean of means
        mu = set_mean_angle(mus);
        % convert back to per-second
        microarray(fly).(s)(behi) = mu * trx(fly).fps^ns;
      else
        % compute mean of each interval
        mus = interval_mean(perframedata',seg(fly).t1,seg(fly).t2-1);
        % compute mean of means
        microarray(fly).(s)(behi) = mean(mus);
      end
      
      % averaging = median
      s = sprintf('%s_during_median',prop);
      
      if isangle,
        % compute median of each interval
        mus = interval_median_angle(perframedata',seg(fly).t1,seg(fly).t2-1);
        % compute mean of medians
        mu = set_mean_angle(mus);
        % convert back to per-second
        microarray(fly).(s)(behi) = mu * trx(fly).fps^ns;
      else
        % compute median of each interval
        mus = interval_median(perframedata',seg(fly).t1,seg(fly).t2-1);
        % compute mean of medians
        microarray(fly).(s)(behi) = mean(mus);
      end
      
      % averaging = start
      s = sprintf('%s_during_start',prop);
            
      % property at start of each interval
      mus = perframedata(seg(fly).t1);
      if isangle,
        % compute mean of starts
        mu = set_mean_angle(mus);
        % convert back to per-second
        microarray(fly).(s)(behi) = mu * trx(fly).fps^ns;
      else
        microarray(fly).(s)(behi) = mean(mus);
      end
      
      % averaging = end
      s = sprintf('%s_during_end',prop);
            
      % property at end of each interval
      mus = perframedata(seg(fly).t2-1);
      if isangle,
        % compute mean of ends
        mu = set_mean_angle(mus);
        % convert back to per-second
        microarray(fly).(s)(behi) = mu * trx(fly).fps^ns;
      else
        microarray(fly).(s)(behi) = mean(mus);
      end
   
      % invert
      
      % averaging = none
      s = sprintf('%s_invert_none',prop);
      
      % averaging angles is annoying
      if isangle,
        % perframedata is converted to per-frame
        mu = set_mean_angle(perframedata(bw_invert));
        % convert back to per-second
        microarray(fly).(s)(behi) = mu * trx(fly).fps^ns;
      else
        microarray(fly).(s)(behi) = mean(perframedata(bw_invert));
      end
      
      % averaging = mean
      s = sprintf('%s_invert_mean',prop);
      
      if isangle,
        % compute mean of each interval
        mus = interval_mean_angle(perframedata',seg_invert.t1,seg_invert.t2-1);
        % compute mean of means
        mu = set_mean_angle(mus);
        % convert back to per-second
        microarray(fly).(s)(behi) = mu * trx(fly).fps^ns;
      else
        % compute mean of each interval
        mus = interval_mean(perframedata',seg_invert.t1,seg_invert.t2-1);
        % compute mean of means
        microarray(fly).(s)(behi) = mean(mus);
      end
      
      % averaging = median
      s = sprintf('%s_invert_median',prop);
      
      if isangle,
        % compute median of each interval
        mus = interval_median_angle(perframedata',seg_invert.t1,seg_invert.t2-1);
        % compute mean of medians
        mu = set_mean_angle(mus);
        % convert back to per-second
        microarray(fly).(s)(behi) = mu * trx(fly).fps^ns;
      else
        % compute median of each interval
        mus = interval_median(perframedata',seg_invert.t1,seg_invert.t2-1);
        % compute mean of medians
        microarray(fly).(s)(behi) = mean(mus);
      end
      
      % averaging = start
      s = sprintf('%s_invert_start',prop);
            
      % property at start of each interval
      mus = perframedata(seg_invert.t1);
      if isangle,
        % compute mean of starts
        mu = set_mean_angle(mus);
        % convert back to per-second
        microarray(fly).(s)(behi) = mu * trx(fly).fps^ns;
      else
        microarray(fly).(s)(behi) = mean(mus);
      end
      
      % averaging = end
      s = sprintf('%s_invert_end',prop);
            
      % property at end of each interval
      mus = perframedata(seg_invert.t2-1);
      if isangle,
        % compute mean of ends
        mu = set_mean_angle(mus);
        % convert back to per-second
        microarray(fly).(s)(behi) = mu * trx(fly).fps^ns;
      else
        microarray(fly).(s)(behi) = mean(mus);
      end
      
    end % end loop over props

    statuscount = statuscount + 1;
    waitbar(statuscount/nstatus,hstatus);
    
  end % end loop over flies
  
end % end loop over behaviors

% ignoring behavior segmentation, all frames
if ~isinput,
  
  for fly = 1:nflies,
    
    for propi = 1:nprops,
      
      prop = props{propi};
      isangle = any(strcmpi(trx(fly).units.(prop).num,'rad'));
      if isangle,
        ns = nnz(strcmpi(trx(fly).units.(prop).den,'s'));
        perframedata = trx(fly).(prop) / trx(fly).fps^ns;
      else
        perframedata = trx(fly).(prop);
      end
      perframedata = removenansfromperframedata(perframedata,prop);
      
      s = sprintf('%s_allframes',prop);
      
      % averaging angles is annoying
      if isangle,
        % perframedata is converted to per-frame
        mu = set_mean_angle(perframedata);
        % convert back to per-second
        microarray(fly).(s) = mu * trx(fly).fps^ns;
      else
        microarray(fly).(s) = mean(perframedata);
      end
      
    end

    statuscount = statuscount + 1;
    waitbar(statuscount/nstatus,hstatus);
    
  end

end

% set types for each fly
if isfield(trx,'type') && isfield(trx,'sex'),
  for i = 1:nflies,
   microarray(i).type = sprintf('%s, %s, %s',trx(i).type,trx(i).sex,basename);
  end
elseif isfield(trx,'sex') && ~isfield(trx,'type'),
  for i = 1:nflies,
    microarray(i).type = sprintf('%s, %s',trx(i).sex,basename);
  end
elseif isfield(trx,'type'),
  for i = 1:nflies,
    microarray(i).type = sprintf('%s, %s',trx(i).type,basename);
  end
else
  for i = 1:nflies,
    microarray(i).type = sprintf('%s',basename);
  end
end

succeeded = true;
delete(hstatus);