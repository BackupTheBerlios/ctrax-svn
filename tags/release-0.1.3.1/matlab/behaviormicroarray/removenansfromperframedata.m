function perframedata = removenansfromperframedata(perframedata,prop)

if any(isnan(perframedata)),
  if all(isnan(perframedata)),
    warning('all data for prop %s is nan\n',prop);
  else
    warning('nan entries for prop %s\n',prop);
    idx = isnan(perframedata);
    idx1 = find(idx);
    for i = idx1,
      i0 = find(~idx(1:i-1),1,'last');
      if isempty(i0),
        i1 = find(~idx(i+1:end),1,'first')+i;
        perframedata(i) = perframedata(i1);
      else
        perframedata(i) = perframedata(i0);
      end
    end
  end
end
