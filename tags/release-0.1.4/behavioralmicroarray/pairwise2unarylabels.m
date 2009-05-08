function pairlabels = pairwise2unarylabels(labels,pairtrx)

nflies = length(pairtrx);
pairlabels = struct('starts',cell(1,nflies),'ends',cell(1,nflies));

for i = 1:length(labels.otherfly),
  fly2 = labels.otherfly(i);
  pairlabels(fly2).starts(end+1) = pairtrx(fly2).f2i(labels.starts(i));
  pairlabels(fly2).ends(end+1) = pairtrx(fly2).f2i(labels.ends(i));
end