function [trx,savename] = cleanup_ctrax_data(matname,moviename,in,ds)

if ~exist('ds','var')
  ds = '';
end

idscurr = unique(in.identity);

in.x_pos = in.x_pos + 1;
in.y_pos = in.y_pos + 1;

% frame number
framenumber = zeros(size(in.x_pos));
j = 0;
for i = 1:length(in.ntargets),
  framenumber(j+(1:in.ntargets(i))) = i;
  j = j + in.ntargets(i);
end;

newidentity = nan(size(in.identity));
for id = idscurr,
  idx = in.identity == id;
  datacurr.x = in.x_pos(idx);
  datacurr.y = in.y_pos(idx);
  datacurr.theta = in.angle(idx);
  datacurr.a = in.maj_ax(idx);
  datacurr.b = in.min_ax(idx);
  datacurr.id = id;
  datacurr.moviename = moviename;
  datacurr.firstframe = framenumber(find(idx,1));
  datacurr.arena.x = nan;
  datacurr.arena.y = nan;
  datacurr.arena.r = nan;
  datacurr.f2i = @(f) f - datacurr.firstframe + 1;
  datacurr.nframes = length(datacurr.x);
  datacurr.endframe = datacurr.nframes + datacurr.firstframe - 1;
  if ~exist('trx','var'),
    trx = datacurr;
  else
    trx(end+1) = datacurr; %#ok<AGROW>
  end;
  newidentity(idx) = length(trx);
end

savename = strrep(matname,'.mat',['trx',ds,'.mat']);
save(savename,'trx');