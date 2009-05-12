function [trx,savename] = cleanup_ctrax_data(matname,moviename,in,ds)

if ~exist('ds','var')
  ds = '';
end

in.x_pos = in.x_pos(:)';
in.y_pos = in.y_pos(:)';
in.maj_ax = in.maj_ax(:)';
in.min_ax = in.min_ax(:)';
in.angle = in.angle(:)';
in.identity = in.identity(:)';

idscurr = unique(in.identity);

in.x_pos = in.x_pos + 1;
in.y_pos = in.y_pos + 1;

isconverted = isfield(in,'pxpermm') && isfield(in,'fps');

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
  if isconverted,
    datacurr.pxpermm = in.pxpermm;
    datacurr.fps = in.fps;
    datacurr.x_mm = in.x_pos_mm(idx);
    datacurr.y_mm = in.y_pos_mm(idx);
    datacurr.a_mm = in.maj_ax_mm(idx);
    datacurr.b_mm = in.min_ax_mm(idx);
  end
  if ~exist('trx','var'),
    trx = datacurr;
  else
    trx(end+1) = datacurr; %#ok<AGROW>
  end;
  newidentity(idx) = length(trx);
end

savename = strrep(matname,'.mat',['trx',ds,'.mat']);
save(savename,'trx');