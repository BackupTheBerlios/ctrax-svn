% [readframe,nframes,fid] = get_readframe_fcn(filename)

function [readframe,nframes,fid] = get_readframe_fcn(filename)

% allow videoio library to be used if it is installed and on the path
CTRAX_ISVIDEOIO = exist('videoReader','file');

[base,ext] = splitext(filename);
if strcmpi(ext,'.fmf'),
  [header_size,version,nr,nc,bytes_per_chunk,nframes,data_format] = fmf_read_header(filename);
  fid = fopen(filename);
  readframe = @(f) fmfreadframe(fid,header_size+(f-1)*bytes_per_chunk,nr,nc,bytes_per_chunk,data_format);
elseif strcmpi(ext,'.sbfmf'),
  [nr,nc,nframes,bgcenter,bgstd,frame2file] = sbfmf_read_header(filename);
  fid = fopen(filename);
  readframe = @(f) sbfmfreadframe(f,fid,frame2file,bgcenter);
else
  fid = -1;
  if CTRAX_ISVIDEOIO,
    readerobj = videoReader(filename,'preciseFrames',30,'frameTimeoutMS',5000);
    info = getinfo(readerobj);
    nframes = info.numFrames;
    seek(readerobj,0);
    seek(readerobj,1);
    readframe = @(f) videoioreadframe(readerobj,f);
  else
    readerobj = mmreader(filename);
    nframes = get(readerobj,'NumberOfFrames');
    if isempty(nframes),
      % approximate nframes from duration
      nframes = get(readerobj,'Duration')*get(readerobj,'FrameRate');
    end
    readframe = @(f) flipdim(read(readerobj,f),1);
  end
end