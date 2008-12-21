% [readframe,nframes,fid] = get_readframe_fcn(filename)

function [readframe,nframes,fid] = get_readframe_fcn(filename)

[base,ext] = splitext(filename);
if strcmpi(ext,'.fmf'),
  [header_size,version,nr,nc,bytes_per_chunk,nframes,data_format] = fmf_read_header(filename);
  fid = fopen(filename);
  readframe = @(f) fmfreadframe(fid,header_size+(f-1)*bytes_per_chunk,nr,nc,bytes_per_chunk,data_format);
elseif strcmpi(ext,'.sbfmf'),
  [nr,nc,nframes,bgcenter,bgstd,frame2file] = sbfmf_read_header(filename);
  fid = fopen(filename);
  readframe = @(f) sbfmfreadframe(f,fid,frame2file,bgcenter);
elseif strcmpi(ext,'avi'),
  fid = -1;
  readerobj = mmreader(filename);
  nframes = get(readerobj,'NumberOfFrames');
  readframe = @(f) read(readerobj,f);
else
  error(['Unrecognized file extension ',ext]);
end