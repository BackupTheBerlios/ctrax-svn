function fmf2avi(infilename,outfilename,nframes,scale)
% fmf2avi(infilename,outfilename,nframes,scale)
%
% reads first NFRAMES from INFILENAME, resizes by a
% factor of SCALE, and writes to OUTFILENAME

% updated JAB 6/9/10
%   added docstring
%   made some input arguments optional
%   added option to set output framerate based on input
% updated JAB 5/27/11
%   padded widths to a multiple of 4 (required for AVI)

if ~exist( 'outfilename', 'var' )
   outfilename = [infilename(1:end-3) 'avi'];
end
if ~exist( 'nframes', 'var' ) || isempty( nframes )
   nframes = inf; % all
end
if ~exist( 'scale', 'var' )
   scale = 1;
end   

if 1 % calculate output framerate
   [max_n_frames, frame_period] = fmf_frame_count( infilename, 1 );
   framerate = 1/frame_period;
else
   framerate = 30;
end

% open input file
[header_size, version, f_height, f_width, bytes_per_chunk, max_n_frames, data_format] = fmf_read_header( infilename );

if isinf( nframes )
   nframes = max_n_frames;
end
fp = fopen(infilename,'r');
fseek(fp,header_size,'bof');

% open output file
aviobj = avifile(outfilename);
aviobj.fps = framerate;
aviobj.compression = 'none';

% read and write
wb = waitbar( 0, 'converting to AVI' );
for i = 1:nframes,
   waitbar( i/nframes, wb, 'converting to AVI' )
  data = fmf_read_frame(fp,f_height,f_width,bytes_per_chunk, data_format);
  data = imresize(uint8(data),scale);
  
  % pad width to a multiple of 4
  if mod( f_width, 4 ) ~= 0
     data = padarray( data, [0 4*ceil( size( data, 2 )/4 ) - size( data, 2 )], 0, 'post' );
  end
  
  aviobj = addframe(aviobj,repmat(data,[1,1,3]));
end

close( wb )
aviobj = close(aviobj);
