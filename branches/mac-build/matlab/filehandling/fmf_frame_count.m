function [f, rate] = fmf_frame_count( filename, force_rate )
% [nframes, framerate] = fmf_frame_count( filename )
%
% Counts the number of frames in a FlyMovieFormat version 1 file.
% Also returns the mean framerate.
%
% JAB 6/30/04

[header_size, version, f_height, f_width, bytes_per_chunk, ...
   max_n_frames, data_format] = fmf_read_header( filename );
if max_n_frames ~= 0
	f = max_n_frames;
   rate = nan;
   if ~exist( 'force_rate', 'var' ) || ~force_rate
      return
   else
      fprintf( 1, 'forced reading frames to calculate framerate\n' )
   end
end

% count frames
fp = fopen( filename, 'r' );
fseek( fp, header_size, 'bof' );
f = 0;
if nargout > 1 % read timestamps
   stamps = [];
end
while 1,
   if exist( 'stamps', 'var' )
      stamps = [stamps fread( fp, 1, 'double' )]; % read timestamp
      fseek( fp, -8, 'cof' ); % go backward to beginning of frame
   end
	err = fseek( fp, bytes_per_chunk, 'cof' );
	if err == -1, break, end
	f = f + 1;
end

fclose( fp );

if exist( 'stamps', 'var' )
   rate = mean( diff( stamps ) );
end

% update frame count in file
return % update not working!!
fp = fopen( filename, 'r+' );
fseek( fp, 20, 'bof' );
fwrite( fp, f, 'long' );
fclose( fp );
