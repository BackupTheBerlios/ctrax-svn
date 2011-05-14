function make_test_data
% make_test_data
%
% makes Python-readable test data for Ctrax test suite
%
% JAB 5/12/11

% find directory for Python code -- same as own directory!
base_name = which( 'make_test_data' );
base_dir = base_name(1:strfind( base_name, 'make_test_data.m' ) - 1);

% use Python to make lists of filenames
eval( ['!python ' base_dir 'dump_mat_names.py'] )

% read in list of original filenames
fp = fopen( 'mat_name_dump.tmp', 'r' );
tracked_filenames = {};
line = fgetl( fp );
while line ~= -1
   tracked_filenames{end+1} = line;
   line = fgetl( fp );
end
fclose( fp );
delete mat_name_dump.tmp

% resave files
for li = 1:length( tracked_filenames )
   % load data
   load( tracked_filenames{li} )
   a = trx.a;
   b = trx.b;
   theta = trx.theta;
   x = trx.x;
   y = trx.y;
   
   % resave data in a Python-friendly MAT-file
   newname = [tracked_filenames{li}(1:strfind( tracked_filenames{li}, '.mat' ) - 1) '_forpython.mat'];
   save( newname, 'x', 'y', 'theta', 'a', 'b' )
end

% read in list of newly tracked filenames
fp = fopen( 'mat_new_name_dump.tmp', 'r' );
new_filenames = {}; % need load_trx run on them
line = fgetl( fp );
while line ~= -1
   new_filenames{end+1} = line;
   line = fgetl( fp );
end
fclose( fp );
delete mat_new_name_dump.tmp

% resave files
for li = 1:length( new_filenames )
   [trx, matname, succeeded] = load_tracks( new_filenames{li} );
   if succeeded
      a = trx.a;
      b = trx.b;
      theta = trx.theta;
      x = trx.x;
      y = trx.y;
      
      % resave
      newname = [new_filenames{li}(1:strfind( new_filenames{li}, '.mat' ) - 1) '_forpython.mat'];
      save( newname, 'x', 'y', 'theta', 'a', 'b' )
   end
end

