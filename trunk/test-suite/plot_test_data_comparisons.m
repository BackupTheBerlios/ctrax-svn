function plot_test_data_comparisons( dump_filename )
% plot_test_data_comparisons( dump_filename )
%
% shows comparisons of groundtruthed and newly tracked data
%
% JAB 5/17/11

if ~exist( 'dump_filename', 'var' )
   dump_filename = '/tmp/Ctrax-test-data/stats_2011-05-17_20-46-20.tmp';
end

% read in list of stat filenames
fp = fopen( dump_filename, 'r' );
stat_filenames = {};
line = fgetl( fp );
while line ~= -1
   stat_filenames{end+1} = line;
   line = fgetl( fp );
end
fclose( fp );

for fi = 1%:length( stat_filenames )
   stat_filenames{fi}
   load( stat_filenames{fi} )

   figure(1); clf; hold on; plot( truedata.x, truedata.y, 'r', 'linewidth', 2 ), plot( newdata.x, newdata.y )
   use_len = min( [length( truedata.x ), length( newdata.x )] );
   use_true.x = truedata.x(1:use_len);
   use_true.y = truedata.y(1:use_len);
   use_true.theta = truedata.theta(1:use_len);
   use_true.a = truedata.a(1:use_len);
   use_true.b = truedata.b(1:use_len);
   use_new.x = newdata.x(1:use_len);
   use_new.y = newdata.y(1:use_len);
   use_new.theta = newdata.theta(1:use_len);
   use_new.a = newdata.a(1:use_len);
   use_new.b = newdata.b(1:use_len);
   
   % 1: how many frames in the track
   comp.length = length( truedata.x ) - use_len;
   
   % 2: position difference
   comp.pos = mean( sqrt( (use_true.x - use_new.x).^2 + (use_true.y - use_new.y).^2 ) );
   
   % 3: angle difference
   comp.ang = mean( unwrap( use_true.theta ) - unwrap( use_new.theta ) );
   
   % 4: size difference
   comp.size = mean( use_true.a.*use_true.b - use_new.a.*use_new.b );

   comp.runtime = runtime;
   save( stat_filenames{fi}, '-append', 'comp' )
end % for each newly tracked file
