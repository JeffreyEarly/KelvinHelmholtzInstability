%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Make a movie from the Kelvin Helmholtz instability example
%	in GLNumericalModelingKit.
%
%	2012 August 7 -- Jeffrey J. Early
%
%	The frames output from this appear are 1556x465 points.
%	Using Image2Movie set the frames to 1920x574 to preserve
%	the aspect ratio. 15 fps seems about right.
%	
%	I ran the simulation at 512x128 grid, 80x20 meters for 200 seconds.


file = '/Users/jearly/Desktop/KelvinHelmholtzInstability.nc';
FramesFolder ='/Users/jearly/Desktop/KelvinHelmholtzInstabilityFrames';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Make the frames folder
%
if exist(FramesFolder,'dir') == 0
	mkdir(FramesFolder);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Read in the problem dimensions
%
x = ncread(file, 'x');
y = ncread(file, 'y');
t = ncread(file, 't');

xDomainLength = max(x)-min(x) + x(2)-x(1);

xunits = ncreadatt(file, 'x', 'units');
yunits = ncreadatt(file, 'y', 'units');
tunits = ncreadatt(file, 't', 'units');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	The stride indicates how many floats we will skip
%

% stride=1, size=3 works well for 512x128
stride = 2;
floatSize = 4;

% Read in the initial position of the floats.
% We will use this information to maintain a constant color on each float.
xposInitial = double(ncread(file, 'x-position', [ceil(stride/2) ceil(stride/2) 1], [length(y)/stride length(x)/stride 1], [stride stride 1]));
xposInitial = reshape(xposInitial, length(y)*length(x)/(stride*stride), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Setup the figure
%
fig = figure('Position', [50 50 1920 585]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

for iTime=1:length(t)
%for iTime=100:100
 	
 	% read in the position of the floats for the given time
	xpos = double(ncread(file, 'x-position', [ceil(stride/2) ceil(stride/2) iTime], [length(y)/stride length(x)/stride 1], [stride stride 1]));
	ypos = double(ncread(file, 'y-position', [ceil(stride/2) ceil(stride/2) iTime], [length(y)/stride length(x)/stride 1], [stride stride 1]));
	
	% make everything a column vector
	xpos = reshape(xpos, length(y)*length(x)/(stride*stride), 1);
	ypos = reshape(ypos, length(y)*length(x)/(stride*stride), 1);
	
	% the x direction is periodic, but the floats don't know they wrapped around. 
	xpos = mod( xpos-min(x), xDomainLength ) + min(x);

	% default color map is only 128 shades---we need more!
	colormap(jet(1024))
	
	% now plot the floats, colored by initial position
	% Scatter works, but is substantially slower than using mesh.
	% scatter(xpos, ypos, floatSize*floatSize, xposInitial, 'filled')	
	mesh([xpos';xpos'],[ypos';ypos'],[xposInitial';xposInitial'],'mesh','column','marker','.','MarkerSize',floatSize*floatSize), view(2)
	grid off
	
	% make the axes look better
	set( gca, 'TickDir', 'out');
	set( gca, 'Linewidth', 1.0);
	axis equal tight
	
	% get rid of the xticks because we're going to put a colorbar below with the same info.
	set( gca, 'xtick', [])
	
	% we need to extended the axis *slightly* beyond the min and max to accommodate the size of the float we've drawn
	delta = abs( x(2)-x(1) )*1.5;
	xlim([min(x)-delta max(x)+delta])
	ylim([min(y)-delta max(y)+delta])
	
	% label everything
	title( sprintf('Floats advected by an unstable shear flow, colored by initial position @ %03d %s', round(t(iTime)), tunits), 'fontsize', 28, 'FontName', 'Helvetica' );
	ylabel( 'distance (meters)', 'FontSize', 24.0, 'FontName', 'Helvetica');
	
	% add a color bar
	cb = colorbar( 'location', 'SouthOutside' );
	set(get(cb,'xlabel'),'String', 'distance (meters)', 'FontSize', 24.0, 'FontName', 'Helvetica');
	set( gca, 'clim', [min(x) max(x)] );
    
	% write everything out
    output = sprintf('%s/t_%03d', FramesFolder,iTime-1);
    export_fig(output,'-r300')
% 	print('-depsc2', output)
end