%% Import and manipulate data
% Load data from Mathematica
load('Epoints.mat');

% Store the histogram without plotting
[N,edges]=histcounts(Expression1,300,'Normalization','pdf');

% Renormalize the counts because we also have the left branch
N = N/2;

% Taking mean values of edges two by two in order to get centers
centers = (edges(1 : end - 1) + edges(2 : end)) / 2;

% Calculating the max of the histogram
[Nmax,Nmaxidx] = max(N);

% Create figure
figure1 = figure('PaperType','A2');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot of the histogram
plot1 = plot([-flip(centers) centers],[flip(N) N],'DisplayName','Density of states','LineWidth',3,...
    'Color',[0.929411768913269 0.694117665290833 0.125490203499794]);


%% Linear fitting part
% Select xdata and ydata for linear fit as half of the data
% before the maximum
xdata1 = centers(1:floor(Nmaxidx/2));
ydata1 = N(1:floor(Nmaxidx/2));

% Make sure data are column vectors
xdata1 = xdata1(:);
ydata1 = ydata1(:);

% Remove NaN values and warn
nanMask1 = isnan(xdata1(:)) | isnan(ydata1(:));
if any(nanMask1)
    warning('GeneratedCode:IgnoringNaNs', ...
        'Data points with NaN coordinates will be ignored.');
    xdata1(nanMask1) = [];
    ydata1(nanMask1) = [];
end

% Setting linear fit plot limits
xplot1 = linspace(min(centers), centers(Nmaxidx));
% Find x values for plotting the fit based on xlim
%axesLimits1 = xlim(axes1);
%xplot1 = linspace(axesLimits1(1), axesLimits1(2));

% Find coefficients for polynomial (order = 1)
fitResults1 = polyfit(xdata1,ydata1,1);
% Evaluate polynomial
yplot1 = polyval(fitResults1,xplot1);
% Plot the fit
fitLine1 = plot(xplot1,yplot1,'DisplayName','Linear Fit','Tag','linear',...
    'Parent',axes1,...
    'LineWidth',3,...
    'LineStyle','--',...
    'Color',[0.850980401039124 0.325490206480026 0.0980392172932625]);


%% Log fitting part
% Select xdata and ydata for log fit as the data
% after the maximum
xdata2 = centers(Nmaxidx:end);
ydata2 = N(Nmaxidx:end);

% Make sure data are column vectors
xdata2 = xdata2(:);
ydata2 = ydata2(:);

% Remove NaN values and warn
nanMask2 = isnan(xdata2(:)) | isnan(ydata2(:));
if any(nanMask2)
    warning('GeneratedCode:IgnoringNaNs', ...
        'Data points with NaN coordinates will be ignored.');
    xdata1(nanMask2) = [];
    ydata1(nanMask2) = [];
end

% Setting log fit plot limits
xplot2 = linspace(centers(floor(4/5*Nmaxidx)),max(centers));

% Set up fittype and options.
ft = fittype( '-a*log(x-b)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0.0001 0 -Inf];
opts.StartPoint = [0.1 0.777 0.3];
opts.Upper = [Inf 1 Inf];

% Fit model to data.
[fitResults2, gof] = fit( xdata2, ydata2, ft, opts );

% Evaluate log
yplot2 = fitResults2(xplot2);

% Plot the fit
fitLine2 = plot(xplot2,yplot2,'DisplayName','Log Fit','Tag','linear',...
    'Parent',axes1,...
    'LineWidth',3,...
    'LineStyle','--',...
    'Color',[0 0.447058826684952 0.74117648601532]);


%% Set new line in proper position
% Get the axes children
hChildren = get(axes1,'Children');
% Remove the new line
hChildren(hChildren==fitLine1) = [];
% Get the index to the associatedLine
lineIndex = find(hChildren==plot1);
% Reorder lines so the new line appears with associated data
hNewChildren = [hChildren(1:lineIndex-1);fitLine1;hChildren(lineIndex:end)];
% Set the children:
set(axes1,'Children',hNewChildren);


%% Adjusting axes and related stuff
% Create xlabel
xlabel('E [in units of t]');

% Create ylabel
ylabel('g(E)');

% Set box
box(axes1,'on');

% Set lower limit
ylim(axes1,[0 inf]);

% Set the remaining axes properties
set(axes1,'FontSize',16);

% Create legend
legend(axes1,'show');
