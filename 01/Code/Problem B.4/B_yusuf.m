% Load data from Mathematica
%load('Bands.mat')
% Decomment to read Yusuf's output instead of mine
Expression1 = dlmread('eigenval.out');

% Specify the number of cells
Ncells = 50;

% Get the discretized kx and the sorted energy values
kx = reshape(Expression1(:,1),4*Ncells,length(Expression1(:,1))/(4*Ncells));
bandpoints = reshape(Expression1(:,2),4*Ncells,length(Expression1(:,2))/(4*Ncells));
%kx = kx(1,:);
bandpoints_sorted = sort(bandpoints);

% Create figure
figure1 = figure('PaperType','A2','Renderer','painters');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
title('Armchair nanotube energy bands','FontSize',18);

for ii = 1:4*Ncells
 plot(kx(ii,:), bandpoints_sorted(ii,:),...
     'LineWidth',1)
end

%% Adjusting axes and related stuff
% Create xlabel
xlabel('$k_x$ [in units of $1/a$]','Interpreter','latex');

% Create ylabel
ylabel('$E-\varepsilon$ [in units of $t$]','Interpreter','latex');

% Axes limits
xlim(axes1,[-pi pi]);
ylim(axes1,[-3.5 3.5]);

% Set box
box(axes1,'on');

% Set the remaining axes properties
set(axes1,'FontSize',16,'XGrid','on','XTick',[-3*pi/2 -pi -pi/2 0 pi/2 pi 3*pi/2],...
    'XTickLabel',{'-3\pi/2','-\pi','-\pi/2','0','\pi/2','\pi','3\pi/2'},...
    'YAxisLocation','left',...
    'YGrid','on');
