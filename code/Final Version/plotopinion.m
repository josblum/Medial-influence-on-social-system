function p = plotopinion(x, timesteps, percent, p4)

p1 = x(:,1);
p2 = x(:,2);
p3 = x(:,3);
p5 = x(:,4);
finaltitle = ['\fontsize{20}Change of mean opinion if ' num2str(percent) char(37) ' get initially influenced' ];

%define figure properties
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.width      = 15;
opts.height     = 8;
opts.fontType   = 'Times';
opts.fontSize   = 9;
%%
% create new figure
fig = figure; clf
plot(p1,'LineWidth', 2);
hold on;
plot(p2,'LineWidth', 2);
plot(p3,'LineWidth', 2);
plot(p5,'LineWidth', 2);
plot(p4,'LineWidth', 2);
% add axis labes and legend
iFontSize = 14;
strFontUnit = 'points'; % [{points} | normalized | inches | centimeters | pixels]
strFontName = 'Times';  % [Times | Courier | ]              TODO complete the list
strFontWeight = 'normal'; % [light | {normal} | demi | bold]
strFontAngle = 'normal'; % [{normal} | italic | oblique]     ps: only for axes
strInterpreter = 'latex';  % [{tex} | latex]
fLineWidth = 1.0;      % width of the line of the axes
set(gca,'XGrid','on','YGrid','on','GridLineStyle',':',...
    'FontName',strFontName,'FontSize',iFontSize,'FontUnits',strFontUnit,'FontWeight',...
    strFontWeight,'FontAngle',strFontAngle,'LineWidth',fLineWidth)
axis tight;
xlabel('\fontsize{14}Timesteps');
ylabel('\fontsize{14}Average Opinion');


axis([50 timesteps -1 0.1]);

title(finaltitle);

% scaling
fig.Units               = 'centimeters';
fig.Position(3)         = opts.width;
fig.Position(4)         = opts.height;

% set text properties
set(fig.Children, ...
     'FontName',     'Times', ...
     'FontSize',     9);

% remove unnecessary white space
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))

fig.PaperPositionMode   = 'auto';

leg = legend('Worst Connected Agents', 'Average Connected Agents', 'Best Connected Agents', 'Randomly chosen Agents');
leg.Title.String = 'Legend:';
leg.Title.FontSize = 16;
leg.FontSize = 14;

end