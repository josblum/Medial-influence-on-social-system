%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File for the Course: MSSSM
% Students:
% 28.11.2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% Loading a Network
tic
network_size = 3e2;                             % total number of agents in a network


reconnection_prob(1,1)=0.8;                     %"Randomness-factor" for erdosRenyi-Model (close FS)
reconnection_prob(2,1)=0.8;                     %"Randomness-factor" for erdosRenyi-Model (FB-friends)
mean_nbr_of_friends_per_agent(1,1)=10;           %Mean number of close friends per agent (even number!)
mean_nbr_of_friends_per_agent(2,1)=50;           %Mean number of FB-friends per agent (even number!)

% Setting parameters
% 1 means to evaluate
% 0 means to do nothing there
NewNetwork = 0; 
Influence_low = 0;
Influence_middle = 1;
Influence_high = 0;


if (NewNetwork)
    [agents, Sim_network, struct_CF, struct_FB] = CreateErdosRenyiNetwork(network_size, reconnection_prob, mean_nbr_of_friends_per_agent);
   
    %plotGraphBasic(struct_CF,2,1)
    %plotGraphBasic(struct_CF,2,1,'FB friends network')
end

load('agents.mat')
load('Sim_network')

% Analysing the network

number_mean_CF = mean(agents(:,1));
number_mean_FB = mean(agents(:,2));

Ctrs_CF = [0:1:mean_nbr_of_friends_per_agent(1,1)*2];
Xcts_CF = hist(agents(:,1), Ctrs_CF);
Ctrs_FB = [0:1:mean_nbr_of_friends_per_agent(2,1)*2];
Xcts_FB = hist(agents(:,2), Ctrs_FB);

% figure 
% subplot(2,1,1)
% bar(Ctrs_CF,Xcts_CF)
% 
% subplot(2,1,2)
% bar(Ctrs_FB,Xcts_FB)

%% Finding the right angents to be influenced

Number_Influenced = 30;

for i = 1:length(agents)
    Friendship_Sum(i) = 2*agents(i,1)+agents(i,2);
end

[Friendship_Sum_order, idx_order] = sort(Friendship_Sum,'ascend');

% Find low-connected and high-connected agents

LowConnedcted = Friendship_Sum_order(1:Number_Influenced);
idx_low = idx_order(1:Number_Influenced);
HighConnected = Friendship_Sum_order(end-Number_Influenced:end);
idx_high = idx_order(end-Number_Influenced+1:end);

% Find medium-connected agents

MeanConnection = floor(sum(Friendship_Sum)/length(Friendship_Sum));
idx_order_mean = find(Friendship_Sum_order == MeanConnection);

missing = Number_Influenced - length(idx_order_mean);

if (missing > 0)
    idx_order_mean_start = idx_order_mean(1) - floor(missing/2);
    idx_order_mean_end = idx_order_mean(end) + floor(missing/2);
    idx_order_mean_full = idx_order(idx_order_mean_start:idx_order_mean_end);    
end

if (missing <= 0)
    idx_order_mean_full = idx_order(idx_order_mean(1:Number_Influenced));   
end


%% Setting the influence parameter of the agents

if(Influence_low)
    for i = 1:Number_Influenced
        agents(idx_low(i),3) = -1;

    end
end

if(Influence_middle)
    for i = 1:Number_Influenced
        agents(idx_order_mean_full(i),3) = -1;

    end
end

if(Influence_high)
    for i = 1:Number_Influenced
        agents(idx_high(i),3) = -1;

    end
end

%% Starting simulation and updating oppinion of every agent

timesteps = 2500;                               % Number of Network Simulations
temp_opinion = zeros(network_size,1);        % Initializing a temporaray opinion vector 
k_CF = 2;                                    % Opinion multiplication factor for CF_friends
k_FB = 1;                                  % Opinion multiplication factor for FB_friends
k_IF = 5;                                   % Opinion multiplication factor for Influence
opinion_development = zeros(network_size,timesteps);

for t = 1:timesteps
    for xx = 1:network_size     % Iterating through all agents in the network
                
            listCF_hor = find(Sim_network(xx,:) == 2 | Sim_network(xx,:) == 3);
            listCF_ver = find(Sim_network(:,xx) == 2 | Sim_network(:,xx) == 3);
            listFB_hor = find(Sim_network(xx,:) == 1 | Sim_network(xx,:) == 3);
            listFB_ver = find(Sim_network(:,xx) == 1 | Sim_network(:,xx) == 3);
            
            
            temp_opinionCF_hor = agents(listCF_hor,5);
            temp_opinionCF_ver = agents(listCF_ver,5);
            temp_opinionCF = [temp_opinionCF_hor; temp_opinionCF_ver];
            participantsCF = length(find(temp_opinionCF ~= 0));
            temp_opinionCF_sum = sum(temp_opinionCF);
            
            temp_opinionFB_hor = agents(listFB_hor,5);
            temp_opinionFB_ver = agents(listFB_ver,5);
            temp_opinionFB = [temp_opinionFB_hor; temp_opinionFB_ver];
            participantsFB = length(find(temp_opinionFB ~= 0));
            temp_opinionFB_sum = sum(temp_opinionFB);
            
            if (participantsCF + participantsFB == 0)           
                enhance_friends_number = 1;
            else
                enhance_friends_number = 0;
            end
                       
                if (agents(xx,3) == 0) 
                    % His own opinion 
                    temp_opinion(xx) = (1-agents(xx,4)) * agents(xx,5) + agents(xx,4) * ((temp_opinionCF_sum * k_CF + temp_opinionFB_sum * k_FB) / (participantsCF*k_CF +participantsFB*k_FB+enhance_friends_number));
                end

                if (agents(xx,3) == -1)        % If agent is influenced by Specialist or Spreader
                    % His own opinion 
                    temp_opinion(xx) = (1-agents(xx,4)) * agents(xx,5) + agents(xx,4) * ((temp_opinionCF_sum * k_CF + temp_opinionFB_sum * k_FB + agents(xx,3)*k_IF) / (participantsCF*k_CF + participantsFB*k_FB + k_IF+enhance_friends_number));
                end
                
            
    end
    
    for xx = 1:network_size                                           
        agents(xx,5) = temp_opinion(xx);                    % Updating opinion for every agent
        opinion_development(xx,t) = temp_opinion(xx);       % Updating opinion_development over time
    end
end

% Find the time at which the overall averaged opinion is -
%% Producing Plots

% define figure properties
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'C:\Users\Dave\Desktop\Schule\7. Semester\Modeling of Social Systems';
opts.width      = 15;
opts.height     = 8;
opts.fontType   = 'Times';
opts.fontSize   = 9;

% create new figure
fig = figure; clf
hist(agents(:,1));

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
axis tight
xlabel('Close Friends');
ylabel('Agents');
title('Distribution of number of Close Friends')
xlim([0 max(agents(:,1))+1]);
% ylim([0 15]);

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


% add a second figure
fig = figure; clf
hist(agents(:,2));

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
axis tight
xlabel('Facebook Friends');
ylabel('Agents');
title('Distribution of number of Facebook Friends')
xlim([0 max(agents(:,2))+1]);
% ylim([0 15]);

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



print([opts.saveFolder 'Facebook Friends'], '-dpng', '-r600')

% add a second figure
fig = figure; clf
    for xx = 1:network_size
        plot(opinion_development(xx,:))
        hold on
    end

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
axis tight
xlabel('Time [Number of opinion updates]');
ylabel('Opinion');
if (Influence_high)
    title('Opinion development "Highly connected targets"')
end
if (Influence_low)
    title('Opinion development "Low connected targets"')
end
if (Influence_middle)
    title('Opinion development "Medium connected targets"')
end
    xlim([0 timesteps]);
% ylim([0 15]);

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
toc
%%
% if (Influence_low)
%     figure
%     for xx = 1:network_size
%         plot(opinion_development(xx,:))
%         hold on
%     end
%     xlabel('Time iterations')
%     ylabel('Opinion')
%     title('Opinion development when influencing low-connected agents')
% end
% 
% if (Influence_middle)
%     figure
%     for xx = 1:network_size
%         plot(opinion_development(xx,:))
%         hold on
%     end
%     xlabel('Time iterations')
%     ylabel('Opinion')
%     title('Opinion development when influencing medium-connected agents')
% end
% 
% if (Influence_high)
%     figure
%     for xx = 1:network_size
%         plot(opinion_development(xx,:))
%         hold on
%     end
%     xlabel('Time iterations')
%     ylabel('Opinion')
%     title('Opinion development when influencing highly connected agents')
% end