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
network_size = 1e2;                            % total number of agents in a network


maxclosefriends = 5;                    %maximum #of closefriends
maxfbfriends = 30;                       %maximum #of fbfriends
minclosefriendships = 0;                 %Set minimal amount of formed close friendships
minfacebookFS = 5;                       %Set minimal amount of formed FB-friendships

friends = [maxclosefriends; maxfbfriends; minclosefriendships; minfacebookFS];

% Setting parameters
% 1 means to evaluate
% 0 means to do nothing there
NewNetwork = 0; 
Influence_low = 0;
Influence_middle = 0;
Influence_high = 1;


if (NewNetwork)
    [agents, Sim_network] = CreateaNetwork(network_size, friends);
end

load('agents.mat')
load('Sim_network')
toc

% Analysing the network

number_mean_CF = mean(agents(:,1));
number_mean_FB = mean(agents(:,2));
%%
Ctrs_CF = [0:1:5];
Xcts_CF = hist(agents(:,1), Ctrs_CF);
Ctrs_FB = [0:1:30];
Xcts_FB = hist(agents(:,2), Ctrs_FB);

figure 
bar(Ctrs_CF,Xcts_CF)

figure
bar(Ctrs_FB,Xcts_FB)

%% Finding the right angents to be influenced

Number_Influenced = 20;

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

timesteps = 2000;                               % Number of Network Simulations
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

if (Influence_low)
    figure
    for xx = 1:network_size
        plot(opinion_development(xx,:))
        hold on
    end
    xlabel('Time iterations')
    ylabel('Opinion')
    title('Opinion development when influencing low-connected agents')
end

if (Influence_middle)
    figure
    for xx = 1:network_size
        plot(opinion_development(xx,:))
        hold on
    end
    xlabel('Time iterations')
    ylabel('Opinion')
    title('Opinion development when influencing medium-connected agents')
end

if (Influence_high)
    figure
    for xx = 1:network_size
        plot(opinion_development(xx,:))
        hold on
    end
    xlabel('Time iterations')
    ylabel('Opinion')
    title('Opinion development when influencing highly connected agents')
end