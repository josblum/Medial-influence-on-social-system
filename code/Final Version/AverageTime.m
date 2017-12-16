%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File for the Course: MSSSM
% Students:
% 28.11.2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Looping a 1000 times
tic
%% Loading a Network
network_size = 3e2;                             % total number of agents in a network


reconnection_prob(1,1)=0.8;                     %"Randomness-factor" for erdosRenyi-Model (close FS)
reconnection_prob(2,1)=0.8;                     %"Randomness-factor" for erdosRenyi-Model (FB-friends)
mean_nbr_of_friends_per_agent(1,1)=10;           %Mean number of close friends per agent (even number!)
mean_nbr_of_friends_per_agent(2,1)=50;           %Mean number of FB-friends per agent (even number!)

% Setting parameters
% 1 means to evaluate
% 0 means to do nothing there
NewNetwork = 0; 
Influence_low = 1;
Influence_middle = 0;
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

for loop = 1:2
    % Setting all influenced agents back to normal
    % Setting all opinions of agents back to 0
    agents(:,3) = zeros(length(agents),1);
    agents(:,5) = zeros(length(agents),1);
%% Finding the right agents to be influenced

Number_Influenced = 30;

% Setting the influence parameter of the agents
    for i = 1:Number_Influenced
        thisagent = randi(network_size);
        if (agents(thisagent,3) == -1)
            --i;       
        else
        agents(thisagent,3) = -1;
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

% Find the time at which the overall averaged opinion is -0.9
opinion_average = zeros(1,timesteps);
for i = 1:timesteps
   opinion_average(i) = mean(opinion_development(:,i));    
end
opinion_average = opinion_average + 0.9;
[opinion_90(loop), t_90(loop)] = min(abs(opinion_average));

number_cf{loop} = agents(:,1);
number_fb{loop} = agents(:,2);
end
toc