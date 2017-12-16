clear all;
clc;
timesteps = 2500;           % Amount of timesteps, default: 2500
percent = 3;                % Percent of initially influenced agents
NewNetwork = 1;             % 1 = generate new network, 0 = keep old network


%%%PARAMETERS DONT CHANGE%%%
network_size = 3e2;                             % total number of agents in a network
reconnection_prob(1,1)=0.8;                     %"Randomness-factor" for erdosRenyi-Model (close FS)
reconnection_prob(2,1)=0.8;                     %"Randomness-factor" for erdosRenyi-Model (FB-friends)
mean_nbr_of_friends_per_agent(1,1)=10;           %Mean number of close friends per agent (even number!)
mean_nbr_of_friends_per_agent(2,1)=50;           %Mean number of FB-friends per agent (even number!)


if (NewNetwork)   % Creates new network
    [agents, Sim_network, struct_CF, struct_FB] = CreateErdosRenyiNetwork(network_size, reconnection_prob, mean_nbr_of_friends_per_agent);
   
    %plotGraphBasic(struct_CF,2,1)
    %plotGraphBasic(struct_CF,2,1,'FB friends network')
end


p4 = [];
for i = 1: timesteps
    p4(i) = -0.9;
end


sol = main3(percent, timesteps, 0);              % 1 Percent with timesteps steps
sol(:,4) = main3(percent, timesteps,1);
plotopinion(sol, timesteps, percent, p4); % sol, timesteps, which plot, title
toc



