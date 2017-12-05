function [agents, Sim_network, struct_CF, struct_FB] = CreateErdosRenyiNetwork(network_size, reconnection_prob, mean_nbr_of_friends_per_agent)

%% defining global parameters



% Agent container is a Nx5 structure, where N is the size of network 

% and each agent has 5 parameters.

%

% Agent parameters:

% 1. # close friends

% 2. # facebook friends

% 3. Status of influence: news source, normal, expert [-1 0 1]

% 4. susceptibility

% 5. opinion: [-1 .. 1]

num_parameters = 5;

agents = zeros(network_size, num_parameters);   % agents vector

Sim_network = zeros(network_size, network_size);       % adjacency matrix of the network

%%% SET THE AGENTS

for i = 1:network_size           
        agents(i, 1) = 0;           % #close friends
        agents(i, 2) = 0;           % #fb friends, double than normal person
        agents(i, 3) = 0;           % Type normal
        agents(i, 4) = rand();      % susceptibility
        agents(i, 5) = 0;           % initial opinion to 0   
end

%Redefine values for function call
       reconnection_prob_CF=reconnection_prob(1);
       mean_nbr_of_friends_per_agent_CF=mean_nbr_of_friends_per_agent(1);
       reconnection_prob_FB=reconnection_prob(2);
       mean_nbr_of_friends_per_agent_FB=mean_nbr_of_friends_per_agent(2);
       
%Call function to create a random network
struct_CF=erdosRenyi(network_size,reconnection_prob_CF,mean_nbr_of_friends_per_agent_CF);
struct_FB=erdosRenyi(network_size,reconnection_prob_FB,mean_nbr_of_friends_per_agent_FB);

%Update parameter of agents
for i = 1:network_size
agents(i,1)=sum(struct_CF.network(:,i))+sum(struct_CF.network(i,:));
agents(i,2)=sum(struct_FB.network(:,i))+sum(struct_FB.network(i,:));
end

Sim_network=struct_CF.network.*2 + struct_FB.network;





save('agents.mat','agents')
save('Sim_network.mat','Sim_network')






