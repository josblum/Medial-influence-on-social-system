close all;

clear all;

clc





network_size = 1e2;                            % total number of agents in a network



connectivity_factor = 0.1;                     % mean number of connections in the network

                                               % normalized to network size



mean_connections = network_size *connectivity_factor;        % Gaussian distribution

dev_connections = mean_connections * 0.1;



num_parameters = 7;



agents = zeros(network_size, num_parameters);   % agents vector

Sim_network = zeros(network_size);       % adjacency matrix of the network


% Agent parameters:

% 1. # close friends
% 2. # facebook friends
% 3. type of agent: news source, normal, expert [-1 0 1]
% 4. susceptibility
% 5. opinion: [-1 .. 1]
% 6. # close friends running variable for network generation
% 7. # facebook friends running variable for network generation

%%%%%%Workaround%%%%%%

maxclosefriends = 10;   %maximum #of closefriends
maxfbfriends = 30;      %maximum #of fbfriends

newsspreader = 3;       %Sets the #fake news speader (Type: -1)
specialists = 5;         %Sets the #specialists (Type: 1)


%%% Set the NEWSSPRADER

for i = 1:newsspreader  %Randomly Distributes the newsspreader
    thisagent = randi(network_size);
    if (agents(thisagent,3) ~= -1)     %assures that newsspreader dont overwrite each other
        
        agents(thisagent, 1) = randi(maxclosefriends);   %#close friends
        agents(thisagent, 2) = maxfbfriends * 2;   %#fb friends, double than normal person
        agents(thisagent, 3) = -1;  %Type newsspreader
        agents(thisagent, 4) = 0;   %susceptibility
        agents(thisagent, 5) = -1;   %opinion
        agents(thisagent, 6) = 0;  %initial close friends to 0
        agents(thisagent, 7) = 0;  %initial fb friends to 0
        
    else
        --i;        
    end
    
end

%%% SET THE SPECIALISTS

for i = 1:specialists  %Randomly Distributes the newsspreader
    thisagent = randi(network_size);
       if (agents(thisagent,3) ~= -1 && agents(thisagent,3) ~= 1)   %assures that specialists dont overwrite each other or newsspreader
        
        agents(thisagent, 1) = maxclosefriends;   %#close friends
        agents(thisagent, 2) = randi(maxfbfriends);   %#fb friends, double than normal person
        agents(thisagent, 3) = 1;  %Type specialist
        agents(thisagent, 4) = 0;   %susceptibility
        agents(thisagent, 5) = 1;   %opinion
        agents(thisagent, 6) = 0;  %initial close friends to 0
        agents(thisagent, 7) = 0;  %initial fb friends to 0
        
    else
        --i;        
       end
end

%%% SET EVERYONE ELSE

for i = 1:network_size  %Randomly Distributes the newsspreader

       if (agents(i,3) ~= -1 && agents(i,3) ~= 1)   %assures that specialists and newsspreader are not overwriten
           
        agents(i, 1) = randi(maxclosefriends);   %#close friends
        agents(i, 2) = randi(maxfbfriends);   %#fb friends, double than normal person
        agents(i, 3) = 0;  %Type normal
        agents(i, 4) = rand();   %susceptibility
        agents(i, 5) = 0;   %initial opinion to 0
        agents(i, 6) = 0;  %initial close friends to 0
        agents(i, 7) = 0;  %initial fb friends to 0
       
              
       end
end
    


