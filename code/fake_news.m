close all
clear all
clc

%% defining global parameters

% Agent container is a Nx5 structure, where N is the size of network 
% and each agent has 5 parameters.
%
% Agent parameters:
% 1. # close friends
% 2. # facebook friends
% 3. type of agent: news source, normal, expert [-1 0 1]
% 4. susceptibility
% 5. opinion: [-1 .. 1]

network_size = 1e2;                            % total number of agents in a network

connectivity_factor = 0.1;                     % mean number of connections in the network
                                               % normalized to network size

mean_connections = network_size * ...
                   connectivity_factor;        % Gaussian distribution
dev_connections = mean_connections * 0.1;

num_parameters = 5;

agents = zeros(network_size, num_parameters);   % agents vector
Sim_network = zeros(network_size);       % adjacency matrix of the network



%% forming a network
% network is represented in a form of an upper triangular adjacency matrix

for xx = 2:length(Sim_network(1,:))
    
    if (agents(xx, 1) + agents(xx, 2) == mean_connections)
        continue;
    end
    
    for yy = 1:xx - 1
        if (agents(yy, 1) + agents(yy, 2) == mean_connections)
            continue;
        end
        r = rand;
        if (r > 0.5)
            Sim_network(yy,xx) = 1;
            agents(xx, 1) = agents(xx, 1) + 1;
            agents(yy, 1) = agents(xx, 1) + 1;
        else
            Sim_network(yy,xx) = 0;
        end
    end
end

G = graph(Sim_network,'upper');
figure
plot(G)
axis square
axis off
title('Network')

