close all;
clear all;
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
% 6. # close friends running variable for network generation
% 7. # facebook friends running variable for network generation

network_size = 1e2;                            % total number of agents in a network

connectivity_factor = 0.1;                     % mean number of connections in the network
                                               % normalized to network size

mean_connections = network_size * ...
                   connectivity_factor;        % Gaussian distribution
dev_connections = mean_connections * 0.1;

num_parameters = 7;

agents = zeros(network_size, num_parameters);   % agents vector
Sim_network = zeros(network_size);       % adjacency matrix of the network

%%%%%%Workaround%%%%%%
agents(:,1)=randi([1 5], network_size, 1); %Set #close friends
agents(:,2)=randi([1 20], network_size, 1); %Set #Facebook friends
%%%%%%%Workaround%%%%%%%

%% forming a network
% network is represented in a form of an upper triangular adjacency matrix
r=0; %initalize random variable

%Forming close friend network
while((sum(agents(:,6))+5) < (sum(agents(:,1)))) %Continues until all friendships are formed  %%%"Bug" need to add maximal number of friends
    for xx = 1:network_size
    
       if ( agents(xx,6) < agents(xx, 1)) %Check if agent needs more friends

        t=0;
         while(t==0)
         r = randi([1 network_size],1,1);  %Random integer between 1 and number of agents
         if(r ~= xx && agents(r,1) > agents(r,6) && Sim_network(xx,r) == 0 && Sim_network(r,xx)== 0)  %Check if friendship is possible
         t=1;
         end
         end

            if(xx<r) 
             Sim_network(xx,r) = 2;
            end
            if(r<xx) %Make sure to have upper triangular matrix
                Sim_network(r,xx) = 2;
            end
            
                agents(xx, 6) = agents(xx, 6) + 1;%increase nbr of friendships
              agents(r, 6) = agents(r, 6) + 1; 
              
       
       end
    end

end

%forming facebook friend network
while((sum(agents(:,7))+20)<(sum(agents(:,2)))) %Continues until all friendships are formed          %%%"Bug" need to add maximal number of friends
    for xx = 1:network_size
    
 if( agents(xx,7) < agents(xx, 2)) %Check if current agent needs more friends
     
            t=0;
            while(t==0)
                r = randi([1 network_size],1,1);  %Random integer between 1 and number of agents
         if(r ~= xx && agents(r,2)>agents(r,7) && Sim_network(xx,r) ~=3 && Sim_network(xx,r) ~=1 && Sim_network(r,xx) ~= 3 && Sim_network(r,xx) ~= 1 )  %Check if friendship is possible
         t=1;
         else
         end
            end
            
            if(xx<r) 
             Sim_network(xx,r) = Sim_network(xx,r)+1;
            end
            if(r<xx) %Make sure to have upper triangular matrix
                Sim_network(r,xx) = Sim_network(r,xx)+1;
            end
            
                agents(xx, 7) = agents(xx, 7) + 1;%increase nbr of friendships done
              agents(r, 7) = agents(r, 7) + 1; 
 
    end
    end
end
FehlerCloseFriends=(sum(agents(:,1))-sum(agents(:,6)));
FehlerFBFriends=(sum(agents(:,2))-sum(agents(:,7)));

    



G = graph(Sim_network,'upper');
figure
plot(G)
axis square
axis off
title('Network')

