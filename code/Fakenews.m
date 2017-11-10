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



mean_connections = network_size * connectivity_factor;        % Gaussian distribution

dev_connections = mean_connections * 0.1;



num_parameters = 7;



agents = zeros(network_size, num_parameters);   % agents vector

Sim_network = zeros(network_size);       % adjacency matrix of the network



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
        agents(thisagent, 2) = randi(maxfbfriends);   %#fb friends, double than normal person
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
        
        agents(thisagent, 1) = randi(maxclosefriends);   %#close friends
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