function [agents, Sim_network] = CreateaNetwork(network_size, friends)
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








connectivity_factor = 0.1;                     % mean number of connections in the network

                                               % normalized to network size



mean_connections = network_size * connectivity_factor;        % Gaussian distribution

dev_connections = mean_connections * 0.1;



num_parameters = 5;



agents = zeros(network_size, num_parameters);   % agents vector

Sim_network = zeros(network_size);       % adjacency matrix of the network



%%%%%%Workaround%%%%%%

maxclosefriends = friends(1);   %maximum #of closefriends
maxfbfriends = friends(2);      %maximum #of fbfriends
minclosefriendships = friends(3); %Set minimal amount of formed close friendships
minfacebookFS = friends(4); %%Set minimal amount of formed FB-friendships

%%% SET THE AGENTS

for i = 1:network_size  %Randomly Distributes the newsspreader

       if (agents(i,3) ~= -1 && agents(i,3) ~= 1)   %assures that specialists and newsspreader are not overwriten
           
        agents(i, 1) = 0;   %#close friends
        agents(i, 2) = 0;   %#fb friends, double than normal person
        agents(i, 3) = 0;  %Type normal
        agents(i, 4) = rand();   %susceptibility
        agents(i, 5) = 0;%rand()*2-1;   %initial opinion to 0
       
              
       end
end

%% forming a network

% network is represented in a form of an upper triangular adjacency matrix

r=0; %initalize random variable



%Forming close friend network
TotalNbrOfCF=randi([minclosefriendships,maxclosefriends*network_size],1,1);%Set random amount of friendships

while((sum(agents(:,1))+1) < TotalNbrOfCF) %Continues until all friendships are formed

    for xx = 1:network_size
       if ( agents(xx,1) < maxclosefriends &&  (sum(agents(:,1))+1<TotalNbrOfCF)) %Check if agent can have more friends
        t=randi(10);
        if(t<3)
            t=0;
            
         while(t==0) %Chose possible partner for friendship
         r = randi([1 network_size],1,1);
         if(r ~= xx && agents(r,1) < maxclosefriends && Sim_network(xx,r) == 0 && Sim_network(r,xx)== 0)  %Check if friendship is possible
         t=1;
         end
         end
        
            if(xx<r) %Make sure to have upper triangular matrix
             Sim_network(xx,r) = 2;
            end
            if(r<xx) 
                Sim_network(r,xx) = 2;
            end
            
              agents(xx, 1) = agents(xx, 1) + 1;%increase nbr of friendships
              agents(r, 1) = agents(r, 1) + 1; 
        end
       end
    end
end



%forming facebook friend network

TotalNbrOfFF=randi([minfacebookFS,maxfbfriends*network_size],1,1);%Set random amount of friendships

while((sum(agents(:,2))+1) < TotalNbrOfFF) %Continues until all friendships are formed

    for xx = 1:network_size
       if ( agents(xx,2) < maxfbfriends && (sum(agents(:,2))+1<TotalNbrOfFF)) %Check if agent can have more friends
        t=randi(10);
        if(t<3)
            t=0;
         while(t==0)
           r = randi([1 network_size],1,1);  %Random integer between 1 and number of agents
             if(r ~= xx && maxfbfriends>agents(r,2) && Sim_network(xx,r) ~=3 && Sim_network(xx,r) ~=1 && Sim_network(r,xx) ~= 3 && Sim_network(r,xx) ~= 1 )  %Check if friendship is possible
             t=1;
             end
         end
         
         if(xx<r)
             Sim_network(xx,r) = Sim_network(xx,r)+1;
         end
         if(r<xx) %Make sure to have upper triangular matrix
             Sim_network(r,xx) = Sim_network(r,xx)+1;
         end
              
        agents(xx, 2) = agents(xx, 2) + 1;%increase nbr of friendships
        agents(r, 2) = agents(r, 2) + 1; 
        end
       end
    end
end

TestNetworkCF=TotalNbrOfCF-sum(agents(:,1));
TestNetworkFB=TotalNbrOfFF-sum(agents(:,2));
if(TestNetworkCF>1 || TestNetworkFB>1)
    print('Error forming network')
end


%%Check if every agents has at least 1 friend
for xx=1:network_size
   if(agents(xx,1)==0)
       while(r==xx)
       r=randi([0,network_size],1,1);
       end
       if(xx<r)
       Sim_network(xx,r)=Sim_network(xx,r)+2;
       agents(xx,1)=agents(xx,1)+1;
       agents(r,1)=agents(r,1)+1;
       end
       if(xx>r)
           Sim_network(r,xx)=Sim_network(r,xx)+2;
                 agents(xx,1)=agents(xx,1)+1;
       agents(r,1)=agents(r,1)+1;
       end
       
   end
   if(agents(xx,2)==0)
       while(r==xx)
       r=randi([0,network_size],1,1);
       end
       if(xx<r)
       Sim_network(xx,r)=Sim_network(xx,r)+1;
              agents(xx,2)=agents(xx,2)+1;
       agents(r,2)=agents(r,2)+1;
       end
       if(xx>r)
           Sim_network(r,xx)=Sim_network(r,xx)+1;
                  agents(xx,2)=agents(xx,2)+1;
       agents(r,2)=agents(r,2)+1;
       end
       
   end
    
end

G = graph(Sim_network,'upper');

figure

plot(G)

axis square

axis off

title('Network')

save('agents.mat','agents')
save('Sim_network.mat','Sim_network')
end

