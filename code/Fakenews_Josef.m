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





network_size = 1e2;                            % total number of agents in a network



connectivity_factor = 0.1;                     % mean number of connections in the network

                                               % normalized to network size



mean_connections = network_size * connectivity_factor;        % Gaussian distribution

dev_connections = mean_connections * 0.1;



num_parameters = 5;



agents = zeros(network_size, num_parameters);   % agents vector

Sim_network = zeros(network_size);       % adjacency matrix of the network



%%%%%%Workaround%%%%%%

maxclosefriends = 5;   %maximum #of closefriends
maxfbfriends = 10;      %maximum #of fbfriends
minclosefriendships=network_size*2; %Set minimal amount of formed close friendships
minfacebookFS=network_size*5; %%Set minimal amount of formed FB-friendships

newsspreader = 3;       %Sets the #fake news speader (Type: -1)
specialists = 0;         %Sets the #specialists (Type: 1)

MaxCFNewsspreader=10; %Maximum #of closefriends for newsspreader
MaxFFNewsspreader=30; %Maximum #of fbfriends for newsspreader
MinCFNewsspreader=5; %Minimal #of closefriends for newsspreader
MinFFNewsspreader=10; %Minimal #of fbfriends for newsspreader

MaxCFSpecialists=10; %Maximum #of closefriends for specialists
MaxFFSpecialists=30; %Maximum #of fbfriends for specialists
MinCFSpecialists=1; %Minimal #of closefriends for specialists
MinFFSpecialists=10; %Minimal #of fbfriends for specialists


%%% Set the NEWSSPRADER

for i = 1:newsspreader  %Randomly Distributes the newsspreader
    thisagent = randi(network_size);
    if (agents(thisagent,3) ~= -1)     %assures that newsspreader dont overwrite each other
        
        agents(thisagent, 1) = 0;   %#close friends
        agents(thisagent, 2) = 0;   %#fb friends, double than normal person
        agents(thisagent, 3) = -1;  %Type newsspreader
        agents(thisagent, 4) = 0;   %susceptibility
        agents(thisagent, 5) = -1;   %opinion
        
    else
        --i;        
    end
    
end

%%% SET THE SPECIALISTS

for i = 1:specialists  %Randomly Distributes the newsspreader
    thisagent = randi(network_size);
       if (agents(thisagent,3) ~= -1 && agents(thisagent,3) ~= 1)   %assures that specialists dont overwrite each other or newsspreader
        
        agents(thisagent, 1) = 0;   %#close friends
        agents(thisagent, 2) = 0;   %#fb friends, double than normal person
        agents(thisagent, 3) = 1;  %Type specialist
        agents(thisagent, 4) = 0;   %susceptibility
        agents(thisagent, 5) = 1;   %opinion
        
    else
        --i;        
       end
end

%%% SET EVERYONE ELSE

for i = 1:network_size  %Randomly Distributes the newsspreader

       if (agents(i,3) ~= -1 && agents(i,3) ~= 1)   %assures that specialists and newsspreader are not overwriten
           
        agents(i, 1) = 0;   %#close friends
        agents(i, 2) = 0;   %#fb friends, double than normal person
        agents(i, 3) = 0;  %Type normal
        agents(i, 4) = rand();   %susceptibility
        agents(i, 5) = 0;%rand()*2-1;   %initial opinion to 0
       
              
       end
end

%%%%%%%Workaround%%%%%%%



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

%%
%%Introduce specialists and newsspreader to network

%%%Cancel all existing friendships of specialists and newsspreader
for xx = 1:network_size
    if(agents(xx,3)==-1 || agents(xx,3) == 1) %Only for specialists and newsspreader
        for yy=1:network_size
        if(Sim_network(xx,yy)~=0)
            if(Sim_network(xx,yy) == 3) %if agent xx and yy are close and fb friends
                Sim_network(xx,yy) = 0;
                agents(yy,1)=agents(yy,1)-1; %Update nbr of friends that agent yy has
                agents(yy,2)=agents(yy,2)-1;
            end
            if(Sim_network(xx,yy) == 2) %if agent xx and yy are close friends
                Sim_network(xx,yy) = 0;
                agents(yy,1)=agents(yy,1)-1; %Update nbr of friends that agent yy has
            end
            if(Sim_network(xx,yy) == 1) %if agent xx and yy are fb friends
                Sim_network(xx,yy) = 0;
                agents(yy,2)=agents(yy,2)-1; %Update nbr of friends that agent yy has
            end
        end
        
                if(Sim_network(yy,xx)~=0)
            if(Sim_network(yy,xx) == 3) %if agent xx and yy are close and fb friends
                Sim_network(yy,xx) = 0;
                agents(yy,1)=agents(yy,1)-1; %Update nbr of friends that agent yy has
                agents(yy,2)=agents(yy,2)-1;
            end
            if(Sim_network(yy,xx) == 2) %if agent xx and yy are close
                Sim_network(yy,xx) = 0;
                agents(yy,1)=agents(yy,1)-1; %Update nbr of friends that agent yy has
            end
            if(Sim_network(yy,xx) == 1) %if agent xx and yy are fb friends
                Sim_network(yy,xx) = 0;
                agents(yy,2)=agents(yy,2)-1; %Update nbr of friends that agent yy has
            end
                end
        end
    end
end

%%
%Test implementation
testCF=0;
testFF=0;
for xx=1:network_size
if(agents(xx,3)==1 || agents(xx,3)==-1)
    testCF=testCF + agents(xx,1);
    testFF=testFF + agents(xx,2);
end
end
TestNetworkCF_2=TotalNbrOfCF-sum(agents(:,1))-testCF;
TestNetworkFB_2=TotalNbrOfFF-sum(agents(:,2))-testFF;
if(TestNetworkCF_2>TestNetworkCF || TestNetworkFB_2>TestNetworkFB)
    print('error erasing friendships of newsspreader & specialists')
end

%%
%%Create close friends for newsspreader & specialists:

for xx=1:network_size
    if(agents(xx,3)==-1 || agents(xx,3)==1) %Check if agent is newsspreader/specialist
        agents(xx,1)=sum(Sim_network(xx,:))+sum(Sim_network(:,xx)); %update number of close friends
        if(agents(xx,3)==-1)
        ss=randi([max(MinCFNewsspreader,agents(xx,1)),max(agents(xx,1),MaxCFNewsspreader)],1,1); %Gives nbr of CF of each Newsspreader; Minimal value is either specified minimum or already existing friendships-->Only to newsspreader or specialists
        end
        if(agents(xx,3)==1)
            ss=randi([max(MinCFSpecialists,agents(xx,1)),max(agents(xx,1),MaxCFSpecialists)],1,1); %Gives nbr of CF of each Specialist
        end
        
        
        while(agents(xx,1)<ss)        
            t=0;
         while(t==0)
         r = randi([1 network_size],1,1);  %Random integer between 1 and number of agents
         if(r ~= xx && Sim_network(xx,r) == 0 && Sim_network(r,xx)== 0)  %Check if friendship is possible
         t=1;
         end
         end
        
            if(xx<r) 
             Sim_network(xx,r) = 2;
            end

            if(r<xx) %Make sure to have upper triangular matrix
                Sim_network(r,xx) = 2;
            end
            
              agents(xx, 1) = agents(xx, 1) + 1;%increase nbr of friendships
              agents(r, 1) = agents(r, 1) + 1; 
        end
        if(agents(xx,1) ~= ss)
            print('Error forming close friendships of Specialists and Newsspreader')
        end
    end
end


        %%%Create FB-friends for newsspreader & specialists:

for xx=1:network_size
    if(agents(xx,3)==-1 || agents(xx,3)==1)
        
        agents(xx,2)=0;%update number of friends
        for yy=1:network_size
            if(Sim_network(xx,yy)==3 || Sim_network(yy,xx)==3 || Sim_network(xx,yy)==1 || Sim_network(yy,xx)==1)
            agents(xx,2)=agents(xx,2)+1;
            end
        end
            
        if(agents(xx,3)==-1)
        ss=randi([max(MinFFNewsspreader,agents(xx,2)),max(agents(xx,2),MaxFFNewsspreader)],1,1); %Gives nbr of CF of each Newsspreader; Minimal value is either specified minimum or already existing friendships-->Only to newsspreader or specialists
        end
        if(agents(xx,3)==1)
            ss=randi([max(MinFFSpecialists,agents(xx,2)),max(agents(xx,2),MaxFFSpecialists)],1,1); %Gives nbr of CF of each Specialist
        end
                
        while(agents(xx,2)<ss)
            t=0;
            while(t==0)
                r = randi([1 network_size],1,1);  %Random integer between 1 and number of agents
                if(r ~= xx && Sim_network(xx,r) ~=3 && Sim_network(xx,r) ~=1 && Sim_network(r,xx) ~= 3 && Sim_network(r,xx) ~= 1 )  %Check if friendship is possible
                    t=1;
                end
            end
            if(xx<r) 
                Sim_network(xx,r) = Sim_network(xx,r)+1;
            end
            if(r<xx) %Make sure to have upper triangular matrix
                Sim_network(r,xx) = Sim_network(xx,r)+1;
            end
            
            agents(xx, 2) = agents(xx, 2) + 1;%increase nbr of friendships
            agents(r, 2) = agents(r, 2) + 1;
        end
        if(agents(xx,2)~=ss)
            print('Error forming FB-network of specialists and newsspreader')
        end
    end
end
%%Check if agents has 1 friend
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

%% Starting Simulation

timesteps = 1000;                               % Number of Network Simulations
temp_opinion = zeros(network_size,1);        % Initializing a temporaray opinion vector 
k_CF = 0.8;                                    % Opinion multiplication factor for CF_friends
k_FB = 1-k_CF;                                  % Opinion multiplication factor for FB_friends
opinion_development = zeros(network_size,timesteps);

for t = 1:timesteps
    for xx = 1:network_size     % Iterating through all agents in the network
        if (agents(xx,3) == -1 || agents(xx,3) == 1)        % If agent is specialist or Newsspreader
            temp_opinion(xx) = agents(xx,5);                    % His opinion stays the same
        end
        
        if (agents(xx,3) == 0)          % If agent is normal agent adjust his temp_opinion
            temp_opinion_CF_hor = 0;        % Initializing for following sum
            temp_opinion_CF_ver = 0;
            temp_opinion_FB_hor = 0;        % Initializing for following sum
            temp_opinion_FB_ver = 0;
            
            listCF_hor = find(Sim_network(xx,:) == 2 | Sim_network(xx,:) == 3);
            listCF_ver = find(Sim_network(:,xx) == 2 | Sim_network(:,xx) == 3);
            listFB_hor = find(Sim_network(xx,:) == 1 | Sim_network(xx,:) == 3);
            listFB_ver = find(Sim_network(:,xx) == 1 | Sim_network(:,xx) == 3);
            
            for zz = 1:length(listCF_hor)
                temp_opinion_CF_hor = temp_opinion_CF_hor + agents(listCF_hor(zz),5);
            end
            
            for zz = 1:length(listCF_ver)
                temp_opinion_CF_ver = temp_opinion_CF_ver + agents(listCF_ver(zz),5);
            end
            
            for zz = 1:length(listFB_hor)
                temp_opinion_FB_hor = temp_opinion_FB_hor + agents(listFB_hor(zz),5);
            end
            
            for zz = 1:length(listFB_ver)
                temp_opinion_FB_ver = temp_opinion_FB_ver + agents(listFB_ver(zz),5);
            end
            % His own opinion 
            temp_opinion(xx) = (1-agents(xx,4)) * agents(xx,5) + agents(xx,4) * ((temp_opinion_CF_hor + temp_opinion_CF_ver)/agents(xx,1) * k_CF + (temp_opinion_FB_hor + temp_opinion_FB_ver)/agents(xx,2) * k_FB) ;
        end
    end
    
    for xx = 1:network_size                                           
        agents(xx,5) = temp_opinion(xx);                    % Updating opinion for every agent
        opinion_development(xx,t) = temp_opinion(xx);       % Updating opinion_development over time
    end
end

figure
for xx = 1:network_size
    plot(opinion_development(xx,:))
    hold on
end
xlabel('Time iterations')
ylabel('Opinion')
title('Opinion development of all agents over time')