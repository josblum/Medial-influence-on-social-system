function [G]=erdosRenyi(network_size,reconnection_prob,number_of_friends_per_agent)
%Funciton [G]=edosRenyi(network_size,connection_prob,Kreg) generates a random graph based on
%the Erdos and Renyi algoritm where all possible pairs of 'nv' nodes are
%connected with probability 'p'. 
%
% Inputs:
%   network_size - number of agents 
%   reconnection_prob  - Probability that friendship to "cell neighbour" is
%       removed and friendship with other node is established instead
%       -->If set to 1 network is "totally random"
%   number_of_friends - initial number of connections (between 1 and
%       network_size/2)
%
% Output:
%   G is a structure inplemented as data structure in this as well as other
%   graph theory algorithms.
%   G.Adj   - is the adjacency matrix (1 for connected nodes, 0 otherwise).
%   G.x and G.y -   are row vectors of size nv wiht the (x,y) coordinates of
%                   each node of G.
%   G.nv    - number of vertices in G
%   G.ne    - number of edges in G
%
%Created by Pablo Blinder. blinderp@bgu.ac.il
%
%Adapted by Josef Blum.
%Last update 05.12.17

%build regular lattice 
network=zeros(network_size,network_size);
number_of_friends_per_agent=abs(number_of_friends_per_agent);

if(number_of_friends_per_agent<1 || number_of_friends_per_agent>network_size)
     print('Error forming network: variable number_of_friends is unvalid')
end
if(mod(number_of_friends_per_agent,2)~=0)
    print('Error forming network: variable number_of_friends has to be an even integer.')
end

for k=1:(number_of_friends_per_agent/2)
    network=(network+diag(ones(1,length(diag(network,k))),k)+diag(ones(1,length(diag(network,network_size-k))),network_size-k));
end
%Find number of non-zero-elements in network
Initial_Number_of_Friendships=nnz(network);
%find connected pairs (non-zero-entries)
[nz_network_line,nz_network_row]=find(network);
%Creates vector of booleans
Removal=(rand(length(nz_network_line),1)<=reconnection_prob);
%Friendships with "true" boolean are removed
network(nz_network_line(Removal),nz_network_row(Removal))=0;

%Find removed friendship pairs in sorted order (vectors)
%vRemove=unique([nz_network_line(Removal),nz_network_row(Removal)]);

%Number of removed frienships
nRemove=Initial_Number_of_Friendships-nnz(network);

%sum(Removal);

%cycle trough disconnected pairs
disconPairs=[nz_network_line(Removal),nz_network_row(Removal)];
for n=1:nRemove
    %choose one of the vertices from the disconnected pair
    i=ceil(rand*size(disconPairs,1));
    j=logical(1+rand>0.5);
    vDisToRec=disconPairs(i,j);
    %find non adjacent vertices and reconnect
    adj=[find(network(:,vDisToRec)) ; find(network(vDisToRec,:))'];
    nonAdj=setdiff(1:network_size,adj);
    vToRec=vDisToRec;
    %make sure connection is not established to node itself
while vToRec==vDisToRec
vToRec=nonAdj(ceil(rand*length(nonAdj)));
end 
    %write new friendships to network
    S=sort([vDisToRec vToRec]);
    network(S(1),S(2))=1;
end
Friendships_at_end = nnz(network);
if(Friendships_at_end~=Initial_Number_of_Friendships)
    print('Error forming network: Bug in erdosRenyi')
end

%Make sure everybody has at least one friend
for i=1:network_size
    if(sum(network(:,i))+sum(network(i,:))==0)
        r=ceil(rand*network_size);
        while(r==i)
         r=ceil(rand*network_size)
        end
        if(r<i)
            network(r,i)=1;
        else
            network(i,r)=1;
        end
    end
end

[x,y]=getNodeCoordinates(network_size);
G=struct('network',network,'x',x','y',y','Nbr_agents',network_size,'Nbr_friendships',nnz(network));