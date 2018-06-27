function [Adjacency, AgentSet, group1, group2] = getAdjacency(numAgents, sensingRange)

% initial network topology
a = (0.9 - 0.1)/(sqrt(numAgents)-1);
i = reshape(0.1: a :0.9,sqrt(numAgents),1); 
%row = [i;i;i;i;i;i];
row = repmat(i,sqrt(numAgents),1);
column = [];
for i = 0.1:a:0.9
    column = [column; repmat(i,sqrt(numAgents),1)];
end
%column = [repmat(0.2,6,1); repmat(0.2+0.12,6,1);repmat(0.2+0.24,6,1);repmat(0.2+0.36,6,1);repmat(0.2+0.48,6,1);repmat(0.2+0.6,6,1)];
AgentSet = [row + 0.02*(2*rand(numAgents,1)-ones(numAgents,1)), column+ 0.02*(2*rand(numAgents,1)-ones(numAgents,1))];
Adjacency = zeros(numAgents,numAgents);
group1 = randperm(numAgents, numAgents/2);
allmember = 1:numAgents;
group2 = setdiff(allmember, group1);

for k = 1:numAgents
    for l = 1:numAgents
        d = norm(AgentSet(k,:)-AgentSet(l,:));
        if d < sensingRange
            Adjacency(l,k) = 1;
        end
    end
end