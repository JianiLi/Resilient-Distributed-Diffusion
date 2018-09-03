function A = UpdateWeight(normalAgents, numAgents, gamma2, Adjacency)

    A = zeros(numAgents,numAgents);
    for k = normalAgents        
        denominator = 0;
        for l = 1:numAgents 
            if Adjacency(l,k) == 1
                denominator = denominator + 1/gamma2(l,k);
            end
        end

        for l = 1:numAgents
            if Adjacency(l,k) == 1
                A(l,k) = (1/gamma2(l,k))/denominator;
            else
                A(l,k) = 0;
            end
        end
    end