function [newAdjacency, J] = removeLargestCostOfRemainingNeighbors(alpha, F, gamma2, Adjacency, newAdjacency, numAgents, attackers, D, U, phi, storedNum, attacker_phi)

        J = zeros(numAgents, numAgents);
        %remove largest 
        for k = 1: numAgents
            f = 0;
                for l = 1:numAgents
                    if Adjacency(l,k) == 1
                        if gamma2(l,k) < alpha
                            
                            f = f + 1;
                        %J(l,k) = (d(n,k) - u(:,k)'*phi(:,l))^2;
                            if ~any(l == attackers)
                                J(l,k) = sum((D(:,k) - U(:,:,k) * phi(:,l)).^2)/storedNum;
                            else
                                if l == attackers(1)
                                    J(l,k) = sum((D(:,k) - U(:,:,k) * attacker_phi(:,k)).^2)/storedNum;
                                else
                                    J(l,k) = sum((D(:,k) - U(:,:,k) * attacker_phi(:,k)).^2)/storedNum;
                                end
                            end
                        end
                    end
                end

                J(k,k) = 0;
                if f - 1 > F
                    for i = 1:F
                        maxR = find(J(:,k)==max(J(:,k))) ;
                        newAdjacency(maxR,k) = 0; 
                        J(maxR,k) = 0;
                    end
                else
                    for i = 1:f - 1
                        maxR = find(J(:,k)==max(J(:,k))) ;
                        newAdjacency(maxR,k) = 0; 
                        J(maxR,k) = 0;
                    end
                end
                    
%                 maxR = find(J(:,k)==max(J(:,k))) ;
%                 newAdjacency(maxR,k) = 0; 
%                 J(maxR,k) = 0;
%                 secondMax = find(J(:,k)==max(J(:,k)));
%                 newAdjacency(secondMax,k) = 0; 
%                 ratio(secondMax,k) = 0;
%                 thirdMax = find(ratio(:,k)==max(ratio(:,k)));
%                 newAdjacency(thirdMax,k) = 0; 
%                 ratio(thirdMax,k) = 0;
%                 fourthMax = find(ratio(:,k)==max(ratio(:,k)));
%                 newAdjacency(fourthMax,k) = 0; 
                newAdjacency(k,k) = 1; 
        end