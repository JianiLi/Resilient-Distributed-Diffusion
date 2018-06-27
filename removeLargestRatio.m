function [newAdjacency,ratio,J] = removeLargestRatio(n, F, newAdjacency, numAgents, Adjacency, attackers, D, U, phi, storedNum, attacker_phi, ...
    Expectation_noco, Expectation_coop, gamma2 )
        
        J = zeros(numAgents, numAgents);
        ratio = zeros(numAgents, numAgents);
        %remove largest 
    
        for k = 1: numAgents
            f = 0;
            %if Expectation_noco(n-1,k) < Expectation_coop(n-1,k)
                for l = 1:numAgents
                    if Adjacency(l,k) == 1
                        %J(l,k) = (d(n,k) - u(:,k)'*phi(:,l))^2;
                        f = f + 1;
                        if ~any(l == attackers)
                            J(l,k) = sum((D(:,k) - U(:,:,k) * phi(:,l)).^2)/storedNum;
                            ratio(l,k) = J(l,k)/gamma2(l,k)^2;
                        else
                            if l == attackers(1)
                                J(l,k) = sum((D(:,k) - U(:,:,k) * attacker_phi(:,k)).^2)/storedNum;
                            else
                                J(l,k) = sum((D(:,k) - U(:,:,k) * attacker_phi(:,k)).^2)/storedNum;
                            end
                            ratio(l,k) = J(l,k)/gamma2(l,k)^2;

                        end
                    end
                end
                
                ratio(k,k) = 0;
                if f - 1 >= F
                    for i = 1:F
                        maxR = find(ratio(:,k)==max(ratio(:,k))) ;
                        newAdjacency(maxR,k) = 0; 
                        ratio(maxR,k) = 0;
                    end
                else
                    for i = 1:f - 1
                        maxR = find(ratio(:,k)==max(ratio(:,k))) ;
                        newAdjacency(maxR,k) = 0; 
                        ratio(maxR,k) = 0;
                    end
                end
                

                
%                 maxR = find(ratio(:,k)==max(ratio(:,k))) ;
%                 newAdjacency(maxR,k) = 0; 
%                 ratio(maxR,k) = 0;
%                 secondMax = find(ratio(:,k)==max(ratio(:,k)));
%                 newAdjacency(secondMax,k) = 0; 
%                 ratio(secondMax,k) = 0;
                
%                 if size(maxR,1) == 1
%                     newAdjacency(maxR,k) = 0; 
%                     ratio(maxR,k) = 0;
%                     secondMax = find(ratio(:,k)==max(ratio(:,k)));
%                     if size(secondMax,1) == 1
%                         newAdjacency(secondMax,k) = 0; 
%                     end
%                 elseif size(maxR,1) == 2
%                     newAdjacency(maxR,k) = 0; 
%                     ratio(maxR,k) = 0;
%                 end
                
                
%                 ratio(secondMax,k) = 0;
%                 thirdMax = find(ratio(:,k)==max(ratio(:,k)));
%                 newAdjacency(thirdMax,k) = 0; 
%                 ratio(thirdMax,k) = 0;
%                 fourthMax = find(ratio(:,k)==max(ratio(:,k)));
%                 newAdjacency(fourthMax,k) = 0; 
                 newAdjacency(k,k) = 1; 
            %end
        end

end

