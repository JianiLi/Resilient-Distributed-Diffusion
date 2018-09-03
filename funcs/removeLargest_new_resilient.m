function [newAdjacency,ratio,J] = removeLargest_new_resilient(n, F, newAdjacency, numAgents, Adjacency, attackers, D, U, phi, storedNum, attacker_phi, ...
    Expectation_noco, Expectation_coop, gamma2 )
        
        J = zeros(numAgents, numAgents);
        ratio = zeros(numAgents, numAgents);
        %remove largest 
    
        for k = 1: numAgents
            if ~any(k == attackers)
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
                    gamma2(k,k) = 0;
                
                if F == 1
                    if f - 1 >= F
                        for i = 1:F
                            minRatio = Inf;
                            for l = 1:numAgents
                                if Adjacency(l,k) == 1 && l ~= k
                                    if l == 1
                                        ratioWithoutl = [ratio(2:numAgents,k)];
                                        gamma2Withoutl = [gamma2(2:numAgents, k)];
                                        gamma2Withoutl(gamma2Withoutl  == 0) = [];
                                    elseif l == numAgents
                                        ratioWithoutl = [ratio(1:numAgents-1, k)];
                                        gamma2Withoutl = [gamma2(1: numAgents-1, k)];
                                        gamma2Withoutl(gamma2Withoutl  == 0) = [];
                                    else
                                        ratioWithoutl = [ratio(1:l-1,k); ratio(l+1:numAgents,k)];
                                        gamma2Withoutl = [gamma2(1:l-1,k); gamma2(l+1:numAgents,k)];
                                        gamma2Withoutl(gamma2Withoutl  == 0) = [];
                                    end
                                    new_ratio = sum(ratioWithoutl) ./ sum(1./gamma2Withoutl).^2;
                                    if new_ratio < minRatio
                                        minRatio = new_ratio;
                                        removeR = l;
                                    end        
                                end
                            end
                            newAdjacency(removeR,k) = 0; 
                            Adjacency(removeR,k) = 0;
                            ratio(removeR,k) = 0;
                            gamma2(removeR,k) = 0;
                        end
                    else
                        newAdjacency(Adjacency(:,k)~=0,k) = 0; 
                        Adjacency(Adjacency(:,k)~=0,k) = 0;
                    end
                    
                else
                    if f - 1 > F
                        minRatio = Inf;
                        neighborsOfk = find(Adjacency(:,k)~=0);
                        neighborsOfk(neighborsOfk==k) = [];
                        removeAllList = nchoosek(neighborsOfk, F);
                        for i = 1:size(removeAllList,1)
                            removeList = removeAllList(i,:);
                            ratioWithoutRemoveList = ratio(:,k);
                            ratioWithoutRemoveList(removeList)=0;
                            gamma2WithoutRemoveList = gamma2(:,k);
                            gamma2WithoutRemoveList(removeList) = 0;
                            gamma2WithoutRemoveList(gamma2WithoutRemoveList  == 0) = [];
                            %new_ratio = sum(ratioWithoutRemoveList) ./ sum(1./gamma2WithoutRemoveList).^2;
                            new_ratio = sum(ratioWithoutRemoveList);
                            if new_ratio <= minRatio
                                minRatio = new_ratio;
                                actualRemoveList = removeList;
                            end  
                        end
                        newAdjacency(actualRemoveList,k) = 0;                         
                    else
                        newAdjacency(Adjacency(:,k)~=0,k) = 0;       
                    end
                end
                
            newAdjacency(k,k) = 1;     
            end
        end
            
end


