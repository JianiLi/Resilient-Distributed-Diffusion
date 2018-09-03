function w = UpdateW(attackers, numAgents, A, w,  phi, attacker_phi, w0_attacker)
    
    if size(size(attacker_phi),2) == 3
        for k = 1:numAgents
            if ~any(k == attackers)
                if ~isempty(attackers)
                    for i = 1:size(attackers,2)
                        if i == 1
                            phi(:,attackers(i)) = attacker_phi(1,:,k);
                        else
                            phi(:,attackers(i)) = attacker_phi(2,:,k);
                        end
                    end
                end
                w(:,k) = phi*A(:,k);
            else
                w(:,k) = w0_attacker;
            end
        end
    else
        for k = 1:numAgents
            if ~any(k == attackers)
                if ~isempty(attackers)
                    for i = 1:size(attackers,2)
                        phi(:,attackers(i)) = attacker_phi(:,k);
                    end
                end
                w(:,k) = phi*A(:,k);
            else
                w(:,k) = w0_attacker;
            end
        end
    end