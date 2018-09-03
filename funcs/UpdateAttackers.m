function attackers_new = UpdateAttackers(numAgents, w, w0_attacker, attackers_new)  

    for k = 1:numAgents
        if ~ismember(k,attackers_new) && norm(w(:,k) - w0_attacker) <= 0.01
            attackers_new = [attackers_new,k];
        elseif ismember(k,attackers_new) && norm(w(:,k) - w0_attacker) > 0.01
            attackers_new(find(attackers_new==k)) = [];
        end
        
    end
