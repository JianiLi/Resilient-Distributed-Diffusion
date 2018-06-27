function attackers_new = UpdateAttackers(numAgents, w, w0_attacker, attackers_new)  

    for k = 1:numAgents
        if norm(w(:,k) - w0_attacker) <= 0.01
            attackers_new = unique([attackers_new,k]);
        end
    end
