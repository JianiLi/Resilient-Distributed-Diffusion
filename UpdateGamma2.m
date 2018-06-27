function gamma2 = UpdateGamma2(normalAgents, attacker, attacker_phi, numAgents, gamma2, Adjacency, w, phi, niu)

    if size(size(attacker_phi),2) == 3
        for k = normalAgents
            for l = 1:numAgents
                if Adjacency(l,k) == 1
                    if any(attacker == l)
                        if l == attacker(1)
                            gamma2(l,k) = (1-niu)*gamma2(l,k) + niu*(norm(w(:,k)-reshape(attacker_phi(1,:,k),2,1))^2);
                        else
                            gamma2(l,k) = (1-niu)*gamma2(l,k) + niu*(norm(w(:,k)-reshape(attacker_phi(2,:,k),2,1))^2);
                        end

                    else
                        gamma2(l,k) = (1-niu)*gamma2(l,k) + niu*(norm(w(:,k)-phi(:,l))^2);
                    end
                end
            end
        end
    else
        for k = normalAgents
            for l = 1:numAgents
                if Adjacency(l,k) == 1
                    if any(attacker == l)
                        gamma2(l,k) = (1-niu)*gamma2(l,k) + niu*(norm(w(:,k)-attacker_phi(:,k))^2);
                    else
                        gamma2(l,k) = (1-niu)*gamma2(l,k) + niu*(norm(w(:,k)-phi(:,l))^2);
                    end
                end
            end
        end
    end
    