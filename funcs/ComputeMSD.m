function [MSD_coop, MSD_ncop] = ComputeMSD(group1, numAgents, attackers, MSD_coop, MSD_ncop, n, w0, w, w_noco, normalAgents)

    if size(size(w0),2) == 2
        for k = 1:numAgents
            if ~any(k==attackers)
                if any(k == group1)
                    MSD_coop(n-1) = MSD_coop(n-1) + norm(w0(:,1) - w(:,k))^2;
                    MSD_ncop(n-1) = MSD_ncop(n-1) + norm(w0(:,1) - w_noco(:,k))^2;
                else
                    MSD_coop(n-1) = MSD_coop(n-1) + norm(w0(:,2) - w(:,k))^2;
                    MSD_ncop(n-1) = MSD_ncop(n-1) + norm(w0(:,2) - w_noco(:,k))^2;
                end
            end
        end
    else
        for k = 1:numAgents
                if ~ any(k == attackers)
                    MSD_coop(n-1) = MSD_coop(n-1) + norm(reshape(w0(n-1,:,k),2,1) - w(:,k))^2;
                    MSD_ncop(n-1) = MSD_ncop(n-1) + norm(reshape(w0(n-1,:,k),2,1) - w_noco(:,k))^2;
                end
        end
    end
    MSD_coop(n-1) = MSD_coop(n-1)/size(normalAgents,2);
    MSD_ncop(n-1) = MSD_ncop(n-1)/size(normalAgents,2);