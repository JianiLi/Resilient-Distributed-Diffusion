function wAverageMatrix = getwAverage(numAgents, wAverageMatrix, attackers, Adjacency, w)
    w_average = 0;
    num = 0;
    for k = 1 : numAgents
        flag = 0;
        for l = attackers
            if Adjacency(l,k) == 1
                flag = 1;
                num = num + 1;
            end
        end
        if flag == 1
            w_average = w_average + w(:,k);
        end
    end
    w_average = w_average / num;

    wAverageMatrix = [wAverageMatrix w_average];