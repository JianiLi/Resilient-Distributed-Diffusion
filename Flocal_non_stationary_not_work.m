clear; close all; clc
numF = 1;
numPoints = 1000;
MSD_coop = zeros(numPoints-1,numF+1); 
MSD_ncop = zeros(numPoints-1,numF+1);
for F = 0:1:numF
    rng('default');

    %% PARAMETERS
    numAgents = 100;
    numTaps = 2;		% channel number
    Mu = 0.01;          % step size
    niu = 0.01;         % forgetting factor
    w = rand(numTaps,numAgents);
    w_noco = w;
    phi = zeros(numTaps,numAgents);
    gamma2 = zeros(numAgents,numAgents);
    A = zeros(numAgents,numAgents);
    sensingRange = 0.16;

    %% DETECTION PARAMETERS
    storedNum = 100;
    D = zeros(storedNum, numAgents);
    U = zeros(storedNum,numTaps, numAgents);
    Expectation_noco = zeros(numPoints,numAgents);
    Expectation_coop = zeros(numPoints,numAgents);
    J = zeros(numAgents, numAgents);
    alpha = 0.05;
    beta = 0.1;


    %% ATTACKER SETTINGS
    %attackers = [24 37 63 88];
    attackers = [24 37 63 88 ];
    attackers_new = [];
    attacker_phi = zeros(2,numTaps, numAgents);
    Agents = 1:numAgents;
    normalAgents = setdiff(Agents,attackers);
    wAverageMatrix = [];

    %% INPUTS (GAUSSIAN)
    mu_x = 0;
    sigma_x2 = 0.8 + 0.4*rand(numAgents,1);
    x = zeros(numPoints,numAgents);
    for k = 1:numAgents
        x(:,1,k) = mvnrnd(mu_x, sigma_x2(k), numPoints);
        x(:,2,k) = mvnrnd(mu_x, sigma_x2(k), numPoints);
    end

    sigma_v2 = 0.15+0.05*rand(numAgents,1);
    v = zeros(numPoints,numAgents);
    for k = 1:numAgents
        v(:,k) = mvnrnd(0, sigma_v2(k), numPoints);
    end

    % NETWORK TOPOLOGY
    [Adjacency, AgentSet, group1, group2] = getAdjacency(numAgents, sensingRange);
    AdjacencyMatrix = Adjacency;
    %plotNetworkTopology(1, AgentSet, w0, AdjacencyMatrix, group1, group2, attackers, attackers_new, numTaps);
    
    omega = 1/2000;
    target = zeros(numPoints,2,numAgents);
    for i = 1:numPoints
        for k = group1
            target(i,:,k) = [0.1*cos(2*pi*omega*i)+0.1,0.1*sin(2*pi*omega*i)+0.1];
        end
        for k = group2
            target(i,:,k) = [0.1*cos(2*pi*omega*i)+0.9,0.1*sin(2*pi*omega*i)+0.9];
        end
    end
    w0 = target;
    
    ra = 0.002;
    omega_a = 1/2000;
    attacker_target = zeros(numPoints,2);
    attacker_delta_theta = zeros(numPoints,2);
    for i = 1:numPoints
        attacker_target(i,:) = [0.1*cos(2*pi*omega_a*i)+0.5,0.1*sin(2*pi*omega_a*i)+0.5];
        attacker_delta_theta(i,:) = [-0.2*pi*omega_a*sin(2*pi*omega_a*i),0.2*pi*omega_a*cos(2*pi*omega_a*i)];
    end

    w0_attacker = attacker_target;

    d = zeros(numPoints,numAgents);
    for k = 1:numAgents
        x_k(:,1) = x(:,1,k);
        x_k(:,2) = x(:,2,k);
        for i = 1:numPoints
            d(i,k) = x_k(i,:)*w0(i,:,k)';
        end
    end
    d = d+v;



    %% DIFFUSION LMS ALGORITHM
    for n = numTaps : numPoints
        newAdjacency = Adjacency;

        if mod(n,1000) == 0 || n == numPoints
            %update adjacency matrix
            AdjacencyMatrix = (A>=0.01) & Adjacency;
            %plotNetworkTopology(n, AgentSet, w0, AdjacencyMatrix, group1, group2, attackers, attackers_new, numTaps)
        end

        if size(attackers,2) ~= 0
            for k = normalAgents
                attacker_phi(:,k) = w(:,k) + ra*(attacker_target(n,:)' + attacker_delta_theta(n,:)'/ra -w(:,k));
            end
        end


        for k = normalAgents
            u(:,k) = x(n:-1:n-numTaps+1,k);     % select part of training input
            % Cooperative state
            h(n,k) = u(:,k)'*w(:,k);            % hypothsis function
            e(n,k) = d(n,k)-h(n,k);             % error
            phi(:,k) = w(:,k) + Mu*( u(:,k)*e(n,k) );
            % non-cooperative state
            h_noco(n,k) = u(:,k)' * w_noco(:,k);    % hypothsis function
            e_noco(n,k) = d(n,k)-h_noco(n,k);    % error
            w_noco_old(:,k) = w_noco(:,k);
            w_noco(:,k) = w_noco(:,k) + Mu*( u(:,k)*e_noco(n,k) );
            phi_noco(:,k) = w_noco(:,k);
        end

        gamma2 = UpdateGamma2(normalAgents, attackers, attacker_phi, numAgents, gamma2, newAdjacency, w, phi, niu);
        %gamma2 = UpdateGamma2(normalAgents, attackers, attacker_phi, numAgents, gamma2, Adjacency, w_noco_old, w_noco, niu);
        %gamma2 = UpdateGamma2(normalAgents, attackers, attacker_phi, numAgents, gamma2, Adjacency, w_noco_old, phi_noco, niu);

        if n >numTaps
                [newAdjacency,ratio,J] = removeLargestRatio(n, F, newAdjacency, numAgents, Adjacency, attackers, D, U, phi, storedNum, attacker_phi, ...
        Expectation_noco, Expectation_coop, gamma2 );
            %[newAdjacency, J] = removeLargestCostOfRemainingNeighbors(0.01, 3, gamma2, AdjacencyMatrix, newAdjacency, numAgents, attackers, D, U, phi, storedNum, attacker_phi);
        end

        A = UpdateWeight(normalAgents, numAgents, gamma2, newAdjacency);
        w = UpdateW(attackers, numAgents, A, w,  phi, attacker_phi, w0_attacker(n,:)');

        % update attackers
        if n > 1000
            attackers_new = UpdateAttackers(numAgents, w, w0_attacker(n,:)', attackers_new);
        end

        [MSD_coop(:,F+1), MSD_ncop(:,F+1)] = ComputeMSD(group1, numAgents, attackers, reshape(MSD_coop(:,F+1),numPoints-1,1), reshape(MSD_ncop(:,F+1), numPoints-1, 1), n, w0, w, w_noco, normalAgents);

        [D, U] = updateStoredStreamingData(d(n,:),u, D, U);
        %Expectation_noco = Expectation(n, numAgents, Expectation_noco, D, U, w_noco, storedNum);
        %Expectation_coop = Expectation(n, numAgents, Expectation_coop, D, U, w, storedNum);
        wAverageMatrix = getwAverage(numAgents, wAverageMatrix, attackers, Adjacency, w);
    end

    %% PLOT
    % MSD
    %figure(2)
    if F == 0
        set (gcf,'Position',[0,0,450,450], 'color','w');
        set(gca,'XTick', [0:1000:5000])
        plot(mag2db(MSD_ncop(:,F+1)), 'linewidth',1);
        hold on;
        mylgd{1} = ['Noncooperative LMS'];
    end
    plot(mag2db(MSD_coop(:,F+1)), 'linewidth',1);
    hold on;
    set(gca,'FontSize',15);
    mylgd{F+2} = ['F-local, F = ', num2str(F)];
    if F == numF
        gca = legend(mylgd,'NorthEastOutside');
        set(gca,'FontSize',12);
        xlabel('Iteration $i$', 'interpreter','latex','fontsize',20);ylabel('MSD(dB)', 'interpreter','latex','fontsize',20);
        box on;
    end
end

% figure(3);
% set (gcf,'Position',[0,0,450,450], 'color','w');
% for i = 1:numPoints-1
%     delta(i) = norm(w0_attacker - wAverageMatrix(:,i));
% end
% plot(delta,'linewidth',2);
% 
% set(gca,'FontSize',15);
% set(gcf,'color','white');
% xlabel('Iteration $i$','Interpreter','LaTex','FontSize',20);
% ylabel('$\|w_k^a - \bar{w}_{k,i}\|$','Interpreter','LaTex','FontSize',20);

% Expectation
% figure(3)
% plot(Expectation_noco(:,14), 'linewidth',1);
% hold on;
% plot(Expectation_coop(:,14), 'linewidth',1);
% set(gcf,'color','white');
% set(gca,'FontSize',15);
% xlabel('Iteration i');ylabel('Expectation');
% box on;
% legend('Noncooperative LMS','Distributed clustering');
% 
% % Expectation difference
% figure(4)
% E = Expectation_noco - Expectation_coop;
% plot(E, 'linewidth',1);
% set(gcf,'color','white');
% set(gca,'FontSize',15);
% xlabel('Iteration i');ylabel('Expectation difference');
% box on;
