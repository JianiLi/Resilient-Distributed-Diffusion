clear; close all; clc
rng('default');
% figure(4)
% plot(b_average(:),'linewidth',2);
% set(gca,'FontSize',18);
% set(gcf,'color','white');
% xlabel('Iteration $i$','Interpreter','LaTex','FontSize',24);
% ylabel('Average $b_i$','Interpreter','LaTex','FontSize',24);

eta_matrix = [];
falseAlarm_matrix = [];
wAverageMatrix = [];

%for eta = 0:0.005:0.2
    numPoints = 5004;
    numTaps = 2;
    Mu = 0.01;
    numAgents = 100;
    niu = 0.01; %forgetting factor
    gamma2 = zeros(numAgents,numAgents);

    % input is guassian
    mu_x = 0;
    sigma_x2 = 0.8 + 0.4*rand(numAgents,1);

    x = zeros(numPoints,numAgents);
    for k = 1:numAgents
        x(:,1,k) = mvnrnd(mu_x, sigma_x2(k), numPoints);
        x(:,2,k) = mvnrnd(mu_x, sigma_x2(k), numPoints);
    end

    % initial network topology
    sensingRange = 0.16;
    [Adjacency, AgentSet, group1, group2] = getAdjacency(numAgents, sensingRange);

    for k = 1:numAgents
        for l = 1:numAgents
            d = norm(AgentSet(k,:)-AgentSet(l,:));
            if d < sensingRange
                Adjacency(l,k) = 1;
            end
        end
    end
    AdjacencyMatrix = Adjacency;

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
    attacker = [24 37 63 88];
    attacker_new = [];
    ra = 0.01;
    omega_a = 1/2000;
    attacker_target = zeros(numPoints,2);
    attacker_delta_theta = zeros(numPoints,2);
    for i = 1:numPoints
        attacker_target(i,:) = [0.1*cos(2*pi*omega_a*i)+0.5,0.1*sin(2*pi*omega_a*i)+0.5];
        attacker_delta_theta(i,:) = [-0.2*pi*omega_a*sin(2*pi*omega_a*i),0.2*pi*omega_a*cos(2*pi*omega_a*i)];
    end


    sigma_v2 = 0.1+0.1*rand(numAgents,1);
    v = zeros(numPoints,numAgents);
    for k = 1:numAgents
        v(:,k) = mvnrnd(0, sigma_v2(k), numPoints);
    end


    d = zeros(numPoints,numAgents);
    for k = 1:numAgents
        x_k(:,1) = x(:,1,k);
        x_k(:,2) = x(:,2,k);
        for i = 1:numPoints
            d(i,k) = x_k(i,:)*w0(i,:,k)';
        end
    end
    d = d+v;


    % initialize variables
    h = [];
    u = [];
    e = []; % error, final result to be computed
    w = rand(numTaps,numAgents);
    w_0 = w;
    w_noco = w;
    w_noco_old = w_noco;
    z = zeros(numPoints, numAgents);
    S = zeros(numPoints, numAgents);
    eta = 0.1008;
    alarmSet = [ ];
    MSD_coop = zeros(numPoints-1,1); MSD_ncop = zeros(numPoints-1,1);

    plotNetworkTopology(1, AgentSet, w0, Adjacency, group1, group2, attacker, attacker_new, numTaps);


    % LMS Adaptation
    for n = numTaps : numPoints
        if mod(n,1000) == 0 || n == numPoints
            %update adjacency matrix
            AdjacencyMatrix = (A>=0.01) & Adjacency;
            plotNetworkTopology(n, AgentSet, w0, AdjacencyMatrix, group1, group2, attacker, attacker_new, numTaps)
        end
       

        if n == numTaps
            for k = attacker
                fprintf('attacker %d has %d neighbors to attack\n', k, size(find(Adjacency(:,k)==1),1)-1)
            end
        end

        for k = 1:numAgents
            if ~any(attacker == k)
                attacker_phi(:,k) = w(:,k) + ra*(attacker_target(n,:)' + attacker_delta_theta(n,:)'/ra -w(:,k));
            end
        end
        
        wAverageMatrix = getwAverage(numAgents, wAverageMatrix, attacker, Adjacency, w);
      
        
        for k = 1 : numAgents
            if ~any(attacker == k)
                u(:,k) = [x(n,1,k),x(n,2,k)];     % select part of training input
                h(n,k) = u(:,k)'*w(:,k);    % hypothsis function
                e(n,k) = d(n,k)-h(n,k);    % error
                w_old(:,k) = w(:,k);
                phi(:,k) = w(:,k) + Mu*( u(:,k)*e(n,k) );
            end
        end

        % update gamma2
        for k = 1:numAgents
            if ~any(k == attacker)
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
        
        if n >numTaps
                [newAdjacency,ratio,J] = removeLargestRatio(n, F, newAdjacency, numAgents, Adjacency, attackers, D, U, phi, storedNum, attacker_phi, ...
        Expectation_noco, Expectation_coop, gamma2 );
            %[newAdjacency, J] = removeLargestCostOfRemainingNeighbors(0.01, 3, gamma2, AdjacencyMatrix, newAdjacency, numAgents, attackers, D, U, phi, storedNum, attacker_phi);
        end

        A = zeros(numAgents,numAgents);
        for k = 1:numAgents
            if ~any(k == attacker)
                denominator = 0;
                for l = 1:numAgents
                    if Adjacency(l,k) == 1
                        denominator = denominator + 1/gamma2(l,k);
                    end
                end

                for l = 1:numAgents
                    if Adjacency(l,k) == 1
                        A(l,k) = (1/gamma2(l,k))/denominator;
                    else
                        A(l,k) = 0;
                    end
                end
            end
        end

        %update adjacency matrix

        if n <= numPoints/2 & mod(n,500)==0
            AdjacencyMatrix = (A>=0.0001) & Adjacency;
        elseif mod(n,500)==0
            AdjacencyMatrix = (A>=0.1) & Adjacency;
        end

        %update w
        for k = 1:numAgents
            if ~any(k == attacker)
                if ~isempty(attacker)
                    for i = 1:size(attacker,2)
                        phi(:,attacker(i)) = attacker_phi(:,k);
                    end
                end
                w(:,k) = phi*A(:,k);
            else
                w(:,k) = attacker_target(n,:)';
            end
        end

        if n > 1000
            for k = 1:numAgents
                if norm(w(:,k) - attacker_target(n,:)') <= 0.1
                    attacker_new = unique([attacker_new,k]);
                end
            end
        end


        % non-cooperative w
        for k = 1 : numAgents
            if ~any(attacker == k)
                h_noco(n,k) = u(:,k)' * w_noco(:,k);    % hypothsis function
                e_noco(n,k) = d(n,k)-h_noco(n,k);    % error
                w_noco_old(:,k) = w_noco(:,k);
                w_noco(:,k) = w_noco(:,k) + Mu*( u(:,k)*e_noco(n,k) );
                error(n-1,k) = norm(w_noco(:,k)-w(:,k));
                % X(n-1,k) = norm(w_noco(:,k) - w_noco_old(:,k));
            end
        end

%         %detector
%         for k = setdiff(allmember,alarmSet)
%             z(n-1,k) = norm(w(:,k) - w_noco(:,k)) - b_average(n-1);
%             if S(n-1,k) + z(n-1,k) > 0
%                 S(n-1,k) = S(n-1,k) +z(n-1,k);
%             else
%                 S(n-1,k) = 0;
%             end
%             if S(n-1,k) > eta
%                 alarmSet = [alarmSet k];
%             end
%         end

            % MSD 

        for k = 1:numAgents
                if ~ any(k == attacker)
                    MSD_coop(n-1) = MSD_coop(n-1) + norm(reshape(w0(n-1,:,k),2,1) - w(:,k))^2;
                    MSD_ncop(n-1) = MSD_ncop(n-1) + norm(reshape(w0(n-1,:,k),2,1) - w_noco(:,k))^2;
                end
        end
        MSD_coop(n-1) = MSD_coop(n-1)/(numAgents-size(attacker,2));
        MSD_ncop(n-1) = MSD_ncop(n-1)/(numAgents-size(attacker,2));



    end

    falseAlarmNum = 0;
    for k = 1:numAgents
        if (~any(attacker_new == k)) & (any(alarmSet==k))
            falseAlarmNum = falseAlarmNum +1;
        end
    end

    falseAlarmRate = falseAlarmNum/numAgents;
    fprintf('eta = %d, false alarm rate = %d\n', eta, falseAlarmRate);
    eta_matrix = [eta_matrix; eta];
    falseAlarm_matrix = [falseAlarm_matrix; falseAlarmRate];

%end

figure(3);
set (gcf,'Position',[0,0,450,450], 'color','w');
drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );
scatter_attacker_target = [];
scatter_wAverageMatrix = [];
for i = 1: numPoints -10 
    if mod(i,10) == 0 & i < 1900
        scatter_attacker_target = [scatter_attacker_target; attacker_target(i,:) ];
        scatter_wAverageMatrix = [scatter_wAverageMatrix wAverageMatrix(:,i)];
    end
end
plot(scatter_attacker_target(:,1),scatter_attacker_target(:,2), 'or','MarkerFaceColor','r','MarkerSize',5);
hold on;
plot(scatter_wAverageMatrix(1,:),scatter_wAverageMatrix(2,:), '+g','MarkerFaceColor','g','MarkerSize',5);
axis('equal');
axis([0.39 0.64 0.39 0.65]);
h=legend('Attacker selected state $w_{k,i}^a$','Average state $\bar{w}_{k,i}$ of compromised nodes'' neighbors');
set(h,'interpreter','Latex','FontSize',11)

%x1 = [0.615 0.615];
%y1 = [0.5 0.52 ];
x1 = [0.6 0.6];
y1 = [0.475 0.495 ];
drawArrow(x1,y1,'MaxHeadSize',1,'linewidth',2,'color','r'); hold on

x2 =  [0.605 0.60 ];
y2 =  [0.59 0.57];
drawArrow(x2,y2,'MaxHeadSize',1,'linewidth',2,'color','g'); hold on
xlim([0.39 0.64]);
ylim([0.39, 0.65]);

h1 = text(0.61, 0.6,'$\bar{w}_{k,0}$');
h2 = text(0.605,0.47,'$w_{k,0}^a$');
set(h1,'interpreter','Latex','FontSize',24);
set(h2,'interpreter','Latex','FontSize',24);
set(gca,'FontSize',15);
set(gcf,'color','white')
xlabel('X','fontsize',20);ylabel('Y','fontsize',20);
% figure(4)
% D = [];
% for k = 1:numPoints-1
%     d = norm(w4Matrix(1,k),w4Matrix(2,k)) - norm(attacker_target(k,1),attacker_target(k,i));
%     D = [D d];
% end
% plot(D);

% figure(5)
% wdist = [];
% for k = 1: numPoints-2
%     wdist = [wdist norm(w4Matrix(:,k+1)-w4Matrix(:,k))];
% end
% plot(wdist);
% hold on;
% thetadist = [];
% for k = 1: numPoints-2
%     thetadist = [thetadist norm(target(k+1,:,4)-target(k,:,4))];
% end
% plot(thetadist);
% 
% figure(6)
% wdist = [];
% for k = 1: numPoints-2
%     wdist = [wdist w4Matrix(:,k+1)-w4Matrix(:,k)];
% end
% ww0dist = [];
% for k = 1:numPoints-2
%     w0 = attacker_target(k,:) - reshape(wdist(:,k),1,2)/0.002;
%     ww0dist = [ww0dist norm(w0 - reshape(w4Matrix(:,k),1,2))];
% end
% plot(ww0dist)

figure(2)
set (gcf,'Position',[0,0,450,450], 'color','w');
set(gcf,'color','white');
axis([0 5000 -60 0])
plot(mag2db(MSD_ncop(1:5000)), 'linewidth',1);
hold on;
plot(mag2db(MSD_coop(1:5000)), 'linewidth',1);
set(gca,'FontSize',15);
gca = legend('Noncooperative LMS','DLMSAW', 'interpreter','latex', 'fontsize',20);
set(gca,'FontSize',20);
xlabel('Iteration $i$', 'interpreter','latex','fontsize',20);ylabel('MSD(dB)', 'interpreter','latex','fontsize',20);
box on;
