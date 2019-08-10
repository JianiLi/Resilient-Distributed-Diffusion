% Test black-box attack: attacker does not know the streaming data
% information, by communication message of the target node and its
% neighbors, the attacker calculates the state

clear; close all; clc
rng('default');

%% PARAMETERS
numAgents = 16;
numTaps = 2;		% channel number
%numPoints = 5050;
numPoints = 5000;
Mu = 0.01;          % step size
niu = 0.02;         % forgetting factor
w = rand(numTaps,numAgents);
w_noco = w;
w0 = [0.1 0.9; 0.1 0.9];
phi = zeros(numTaps,numAgents);
gamma2 = zeros(numAgents,numAgents);
A = zeros(numAgents,numAgents);
sensingRange = 0.6;

%% DETECTION PARAMETERS 
MSD_coop = zeros(numPoints-1,1); 
MSD_ncop = zeros(numPoints-1,1);
W1 = zeros(numPoints-1,numAgents);
W2 = zeros(numPoints-1,numAgents);
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
attackers = [  6 11];
%attackers = [6 11];
%attackers = [12 14 17 32 36 39 63 66 68 82 85 88 ];
attackers_new = [];
ra = 0.0005;
attacker_phi = zeros(numTaps, numAgents);
w0_attacker = [0.2;0.2];          % attacker's goal state
Agents = 1:numAgents;
normalAgents = setdiff(Agents,attackers);
wAverageMatrix = [];

%% INPUTS (GAUSSIAN)
mu_x = 0;
sigma_x2 = 0.3 + 0.15*rand(numAgents,1);
x = zeros(numPoints,numAgents);
for k = 1:numAgents
    x(:,k) = mvnrnd(mu_x, sigma_x2(k), numPoints);
end

sigma_v2 = 0.5+0.05*rand(numAgents,1);
v = zeros(numPoints,numAgents);
for k = 1:numAgents
    v(:,k) = mvnrnd(0, sigma_v2(k), numPoints);
end

% NETWORK TOPOLOGY
[Adjacency, AgentSet, group1, group2] = getAdjacency(numAgents, sensingRange);
AdjacencyMatrix = Adjacency;
plotNetworkTopology(1, AgentSet, w0, AdjacencyMatrix, group1, group2, attackers, attackers_new, numTaps);

d = zeros(numPoints,numAgents);
for k = group1
    d(:,k) = filter(w0(:,1), 1, x(:,k));
end
for k = group2
    d(:,k) = filter(w0(:,2), 1, x(:,k));
end
d = d+v;

%estimated_A = zeros(numAgents,numAgents);
%estimated_A_bef_norm = estimated_A;
%error_A = zeros(numTaps,numAgents);

estimated_A = cell(1,numAgents);
sz = zeros(numAgents,1);
for k = normalAgents
    sz(k) = size(find(Adjacency(:,k)==1),1);
    estimated_A(k) = mat2cell(1/sz(k) * rand(sz(k),1),sz(k),1);
end
error_A = zeros(numPoints,numAgents,numTaps);
error = zeros(numPoints,numAgents,numTaps);
phi_last_iter = phi;
difference_error_A_with_error = zeros(numPoints,numAgents,numTaps);
estimated_w = w;


%% DIFFUSION LMS ALGORITHM
for n = numTaps : numPoints
    
    newAdjacency = Adjacency;

    if mod(n,1000) == 0 || n == numPoints
        %update adjacency matrix
        AdjacencyMatrix = (A>=0.01) & Adjacency;
        plotNetworkTopology(n, AgentSet, w0, AdjacencyMatrix, group1, group2, attackers, attackers_new, numTaps)
        pause(1)
    end
    
%     if size(attackers,2) ~= 0
%        for k = normalAgents
%            attacker_phi(1,:,k) = w(:,k) + ra(1)*(w0_attacker-w(:,k));
%            attacker_phi(2,:,k) = w(:,k) + ra(2)*(w0_attacker-w(:,k));
%        end
%     end

    phi_last_iter = phi;
    
    for k = normalAgents
        u(:,k) = x(n:-1:n-numTaps+1,k);     % select part of training input
        % Cooperative state
        h(n,k) = u(:,k)'*w(:,k);            % hypothsis function
        e(n,k) = d(n,k)-h(n,k);             % error
        phi(:,k) = w(:,k) + Mu*( u(:,k)*e(n,k) );
        error(n,k,:) = Mu*( u(:,k)*e(n,k) );
        
        % non-cooperative state
        h_noco(n,k) = u(:,k)' * w_noco(:,k);    % hypothsis function
        e_noco(n,k) = d(n,k)-h_noco(n,k);    % error
        w_noco_old(:,k) = w_noco(:,k);
        w_noco(:,k) = w_noco(:,k) + Mu*( u(:,k)*e_noco(n,k) );
        phi_noco(:,k) = w_noco(:,k);
    end
    
    % attack on node 1. neighbor of node 1: 1, 2, 3, 11, 12
    %error_A = phi(:,1) - phi(:,[1, 2, 3, 11, 12]) * estimated_A;
    %estimated_A_bef_norm = estimated_A + 0.02*(phi(:,[1, 2, 3, 11, 12])'*error_A);
    %estimated_A = estimated_A_bef_norm / sum(estimated_A_bef_norm);
    
    if n > 500
        for k = normalAgents
            for a = attackers
                phi_last_iter(:,a) = attacker_phi(:,k);
            end
            error_A(n,k,:) = phi(:,k) - phi_last_iter(:, [find(Adjacency(:,k)==1)]) * cell2mat(estimated_A(:,k));
            %error_A(n,k,:) = phi(:,k) - phi_last_iter(:, [find(Adjacency(:,k)==1)]) * cell2mat(estimated_A(:,k)) - (phi(:,k) - phi_last_iter(:,a));
            difference_error_A_with_error(n,k,:) = error(n,k,:) - error_A(n,k,:);
            estimated_A_bef_norm = cell2mat(estimated_A(:,k)) + 0.01*phi_last_iter(:, [find(Adjacency(:,k)==1)])'*reshape(error_A(n,k,:),2,1);
            estimated_A_bef_norm([find(estimated_A_bef_norm < 0)]) = 0;
            estimated_A(:,k) = mat2cell( estimated_A_bef_norm/sum(estimated_A_bef_norm), sz(k), 1);
        end
    end
    
    if size(attackers,2) ~= 0
       for k = normalAgents
           % attacker_phi(:,k) = w(:,k) + ra*(w0_attacker-w(:,k));
           estimated_w(:,k) = phi(:,k) - reshape(error_A(n,k,:),2,1);
           %estimated_w(:,k) = phi_last_iter(:, [find(Adjacency(:,k)==1)]) * cell2mat(estimated_A(:,k));
           %estimated_w = phi(:,k) - Mu*( u(:,k)*e(n,k) );
           attacker_phi(:,k) = estimated_w(:,k) + ra*(w0_attacker-estimated_w(:,k));
       end
    end

    gamma2 = UpdateGamma2(normalAgents, attackers, attacker_phi, numAgents, gamma2, newAdjacency, w, phi, niu);

    if n >numTaps
            %[newAdjacency,ratio,J] = removeLargestRatio(n, 2, newAdjacency, numAgents, Adjacency, attackers, D, U, phi, storedNum, attacker_phi, ...
   % Expectation_noco, Expectation_coop, gamma2 );
   % [newAdjacency,ratio,J] = removeLargest_new_resilient(n, 2, newAdjacency, numAgents, Adjacency, attackers, D, U, phi, storedNum, attacker_phi, ...
 %   Expectation_noco, Expectation_coop, gamma2 );
    end
    
    A = UpdateWeight(normalAgents, numAgents, gamma2, newAdjacency);
    w = UpdateW(attackers, numAgents, A, w,  phi, attacker_phi, w0_attacker);
   
    
    % update attackers
    if n > 1000
        attackers_new = UpdateAttackers(numAgents, w, w0_attacker, attackers_new);
    end
    
    [MSD_coop, MSD_ncop] = ComputeMSD(group1, numAgents, attackers, MSD_coop, MSD_ncop, n, w0, w, w_noco, normalAgents);
    W1(n,:) = w(1,:);
    W2(n,:) = w(2,:);
    [D, U] = updateStoredStreamingData(d(n,:),u, D, U);
    wAverageMatrix = getwAverage(numAgents, wAverageMatrix, attackers, Adjacency, w);    
end


%% PLOT
% MSD
figure(2)
set (gcf,'Position',[0,0,450,450], 'color','w');
set(gcf,'color','white');
set(gca,'XTick', [0:1000:5000])
plot(mag2db(MSD_ncop), 'linewidth',1);
hold on;
plot(mag2db(MSD_coop), 'linewidth',1);
hold on;
set(gca,'FontSize',15);
gca = legend('Noncooperative LMS','DLMSAW', 'interpreter','latex', 'fontsize',20);
set(gca,'FontSize',20);
xlabel('Iteration $i$', 'interpreter','latex','fontsize',20);ylabel('MSD(dB)', 'interpreter','latex','fontsize',20);
box on;

figure(3)
%plot(error_node1(1,:))
%hold on;
set (gcf,'Position',[0,0,450,450], 'color','w');
set(gcf,'color','white');
line_color = hsv(numAgents);
for k = setdiff(group1,attackers)
    g1 = plot(difference_error_A_with_error(:,k,1),'Color', 'b','linewidth',1);
    hold on;
end
for k = setdiff(group2,attackers)
    g2 = plot(difference_error_A_with_error(:,k,1),'Color', 'g','linewidth',1);
    hold on;
end
legend([g1(1),g2(1)],'Blue nodes','Green nodes', 'interpreter','latex', 'fontsize',20);
xlim([500 5000])
xlabel('Iteration $i$', 'interpreter','latex','fontsize',20);ylabel('$$|\hat{\boldmath{w}}_{k,i} - \boldmath{w}_{k,i}|$$', 'interpreter','latex','fontsize',22);
box on;

figure(4)
set (gcf,'Position',[0,0,450,450], 'color','w');
set(gcf,'color','white');
%set(gca,'XTick', [0:1000:5000])
bl = plot(W1(:,setdiff(group1,attackers)), 'color','b','linewidth',1);
hold on;
gr = plot(W1(:,setdiff(group2,attackers)), 'color','g','linewidth',1);
hold on;
line([0,5000],[0.9,0.9],'linestyle',':', 'color', 'm','linewidth',2);
hold on;
if ~isempty(attackers)
    line([0,5000],[0.2,0.2],'linestyle',':', 'color', 'm','linewidth',2);
else
    line([0,5000],[0.1,0.1],'linestyle',':', 'color', 'm','linewidth',2);
end
gca = legend([bl(1), gr(1)],'Blue nodes','Green nodes', 'interpreter','latex', 'fontsize',20);
set(gca,'FontSize',15);
ylim([0,1]);
xlabel('Iteration $i$', 'interpreter','latex','fontsize',20);ylabel('State $w_{k,i}(1)$', 'interpreter','latex','fontsize',20);
box on;

figure(5)
set (gcf,'Position',[0,0,450,450], 'color','w');
set(gcf,'color','white');
%set(gca,'XTick', [0:1000:5000])
bl = plot(W2(:,setdiff(group1,attackers)), 'color','b','linewidth',1);
hold on;
gr = plot(W2(:,setdiff(group2,attackers)), 'color','g','linewidth',1);
hold on;
line([0,5000],[0.9,0.9],'linestyle',':', 'color', 'm','linewidth',2);
hold on;
if ~isempty(attackers)
    line([0,5000],[0.2,0.2],'linestyle',':', 'color', 'm','linewidth',2);
else
    line([0,5000],[0.1,0.1],'linestyle',':', 'color', 'm','linewidth',2);
end
gca = legend([bl(1), gr(1)],'Blue nodes','Green nodes', 'interpreter','latex', 'fontsize',20);
set(gca,'FontSize',15);
xlabel('Iteration $i$', 'interpreter','latex','fontsize',20);ylabel('State $w_{k,i}(2)$', 'interpreter','latex','fontsize',20);
ylim([0,1]);
box on;
