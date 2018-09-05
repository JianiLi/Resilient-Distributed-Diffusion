% Test black-box attack: attacker does not know the streaming data
% information, by communication message of the target node and its
% neighbors, the attacker calculates the state

clear; close all; clc
rng('default');

%% PARAMETERS
numAgents = 100;
numTaps = 2;		% channel number
%numPoints = 5050;
numPoints = 5000;
Mu = 0.01;          % step size
niu = 0.01;         % forgetting factor
w = rand(numTaps,numAgents);
w_noco = w;
w0 = [0.1 0.9; 0.1 0.9];
phi = zeros(numTaps,numAgents);
gamma2 = zeros(numAgents,numAgents);
A = zeros(numAgents,numAgents);
sensingRange = 0.16;

%% DETECTION PARAMETERS 
MSD_coop = zeros(numPoints-1,1); 
MSD_ncop = zeros(numPoints-1,1);
difference_w = zeros(numPoints-1,numAgents);
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
%attackers = [];
%attackers = [2 24 33  37 42 63 88];
attackers = [24 37 33 42 63 88];
%attackers = [2 11 24 33 35 37 53 55 63 88];
%attackers = [12 14 17 32 36 39 63 66 68 82 85 88 ];
attackers_new = [];
ra = 0.001;
attacker_phi = zeros(numTaps, numAgents);
w0_attacker = [0.5;0.5];          % attacker's goal state
Agents = 1:numAgents;
normalAgents = setdiff(Agents,attackers);
wAverageMatrix = [];

%% INPUTS (GAUSSIAN)
mu_x = 0;
sigma_x2 = 0.8 + 0.4*rand(numAgents,1);
x = zeros(numPoints,numAgents);
for k = 1:numAgents
    x(:,k) = mvnrnd(mu_x, sigma_x2(k), numPoints);
end

sigma_v2 = 0.00+0.05*rand(numAgents,1);
v = zeros(numPoints,numAgents);
for k = 1:numAgents
    v(:,k) = mvnrnd(0, sigma_v2(k), numPoints);
end

% NETWORK TOPOLOGY
[Adjacency, AgentSet, group1, group2] = getAdjacency(numAgents, sensingRange);
AdjacencyMatrix = Adjacency;
plotNetworkTopology(1, AgentSet, w0, AdjacencyMatrix, group1, group2, attackers, attackers_new, numTaps);
pause(1)

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


%% DIFFUSION LMS ALGORITHM
for n = numTaps : numPoints
    n
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
    
    for k = normalAgents
        for a = attackers
            phi_last_iter(:,a) = attacker_phi(:,k);
        end
        error_A(n,k,:) = phi(:,k) - phi_last_iter(:, [find(Adjacency(:,k)==1)]) * cell2mat(estimated_A(:,k));
        difference_error_A_with_error(n,k,:) = error(n,k,: ) - error_A(n,k,:);
        estimated_A_bef_norm = cell2mat(estimated_A(:,k)) + 0.01*phi_last_iter(:, [find(Adjacency(:,k)==1)])'*reshape(error_A(n,k,:),2,1);
        estimated_A_bef_norm([find(estimated_A_bef_norm < 0)]) = 0;
        estimated_A(:,k) = mat2cell( estimated_A_bef_norm/sum(estimated_A_bef_norm), sz(k), 1);
    end
    
    if size(attackers,2) ~= 0
       for k = normalAgents
           % attacker_phi(:,k) = w(:,k) + ra*(w0_attacker-w(:,k));
           estimated_w_k = phi(:,k) - reshape(error_A(n,k,:),2,1);
           difference_w(n,k) = norm(estimated_w_k - w(:,k));
           %estimated_w = phi(:,k) - Mu*( u(:,k)*e(n,k) );
           attacker_phi(:,k) = estimated_w_k + ra*(w0_attacker-estimated_w_k);
       end
    end

    gamma2 = UpdateGamma2(normalAgents, attackers, attacker_phi, numAgents, gamma2, newAdjacency, w, phi, niu);

     if n >numTaps
       % [newAdjacency,ratio,J] = removeLargestRatio(n, 2, newAdjacency, numAgents, Adjacency, attackers, D, U, phi, storedNum, attacker_phi, ...
   % Expectation_noco, Expectation_coop, gamma2 );
        %    [newAdjacency,ratio,J] = removeLargest_new_resilient(n, 2, newAdjacency, numAgents, Adjacency, attackers, D, U, phi, storedNum, attacker_phi, ...
  %Expectation_noco, Expectation_coop, gamma2 );
     end
    
    A = UpdateWeight(normalAgents, numAgents, gamma2, newAdjacency);
    w = UpdateW(attackers, numAgents, A, w,  phi, attacker_phi, w0_attacker);
   
    
    % update attackers
    if n > 1000
        attackers_new = UpdateAttackers(numAgents, w, w0_attacker, attackers_new);
    end
    
    [MSD_coop, MSD_ncop] = ComputeMSD(group1, numAgents, attackers, MSD_coop, MSD_ncop, n, w0, w, w_noco, normalAgents);
    
    [D, U] = updateStoredStreamingData(d(n,:),u, D, U);
    wAverageMatrix = getwAverage(numAgents, wAverageMatrix, attackers, Adjacency, w); 
    
end


%% PLOT
% MSD
figure(2)
set(gcf,'Position',[0,0,450,450], 'color','w');
set(gcf,'color','white');
set(gca,'XTick', [0:1000:5000])
plot(mag2db(MSD_ncop), 'linewidth',1);
hold on;
plot(mag2db(MSD_coop), 'linewidth',1);
hold on;
set(gca,'FontSize',15);
L = legend('Noncooperative LMS','DLMSAW');
set(L,'Interpreter','latex')
set(L,'FontSize',20);
xlabel('Iteration $i$', 'interpreter','latex','fontsize',20);ylabel('MSD(dB)', 'interpreter','latex','fontsize',20);
box on;

% figure(3)
% %plot(error_node1(1,:))
% %hold on;
% for k = Agents
%     plot(difference_error_A_with_error(:,k,1))
%     hold on;
%     mylgd{k} = ['error', num2str(k)];
% end
% legend(mylgd);

figure(4)
set(gcf,'Position',[0,0,450,450], 'color','w');
set(gcf,'color','white');
set(gca,'XTick', [0:1000:5000])
H = plot(difference_w, 'linewidth',1);
hold on;
set(gca,'FontSize',15);
%L = legend(H,'all normal agents');
set(L,'Interpreter','latex')
set(L,'FontSize',20);
xlabel('Iteration $i$', 'interpreter','latex','fontsize',20);ylabel('$$|\hat{\boldmath{w}}_{k,i} - \boldmath{w}_{k,i}|$$', 'interpreter','latex','fontsize',22);
box on;