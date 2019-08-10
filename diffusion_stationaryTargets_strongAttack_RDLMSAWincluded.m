clear; close all; clc
rng('default');

%% PARAMETERS
numAgents = 100;
numTaps = 2;		% channel number
%numPoints = 2005;
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
%sensingRange = 0.25;

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
%attackers = [24 37 53 63 65 67 88];
attackers = [24 37 63 88];
%attackers = [];
attackers_new = [];
ra = [0.002, 0.002];
attacker_phi = zeros(2,numTaps, numAgents);
w0_attacker = [0.5;0.5];          % attacker's goal state
Agents = 1:numAgents;
normalAgents = setdiff(Agents,attackers);
wAverageMatrix = [];

%% INPUTS (GAUSSIAN)
mu_x = 0;
%sigma_x2 = 0.8 + 0.4*rand(numAgents,1);
sigma_x2 = 0.8 + 0.4*rand(numAgents,1);
sigma_v2 = 0.15+0.05*rand(numAgents,1);
x = zeros(numPoints,numAgents);
for k = 1:numAgents
    x(:,k) = mvnrnd(mu_x, sigma_x2(k), numPoints);
end

%sigma_v2 = 0.15+0.05*rand(numAgents,1);
v = zeros(numPoints,numAgents);
for k = 1:numAgents
    v(:,k) = mvnrnd(0, sigma_v2(k), numPoints);
end

% NETWORK TOPOLOGY
[Adjacency, AgentSet, group1, group2] = getAdjacency(numAgents, sensingRange);
AdjacencyMatrix = Adjacency;
fig=figure(1);
filename = './figures/diffusionAlgorithm.gif';
plotNetworkTopology(1, AgentSet, w0, AdjacencyMatrix, group1, group2, attackers, attackers_new, numTaps);
%title('F = 5');
drawnow
frame = getframe(fig);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256); 
% Write to the GIF File 
imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
imwrite(imind,cm,filename,'gif','WriteMode','append'); 
imwrite(imind,cm,filename,'gif','WriteMode','append'); 
imwrite(imind,cm,filename,'gif','WriteMode','append'); 

d = zeros(numPoints,numAgents);
for k = group1
    d(:,k) = filter(w0(:,1), 1, x(:,k));
end
for k = group2
    d(:,k) = filter(w0(:,2), 1, x(:,k));
end
d = d+v;


%% DIFFUSION LMS ALGORITHM
for n = numTaps : numPoints
    
    %if n > 2000
    %    attackers = [24 37 63 68];
    %    normalAgents = setdiff(Agents,attackers);
    %end
    
    newAdjacency = Adjacency;

    if  mod(n,1000) == 0 || n == numPoints
        %update adjacency matrix
        AdjacencyMatrix = (A>=0.01) & Adjacency;
        
        plotNetworkTopology(n, AgentSet, w0, AdjacencyMatrix, group1, group2, attackers, attackers_new, numTaps)
        %title('F = 5');
        drawnow
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end
    
    if size(attackers,2) ~= 0
        for k = normalAgents
            attacker_phi(1,:,k) = w(:,k) + ra(1)*(w0_attacker-w(:,k));
            attacker_phi(2,:,k) = w(:,k) + ra(2)*(w0_attacker-w(:,k));
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
    
    %% R-DLMSAW starts, select F
    F = 1;
    if n >numTaps
      %  [newAdjacency,ratio,J] = removeLargestRatio(n, 1, newAdjacency, numAgents, Adjacency, attackers, D, U, phi, storedNum, attacker_phi, ...
  %  Expectation_noco, Expectation_coop, gamma2 );
        %[newAdjacency, J] = removeLargestCostOfRemainingNeighbors(0.01, 3, gamma2, AdjacencyMatrix, newAdjacency, numAgents, attackers, D, U, phi, storedNum, attacker_phi);
    [newAdjacency,ratio,J] = removeLargest_new_resilient(n, F, newAdjacency, numAgents, Adjacency, attackers, D, U, phi, storedNum, attacker_phi, ...
                Expectation_noco, Expectation_coop, gamma2 );
    end
    
    %% R-DLMSAW ends
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
    %Expectation_noco = Expectation(n, numAgents, Expectation_noco, D, U, w_noco, storedNum);
    %Expectation_coop = Expectation(n, numAgents, Expectation_coop, D, U, w, storedNum);
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
L = legend('Noncooperative LMS','DLMSAW', 'interpreter','latex', 'fontsize',20);
set(L,'FontSize',20);
xlabel('Iteration $i$', 'interpreter','latex','fontsize',20);ylabel('MSD(dB)', 'interpreter','latex','fontsize',20);
box on;

figure(3)
set (gcf,'Position',[0,0,450,450], 'color','w');
set(gcf,'color','white');
%set(gca,'XTick', [0:1000:5000])
bl = plot(W1(:,setdiff(group1,attackers)), 'color','b','linewidth',1);
hold on;
gr = plot(W1(:,setdiff(group2,attackers)), 'color','g','linewidth',1);
hold on;
line([0,5000],[0.9,0.9],'linestyle',':', 'color', 'm','linewidth',2);
hold on;
line([0,5000],[0.1,0.1],'linestyle',':', 'color', 'm','linewidth',2);
hold on;
set(gca,'TickLabelInterpreter','latex','FontSize',15);
if ~isempty(attackers)
    neighbor_attacker = [78,97,99,72,52,98,14,89,54,26,47,23,25,13,53,34,15,27,36,35,64,79,73,38,68,74,62,87,48,77,61,46,33];
    re = plot(W1(:,neighbor_attacker), 'color','r','linewidth',1);
    L = legend([bl(1), gr(1), re(1)],'Blue nodes','Green nodes', 'Attacked nodes', 'interpreter','latex', 'fontsize',22, 'Location', 'Best');
else
    L = legend([bl(1), gr(1)],'Blue nodes','Green nodes', 'interpreter','latex', 'fontsize',22, 'Location', 'Best');
end
ylim([0,1]);
xlabel('Iteration $i$', 'interpreter','latex','fontsize',25);ylabel('State $w_{k,i}(1)$', 'interpreter','latex','fontsize',25);
box on;

figure(4)
set (gcf,'Position',[0,0,450,450], 'color','w');
set(gcf,'color','white');
%set(gca,'XTick', [0:1000:5000])
bl = plot(W2(:,setdiff(group1,attackers)), 'color','b','linewidth',1);
hold on;
gr = plot(W2(:,setdiff(group2,attackers)), 'color','g','linewidth',1);
hold on;
line([0,5000],[0.9,0.9],'linestyle',':', 'color', 'm','linewidth',2);
hold on;
line([0,5000],[0.1,0.1],'linestyle',':', 'color', 'm','linewidth',2);
hold on;
set(gca,'TickLabelInterpreter','latex', 'FontSize',15);
if ~isempty(attackers)
    neighbor_attacker = [78,97,99,72,52,98,14,89,54,26,47,23,25,13,53,34,15,27,36,35,64,79,73,38,68,74,62,87,48,77,61,46,33];
    re = plot(W2(:,neighbor_attacker), 'color','r','linewidth',1);
    L = legend([bl(1), gr(1), re(1)],'Blue nodes','Green nodes', 'Attacked nodes', 'interpreter','latex', 'fontsize',22, 'Location', 'Best');
else
    L = legend([bl(1), gr(1)],'Blue nodes','Green nodes', 'interpreter','latex', 'fontsize',22, 'Location', 'Best');
end
xlabel('Iteration $i$', 'interpreter','latex','fontsize',25);ylabel('State $w_{k,i}(2)$', 'interpreter','latex','fontsize',25);
ylim([0,1]);
box on;

figure(5);
set (gcf,'Position',[0,0,450,450], 'color','w');
set(gcf,'color','white');
for i = 1:numPoints-1
    delta(i) = norm(w0_attacker - wAverageMatrix(:,i));
end
plot(delta,'linewidth',2);
set(gca,'TickLabelInterpreter','latex', 'FontSize',15);
xlabel('Iteration $i$','Interpreter','LaTex','FontSize',25);
ylabel('$\|w_k^a - \bar{w}_{k,i}\|$','Interpreter','LaTex','FontSize',25);

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
