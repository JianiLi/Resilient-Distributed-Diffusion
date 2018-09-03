    
function  plotNetworkTopology(n, AgentSet, w0, AdjacencyMatrix, group1, group2, attacker, attacker_new, numTaps)
    target = w0;
    pause(1);
    clf;
    H = figure(1);
    set (gcf,'Position',[0,0,500,500], 'color','w')
    axis('equal')
    axis([0 1 0 1]);
    %axis normal;
    %set(gca,'xtick',[0 1]);
    %set(gca,'ytick',[0 1]);
    set(gca,'xtick',[],'xticklabel',[])
    set(gca,'ytick',[],'yticklabel',[])
    set(gca,'FontSize',16);
    set(gcf,'color','white')
    %xlabel('X');ylabel('Y');
    box on;
    %title('Resilient Cooperative Multi-Target Localization')
    hold on;
    gplot(AdjacencyMatrix,AgentSet);
    hg=findobj('type','line');
    set(hg,'linewidth',0.5,'Color',[0.6 0.6 0.6])

    for k = group1
        plot(AgentSet(k,1), AgentSet(k,2), 'ob', 'MarkerFaceColor','b', 'MarkerSize',8);
    end
    for k = group2
        plot(AgentSet(k,1), AgentSet(k,2), 'og', 'MarkerFaceColor','g', 'MarkerSize',8);
    end

    for k = attacker_new
        plot(AgentSet(k,1), AgentSet(k,2), 'or', 'MarkerFaceColor','r', 'MarkerSize',8);
    end
    for k = attacker
        plot(AgentSet(k,1), AgentSet(k,2), 'or', 'MarkerFaceColor','r', 'MarkerSize',8);
    end
    for k = attacker
        plot(AgentSet(k,1), AgentSet(k,2), 'oy', 'MarkerFaceColor','y', 'MarkerSize',4);
    end
    %plot(target(1,1),target(2,1),'ob', 'MarkerSize',3, 'LineWidth',2);
    %plot(target(1,2),target(2,2),'og', 'MarkerSize',3, 'LineWidth',2);
    %plot(target(1,1),target(2,1),'ob', 'MarkerSize',10, 'LineWidth',2);
    %plot(target(1,2),target(2,2),'og', 'MarkerSize',10, 'LineWidth',2);
    %plot(attacker_target(1,1),attacker_target(2,1),'or', 'MarkerSize',3, 'LineWidth',2);
    %plot(attacker_target(1,1),attacker_target(2,1),'or', 'MarkerSize',10, 'LineWidth',2);
    if n == numTaps+2
        saveas(gcf,'network.eps','psc2');
    end