for F = 0:5
    x = 0:0.01:1;
    y = F*x+0.1;
    if F == 1
        plot(x,y, 'linewidth',1, 'color', [1 1 0]);
    elseif F == 2
        plot(x,y, 'linewidth',1, 'color', [0 0 1]);
    elseif F == 3
        plot(x,y, 'linewidth',1, 'color', [0 1 1]);
    elseif F == 4
        plot(x,y, 'linewidth',1, 'color', [1 0 0]);
    elseif F == 5
        plot(x,y, 'linewidth',1, 'color', [0 1 0]);
    end
    hold on;
end