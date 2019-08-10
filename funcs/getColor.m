function color = getColor(value, min, max)
set(0,'DefaultFigureColormap',feval('prism'));
 cmap = get(groot,'defaultfigurecolormap');
 [m, ~] = size(cmap);
 row = round((value/(max-min))*(m-1)) + 1;
 color = cmap(row, :);  
end