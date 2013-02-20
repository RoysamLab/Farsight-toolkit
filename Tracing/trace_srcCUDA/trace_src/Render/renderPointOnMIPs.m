function renderPointOnMIPs(MIPfig, x, d1,d2,d3)
% x is in ColumnFirst(ITK) format
%d1 = maxRows = ymax
%d2 = maxCols = xmax
%d3 = maxPages = zmax
%plot works correctly because it follows columnFirst(ITK) format

figure(MIPfig);
hold on;
plot(x(1), x(2),'r.');      % YX plot
plot(x(3)+d2, x(2),'y.');   % YZ plot
plot(x(1), x(3)+d1, 'g.');  % XZ plot

hold off;