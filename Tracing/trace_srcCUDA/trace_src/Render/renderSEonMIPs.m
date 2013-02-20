function renderSEonMIPs(MIPfig, Segment, d1,d2,d3)
%Segment is in ColumnFirst(ITK) format
%d1 = maxRows = ymax
%d2 = maxCols = xmax
%d3 = maxPages = zmax
%plot works correctly because it follows columnFirst(ITK) format

NN = 17;
figure(MIPfig);
hold on;
W = generateShell(Segment, NN);
 
% x = Segment.mu;
% plot(x(1), x(2),'r.');      % YX plot
% plot(x(3)+d2, x(2),'y.');   % YZ plot
% plot(x(1), x(3)+d1, 'g.');  % XZ plot

% plot on x y pane
plot(W(1,1:NN+1), W(2,1:NN+1) ,'r');
plot(W(1,NN+2:2*NN+2), W(2,NN+2:2*NN+2) ,'b');
plot(W(1,2*NN+3:3*NN+3), W(2,2*NN+3:3*NN+3) ,'g');

% plot on y z pane
plot(W(3,1:NN+1)+d2, W(2,1:NN+1),'r');
plot(W(3,NN+2:2*NN+2)+d2, W(2,NN+2:2*NN+2),'b');
plot(W(3,2*NN+3:3*NN+3)+d2, W(2,2*NN+3:3*NN+3),'g');

% plot on x z pane
plot(W(1,1:NN+1), W(3,1:NN+1)+d1,'r');
plot(W(1,NN+2:2*NN+2), W(3,NN+2:2*NN+2)+d1,'b');
plot(W(1,2*NN+3:3*NN+3), W(3,2*NN+3:3*NN+3)+d1,'g');
hold off;