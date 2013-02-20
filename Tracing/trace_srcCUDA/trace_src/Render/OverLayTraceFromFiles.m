function OverLayTraceFromFiles(Vol, nodefile, linkfile)

MIPshow(Vol, 'min');
[d1,d2,d3] = size(Vol);
Seg = ReadSegment(nodefile);
mu = [Seg.mu];
numTrace = double(max([Seg.TraceID]));
hold on;
plot(mu(1,:), mu(2,:),'r.');      % YX plot
plot(mu(3,:)+d2, mu(2,:),'y.');   % YZ plot
plot(mu(1,:), mu(3,:)+d1, 'g.');  % XZ plot


col = hsv(numTrace);
col = col(randperm(numTrace),:);

fid = fopen(linkfile,'r');

% % while(~feof(fid))
% %     w = fgetl(fid);
% %     c = str2num(w)+1;
% %     a = Seg(c(1));
% %     
% %     tID = Seg(c(1)).TraceID+1;
% %     for i=2:length(c)
% %         b = Seg(c(i));
% %         line([a.mu(1) b.mu(1)], [a.mu(2) b.mu(2)], 'color', 'red', 'linewidth', 2);
% %         line([a.mu(3) b.mu(3)]+d2, [a.mu(2) b.mu(2)], 'color', 'green', 'linewidth', 2);
% %         line([a.mu(1) b.mu(1)], [a.mu(3) b.mu(3)]+d1, 'color', 'blue', 'linewidth', 2);
% %     end
% % end

hold off;