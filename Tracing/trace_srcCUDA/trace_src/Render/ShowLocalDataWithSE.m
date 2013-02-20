function ShowLocalDataWithSE(Vol,s, tle, b)
if nargin == 3
    b = 25;
end
figure(1);

[d1,d2,d3] = size(Vol);
Bl = round(s.mu - b);
Bh = round(s.mu + b);
v = 255*ones(2*b+1,2*b+1,2*b+1);
for x=Bl(1):Bh(1)
    for y=Bl(2):Bh(2)
        for z=Bl(3):Bh(3)
            if (x>0)&&(y>0)&&(z>0)&&(x<=d2)&&(y<=d1)&&(z<=d3)
                v(y-Bl(2)+1,x-Bl(1)+1,z-Bl(3)+1) = Vol(y,x,z);
            end
        end
    end
end

MIPshow(v,'min');title(tle);
s.mu = s.mu - [Bl(1) Bl(2) Bl(3)]';
renderSEonMIPs(1,s,2*b+1,2*b+1,2*b+1);