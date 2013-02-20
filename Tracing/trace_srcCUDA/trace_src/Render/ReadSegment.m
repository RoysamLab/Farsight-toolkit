function Segment = ReadSegment(fname)
fid = fopen(fname);
if fid>0
    C = textscan(fid, 'ID:%d TraceID:%d mu:[%f,%f,%f] A:[%f,%f,%f] Q:[%f,%f,%f,%f] e1:[%f] R1:[%f,%f,%f] R2:[%f,%f,%f] R3:[%f,%f,%f] Foregd:%f Backgd:%f Lhood:%f MAD:%f NumNeighbors:%d');
end

for i=1:length(C{1})
    Segment(i).ID = C{1}(i);
    Segment(i).TraceID = C{2}(i);
    Segment(i).mu = [C{3}(i), C{4}(i), C{5}(i)]'+1;
    Segment(i).a1 = C{6}(i);
    Segment(i).a2 = C{7}(i);
    Segment(i).a3 = C{8}(i);
    Segment(i).q = [C{9}(i), C{10}(i),C{11}(i) C{12}(i)];
    Segment(i).e1 = C{13}(i);
    Segment(i).e2 = 1;
    Segment(i).R1 = [C{14}(i), C{15}(i),C{16}(i)]';
    Segment(i).R2 = [C{17}(i), C{18}(i),C{19}(i)]';
    Segment(i).R3 = [C{20}(i), C{21}(i),C{22}(i)]';    
    Segment(i).f = C{23}(i);
    Segment(i).b = C{24}(i);
    Segment(i).L = C{25}(i);
    Segment(i).MAD = C{26}(i);
    Segment(i).numNbr = C{27}(i);
end