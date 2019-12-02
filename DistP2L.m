function [h] = DistP2L(LinePt,Pt)
% calculate the distance between
% the point, Pt
% and the line determined by points LinePt = [LinePt1 Linept2];
% all points are 3 dimensional
a = Pt - LinePt(:,1);
b = Pt - LinePt(:,2);
c = LinePt(:,1) - LinePt(:,2);
if c == [0;0;0];
    disp('LinePt1 equals LinePt2!');
end
h = norm(cross(a,b))/norm(c);
end

