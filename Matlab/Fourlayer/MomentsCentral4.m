function [M,V] = MomentsCentral(M1,M2)

    M = M1;
    V = @(x0,p1,p2,p3,p4) M2(x0,p1,p2,p3,p4) - M1(x0,p1,p2,p3,p4)^2;

end