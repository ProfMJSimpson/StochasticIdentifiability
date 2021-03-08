function [M,V] = MomentsCentral(M1,M2)

    M = M1;
    V = @(x0,p1,p2) M2(x0,p1,p2) - M1(x0,p1,p2)^2;

end