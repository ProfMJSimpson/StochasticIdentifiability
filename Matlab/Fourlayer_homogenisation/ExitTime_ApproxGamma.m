function f = ExitTime_ApproxGamma(l,x0,p)

    % Inputs:
    %   l       = [l(1),l(2)]
    %   x0      = initial location
    %   p       = [p(1),p(2)] = [D1,D2]
    %
    % Outputs:
    %   f  = @(t) ....
    
    %% Get raw moments
    [M1,M2] = MomentsRaw2(l(1),l(2));
    
    %% Get central moments
    [M,V] = MomentsCentral2(M1,M2);

    %m = M(x0,p(1),p(2),p(3));
    %v = V(x0,p(1),p(2),p(3));
    
    %% Get Gamma parameters
    %a = m^2 / v;
    %b = v / m;
    [a,b] = GammaParams2(M,V);
    
    
    %% Return interpolating function
    f  = @(t) gampdf(t,a(x0,p(1),p(2)),b(x0,p(1),p(2)));
    
end