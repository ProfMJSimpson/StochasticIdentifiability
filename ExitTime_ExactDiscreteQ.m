function f = ExitTime_ExactDiscreteQ(l,x0,p)

    % Inputs:
    %   l       = [l(1),l(2),l(3)]
    %   x0      = initial location
    %   p       = [p(1),p(2),p(3)] = [D1,D2,D3]
    %   tmax    = 
    %
    % Outputs:
    %   f  = @(t) ....
    
    
    %% Process parameters
    L = sum(l);
    
    
    %% Set up transition matrix
    
    % Probability of moving left
    PL = zeros(1,L);
    PL(1:l(1)) = p(1);
    PL(l(1)+1:l(1)+l(2)) = p(2);
    PL(l(1)+l(2)+1:l(1)+l(2)+l(3)) = p(3);

    % Probability of moving right
    PR = zeros(1,L);
    PR(1:l(1)-1) = p(1);
    PR(l(1):l(1)+l(2)-1) = p(2);
    PR(l(1)+l(2):l(1)+l(2)+l(3)-1) = p(3);

    % Probability of not moving
    PC = 1 - PL - PR;
    
    % Transition matrix
    u  = PL;
    d  = [1,PC];
    l  = [0,PR(1:end-1)];
    P  = diag(d) + diag(u,1) + diag(l,-1);
    
    % Diagonalise
    [V,D] = eig(P);
        
    % CDF
    % F(t) = [1,0,0,...] * P^t * X0
    %      = [1,0,0,...] * V * D^t * iV * X0
    %      = ([1,0,0,...] * V) * D^t * (iV * X0)
    %      = Vl * D^t * Vr
    D  = diag(D);
    Vl = [1,zeros(1,L)] * V;
    X0 = zeros(L+1,1);
    X0(x0+1) = 1;
    Vr = V \ X0; 
    F_s = @(t) Vl * diag(D.^t) * Vr;    % CDF for single input time
    
    % PDF
    f_s = @(t) max(1e-10,F_s(t) - F_s(t-1));   
    f = @(t) vec_call(f_s,t);
        
end

function res = vec_call(f,x)
        
    res = zeros(size(x));
    for i = 1:length(x)
        res(i) = f(x(i));
    end
    
end