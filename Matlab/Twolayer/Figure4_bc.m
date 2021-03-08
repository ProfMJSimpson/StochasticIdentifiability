clearvars
p1=0.2;
p2=0.4;
P0=[0.5, 0.5];
lb=[0.01, 0.01];
ub=[0.5, 0.5];
MC=50;
x01=70;
x02=70;
l1=30;
l2=40;
npts=20;
[T1] =Stochastic2(l1, l2, x01, p1, p2, MC,1);
[T2] =Stochastic2(l1, l2, x02, p1, p2, MC,1);
fune=@(p) -loglikelihood(T1, [l1, l2], x01, p, 'exact')...
          -loglikelihood(T2, [l1, l2], x02, p, 'exact');
      
[p_mle_e,nLLmin] = fmincon(fune,P0,[],[],[],[],lb,ub);
LLmaxe = -nLLmin;

funa=@(p) -loglikelihood(T1, [l1, l2], x01, p, 'approx')...
          -loglikelihood(T2, [l1, l2], x02, p, 'approx');
[p_mle_a,nLLmin] = fmincon(funa,P0,[],[],[],[],lb,ub);
LLmaxa = -nLLmin;

% Profile P1 using the exact likelihood
P1_lower   = linspace(lb(1),p_mle_e(1),npts);
P1_upper = linspace(p_mle_e(1),ub(1),npts);
P1 = [P1_lower P1_upper];
PP1  = zeros(size(P1));
p0 = P0(2);
for i = 1:length(P1)
    % Minimize out P2 and P3
    fun = @(p) fune([P1(i),p(1)]);
    [p0,PP1(i)] = fmincon(fun,p0,[],[],[],[],lb(2),ub(2));
    
end
PP1 = -PP1 - LLmaxe;
subplot(1,2,1);
plot(P1,PP1,"r","LineWidth",2); hold on;
plot([0,0.5],[-1.92,-1.92],"k","LineWidth",2);
plot([p1,p1],[-3,1],"k","LineWidth",2);
ylim([-3,0]);
xlim([0, 0.5]);
xlabel("p1");
ylabel("Profile likelihood");
%legend("Profile","-1.92","True Value");

% Profile P2 using the exact likelihood
P2_lower   = linspace(lb(2),p_mle_e(2),npts);
P2_upper = linspace(p_mle_e(2),ub(2),npts);
P2 = [P2_lower P2_upper];
PP2  = zeros(size(P2));
p0 = P0(1);
for i = 1:length(P2)
    % Minimize out P1 and P3
    fun = @(p) fune([p(1),P2(i)]);
    [p0,PP2(i)] = fmincon(fun,p0,[],[],[],[],lb(1),ub(1));
end
PP2 = -PP2 - LLmaxe;
subplot(1,2,2);
plot(P2,PP2,"r","LineWidth",2); hold on;
plot([0,0.5],[-1.92,-1.92],"k","LineWidth",2);
plot([p2,p2],[-3,1],"k","LineWidth",2);
ylim([-3,0]);
xlim([0, 0.5]);
xlabel("p2");
ylabel("Profile likelihood");
%legend("Profile","-1.92","True Value");



% Profile P1 using the approximate Gamma likelihood
P1_lower   = linspace(lb(1),p_mle_a(1),npts);
P1_upper = linspace(p_mle_a(1),ub(1),npts);
P1 = [P1_lower P1_upper];
PP1  = zeros(size(P1));
p0 = P0(2);
for i = 1:length(P1)
    % Minimize out P2 
    fun = @(p) funa([P1(i),p(1)]);
    [p0,PP1(i)] = fmincon(fun,p0,[],[],[],[],lb(2),ub(2));
    
end
PP1 = -PP1 - LLmaxa;
subplot(1,2,1);
plot(P1,PP1,"g","LineWidth",2); hold on;
ylim([-3,0]);
xlim([0, 0.5]);



% Profile P2 using the approximate Gamma likelihood
P2_lower   = linspace(lb(2),p_mle_a(2),npts);
P2_upper = linspace(p_mle_a(2),ub(2),npts);
P2 = [P2_lower P2_upper];
PP2  = zeros(size(P2));
p0 = P0(1);
for i = 1:length(P2)
    % Minimize out P1
    fun = @(p) funa([p(1),P2(i)]);
    [p0,PP2(i)] = fmincon(fun,p0,[],[],[],[],lb(1),ub(1));
end
PP2 = -PP2 - LLmaxa;
subplot(1,2,2);
plot(P2,PP2,"g","LineWidth",2); hold on;
ylim([-3,0]);
xlim([0, 0.5]);

