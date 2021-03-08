clearvars 
p1=0.2;
p2=0.3;
p3=0.2;
p4=0.4;
P0=[0.3, 0.3, 0.3, 0.3];
lb=[0.05, 0.05, 0.05, 0.05];
ub=[0.5, 0.5, 0.5, 0.5];
MC=500;
x01=100;
x02=75;
x03=50;
x04=25;
l1=25;
l2=25;
l3=25;
l4=25;
npts=20;

[T1] =Stochastic4(l1, l2, l3, l4, x01, p1, p2, p3, p4, MC,1);
[T2] =Stochastic4(l1, l2, l3, l4, x02, p1, p2, p3, p4, MC,1);
[T3] =Stochastic4(l1, l2, l3, l4, x03, p1, p2, p3, p4, MC,1);
[T4] =Stochastic4(l1, l2, l3, l4, x04, p1, p2, p3, p4, MC,1);

fune=@(p) -loglikelihood(T1, [l1, l2, l3, l4], x01, p, 'exact')... 
          -loglikelihood(T2, [l1, l2, l3, l4], x02, p, 'exact')...
          -loglikelihood(T3, [l1, l2, l3, l4], x03, p, 'exact')...
          -loglikelihood(T4, [l1, l2, l3, l4], x04, p, 'exact');
      
[p_mle_e,nLLmin] = fmincon(fune,P0,[],[],[],[],lb,ub);
LLmaxe = -nLLmin;

funa=@(p) -loglikelihood(T1, [l1, l2, l3, l4], x01, p, 'approx')... 
          -loglikelihood(T2, [l1, l2, l3, l4], x02, p, 'approx')...
          -loglikelihood(T3, [l1, l2, l3, l4], x03, p, 'approx')...
          -loglikelihood(T4, [l1, l2, l3, l4], x04, p, 'approx'); 
[p_mle_a,nLLmin] = fmincon(funa,P0,[],[],[],[],lb, ub);
LLmaxa = -nLLmin;

% Profile P1 using the exact likelihood
P1_lower   = linspace(lb(1),p_mle_e(1),npts);
P1_upper = linspace(p_mle_e(1),ub(1),npts);
P1 = [P1_lower P1_upper];
PP1  = zeros(size(P1));
p0 = P0([2,3,4]);
for i = 1:length(P1)
    i
    % Minimize out P2,P3 and P4
    fun = @(p) fune([P1(i),p(1),p(2),p(3)]);
    [p0,PP1(i)] = fmincon(fun,p0,[],[],[],[],lb([2,3,4]),ub([2,3,4]));
    
end
figure
PP1 = -PP1 - LLmaxe;
subplot(1,4,1);
plot(P1,PP1,"LineWidth",2); hold on;
plot([0,0.5],[-1.92,-1.92],"r:","LineWidth",2);
plot([p1,p1],[-3,1],"k:","LineWidth",2);
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
p0 = P0([1,3,4]);
for i = 1:length(P2)
    i
    % Minimize out P1, P3 and P4
    fun = @(p) fune([p(1),P2(i),p(2),p(3)]);
    [p0,PP2(i)] = fmincon(fun,p0,[],[],[],[],lb([1,3,4]),ub([1,3,4]));
end
PP2 = -PP2 - LLmaxe;
subplot(1,4,2);
plot(P2,PP2,"LineWidth",2); hold on;
plot([0,0.5],[-1.92,-1.92],"r:","LineWidth",2);
plot([p2,p2],[-3,1],"k:","LineWidth",2);
ylim([-3,0]);
xlim([0, 0.5]);
xlabel("p2");
ylabel("Profile likelihood");
%legend("Profile","-1.92","True Value");

% Profile P3 using the exact likelihood
P3_lower   = linspace(lb(3),p_mle_e(3),npts);
P3_upper = linspace(p_mle_e(3),ub(3),npts);
P3 = [P3_lower P3_upper];
PP3  = zeros(size(P3));
p0 = P0([1,2,4]);
for i = 1:length(P3)
    % Minimize out P1, P2 and P4
    fun = @(p) fune([p(1),p(2),P3(i),p(3)]);
    [p0,PP3(i)] = fmincon(fun,p0,[],[],[],[],lb([1,2,4]),ub([1,2,4]));
end
PP3 = -PP3 - LLmaxe;
subplot(1,4,3);
plot(P3,PP3,"LineWidth",2); hold on;
plot([0,0.5],[-1.92,-1.92],"r:","LineWidth",2);
plot([p3,p3],[-3,1],"k:","LineWidth",2);
ylim([-3,0]);
xlim([0, 0.5]);
xlabel("p3");
ylabel("Profile likelihood");
%legend("Profile","-1.92","True Value");



% Profile P4 using the exact likelihood
P4_lower   = linspace(lb(4),p_mle_e(4),npts);
P4_upper = linspace(p_mle_e(4),ub(4),npts);
P4 = [P4_lower P4_upper];
PP4  = zeros(size(P4));
p0 = P0([1,2,3]);
for i = 1:length(P4)
    % Minimize out P1, P2 and P3
    fun = @(p) fune([p(1),p(2),p(3),P4(i)]);
    [p0,PP4(i)] = fmincon(fun,p0,[],[],[],[],lb([1,2,3]),ub([1,2,3]));
end
PP4 = -PP4 - LLmaxe;
subplot(1,4,4);
plot(P4,PP4,"LineWidth",2); hold on;
plot([0,0.5],[-1.92,-1.92],"r:","LineWidth",2);
plot([p4,p4],[-3,1],"k:","LineWidth",2);
ylim([-3,0]);
xlim([0, 0.5]);
xlabel("p4");
ylabel("Profile likelihood");
%legend("Profile","-1.92","True Value");


















% Profile P1 using the approximate Gamma likelihood
P1_lower   = linspace(lb(1),p_mle_a(1),npts);
P1_upper = linspace(p_mle_a(1),ub(1),npts);
P1 = [P1_lower P1_upper];
PP1  = zeros(size(P1));
p0 = P0([2,3,4]);
for i = 1:length(P1)
    i
    % Minimize out P2, P3 and P4
    fun = @(p) funa([P1(i),p(1),p(2),p(3)]);
    [p0,PP1(i)] = fmincon(fun,p0,[],[],[],[],lb([2,3,4]),ub([2,3,4]));
    
end
PP1 = -PP1 - LLmaxa;
subplot(1,4,1);
plot(P1,PP1,"g","LineWidth",2); hold on;
ylim([-3,0]);
xlim([0, 0.5]);


% Profile P2 using the approximate Gamma likelihood
P2_lower   = linspace(lb(2),p_mle_a(2),npts);
P2_upper = linspace(p_mle_a(2),ub(2),npts);
P2 = [P2_lower P2_upper];
PP2  = zeros(size(P2));
p0 = P0([1,3,4]);
for i = 1:length(P2)
    % Minimize out P1, P3 and P4
    fun = @(p) funa([p(1),P2(i),p(2),p(3)]);
    [p0,PP2(i)] = fmincon(fun,p0,[],[],[],[],lb([1,3,4]),ub([1,3,4]));
end
PP2 = -PP2 - LLmaxa;
subplot(1,4,2);
plot(P2,PP2,"g","LineWidth",2); hold on;


% Profile P3 using the approximate Gamma likelihood
P3_lower   = linspace(lb(3),p_mle_a(3),npts);
P3_upper = linspace(p_mle_a(3),ub(3),npts);
P3 = [P3_lower P3_upper];
PP3  = zeros(size(P3));
p0 = P0([1,2,4]);
for i = 1:length(P3)
    % Minimize out P1, P2 and P4
    fun = @(p) funa([p(1),p(2),P3(i),p(3)]);
    [p0,PP3(i)] = fmincon(fun,p0,[],[],[],[],lb([1,2,4]),ub([1,2,4]));
end
PP3 = -PP3 - LLmaxa;
subplot(1,4,3);
plot(P3,PP3,"g","LineWidth",2); hold on;

% Profile P4 using the approximate Gamma likelihood
P4_lower   = linspace(lb(4),p_mle_a(4),npts);
P4_upper = linspace(p_mle_a(4),ub(4),npts);
P4 = [P4_lower P4_upper];
PP4  = zeros(size(P4));
p0 = P0([1,2,3]);
for i = 1:length(P4)
    % Minimize out P1, P2 and P3
    fun = @(p) funa([p(1),p(2),p(3),P4(i)]);
    [p0,PP4(i)] = fmincon(fun,p0,[],[],[],[],lb([1,2,3]),ub([1,2,3]));
end
PP4 = -PP4 - LLmaxa;
subplot(1,4,4);
plot(P4,PP4,"g","LineWidth",2); hold on;




