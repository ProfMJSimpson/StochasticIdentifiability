clearvars 
p1=0.2;
p2=0.3;
p3=0.4;
P0=[0.2, 0.2, 0.2];
lb=[0.05, 0.05, 0.05];
ub=[0.5, 0.5, 0.5];
MC=50;
x01=100;
x02=60;
x03=30;
l1=30;
l2=30;
l3=40;
npts=20;
[T1] =Stochastic3(l1, l2, l3, x01, p1, p2, p3, MC,1);
[T2] =Stochastic3(l1, l2, l3, x02, p1, p2, p3, MC,1);
[T3] =Stochastic3(l1, l2, l3, x03, p1, p2, p3, MC,1);
fune=@(p) -loglikelihood(T1, [l1, l2, l3], x01, p, 'exact')...
          -loglikelihood(T2, [l1, l2, l3], x02, p, 'exact')...
          -loglikelihood(T3, [l1, l2, l3], x03, p, 'exact');
      
[p_mle_e,nLLmin] = fmincon(fune,P0,[],[],[],[],lb,ub);
LLmaxe = -nLLmin;

funa=@(p) -loglikelihood(T1, [l1, l2, l3], x01, p, 'approx')...
          -loglikelihood(T2, [l1, l2, l3], x02, p, 'approx')...
          -loglikelihood(T3, [l1, l2, l3], x03, p, 'approx'); 
[p_mle_a,nLLmin] = fmincon(funa,P0,[],[],[],[],lb,ub);
LLmaxa = -nLLmin;

% Profile P1 using the exact likelihood
P1_lower   = linspace(lb(1),p_mle_e(1),npts);
P1_upper = linspace(p_mle_e(1),ub(1),npts);
P1 = [P1_lower P1_upper];
PP1  = zeros(size(P1));
p0 = P0([2,3]);
for i = 1:length(P1)
    % Minimize out P2 and P3
    fun = @(p) fune([P1(i),p(1),p(2)]);
    [p0,PP1(i)] = fmincon(fun,p0,[],[],[],[],lb([2,3]),ub([2,3]));
    
end
figure
PP1 = -PP1 - LLmaxe;
subplot(1,3,1);
plot(P1,PP1,"LineWidth",2); hold on;
plot([0,0.5],[-1.92,-1.92],"r:","LineWidth",2);
plot([p1,p1],[-3,1],"k:","LineWidth",2);
ylim([-3,0]);
xlim([0, 0.5]);
xlabel("p1");
ylabel("Profile likelihood");
legend("Profile","-1.92","True Value");

% Profile P2 using the exact likelihood
P2_lower   = linspace(lb(2),p_mle_e(2),npts);
P2_upper = linspace(p_mle_e(2),ub(2),npts);
P2 = [P2_lower P2_upper];
PP2  = zeros(size(P2));
p0 = P0([1,3]);
for i = 1:length(P2)
    % Minimize out P1 and P3
    fun = @(p) fune([p(1),P2(i),p(2)]);
    [p0,PP2(i)] = fmincon(fun,p0,[],[],[],[],lb([1,3]),ub([1,3]));
end
PP2 = -PP2 - LLmaxe;
subplot(1,3,2);
plot(P2,PP2,"LineWidth",2); hold on;
plot([0,0.5],[-1.92,-1.92],"r:","LineWidth",2);
plot([p2,p2],[-3,1],"k:","LineWidth",2);
ylim([-3,0]);
xlim([0, 0.5]);
xlabel("p2");
ylabel("Profile likelihood");
legend("Profile","-1.92","True Value");

% Profile P3 using the exact likelihood
P3_lower   = linspace(lb(3),p_mle_e(3),npts);
P3_upper = linspace(p_mle_e(3),ub(3),npts);
P3 = [P3_lower P3_upper];
PP3  = zeros(size(P3));
p0 = P0([1,2]);
for i = 1:length(P3)
    % Minimize out P1 and P2
    fun = @(p) fune([p(1),p(2),P3(i)]);
    [p0,PP3(i)] = fmincon(fun,p0,[],[],[],[],lb([1,2]),ub([1,2]));
end
PP3 = -PP3 - LLmaxe;
subplot(1,3,3);
plot(P3,PP3,"LineWidth",2); hold on;
plot([0,0.5],[-1.92,-1.92],"r:","LineWidth",2);
plot([p3,p3],[-3,1],"k:","LineWidth",2);
ylim([-3,0]);
xlim([0, 0.5]);
xlabel("p3");
ylabel("Profile likelihood");
legend("Profile","-1.92","True Value");



% Profile P1 using the approximate Gamma likelihood
P1_lower   = linspace(lb(1),p_mle_a(1),npts);
P1_upper = linspace(p_mle_a(1),ub(1),npts);
P1 = [P1_lower P1_upper];
PP1  = zeros(size(P1));
p0 = P0([2,3]);
for i = 1:length(P1)
    % Minimize out P2 and P3
    fun = @(p) funa([P1(i),p(1),p(2)]);
    [p0,PP1(i)] = fmincon(fun,p0,[],[],[],[],lb([2,3]),ub([2,3]));
    
end
PP1 = -PP1 - LLmaxa;
%figure;
subplot(1,3,1);
plot(P1,PP1,"g","LineWidth",2); hold on;
%plot([0,0.5],[-1.92,-1.92],"r:","LineWidth",2);
%plot([p1,p1],[-3,1],"k:","LineWidth",2);
ylim([-3,0]);
xlim([0, 0.5]);
%xlabel("p1");
%ylabel("Profile likelihood");
%legend("Profile","-1.92","True Value");

% Profile P2 using the approximate Gamma likelihood
P2_lower   = linspace(lb(2),p_mle_a(2),npts);
P2_upper = linspace(p_mle_a(2),ub(2),npts);
P2 = [P2_lower P2_upper];
PP2  = zeros(size(P2));
p0 = P0([1,3]);
for i = 1:length(P2)
    % Minimize out P1 and P3
    fun = @(p) funa([p(1),P2(i),p(2)]);
    [p0,PP2(i)] = fmincon(fun,p0,[],[],[],[],lb([1,3]),ub([1,3]));
end
PP2 = -PP2 - LLmaxa;
subplot(1,3,2);
plot(P1,PP2,"g","LineWidth",2); hold on;
%plot([0,0.5],[-1.92,-1.92],"r:","LineWidth",2);
%plot([p2,p2],[-3,1],"k:","LineWidth",2);
%ylim([-3,0]);
%xlim([0, 0.5]);
%xlabel("p2");
%ylabel("Profile likelihood");
%legend("Profile","-1.92","True Value");

% Profile P3 using the approximate Gamma likelihood
P3_lower   = linspace(lb(3),p_mle_a(3),npts);
P3_upper = linspace(p_mle_a(3),ub(3),npts);
P3 = [P3_lower P3_upper];
PP3  = zeros(size(P3));
p0 = P0([1,2]);
for i = 1:length(P3)
    % Minimize out P1 and P2
    fun = @(p) funa([p(1),p(2),P3(i)]);
    [p0,PP3(i)] = fmincon(fun,p0,[],[],[],[],lb([1,2]),ub([1,2]));
end
PP3 = -PP3 - LLmaxa;
subplot(1,3,3);
plot(P1,PP3,"g","LineWidth",2); hold on;
%plot([0,0.5],[-1.92,-1.92],"r:","LineWidth",2);
%plot([p3,p3],[-3,1],"k:","LineWidth",2);
%ylim([-3,0]);
%xlim([0, 0.5]);
%xlabel("p3");
%ylabel("Profile likelihood");
%legend("Profile","-1.92","True Value");

% bivariate profile P1/P2
P1   = linspace(lb(1),ub(1),npts);
P2   = linspace(lb(2),ub(2),npts);

PP1  = zeros(length(P1),length(P2));
p0 = P0(1);
for i = 1:length(P1)
    for j = 1:length(P2)
    
    % Minimize out P3
    fun = @(p) funa([P1(i),P2(j),p(1)]);
    [p0,PP1(i,j)] = fmincon(fun,p0,[],[],[],[],lb(3),ub(3));
    end
end
PP1 = -PP1 - LLmaxa;
figure;
subplot(1,3,1)
contourf(P1,P2,PP1,'LevelList',[-3.00],'Linewidth',2);hold on
plot(p2,p1,'go','MarkerFaceColor','g')
axis image
xlabel("p2");
ylabel("p1");

 % Profile P1/P3
 P1   = linspace(lb(1),ub(1),npts);
 P3   = linspace(lb(3),ub(3),npts);
 
 PP2  = zeros(length(P1),length(P3));
 p0 = P0(1);
 for i = 1:length(P1)
     for j = 1:length(P3)
     
     % Minimize out P2
     fun = @(p) funa([P1(i),p(1),P3(j)]);
     [p0,PP2(i,j)] = fmincon(fun,p0,[],[],[],[],lb(3),ub(3));
     end
 end
 PP2 = -PP2 - LLmaxa;
 subplot(1,3,2)
 contourf(P2,P3,PP2,'LevelList',[-3.00],'Linewidth',2);hold on
 plot(p3,p1,'go','MarkerFaceColor','g')
 axis image
 xlabel("p3");
 ylabel("p1");
 
% Profile P2/P3
P2   = linspace(lb(2),ub(2),npts);
P3   = linspace(lb(3),ub(3),npts);

PP3  = zeros(length(P2),length(P3));
p0 = P0(1);
for i = 1:length(P2)
    for j = 1:length(P3)
    
    % Minimize out P1
    fun = @(p) funa([p(1),P2(i),P3(j)]);
    [p0,PP3(i,j)] = fmincon(fun,p0,[],[],[],[],lb(3),ub(3));
    end
end
PP3 = -PP3 - LLmaxa;
subplot(1,3,3)
contourf(P2,P3,PP3,'LevelList',[-3.00],'Linewidth',2);hold on
plot(p3,p2,'go','MarkerFaceColor','g')
axis image
xlabel("p3");
ylabel("p2");

% exact bivariate profile P1/P2
P1   = linspace(lb(1),ub(1),npts);
P2   = linspace(lb(2),ub(2),npts);

PP1  = zeros(length(P1),length(P2));
p0 = P0(1);
for i = 1:length(P1)
    for j = 1:length(P2)
    
    % Minimize out P3
    fun = @(p) fune([P1(i),P2(j),p(1)]);
    [p0,PP1(i,j)] = fmincon(fun,p0,[],[],[],[],lb(1),ub(1));
    end
end
PP1 = -PP1 - LLmaxe;
figure;
subplot(1,3,1)
contourf(P1,P2,PP1,'LevelList',[-3.00],'Linewidth',2);hold on
plot(p2,p1,'ro','MarkerFaceColor','r')
axis image
xlabel("p2");
ylabel("p1");

 % exact bivariate profile P1/P3

 P1   = linspace(lb(1),ub(1),npts);
 P3   = linspace(lb(3),ub(3),npts);
 
 PP2  = zeros(length(P1),length(P3));
 p0 = [P0([1])];
 for i = 1:length(P1)
     for j = 1:length(P3)
     
     % Minimize out P2
     fun = @(p) fune([P1(i),p(1),P3(j)]);
     [p0,PP2(i,j)] = fmincon(fun,p0,[],[],[],[],lb(2),ub(2));
     end
 end
 PP2 = -PP2 - LLmaxe;
 subplot(1,3,2)
 contourf(P2,P3,PP2,'LevelList',[-3.00],'Linewidth',2);hold on
 plot(p3,p1,'ro','MarkerFaceColor','r')
 axis image
 xlabel("p3");
 ylabel("p1");
 
%exact bivariate profile P2/P3

P2   = linspace(lb(2),ub(2),npts);
P3   = linspace(lb(3),ub(3),npts);

PP3  = zeros(length(P2),length(P3));
p0 = P0(1);
for i = 1:length(P2)
    for j = 1:length(P3)
    
    % Minimize out P1
    fun = @(p) fune([p(1),P2(i),P3(j)]);
    [p0,PP3(i,j)] = fmincon(fun,p0,[],[],[],[],lb(3),ub(3));
    end
end
PP3 = -PP3 - LLmaxe;
subplot(1,3,3)
contourf(P2,P3,PP3,'LevelList',[-3.00],'Linewidth',2);hold on
plot(p3,p2,'go','MarkerFaceColor','r')
axis image
xlabel("p3");
ylabel("p2");



