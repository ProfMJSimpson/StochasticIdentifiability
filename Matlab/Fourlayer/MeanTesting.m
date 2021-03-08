
p1=0.5;
p2=0.5;
p3=0.5;
p4=0.5;

x0 = 40;

l1=25;
l2=25;
l3=25;
l4=25;

% Create raw moments function
[M1,M2] = MomentsRaw4(l1,l2,l3,l4);

 MM1 = M1(x0,p1,p2,p3,p4)
 MM2 = M2(x0,p1,p2,p3,p4)