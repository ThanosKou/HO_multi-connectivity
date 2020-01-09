syms a w u
B = [0;0;0;0;0;0;0;1];
A = [-3*u,      a,      a,        0,      0,      0,      0; %0,0
      3*u,-(w+a+2*u),   0,        2*a,     a,      0,      0; %1,0
      0,        w, -(2*u+a),      0,       a,      0,      0; %1,1
      0,        2*u,    0,   -(u+w+2*a),   0,     3*a,      a; %2,0
      0,        0,      2*u,       w,   -(u+2*a),  0,      2*a; %2,1
      0,        0,      0,         u,      0,   -(w+3*a),  0;  %3,0
      0,        0,      0,         0,      u,      w,     -3*a;  %3,1
      1,        1,      1,         1,      1,      1,      1];     %sum_PROB
 
 
X = linsolve(A,B);

X(1)