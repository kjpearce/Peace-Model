%%%% Parameter sensitivities with complex-step approximation

%%%%% Sustainable Peace ODE Model 

%%%%% Author: Kate Pearce


%%%%% From causal loop diagram in Liebovitch et al. "Modeling the Dynamics of
%%%%% Sustainable Peace" Springer 2018


clear
close all

%%%% state vector: x = [x1, x2, x3, x4, x5, x6];
%%%% x1, x3, x5 = positive peace factors
%%%% x2, x4, x6 = negative peace factors

%%%% initial conditions
x1_0 = 1;
x2_0 = 1;
x3_0 = 1;
x4_0 = 1;
x5_0 = 1;
x6_0 = 1;

X0 = [x1_0; x2_0; x3_0; x4_0; x5_0; x6_0];

%%%% memory parameters
mpos = 0.2;
gamma = 4.5;

mems = [mpos; gamma];

%%%% self-reinforcement parameters
b1 = 1;
b2 = 1;
b3 = 1;
b4 = 1;
b5 = 1;
b6 = 1;
selfr = [b1; b2; b3; b4; b5; b6];

%%%% nonzero strength parameters
c1 = 1.5;
c2 = 5;
c3 = [0.3; 1.5];
c4 = [5; 3];
c5 = [3; 3; -5];
c6 = [5; 0.3; -0.3];

strength = [c1; c2; c3; c4; c5; c6];


%%% time interval and solver options
tfinal = 30;
tspan = 0 : 0.01 : tfinal; 
odeoptions = odeset('AbsTol',1e-10, 'RelTol', 1e-10);


[~,Y] = ode15s(@peace_ddt,tspan,X0,odeoptions,mems,selfr,strength);
plot(tspan,Y)