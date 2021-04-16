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


[~,Z] = ode15s(@peace_ddt,tspan,X0,odeoptions,mems,selfr,strength);

figure()
plot(tspan,Z,'LineWidth',3)
set(gca, 'FontSize', 18)
legend('x1: + Hist. Mem.','x2: - Hist. Mem.','x3: + Fut. Exp.', 'x4: - Fut. Exp','x5: PIR','x6: NIR','Location','SouthEast');


%%% local parameter sensitivity analysis

%%% compute parameter sensitivities with complex-step

%%% step size
h = 1e-10;

paramvec = [mems; selfr; strength];
%%% indices of each
nummems = length(mems);
numselfr = length(selfr);
numstrength = length(strength);

numpar = length(paramvec);
numpts = length(tspan);

Chi = zeros(numpts, numpar);
Chi_s = Chi; %%% for scaled sensitivity matrix

for parind = 1:numpar
    %%% temp parameter vector to perturb
    partemp = paramvec;

    %%% parameter to perturb
    par_real = partemp(parind); 
    %%% complex perturbation
    par_pert = complex(par_real,h); 

    %%% full parameter vector with complex pert
    partemp(parind) = par_pert; 
    
    memspert = partemp(1:nummems);
    selfrpert = partemp(nummems+1:nummems+numselfr);
    strengthpert = partemp(nummems+numselfr+1:end);
    
    [~,Y] = ode15s(@peace_ddt,tspan,X0,odeoptions,memspert,selfrpert,strengthpert);

    %%% Y(:,1) model soln for reaction product 
    sens = imag(Y(:,1))/h;

    %%% store parameter sensitivity in corresp column
    Chi(:,parind) = sens;
    
    %%% scaled sensitivity
    Chi_s(:,parind) = Chi(:,parind)*0.2*par_real;
end

paramnames = {'m^{+}','{\bf \gamma}','b_1','b_2','b_3','b_4','b_5','b_6','c_{15}','c_{26}','c_{31}','c_{35}','c_{42}','c_{46}','c_{51}','c_{53}','c_{56}','c_{62}','c_{64}','c_{65}'};

%%%% plot sensitivities (of X1) wrt each parameter
figure()
for parind = 1:numpar
    subplot(4,5,parind)
    plot(tspan, Chi(:,parind),'LineWidth',3)
    title(paramnames{parind})
end

%%%% scaled sensitivities
figure()
for parind = 1:numpar
    subplot(4,5,parind)
    plot(tspan, Chi_s(:,parind),'LineWidth',3)
    title(paramnames{parind})
end

%%%% compute scaled information matrix
IM = Chi_s'*Chi_s;

%%%% extract eigenvalues in increasing order
[eigvecs, lambda] = eig(IM);

[eigval,ind] = sort(diag(lambda)); 

logEvals = zeros(1,numpar);
for par=1:numpar
    logEvals(par)=log10(abs(eigval(par)));
end

%%% create lambda names for plot
eigvalname = cell(1,numpar);
for par=1:numpar
    lam = "\lambda";
    str = sprintf('%s_{%d}',lam,par);
    eigvalname{par} = str;
end

figure()
bar(logEvals)
ylabel('|Eigenvalue| (log scale)')
xticks(1:1:numpar)
xticklabels(eigvalname)
title('Eigenvalues of Info Mat')
ax = gca;
ax.FontSize = 18;

%%% plot normalized eigenvectors corresp to smallest eigenvalues
for vals = 1:6
    subplot(2,3,vals)
    bar(eigvecs(:,vals))
    title(sprintf('v_{%d}',vals))
    ylim([-1 1])
    yticks(-1:.25:1)
    xticks(1:3:numpar)
end

clearvars -except mems selfr strength tspan Chi Chi_s Z h