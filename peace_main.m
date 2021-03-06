%%%% Parameter sensitivities with complex-step approximation

%%%%% Sustainable Peace ODE Model 

%%%%% Author: Kate Pearce


%%%%% From causal loop diagram in Liebovitch et al. "Modeling the Dynamics of
%%%%% Sustainable Peace" Springer 2018


clear
close all


num_states = 6;
%%%% alternating: e.g. num_states = 6
%%%% state vector of peace factors: x = [x1, x2, x3, x4, x5, x6];
%%%% x1, x3, x5 = positive peace factors
%%%% x2, x4, x6 = negative peace factors

%%%% initial conditions
X0 = ones(num_states,1);

%%%% strength parameters C_ij
%%%% indices of nonzero entries as cell (go row by row, inc by col entry)
stren_inds = {[1,5],[2,6],[3,1],[3,5],[4,2],[4,6],[5,1],[5,3],[5,6],[6,2],[6,4],[6,5]};
stren_weights =[1.5, 5,    0.3,  1.5,   5,    3,    3,    3,    -5,   5,   0.3, -0.3];
num_nonz_stren = length(stren_inds);


%%%% strength parameter matrix
C = zeros(num_states,num_states);
%%%% strenght 
stren_pars = cell(num_nonz_stren,1);

for entry = 1:num_nonz_stren
    entry_rowcol = stren_inds{entry};
    row = entry_rowcol(1);
    col = entry_rowcol(2);
    weight = stren_weights(entry);
    C(row,col) = weight;
    stren_pars{entry}=sprintf('c%d%d',row,col);
end

%%%% memory parameters
mem_pars = {'mpos', 'gamma'};
%%% positive memory parameter
mpos = 0.2;
%%% scaling factor
gamma = 4.5;
mem_vals = [mpos, gamma];

mem_pars_cell = [mem_pars; num2cell(mem_vals)];

mems = struct(mem_pars_cell{:});

%%%% self-reinforcement parameters: equal to total number of peace factors
selfr_vals = ones(1,num_states); %%% initial guess is all ones
selfr_pars = cell(1,num_states);
for par = 1:num_states
    selfr_pars{par} = sprintf('b%d',par);
end
self_cell = [selfr_pars; num2cell(selfr_vals)];
selfr = struct(self_cell{:});


%%% time interval and solver options
tfinal = 30;
tspan = 0 : 0.01 : tfinal; 
odeoptions = odeset('AbsTol',1e-10, 'RelTol', 1e-10);


[~,Z] = ode15s(@peace_ddt,tspan,X0,odeoptions,mems,selfr,C);

figure()
plot(tspan,Z,'LineWidth',3)
set(gca, 'FontSize', 18)
legend('x1: + Hist. Mem.','x2: - Hist. Mem.','x3: + Fut. Exp.', 'x4: - Fut. Exp','x5: PIR','x6: NIR','Location','SouthEast');


%%% local parameter sensitivity analysis

%%% compute parameter sensitivities with complex-step

%%% step size
h = 1e-10;

paramvec = [mem_vals, selfr_vals, stren_weights];
%%% indices of each
nummems = length(mem_vals);
numselfr = length(selfr_vals);
numstrength = length(stren_weights);

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
    
    memspert_vals = partemp(1:nummems);
    memspert_cell = [mem_pars; num2cell(memspert_vals)];
    memspert = struct(memspert_cell{:});
    
    selfrpert_vals = partemp(nummems+1:nummems+numselfr);
    selfrpert_cell = [selfr_pars; num2cell(selfrpert_vals)];
    selfrpert = struct(selfrpert_cell{:});
    
    strengthpert_vals = partemp(nummems+numselfr+1:end);
    %%% put back into matrix
    Cpert = zeros(num_states,num_states);
    for entry = 1:length(strengthpert_vals)
        entry_rowcol = stren_inds{entry};
        row = entry_rowcol(1);
        col = entry_rowcol(2);
        weight = strengthpert_vals(entry);
        Cpert(row,col) = weight;
    end

    
    [~,Y] = ode15s(@peace_ddt,tspan,X0,odeoptions,memspert,selfrpert,Cpert);

    %%% Y(:,1) model soln for reaction product 
    sens = imag(Y(:,1))/h;

    %%% store parameter sensitivity in corresp column
    Chi(:,parind) = sens;
    
    %%% scaled sensitivity
    Chi_s(:,parind) = Chi(:,parind)*0.2*par_real;
end

paramnames = {'m^{+}','\gamma','b_1','b_2','b_3','b_4','b_5','b_6','c_{15}','c_{26}','c_{31}','c_{35}','c_{42}','c_{46}','c_{51}','c_{53}','c_{56}','c_{62}','c_{64}','c_{65}'};

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

clearvars -except mems selfr C tspan Chi Chi_s Z h