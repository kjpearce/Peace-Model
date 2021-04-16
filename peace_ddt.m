%%%%% Sustainable Peace ODE Model RHS

%%%%% Author: Kate Pearce


%%%%% From causal loop diagram in Liebovitch et al. "Modeling the Dynamics of
%%%%% Sustainable Peace" Springer 2018



%%%%% 

  function dy = peace_ddt(t,y,memory,selfreinf,strength)

%%%% model parameters 

%%% m: degree of memory for system (neg mem stronger influence than pos) 
m_pos = memory(1);
gamma = memory(2);

m_neg = m_pos * gamma;

%%% b: self-reinforcement of peace factor (unknown values)
b = selfreinf;


%%% nonzero strengths between factors 

c15 = strength(1);
c26 = strength(2);
c31 = strength(3);
c35 = strength(4);
c42 = strength(5);
c46 = strength(6);
c51 = strength(7);
c53 = strength(8);
c56 = strength(9);
c62 = strength(10);
c64 = strength(11);
c65 = strength(12);



%%%% system species
x1 = y(1);
x2 = y(2);
x3 = y(3);
x4 = y(4);
x5 = y(5);
x6 = y(6);


%
  dy = [-m_pos*x1 + b(1) + c15*tanh(x5);
        -m_neg*x2 + b(2) + c26*tanh(x6);
        -m_pos*x3 + b(3) + c31*tanh(x1) + c35*tanh(x5); 
        -m_neg*x4 + b(4) + c42*tanh(x2) + c46*tanh(x6);
        -m_pos*x5 + b(5) + c51*tanh(x1) + c53*tanh(x3) + c56*tanh(x6);
        -m_neg*x6 + b(6) + c62*tanh(x2) + c64*tanh(x4) + c65*tanh(x5)];
  end




