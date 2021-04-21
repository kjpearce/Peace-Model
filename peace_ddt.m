%%%%% Sustainable Peace ODE Model RHS

%%%%% Author: Kate Pearce


%%%%% From causal loop diagram in Liebovitch et al. "Modeling the Dynamics of
%%%%% Sustainable Peace" Springer 2018



  function dy = peace_ddt(t,y,m,b,C)

%%% C is the weighted adjacency matrix for corresp graph of ''strengths''
num_states = length(y);
%%% m: degree of memory for system (neg mem stronger influence than pos) 
m_pos = m.mpos;
m_neg = m.mpos * m.gamma;

%%% create vector of m params:
mvec = zeros(num_states,1);
for state = 1:num_states
    if mod(state,2) == 0
        %%% if the state is even
        mvec(state) = m_neg;
    else
        mvec(state) = m_pos;
    end
end

%%% create b vector (self-reinforcement params)
bvec = zeros(num_states,1);
%%% grab field values from b (in inc order of states)
bcell = struct2cell(b); 
for state = 1:num_states
    bvec(state) = bcell{state};
end

lin = -mvec.*y+bvec;

Cterm = C*tanh(y);

dy = lin + Cterm;

  end




