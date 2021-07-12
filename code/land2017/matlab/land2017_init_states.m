function [states, varargout] = land2017_init_states()
  % % Default state values for ODE model: land2017
  % % --------------------------------------------
  % %
  % % states = land2017_init_states();
  % % [states, states_names] = land2017_init_states();

  % --- Default initial state values --- 
  states = zeros(7, 1);
  states(1) = 0; % XS;
  states(2) = 0; % XW;
  states(3) = 0.01; % CaTrpn;
  states(4) = 1; % TmB;
  states(5) = 0; % Zetas;
  states(6) = 0; % Zetaw;
  states(7) = 0; % Cd;

  if nargout == 2

    % --- State names --- 
    state_names = cell(7, 1);
    state_names{1} = 'XS';
    state_names{2} = 'XW';
    state_names{3} = 'CaTrpn';
    state_names{4} = 'TmB';
    state_names{5} = 'Zetas';
    state_names{6} = 'Zetaw';
    state_names{7} = 'Cd';
    varargout(1) = {state_names};
  end
end