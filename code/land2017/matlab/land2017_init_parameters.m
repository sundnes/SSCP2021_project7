function [parameters, varargout] = land2017_init_parameters()
  % % Default parameter values for ODE model: land2017
  % % ------------------------------------------------
  % %
  % % parameters = land2017_init_parameters();
  % % [parameters, parameters_names] = land2017_init_parameter();

  if nargout < 1 || nargout > 2
    error('Expected 1-2 output arguments.');
  end

  % --- Default parameters values --- 
  parameters = zeros(29, 1);
  parameters(1) = 2.3; % Beta0;
  parameters(2) = -2.4; % Beta1;
  parameters(3) = 25; % Tot_A;
  parameters(4) = 120; % Tref;
  parameters(5) = 0.35; % Trpn50;
  parameters(6) = 0.805; % cat50_ref;
  parameters(7) = 0; % dLambda;
  parameters(8) = 200; % etal;
  parameters(9) = 20; % etas;
  parameters(10) = 0.0085; % gammas;
  parameters(11) = 0.615; % gammaw;
  parameters(12) = 0.1; % ktrpn;
  parameters(13) = 0.04; % ku;
  parameters(14) = 0.182; % kuw;
  parameters(15) = 0.012; % kws;
  parameters(16) = 1; % lmbda;
  parameters(17) = 2.4; % ntm;
  parameters(18) = 2; % ntrpn;
  parameters(19) = 2.1; % p_a;
  parameters(20) = 9.1; % p_b;
  parameters(21) = 7; % p_k;
  parameters(22) = 2.23; % phi;
  parameters(23) = 0.25; % rs;
  parameters(24) = 0.5; % rw;
  parameters(25) = 1.45; % Ca_amplitude;
  parameters(26) = 0.09; % Ca_diastolic;
  parameters(27) = 5; % start_time;
  parameters(28) = 20; % tau1;
  parameters(29) = 110; % tau2;

  if nargout == 2

    % --- Parameter names --- 
    parameter_names = cell(29, 1);
    parameter_names{1} = 'Beta0';
    parameter_names{2} = 'Beta1';
    parameter_names{3} = 'Tot_A';
    parameter_names{4} = 'Tref';
    parameter_names{5} = 'Trpn50';
    parameter_names{6} = 'cat50_ref';
    parameter_names{7} = 'dLambda';
    parameter_names{8} = 'etal';
    parameter_names{9} = 'etas';
    parameter_names{10} = 'gammas';
    parameter_names{11} = 'gammaw';
    parameter_names{12} = 'ktrpn';
    parameter_names{13} = 'ku';
    parameter_names{14} = 'kuw';
    parameter_names{15} = 'kws';
    parameter_names{16} = 'lmbda';
    parameter_names{17} = 'ntm';
    parameter_names{18} = 'ntrpn';
    parameter_names{19} = 'p_a';
    parameter_names{20} = 'p_b';
    parameter_names{21} = 'p_k';
    parameter_names{22} = 'phi';
    parameter_names{23} = 'rs';
    parameter_names{24} = 'rw';
    parameter_names{25} = 'Ca_amplitude';
    parameter_names{26} = 'Ca_diastolic';
    parameter_names{27} = 'start_time';
    parameter_names{28} = 'tau1';
    parameter_names{29} = 'tau2';
    varargout(1) = {parameter_names};
  end
end