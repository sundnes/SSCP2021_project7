function [monitored_names] = land2017_monitored_names()
  % % Monitored value names for ODE model: land2017
  % % ---------------------- ----------------------
  % %
  % % monitored_names = land2017_monitored_names();

  % --- Monitored names --- 
  monitored_names = cell(35, 1);
  monitored_names{1} = 'XS_max';
  monitored_names{2} = 'XW_max';
  monitored_names{3} = 'CaTrpn_max';
  monitored_names{4} = 'kwu';
  monitored_names{5} = 'ksu';
  monitored_names{6} = 'Aw';
  monitored_names{7} = 'As';
  monitored_names{8} = 'cw';
  monitored_names{9} = 'cs';
  monitored_names{10} = 'lambda_min12';
  monitored_names{11} = 'lambda_min087';
  monitored_names{12} = 'h_lambda_prima';
  monitored_names{13} = 'h_lambda';
  monitored_names{14} = 'XU';
  monitored_names{15} = 'gammawu';
  monitored_names{16} = 'gammasu';
  monitored_names{17} = 'cat50';
  monitored_names{18} = 'kb';
  monitored_names{19} = 'Ta';
  monitored_names{20} = 'C';
  monitored_names{21} = 'dCd';
  monitored_names{22} = 'eta';
  monitored_names{23} = 'Fd';
  monitored_names{24} = 'F1';
  monitored_names{25} = 'Tp';
  monitored_names{26} = 'Ttot';
  monitored_names{27} = 'beta';
  monitored_names{28} = 'cai';
  monitored_names{29} = 'dXS_dt';
  monitored_names{30} = 'dXW_dt';
  monitored_names{31} = 'dCaTrpn_dt';
  monitored_names{32} = 'dTmB_dt';
  monitored_names{33} = 'dZetas_dt';
  monitored_names{34} = 'dZetaw_dt';
  monitored_names{35} = 'dCd_dt';
end