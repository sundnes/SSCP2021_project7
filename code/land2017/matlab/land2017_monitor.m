function [monitored] = land2017_monitor(states, t, parameters)
  % Computes monitored expressions of the land2017 ODE

  % Assign states
  if length(states)~=7
    error('Expected the states array to be of size 7.');
  end
  XS=states(1); XW=states(2); CaTrpn=states(3); TmB=states(4);...
    Zetas=states(5); Zetaw=states(6); Cd=states(7);

  % Assign parameters
  if length(parameters)~=29
    error('Expected the parameters array to be of size 29.');
  end
  Beta0=parameters(1); Beta1=parameters(2); Tot_A=parameters(3);...
    Tref=parameters(4); Trpn50=parameters(5); cat50_ref=parameters(6);...
    dLambda=parameters(7); etal=parameters(8); etas=parameters(9);...
    gammas=parameters(10); gammaw=parameters(11); ktrpn=parameters(12);...
    ku=parameters(13); kuw=parameters(14); kws=parameters(15);...
    lmbda=parameters(16); ntm=parameters(17); ntrpn=parameters(18);...
    p_a=parameters(19); p_b=parameters(20); p_k=parameters(21);...
    phi=parameters(22); rs=parameters(23); rw=parameters(24);...
    Ca_amplitude=parameters(25); Ca_diastolic=parameters(26);...
    start_time=parameters(27); tau1=parameters(28); tau2=parameters(29);

  % Init return args
  monitored = zeros(35, 1);

  % Expressions for the Calcium transient (from Rice et al (2008)) component
  monitored(27) = (tau1/tau2)^(-1/(-1 + tau1/tau2)) - (tau1/tau2)^(-1/(1 -...
    tau2/tau1));
  monitored(28) = ((t > start_time)*(Ca_diastolic + (Ca_amplitude -...
    Ca_diastolic)*(-exp((start_time - t)/tau2) + exp((start_time -...
    t)/tau1))/monitored(27)) + ~(t > start_time)*(Ca_diastolic));

  % Expressions for the mechanics component
  monitored(1) = ((XS > 0)*(XS) + ~(XS > 0)*(0));
  monitored(2) = ((XW > 0)*(XW) + ~(XW > 0)*(0));
  monitored(3) = ((CaTrpn > 0)*(CaTrpn) + ~(CaTrpn > 0)*(0));
  monitored(4) = -kws + kuw*(-1 + 1.0/rw);
  monitored(5) = kws*rw*(-1 + 1.0/rs);
  monitored(6) = Tot_A*rs/(rs + rw*(1 - rs));
  monitored(7) = monitored(6);
  monitored(8) = kuw*phi*(1 - rw)/rw;
  monitored(9) = kws*phi*rw*(1 - rs)/rs;
  monitored(10) = ((lmbda < 1.2)*(lmbda) + ~(lmbda < 1.2)*(1.2));
  monitored(11) = ((monitored(10) < 0.87)*(monitored(10)) + ~(monitored(10) <...
    0.87)*(0.87));
  monitored(12) = 1 + Beta0*(-1.87 + monitored(10) + monitored(11));
  monitored(13) = ((monitored(12) > 0)*(monitored(12)) + ~(monitored(12) >...
    0)*(0));
  monitored(14) = 1 - TmB - XS - XW;
  monitored(15) = gammaw*abs(Zetaw);
  monitored(16) = gammas*((Zetas*(Zetas > 0) > (-1 - Zetas)*(Zetas <...
    -1))*(Zetas*(Zetas > 0)) + ~(Zetas*(Zetas > 0) > (-1 - Zetas)*(Zetas <...
    -1))*((-1 - Zetas)*(Zetas < -1)));
  monitored(29) = kws*XW - XS*monitored(16) - XS*monitored(5);
  monitored(30) = kuw*monitored(14) - kws*XW - XW*monitored(15) -...
    XW*monitored(4);
  monitored(17) = cat50_ref + Beta1*(-1 + monitored(10));
  monitored(31) = ktrpn*(-CaTrpn + (monitored(28)/monitored(17))^ntrpn*(1 -...
    CaTrpn));
  monitored(18) = ku*Trpn50^ntm/(1 - rs - rw*(1 - rs));
  monitored(32) = ((CaTrpn^(-ntm/2) < 100)*(CaTrpn^(-ntm/2)) +...
    ~(CaTrpn^(-ntm/2) < 100)*(100))*monitored(14)*monitored(18) -...
    ku*CaTrpn^(ntm/2)*TmB;
  monitored(33) = dLambda*monitored(7) - Zetas*monitored(9);
  monitored(34) = dLambda*monitored(6) - Zetaw*monitored(8);
  monitored(19) = Tref*((1 + Zetas)*XS + XW*Zetaw)*monitored(13)/rs;
  monitored(20) = -1 + monitored(10);
  monitored(21) = -Cd + monitored(20);
  monitored(22) = ((monitored(21) < 0)*(etas) + ~(monitored(21) < 0)*(etal));
  monitored(35) = p_k*(-Cd + monitored(20))/monitored(22);
  monitored(23) = monitored(21)*monitored(22);
  monitored(24) = -1 + exp(p_b*monitored(20));
  monitored(25) = p_a*(monitored(23) + monitored(24));
  monitored(26) = monitored(19) + monitored(25);
end