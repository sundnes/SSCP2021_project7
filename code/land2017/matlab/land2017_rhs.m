function [values] = land2017_rhs(states, t, parameters)
  % Compute the right hand side of the land2017 ODE

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
  Beta1=parameters(2); Tot_A=parameters(3); Trpn50=parameters(5);...
    cat50_ref=parameters(6); dLambda=parameters(7); etal=parameters(8);...
    etas=parameters(9); gammas=parameters(10); gammaw=parameters(11);...
    ktrpn=parameters(12); ku=parameters(13); kuw=parameters(14);...
    kws=parameters(15); lmbda=parameters(16); ntm=parameters(17);...
    ntrpn=parameters(18); p_k=parameters(21); phi=parameters(22);...
    rs=parameters(23); rw=parameters(24); Ca_amplitude=parameters(25);...
    Ca_diastolic=parameters(26); start_time=parameters(27);...
    tau1=parameters(28); tau2=parameters(29);

  % Init return args
  values = zeros(7, 1);

  % Expressions for the Calcium transient (from Rice et al (2008)) component
  beta = (tau1/tau2)^(-1/(-1 + tau1/tau2)) - (tau1/tau2)^(-1/(1 - tau2/tau1));
  cai = ((t > start_time)*(Ca_diastolic + (Ca_amplitude -...
    Ca_diastolic)*(-exp((start_time - t)/tau2) + exp((start_time -...
    t)/tau1))/beta) + ~(t > start_time)*(Ca_diastolic));

  % Expressions for the mechanics component
  kwu = -kws + kuw*(-1 + 1.0/rw);
  ksu = kws*rw*(-1 + 1.0/rs);
  Aw = Tot_A*rs/(rs + rw*(1 - rs));
  As = Aw;
  cw = kuw*phi*(1 - rw)/rw;
  cs = kws*phi*rw*(1 - rs)/rs;
  lambda_min12 = ((lmbda < 1.2)*(lmbda) + ~(lmbda < 1.2)*(1.2));
  XU = 1 - TmB - XS - XW;
  gammawu = gammaw*abs(Zetaw);
  gammasu = gammas*((Zetas*(Zetas > 0) > (-1 - Zetas)*(Zetas <...
    -1))*(Zetas*(Zetas > 0)) + ~(Zetas*(Zetas > 0) > (-1 - Zetas)*(Zetas <...
    -1))*((-1 - Zetas)*(Zetas < -1)));
  values(1) = kws*XW - XS*gammasu - XS*ksu;
  values(2) = kuw*XU - kws*XW - XW*gammawu - XW*kwu;
  cat50 = cat50_ref + Beta1*(-1 + lambda_min12);
  values(3) = ktrpn*(-CaTrpn + (cai/cat50)^ntrpn*(1 - CaTrpn));
  kb = ku*Trpn50^ntm/(1 - rs - rw*(1 - rs));
  values(4) = ((CaTrpn^(-ntm/2) < 100)*(CaTrpn^(-ntm/2)) + ~(CaTrpn^(-ntm/2)...
    < 100)*(100))*XU*kb - ku*CaTrpn^(ntm/2)*TmB;
  values(5) = dLambda*As - Zetas*cs;
  values(6) = dLambda*Aw - Zetaw*cw;
  C = -1 + lambda_min12;
  dCd = -Cd + C;
  eta = ((dCd < 0)*(etas) + ~(dCd < 0)*(etal));
  values(7) = p_k*(-Cd + C)/eta;
end