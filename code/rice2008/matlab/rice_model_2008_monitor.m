function [monitored] = rice_model_2008_monitor(states, t, parameters)
  % Computes monitored expressions of the rice_model_2008 ODE

  % Assign states
  if length(states)~=11
    error('Expected the states array to be of size 11.');
  end
  intf=states(1); SL=states(2); TRPNCaL=states(3); TRPNCaH=states(4);...
    N_NoXB=states(5); P_NoXB=states(6); N=states(7); XBpostr=states(8);...
    XBprer=states(9); xXBprer=states(10); xXBpostr=states(11);

  % Assign parameters
  if length(parameters)~=56
    error('Expected the parameters array to be of size 56.');
  end
  Qfapp=parameters(1); Qgapp=parameters(2); Qgxb=parameters(3);...
    Qhb=parameters(4); Qhf=parameters(5); fapp=parameters(6);...
    gapp=parameters(7); gslmod=parameters(8); gxb=parameters(9);...
    hb=parameters(10); hbmdc=parameters(11); hf=parameters(12);...
    hfmdc=parameters(13); sigman=parameters(14); sigmap=parameters(15);...
    xbmodsp=parameters(16); KSE=parameters(17); PCon_c=parameters(18);...
    PCon_t=parameters(19); PExp_c=parameters(20); PExp_t=parameters(21);...
    SEon=parameters(22); SL_c=parameters(23); SLmax=parameters(24);...
    SLmin=parameters(25); SLrest=parameters(26); SLset=parameters(27);...
    fixed_afterload=parameters(28); kxb_normalised=parameters(29);...
    massf=parameters(30); visc=parameters(31); Ca_amplitude=parameters(32);...
    Ca_diastolic=parameters(33); start_time=parameters(34);...
    tau1=parameters(35); tau2=parameters(36); TmpC=parameters(37);...
    len_hbare=parameters(38); len_thick=parameters(39);...
    len_thin=parameters(40); x_0=parameters(41); Qkn_p=parameters(42);...
    Qkoff=parameters(43); Qkon=parameters(44); Qkp_n=parameters(45);...
    kn_p=parameters(46); koffH=parameters(47); koffL=parameters(48);...
    koffmod=parameters(49); kon=parameters(50); kp_n=parameters(51);...
    nperm=parameters(52); perm50=parameters(53); xPsi=parameters(54);...
    Trop_conc=parameters(55); kxb=parameters(56);

  % Init return args
  monitored = zeros(65, 1);

  % Expressions for the Sarcomere geometry component
  monitored(50) = ((len_thick/2 < SL/2)*(len_thick/2) + ~(len_thick/2 <...
    SL/2)*(SL/2));
  monitored(51) = ((len_thin - SL/2 > len_hbare/2)*(len_thin - SL/2) +...
    ~(len_thin - SL/2 > len_hbare/2)*(len_hbare/2));
  monitored(52) = -monitored(51) + monitored(50);
  monitored(53) = 2*monitored(52)/(len_thick - len_hbare);
  monitored(54) = monitored(52)/len_thin;

  % Expressions for the Thin filament regulation and crossbridge cycling
  % rates component
  monitored(1) = fapp*xbmodsp*Qfapp^(-37/10 + TmpC/10);
  monitored(2) = 1 + gslmod*(1 - monitored(53));
  monitored(3) = gapp*xbmodsp*Qgapp^(-37/10 + TmpC/10)*monitored(2);
  monitored(4) = exp(-hfmdc*xXBprer^2*sign(xXBprer)/x_0^2);
  monitored(5) = exp(hbmdc*(-x_0 + xXBpostr)^2*sign(-x_0 + xXBpostr)/x_0^2);
  monitored(6) = hf*xbmodsp*Qhf^(-37/10 + TmpC/10)*monitored(4);
  monitored(7) = hb*xbmodsp*Qhb^(-37/10 + TmpC/10)*monitored(5);
  monitored(8) = ((xXBpostr < x_0)*(exp(sigmap*(x_0 - xXBpostr)^2/x_0^2)) +...
    ~(xXBpostr < x_0)*(exp(sigman*(-x_0 + xXBpostr)^2/x_0^2)));
  monitored(9) = gxb*xbmodsp*Qgxb^(-37/10 + TmpC/10)*monitored(8);

  % Expressions for the Normalised active and passive force component
  monitored(10) = (fapp*gxb + fapp*hb)/(fapp*gxb + fapp*hb + fapp*hf +...
    gapp*gxb + gapp*hb + gxb*hf);
  monitored(11) = fapp*hf/(fapp*gxb + fapp*hb + fapp*hf + gapp*gxb + gapp*hb...
    + gxb*hf);
  monitored(12) = kxb_normalised*x_0*monitored(11);
  monitored(13) = kxb_normalised*(XBpostr*xXBpostr +...
    XBprer*xXBprer)*monitored(53);
  monitored(14) = monitored(13)/monitored(12);
  monitored(15) = PCon_t*(-1 + exp(PExp_t*abs(SLrest - SL)))*sign(-SLrest +...
    SL);
  monitored(16) = ((SL > SL_c)*(PCon_c*(-1 + exp(PExp_c*abs(SL_c - SL)))) +...
    ~(SL > SL_c)*(0));
  monitored(17) = monitored(15) + monitored(16);
  monitored(18) = PCon_t*(-1 + exp(PExp_t*abs(SLrest - SLset)))*sign(SLset -...
    SLrest);
  monitored(19) = ((SEon == 1)*(KSE*(SLset - SL)) + ~(SEon ==...
    1)*(fixed_afterload));
  monitored(20) = ((SL > SLmin & SL <= SLmax)*((visc*(SLset - SL) +...
    intf)/massf) + ~(SL > SLmin & SL <= SLmax)*(0));
  monitored(55) = -monitored(14) - monitored(17) + monitored(18) +...
    monitored(19);
  monitored(56) = monitored(20);

  % Expressions for the Equation for simulated calcium transient component
  monitored(21) = (tau1/tau2)^(-1/(-1 + tau1/tau2)) - (tau1/tau2)^(-1/(1 -...
    tau2/tau1));
  monitored(22) = ((t > start_time)*(Ca_diastolic + (Ca_amplitude -...
    Ca_diastolic)*(-exp((start_time - t)/tau2) + exp((start_time -...
    t)/tau1))/monitored(21)) + ~(t > start_time)*(Ca_diastolic));

  % Expressions for the Ca binding to troponin to thin filament regulation
  % component
  monitored(23) = kon*Qkon^(-37/10 + TmpC/10);
  monitored(24) = koffL*koffmod*Qkoff^(-37/10 + TmpC/10);
  monitored(25) = koffH*koffmod*Qkoff^(-37/10 + TmpC/10);
  monitored(26) = -TRPNCaL*monitored(24) + (1 -...
    TRPNCaL)*monitored(22)*monitored(23);
  monitored(27) = -TRPNCaH*monitored(25) + (1 -...
    TRPNCaH)*monitored(22)*monitored(23);
  monitored(28) = (1 - monitored(54))*TRPNCaL + TRPNCaH*monitored(54);
  monitored(29) = sqrt(abs(1.0/(1 + (perm50/monitored(28))^nperm)));
  monitored(30) = ((1.0/monitored(29) < 100)*(1.0/monitored(29)) +...
    ~(1.0/monitored(29) < 100)*(100));
  monitored(57) = monitored(26);
  monitored(58) = monitored(27);
  monitored(31) = kn_p*Qkn_p^(-37/10 + TmpC/10)*monitored(29);
  monitored(32) = kp_n*Qkp_n^(-37/10 + TmpC/10)*monitored(30);

  % Expressions for the Regulation and crossbridge cycling state equations
  % component
  monitored(59) = P_NoXB*monitored(32) - N_NoXB*monitored(31);
  monitored(60) = N_NoXB*monitored(31) - P_NoXB*monitored(32);
  monitored(33) = XBprer*monitored(6) - XBpostr*monitored(7) -...
    XBpostr*monitored(9);
  monitored(34) = 1 - N - XBpostr - XBprer;
  monitored(61) = monitored(32)*monitored(34) - N*monitored(31);
  monitored(35) = XBpostr*monitored(7) + monitored(1)*monitored(34) -...
    XBprer*monitored(3) - XBprer*monitored(6);
  monitored(62) = monitored(33);
  monitored(63) = monitored(35);

  % Expressions for the Mean strain of strongly bound states component
  monitored(36) = (monitored(1)*monitored(7) +...
    monitored(1)*monitored(9))/(monitored(1)*monitored(6) +...
    monitored(1)*monitored(7) + monitored(1)*monitored(9) +...
    monitored(3)*monitored(7) + monitored(3)*monitored(9) +...
    monitored(6)*monitored(9));
  monitored(37) = monitored(1)*monitored(6)/(monitored(1)*monitored(6) +...
    monitored(1)*monitored(7) + monitored(1)*monitored(9) +...
    monitored(3)*monitored(7) + monitored(3)*monitored(9) +...
    monitored(6)*monitored(9));
  monitored(38) = monitored(20)/2 + xPsi*((-x_0 - xXBprer +...
    xXBpostr)*monitored(7) - monitored(1)*xXBprer)/monitored(36);
  monitored(39) = monitored(20)/2 + xPsi*(x_0 - xXBpostr +...
    xXBprer)*monitored(6)/monitored(37);
  monitored(64) = monitored(38);
  monitored(65) = monitored(39);

  % Expressions for the Calculation of micromolar per millisecondes of Ca for
  % apparent Ca binding component
  monitored(40) = (XBpostr + XBprer)/(monitored(10) + monitored(11));
  monitored(41) = (monitored(33) + monitored(35))/(monitored(10) +...
    monitored(11));
  monitored(42) = ((SL < len_thick)*(-0.5*monitored(20)) + ~(SL <...
    len_thick)*(0));
  monitored(43) = ((-SL + 2*len_thin > len_hbare)*(-0.5*monitored(20)) +...
    ~(-SL + 2*len_thin > len_hbare)*(0));
  monitored(44) = -monitored(43) + monitored(42);
  monitored(45) = monitored(44)/len_thin;
  monitored(46) = 2*monitored(44)/(len_thick - len_hbare);
  monitored(47) = Trop_conc*((1 - monitored(54))*TRPNCaL + ((1 -...
    monitored(40))*TRPNCaL + TRPNCaH*monitored(40))*monitored(54));
  monitored(48) = Trop_conc*((1 - monitored(54))*monitored(26) + ((1 -...
    monitored(40))*TRPNCaL + TRPNCaH*monitored(40))*monitored(45) + ((1 -...
    monitored(40))*monitored(26) + TRPNCaH*monitored(41) +...
    monitored(27)*monitored(40) - TRPNCaL*monitored(41))*monitored(54) -...
    TRPNCaL*monitored(45));
  monitored(49) = kxb*(XBpostr*xXBpostr + XBprer*xXBprer)*monitored(46) +...
    kxb*(XBpostr*monitored(39) + XBprer*monitored(38) +...
    monitored(33)*xXBpostr + monitored(35)*xXBprer)*monitored(53);
end