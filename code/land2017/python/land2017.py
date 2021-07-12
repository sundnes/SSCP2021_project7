# Gotran generated code for the "land2017" model
from __future__ import division

import numpy as np
import math


def init_state_values(**values):
    """
    Initialize state values
    """
    # Init values
    # XS=0, XW=0, CaTrpn=0.01, TmB=1, Zetas=0, Zetaw=0, Cd=0
    init_values = np.array([0, 0, 0.01, 1, 0, 0, 0], dtype=np.float_)

    # State indices and limit checker
    state_ind = dict([("XS", 0), ("XW", 1), ("CaTrpn", 2), ("TmB", 3),\
        ("Zetas", 4), ("Zetaw", 5), ("Cd", 6)])

    for state_name, value in values.items():
        if state_name not in state_ind:
            raise ValueError("{0} is not a state.".format(state_name))
        ind = state_ind[state_name]

        # Assign value
        init_values[ind] = value

    return init_values

def init_parameter_values(**values):
    """
    Initialize parameter values
    """
    # Param values
    # Beta0=2.3, Beta1=-2.4, Tot_A=25, Tref=120, Trpn50=0.35,
    # cat50_ref=0.805, dLambda=0, etal=200, etas=20, gammas=0.0085,
    # gammaw=0.615, ktrpn=0.1, ku=0.04, kuw=0.182, kws=0.012,
    # lmbda=1, ntm=2.4, ntrpn=2, p_a=2.1, p_b=9.1, p_k=7, phi=2.23,
    # rs=0.25, rw=0.5, Ca_amplitude=1.45, Ca_diastolic=0.09,
    # start_time=5, tau1=20, tau2=110
    init_values = np.array([2.3, -2.4, 25, 120, 0.35, 0.805, 0, 200, 20,\
        0.0085, 0.615, 0.1, 0.04, 0.182, 0.012, 1, 2.4, 2, 2.1, 9.1, 7, 2.23,\
        0.25, 0.5, 1.45, 0.09, 5, 20, 110], dtype=np.float_)

    # Parameter indices and limit checker
    param_ind = dict([("Beta0", 0), ("Beta1", 1), ("Tot_A", 2), ("Tref", 3),\
        ("Trpn50", 4), ("cat50_ref", 5), ("dLambda", 6), ("etal", 7),\
        ("etas", 8), ("gammas", 9), ("gammaw", 10), ("ktrpn", 11), ("ku",\
        12), ("kuw", 13), ("kws", 14), ("lmbda", 15), ("ntm", 16), ("ntrpn",\
        17), ("p_a", 18), ("p_b", 19), ("p_k", 20), ("phi", 21), ("rs", 22),\
        ("rw", 23), ("Ca_amplitude", 24), ("Ca_diastolic", 25),\
        ("start_time", 26), ("tau1", 27), ("tau2", 28)])

    for param_name, value in values.items():
        if param_name not in param_ind:
            raise ValueError("{0} is not a parameter.".format(param_name))
        ind = param_ind[param_name]

        # Assign value
        init_values[ind] = value

    return init_values

def state_indices(*states):
    """
    State indices
    """
    state_inds = dict([("XS", 0), ("XW", 1), ("CaTrpn", 2), ("TmB", 3),\
        ("Zetas", 4), ("Zetaw", 5), ("Cd", 6)])

    indices = []
    for state in states:
        if state not in state_inds:
            raise ValueError("Unknown state: '{0}'".format(state))
        indices.append(state_inds[state])
    if len(indices)>1:
        return indices
    else:
        return indices[0]

def parameter_indices(*params):
    """
    Parameter indices
    """
    param_inds = dict([("Beta0", 0), ("Beta1", 1), ("Tot_A", 2), ("Tref", 3),\
        ("Trpn50", 4), ("cat50_ref", 5), ("dLambda", 6), ("etal", 7),\
        ("etas", 8), ("gammas", 9), ("gammaw", 10), ("ktrpn", 11), ("ku",\
        12), ("kuw", 13), ("kws", 14), ("lmbda", 15), ("ntm", 16), ("ntrpn",\
        17), ("p_a", 18), ("p_b", 19), ("p_k", 20), ("phi", 21), ("rs", 22),\
        ("rw", 23), ("Ca_amplitude", 24), ("Ca_diastolic", 25),\
        ("start_time", 26), ("tau1", 27), ("tau2", 28)])

    indices = []
    for param in params:
        if param not in param_inds:
            raise ValueError("Unknown param: '{0}'".format(param))
        indices.append(param_inds[param])
    if len(indices)>1:
        return indices
    else:
        return indices[0]

def monitor_indices(*monitored):
    """
    Monitor indices
    """
    monitor_inds = dict([("XS_max", 0), ("XW_max", 1), ("CaTrpn_max", 2),\
        ("kwu", 3), ("ksu", 4), ("Aw", 5), ("As", 6), ("cw", 7), ("cs", 8),\
        ("lambda_min12", 9), ("lambda_min087", 10), ("h_lambda_prima", 11),\
        ("h_lambda", 12), ("XU", 13), ("gammawu", 14), ("gammasu", 15),\
        ("cat50", 16), ("kb", 17), ("Ta", 18), ("C", 19), ("dCd", 20),\
        ("eta", 21), ("Fd", 22), ("F1", 23), ("Tp", 24), ("Ttot", 25),\
        ("beta", 26), ("cai", 27), ("dXS_dt", 28), ("dXW_dt", 29),\
        ("dCaTrpn_dt", 30), ("dTmB_dt", 31), ("dZetas_dt", 32), ("dZetaw_dt",\
        33), ("dCd_dt", 34)])

    indices = []
    for monitor in monitored:
        if monitor not in monitor_inds:
            raise ValueError("Unknown monitored: '{0}'".format(monitor))
        indices.append(monitor_inds[monitor])
    if len(indices)>1:
        return indices
    else:
        return indices[0]

def rhs(states, t, parameters, values=None):
    """
    Compute the right hand side of the land2017 ODE
    """

    # Assign states
    assert(len(states) == 7)
    XS, XW, CaTrpn, TmB, Zetas, Zetaw, Cd = states

    # Assign parameters
    assert(len(parameters) == 29)
    Beta1=parameters[1]; Tot_A=parameters[2]; Trpn50=parameters[4];\
        cat50_ref=parameters[5]; dLambda=parameters[6]; etal=parameters[7];\
        etas=parameters[8]; gammas=parameters[9]; gammaw=parameters[10];\
        ktrpn=parameters[11]; ku=parameters[12]; kuw=parameters[13];\
        kws=parameters[14]; lmbda=parameters[15]; ntm=parameters[16];\
        ntrpn=parameters[17]; p_k=parameters[20]; phi=parameters[21];\
        rs=parameters[22]; rw=parameters[23]; Ca_amplitude=parameters[24];\
        Ca_diastolic=parameters[25]; start_time=parameters[26];\
        tau1=parameters[27]; tau2=parameters[28]

    # Init return args
    if values is None:
        values = np.zeros((7,), dtype=np.float_)
    else:
        assert isinstance(values, np.ndarray) and values.shape == (7,)

    # Expressions for the Calcium transient (from Rice et al (2008)) component
    beta = math.pow(tau1/tau2, -1/(-1 + tau1/tau2)) - math.pow(tau1/tau2,\
        -1/(1 - tau2/tau1))
    cai = (Ca_diastolic + (Ca_amplitude -\
        Ca_diastolic)*(-math.exp((start_time - t)/tau2) +\
        math.exp((start_time - t)/tau1))/beta if t > start_time else\
        Ca_diastolic)

    # Expressions for the mechanics component
    kwu = -kws + kuw*(-1 + 1.0/rw)
    ksu = kws*rw*(-1 + 1.0/rs)
    Aw = Tot_A*rs/(rs + rw*(1 - rs))
    As = Aw
    cw = kuw*phi*(1 - rw)/rw
    cs = kws*phi*rw*(1 - rs)/rs
    lambda_min12 = (lmbda if lmbda < 1.2 else 1.2)
    XU = 1 - TmB - XS - XW
    gammawu = gammaw*math.fabs(Zetaw)
    gammasu = gammas*(Zetas*(Zetas > 0) if Zetas*(Zetas > 0) > (-1 -\
        Zetas)*(Zetas < -1) else (-1 - Zetas)*(Zetas < -1))
    values[0] = kws*XW - XS*gammasu - XS*ksu
    values[1] = kuw*XU - kws*XW - XW*gammawu - XW*kwu
    cat50 = cat50_ref + Beta1*(-1 + lambda_min12)
    values[2] = ktrpn*(-CaTrpn + math.pow(cai/cat50, ntrpn)*(1 - CaTrpn))
    kb = ku*math.pow(Trpn50, ntm)/(1 - rs - rw*(1 - rs))
    values[3] = (math.pow(CaTrpn, -ntm/2) if math.pow(CaTrpn, -ntm/2) < 100 else\
        100)*XU*kb - ku*math.pow(CaTrpn, ntm/2)*TmB
    values[4] = dLambda*As - Zetas*cs
    values[5] = dLambda*Aw - Zetaw*cw
    C = -1 + lambda_min12
    dCd = -Cd + C
    eta = (etas if dCd < 0 else etal)
    values[6] = p_k*(-Cd + C)/eta

    # Return results
    return values

def monitor(states, t, parameters, monitored=None):
    """
    Computes monitored expressions of the land2017 ODE
    """

    # Assign states
    assert(len(states) == 7)
    XS, XW, CaTrpn, TmB, Zetas, Zetaw, Cd = states

    # Assign parameters
    assert(len(parameters) == 29)
    Beta0, Beta1, Tot_A, Tref, Trpn50, cat50_ref, dLambda, etal, etas,\
        gammas, gammaw, ktrpn, ku, kuw, kws, lmbda, ntm, ntrpn, p_a, p_b,\
        p_k, phi, rs, rw, Ca_amplitude, Ca_diastolic, start_time, tau1, tau2 =\
        parameters

    # Init return args
    if monitored is None:
        monitored = np.zeros((35,), dtype=np.float_)
    else:
        assert isinstance(monitored, np.ndarray) and monitored.shape == (35,)

    # Expressions for the Calcium transient (from Rice et al (2008)) component
    monitored[26] = math.pow(tau1/tau2, -1/(-1 + tau1/tau2)) -\
        math.pow(tau1/tau2, -1/(1 - tau2/tau1))
    monitored[27] = (Ca_diastolic + (Ca_amplitude -\
        Ca_diastolic)*(-math.exp((start_time - t)/tau2) +\
        math.exp((start_time - t)/tau1))/monitored[26] if t > start_time else\
        Ca_diastolic)

    # Expressions for the mechanics component
    monitored[0] = (XS if XS > 0 else 0)
    monitored[1] = (XW if XW > 0 else 0)
    monitored[2] = (CaTrpn if CaTrpn > 0 else 0)
    monitored[3] = -kws + kuw*(-1 + 1.0/rw)
    monitored[4] = kws*rw*(-1 + 1.0/rs)
    monitored[5] = Tot_A*rs/(rs + rw*(1 - rs))
    monitored[6] = monitored[5]
    monitored[7] = kuw*phi*(1 - rw)/rw
    monitored[8] = kws*phi*rw*(1 - rs)/rs
    monitored[9] = (lmbda if lmbda < 1.2 else 1.2)
    monitored[10] = (monitored[9] if monitored[9] < 0.87 else 0.87)
    monitored[11] = 1 + Beta0*(-1.87 + monitored[10] + monitored[9])
    monitored[12] = (monitored[11] if monitored[11] > 0 else 0)
    monitored[13] = 1 - TmB - XS - XW
    monitored[14] = gammaw*math.fabs(Zetaw)
    monitored[15] = gammas*(Zetas*(Zetas > 0) if Zetas*(Zetas > 0) > (-1 -\
        Zetas)*(Zetas < -1) else (-1 - Zetas)*(Zetas < -1))
    monitored[28] = kws*XW - XS*monitored[15] - XS*monitored[4]
    monitored[29] = kuw*monitored[13] - kws*XW - XW*monitored[14] -\
        XW*monitored[3]
    monitored[16] = cat50_ref + Beta1*(-1 + monitored[9])
    monitored[30] = ktrpn*(-CaTrpn + math.pow(monitored[27]/monitored[16],\
        ntrpn)*(1 - CaTrpn))
    monitored[17] = ku*math.pow(Trpn50, ntm)/(1 - rs - rw*(1 - rs))
    monitored[31] = (math.pow(CaTrpn, -ntm/2) if math.pow(CaTrpn, -ntm/2) <\
        100 else 100)*monitored[13]*monitored[17] - ku*math.pow(CaTrpn,\
        ntm/2)*TmB
    monitored[32] = dLambda*monitored[6] - Zetas*monitored[8]
    monitored[33] = dLambda*monitored[5] - Zetaw*monitored[7]
    monitored[18] = Tref*((1 + Zetas)*XS + XW*Zetaw)*monitored[12]/rs
    monitored[19] = -1 + monitored[9]
    monitored[20] = -Cd + monitored[19]
    monitored[21] = (etas if monitored[20] < 0 else etal)
    monitored[34] = p_k*(-Cd + monitored[19])/monitored[21]
    monitored[22] = monitored[20]*monitored[21]
    monitored[23] = -1 + math.exp(p_b*monitored[19])
    monitored[24] = p_a*(monitored[22] + monitored[23])
    monitored[25] = monitored[18] + monitored[24]

    # Return results
    return monitored
