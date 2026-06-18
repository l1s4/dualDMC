import os
if "KERAS_BACKEND" not in os.environ:
    os.environ["KERAS_BACKEND"] = "torch"

import math
import numpy as np
import seaborn as sns
import pandas as pd
import keras
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import truncnorm
from numba import njit
import bayesflow as bf


@njit
def ddmc_trial(t0, muc, b, A1, A2, tau1, tau2, sigma = 4.0, dt = 1.0, max_time = 7000): 
    """
    Simulate a single DDMC trial with diffusion constant sigma = 4.0. 
    Parameters: 
        - t0: non-decision time
        - muc: drift rate controlled process
        - b: upper boundary (lower boundary: -b)
        - A1, A2: amplitudes of automatic processes
        - tau1, tau2: timepoint of max automatic activation
    Return: list with result
        - choice (correct = 1.0, incorrect = 0.0, no decision = -1.0)
        - rt (ms)
    """

    dX = np.random.beta(3, 3) * (2 * b) - b    # initial position
    const = sigma * np.sqrt(dt)
    t = dt

    # Loop through process and check boundary conditions
    while (b > dX > -b) and t <= max_time:
        
        noise = np.random.randn()
        # with a2 = a2 = 2
        mu_a1 = A1 / tau1 * np.exp(1 - t / tau1) * (1 - t / tau1)
        mu_a2 = A2 / tau2 * np.exp(1 - t / tau2) * (1 - t / tau2)
        mu_t = mu_a1 + mu_a2 + muc
        
        dX += mu_t * dt + const * noise

        t += dt

    rt = t + t0     # add ndt
    
    # Get decision
    if dX >= b:     # correct
        c = 1.0
    elif dX <= -b:  # incorrect
        c = 0.0
    else:           # no decision within max_time
        c = -1.0
        rt = -1.0

    rt = rt / 1000  # convert to seconds

    return [c, rt]

def simulate_ddmc(t0, muc, b, A1, A2, tau1, tau2, n_subjects=1, n_trials=1):
    """
    Simulate a DDMC experiment for n_subjects with n_trials each
    Results are contaminated
    Parameters: 
        - t0: non-decision time
        - muc: drift rate controlled process
        - b: upper boundary (lower boundary: -b)
        - A1, A2: amplitudes of automatic processes
        - tau1, tau2: timepoint of max automatic activation
    Return: dict containing 'sim_data', an array containing the results
        - [choice, rt, condition]
    """

    if isinstance(muc, (float, int)):
        t0 = np.ones((n_subjects, )) * t0
        muc = np.ones((n_subjects, )) * muc
        b = np.ones((n_subjects, )) * b
        A1 = np.ones((n_subjects, )) * A1
        A2 = np.ones((n_subjects, )) * A2
        tau1 = np.ones((n_subjects, )) * tau1
        tau2 = np.ones((n_subjects, )) * tau2

    data = np.zeros((n_subjects, n_trials, 3))

    for i in range(n_subjects):
        quarter = int(np.ceil(n_trials / 4)) 
        conditions = np.repeat(np.arange(4), quarter)[:n_trials]
        for j in range(n_trials):
            auto1 = 1 if (conditions[j] == 0 or conditions[j] == 1) else -1 # CC, CI
            auto2 = 1 if (conditions[j] == 0 or conditions[j] == 2) else -1 # CC, IC
            res = ddmc_trial(
                t0[i], muc[i], b[i], auto1*A1[i], auto2*A2[i], tau1[i], tau2[i]
            )
            res.append(conditions[j])
            data[i, j] = res

    
    # contaminate
    cont_rt = np.abs(
        np.random.standard_t(df = 1, size = n_subjects * n_trials)
    ).reshape(n_subjects, n_trials)
    cont_resp = np.random.binomial(
        n = 1, p = 0.5, size = n_subjects * n_trials
    ).reshape(n_subjects, n_trials)

    cont_prob = 0.1
    replace = np.random.binomial(
        n = 1, p = cont_prob, size = n_subjects * n_trials
    ).reshape(n_subjects, n_trials)

    data[:, :, 1] = (1-replace)*data[:, :, 1] + [replace*cont_rt]
    data[:, :, 0] = (1-replace)*data[:, :, 0] + [replace*cont_resp]

    if n_subjects == 1 and n_trials == 1:
        data = data[0, 0]
    elif n_subjects == 1:
        data = data[0]
    elif n_trials == 1:
        data = data[:, 0]

    return dict(sim_data=data)


# Set group-level priors
global_prior = {
    'mu_t0'     : (360.0, 50.0),
    'mu_muc'    : (0.4, 0.1),
    'mu_b'      : (55, 15),
    'mu_A1'     : (20, 8),
    'mu_A2'     : (20, 8),
    'mu_tau1' : ([30, 180], 15),
    'mu_tau2' : ([30, 180], 15),

    'log_sigma_t0'  : (3.5, 0.3),
    'log_sigma_muc' : (-1.7, 0.5),
    'log_sigma_b'   : (2.3, 0.4),
    'log_sigma_A1'  : (1.6, 0.5),
    'log_sigma_A2'  : (1.6, 0.5),
    'log_sigma_tau1': (3.0, 0.5),
    'log_sigma_tau2': (3.0, 0.5)
}

def sample_hierarchical_priors(rng=None, n_subjects = 1):
    """
    Sample priors
    Group level: t0, muc, b, A1, A2: normal distribution; 
                 tau1, tau2: mixture distribution of two normal distributions
    Subject level: normal distributions for all parameters
    """
    if rng == None: 
        rng = np.random.default_rng()

    ax = 0
    bx = np.inf

    # Group level
    mu_t0   = truncnorm.rvs(
        (ax - global_prior['mu_t0'][0]) / global_prior['mu_t0'][1], 
        (bx - global_prior['mu_t0'][0]) / global_prior['mu_t0'][1], 
        global_prior['mu_t0'][0], global_prior['mu_t0'][1], random_state = rng
    )
    mu_muc  = truncnorm.rvs(
        (ax - global_prior['mu_muc'][0]) / global_prior['mu_muc'][1], 
        (bx - global_prior['mu_muc'][0]) / global_prior['mu_muc'][1], 
        global_prior['mu_muc'][0],  global_prior['mu_muc'][1], random_state = rng
    )
    mu_b    = truncnorm.rvs(
        (ax - global_prior['mu_b'][0]) /  global_prior['mu_b'][1], 
        (bx - global_prior['mu_b'][0]) /  global_prior['mu_b'][1], 
        global_prior['mu_b'][0],  global_prior['mu_b'][1], random_state = rng
    )
    mu_A1   = truncnorm.rvs(
        (ax - global_prior['mu_A1'][0]) / global_prior['mu_A1'][1], 
        (bx - global_prior['mu_A1'][0]) / global_prior['mu_A1'][1], 
        global_prior['mu_A1'][0], global_prior['mu_A1'][1], random_state = rng
    )
    mu_A2   = truncnorm.rvs(
        (ax - global_prior['mu_A2'][0]) / global_prior['mu_A2'][1], 
        (bx - global_prior['mu_A2'][0]) / global_prior['mu_A2'][1], 
        global_prior['mu_A2'][0], global_prior['mu_A2'][1], random_state = rng
    )
    
    # Mixture distribution for tau1, tau2
    c1 = rng.binomial(1, 0.5)
    c2 = rng.binomial(1, 0.5)
    if c1 == 0: 
        mu_tau1 = truncnorm.rvs(
            (ax - global_prior['mu_tau1'][0][0]) / global_prior['mu_tau1'][1], 
            (bx - global_prior['mu_tau1'][0][0]) / global_prior['mu_tau1'][1], 
            global_prior['mu_tau1'][0][0], global_prior['mu_tau1'][1], 
            random_state = rng
        )
    else: 
        mu_tau1 = truncnorm.rvs(
            (ax - global_prior['mu_tau1'][0][1]) / global_prior['mu_tau1'][1], 
            (bx - global_prior['mu_tau1'][0][1]) / global_prior['mu_tau1'][1], 
            global_prior['mu_tau1'][0][1], global_prior['mu_tau1'][1], 
            random_state = rng
        )
       
    if c2 == 0: 
        mu_tau2 = truncnorm.rvs(
            (ax - global_prior['mu_tau2'][0][0]) / global_prior['mu_tau2'][1], 
            (bx - global_prior['mu_tau2'][0][0]) / global_prior['mu_tau2'][1], 
            global_prior['mu_tau2'][0][0], global_prior['mu_tau2'][1], 
            random_state = rng
        )
    else: 
        mu_tau2 = truncnorm.rvs(
            (ax - global_prior['mu_tau2'][0][1]) / global_prior['mu_tau2'][1], 
            (bx - global_prior['mu_tau2'][0][1]) / global_prior['mu_tau2'][1], 
            global_prior['mu_tau2'][0][1], global_prior['mu_tau2'][1], 
            random_state = rng
        )

    log_sigma_t0 = rng.normal(*global_prior['log_sigma_t0'])
    log_sigma_b = rng.normal(*global_prior['log_sigma_b'])
    log_sigma_muc = rng.normal(*global_prior['log_sigma_muc'])
    log_sigma_A1 = rng.normal(*global_prior['log_sigma_A1'])
    log_sigma_A2 = rng.normal(*global_prior['log_sigma_A2'])
    log_sigma_tau1 = rng.normal(*global_prior['log_sigma_tau1'])
    log_sigma_tau2 = rng.normal(*global_prior['log_sigma_tau2'])

    # Subject level
    t0 = truncnorm.rvs(
        (ax - mu_t0) / np.exp(log_sigma_t0), (bx - mu_t0) / np.exp(log_sigma_t0), 
        mu_t0, np.exp(log_sigma_t0), size=n_subjects, random_state=rng
    )
    muc = truncnorm.rvs(
        (ax - mu_muc) / np.exp(log_sigma_muc), (bx - mu_muc) / np.exp(log_sigma_muc), 
        mu_muc, np.exp(log_sigma_muc), size=n_subjects, random_state = rng
    )
    b = truncnorm.rvs(
        (ax - mu_b) / np.exp(log_sigma_b), (bx - mu_b) / np.exp(log_sigma_b), 
        mu_b, np.exp(log_sigma_b), size=n_subjects, random_state = rng
    )
    A1 = truncnorm.rvs(
        (ax - mu_A1) / np.exp(log_sigma_A1), (bx - mu_A1) / np.exp(log_sigma_A1), 
        mu_A1, np.exp(log_sigma_A1), size=n_subjects, random_state = rng
    )
    A2 = truncnorm.rvs(
        (ax - mu_A2) / np.exp(log_sigma_A2), (bx - mu_A2) / np.exp(log_sigma_A2), 
        mu_A2, np.exp(log_sigma_A2), size=n_subjects, random_state = rng
    )
    tau1 = truncnorm.rvs(
        (ax - mu_tau1) / np.exp(log_sigma_tau1), 
        (bx - mu_tau1) / np.exp(log_sigma_tau1), 
        mu_tau1, np.exp(log_sigma_tau1), size=n_subjects, random_state = rng
    )
    tau2 = truncnorm.rvs(
        (ax - mu_tau2) / np.exp(log_sigma_tau2), 
        (bx - mu_tau2) / np.exp(log_sigma_tau2), 
        mu_tau2, np.exp(log_sigma_tau2), size=n_subjects, random_state = rng
    )

    return {
        # group 
        'mu_t0': mu_t0, 'mu_muc': mu_muc, 'mu_b': mu_b, 'mu_A1': mu_A1, 
        'mu_A2': mu_A2, 'mu_tau1': mu_tau1, 'mu_tau2' : mu_tau2, 
        'log_sigma_t0': log_sigma_t0, 'log_sigma_muc': log_sigma_muc, 
        'log_sigma_b': log_sigma_b, 'log_sigma_A1': log_sigma_A1, 
        'log_sigma_A2'  : log_sigma_A2, 'log_sigma_tau1': log_sigma_tau1, 
        'log_sigma_tau2': log_sigma_tau2, 

        # subject
        't0': t0, 'muc': muc, 'b': b, 'A1': A1, 'A2': A2, 'tau1': tau1, 
        'tau2': tau2
    }

def score_log_norm(x, m, s):
    return -(x-m) / s**2

def score_log_mixture(x, m1, s1, m2, s2):
    n1 = np.exp(-0.5 * ((x - m1) / s1) ** 2) / (np.sqrt(2*np.pi) * s1)
    n2 = np.exp(-0.5 * ((x - m2) / s2) ** 2) / (np.sqrt(2*np.pi) * s2)

    p = 0.5 * n1 + 0.5 * n2

    return (
        0.5 * n1 * (-(x - m1) / s1**2)
        + 0.5 * n2 * (-(x - m2) / s2**2)
    ) / p

def prior_global_score(x: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    mu_t0   = x["mu_t0"]
    mu_muc  = x["mu_muc"]
    mu_b    = x["mu_b"]
    mu_A1   = x["mu_A1"]
    mu_A2   = x["mu_A2"]
    mu_tau1 = x["mu_tau1"]
    mu_tau2 = x["mu_tau2"]
    log_sigma_t0    = x["log_sigma_t0"]
    log_sigma_muc   = x["log_sigma_muc"]
    log_sigma_b     = x["log_sigma_b"]
    log_sigma_A1    = x["log_sigma_A1"]
    log_sigma_A2    = x["log_sigma_A2"]
    log_sigma_tau1  = x["log_sigma_tau1"]
    log_sigma_tau2  = x["log_sigma_tau2"]

    parts = {
        "mu_t0": score_log_norm(
            mu_t0, m=global_prior["mu_t0"][0], s=global_prior["mu_t0"][1]
        ),
        "mu_muc": score_log_norm(
            mu_muc, m=global_prior["mu_muc"][0], s=global_prior["mu_muc"][1]
        ),
        "mu_b": score_log_norm(
            mu_b, m=global_prior["mu_b"][0], s=global_prior["mu_b"][1]
        ),
        "mu_A1": score_log_norm(
            mu_A1, m=global_prior["mu_A1"][0], s=global_prior["mu_A1"][1]
        ),
        "mu_A2": score_log_norm(
            mu_A2, m=global_prior["mu_A2"][0], s=global_prior["mu_A2"][1]
        ),
        "mu_tau1": score_log_mixture(
            mu_tau1, 
            m1=global_prior["mu_tau1"][0][0],m2=global_prior["mu_tau1"][0][1], 
            s1=global_prior["mu_tau1"][1], s2=global_prior["mu_tau1"][1]
        ),
        "mu_tau2": score_log_mixture(
            mu_tau2, 
            m1=global_prior["mu_tau2"][0][0],m2=global_prior["mu_tau2"][0][1], 
            s1=global_prior["mu_tau2"][1], s2=global_prior["mu_tau2"][1]
        ),

        "log_sigma_t0": score_log_norm(
            log_sigma_t0, 
            m=global_prior["log_sigma_t0"][0], 
            s=global_prior["log_sigma_t0"][1]
        ),
        "log_sigma_muc": score_log_norm(
            log_sigma_muc, 
            m=global_prior["log_sigma_muc"][0], 
            s=global_prior["log_sigma_muc"][1]
        ),
        "log_sigma_b": score_log_norm(
            log_sigma_b, 
            m=global_prior["log_sigma_b"][0], s=global_prior["log_sigma_b"][1]
        ),
        "log_sigma_A1": score_log_norm(
            log_sigma_A1, 
            m=global_prior["log_sigma_A1"][0], s=global_prior["log_sigma_A1"][1]
        ),
        "log_sigma_A2": score_log_norm(
            log_sigma_A2, 
            m=global_prior["log_sigma_A2"][0], s=global_prior["log_sigma_A2"][1]
        ),
        "log_sigma_tau1": score_log_norm(
            log_sigma_tau1, 
            m=global_prior["log_sigma_tau1"][0], 
            s=global_prior["log_sigma_tau1"][1]
        ),
        "log_sigma_tau2": score_log_norm(
            log_sigma_tau2, 
            m=global_prior["log_sigma_tau2"][0], 
            s=global_prior["log_sigma_tau2"][1]
        ),
    }
    return parts


simulator_hierarchical = bf.make_simulator([sample_hierarchical_priors, simulate_ddmc])


param_names_global = list(global_prior.keys())
param_names_local = ['t0', 'muc', 'b', 'A1', 'A2', 'tau1', 'tau2']

# Global (group-level) workflow
adapter_global = (
    bf.adapters.Adapter()
    .to_array()
    .convert_dtype("float64", "float32")
    .concatenate(param_names_global, into="inference_variables")
    .rename("sim_data", "summary_variables")
)

workflow_global = bf.CompositionalWorkflow(
    adapter=adapter_global,
    simulator=simulator_hierarchical,
    initial_learning_rate=5e-4,
    summary_network=bf.networks.DeepSet(
        dropout=0.01070354852467715,
        summary_dim=24             # 4 * dim of inference variables (4 * 1 * 8)
    ),
    inference_network=bf.networks.DiffusionModel(
        dropout=0.01070354852467715
    ),
    checkpoint_filepath = Path(os.getcwd()).resolve(),  # save in cwd
    checkpoint_name = "HDDM_GLOBAL"        # file name
)


# Local (subject-level) workflow
adapter_local = (
    bf.adapters.Adapter()
    .to_array()
    .convert_dtype("float64", "float32")
    .concatenate(param_names_local, into="inference_variables")
    .concatenate(param_names_global, into="inference_conditions")
    .rename("sim_data", "summary_variables")
)

workflow_local = bf.BasicWorkflow(
    adapter=adapter_local,
    initial_learning_rate=5e-4,
    summary_network=bf.networks.SetTransformer(
        dropout=0.01070354852467715, 
        summary_dim=24             # 4 * dim of inference variables (4 * 1 * 8)
    ),
    inference_network=bf.networks.StableConsistencyModel(
        dropout=0.01070354852467715
    ), 
    checkpoint_filepath = Path(os.getcwd()).resolve(),  # save in cwd
    checkpoint_name = "HDDM_LOCAL"         # file name
)

# Configuration
N_TRAINING_BATCHES = 256
BATCH_SIZE = 64
EPOCHS = 500
N_LOCAL_SUBJECTS = 6    # number of subjects in a single training batch
N_TRIALS = 200          # trials per subject
N_TEST = 100            # number of test data sets
N_SAMPLES = 100         # posterior samples per test dataset




training_data = simulator_hierarchical.sample_parallel(
    (N_TRAINING_BATCHES * BATCH_SIZE),
    n_subjects=N_LOCAL_SUBJECTS, n_trials=N_TRIALS
)

training_data['sim_data'] = training_data['sim_data'].reshape(
    (N_TRAINING_BATCHES * BATCH_SIZE, N_LOCAL_SUBJECTS * N_TRIALS, 3)
)

history = workflow_global.fit_offline(
    training_data,
    epochs=EPOCHS,
    batch_size=BATCH_SIZE
)
# save loss figure
fig_loss = bf.diagnostics.loss(history=history)
fig_loss.savefig('loss_global.pdf')


training_data = simulator_hierarchical.sample_parallel(
    (N_TRAINING_BATCHES * BATCH_SIZE), n_trials=N_TRIALS
)

history = workflow_local.fit_offline(
    training_data,
    epochs=EPOCHS, 
    batch_size=BATCH_SIZE,
)
# save loss figure
fig_loss = bf.diagnostics.loss(history=history)
fig_loss.savefig('loss_local.pdf')
