# note: script is based on https://github.com/simschaefer/amortized-dmc
# note: a1 = a2 = 2

import os
import math
import numpy as np  
import pandas as pd
#import keras
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import truncnorm
if "KERAS_BACKEND" not in os.environ:
    os.environ["KERAS_BACKEND"] = "torch"

import bayesflow as bf


# Priors
def ddmc_prior(rng = None): 
    if rng is None:
        rng = np.random.default_rng()
    
    mu_r    = truncnorm.rvs((0 - 400) / 30, (np.inf - 400) / 30, 400, 30, random_state = rng)
    sd_r    = truncnorm.rvs((0 - 30) / 10, (np.inf - 30) / 10, 30, 10, random_state = rng)
    b       = truncnorm.rvs((0 - 80) / 20, (np.inf - 80) / 20, 80, 20, random_state = rng)
    muc     = rng.beta(2, 2)
    A1      = truncnorm.rvs((0 - 20) / 8, (np.inf - 20) / 8, 20, 8, random_state = rng)
    A2      = truncnorm.rvs((0 - 20) / 8, (np.inf - 20) / 8, 20, 8, random_state = rng)
    tau1    = rng.uniform(20, 180)
    tau2    = rng.uniform(20, 180)
    

    return dict(mu_r=mu_r, sd_r=sd_r, b=b, muc=muc, A1=A1, A2=A2, tau1=tau1, 
                tau2=tau2)

# Simulate DDMC trials
def ddmc_trial(muc, A1, A2, tau1, tau2, b, ndts, noise, t, sigma = 4.0, dt = 1): 

    # prepare output
    num_trials, _ = noise.shape
    rts = np.full(num_trials, -1.0)
    resps = np.full(num_trials, -1)

    # initial position for all trials
    X0 = np.random.beta(3, 3, size=num_trials) * (2 * b) - b

    # with a1 = a2 = 2
    mu_a1 = A1 / tau1 * np.exp(1 - t / tau1) * (1 - t / tau1)
    mu_a2 = A2 / tau2 * np.exp(1 - t / tau2) * (1 - t / tau2)
    mu_t = mu_a1 + mu_a2 + muc
    dX = mu_t[None, :] * dt + sigma * np.sqrt(dt) * noise

    X = np.cumsum(dX, axis=1) + X0[:, None]

    crossed_upper = X >= b
    crossed_lower = X <= -b
    crossed_any = crossed_upper | crossed_lower

    # first crossing index for each trial
    first_crossing = np.argmax(crossed_any, axis=1)
    has_crossed = np.any(crossed_any, axis=1)

    # fill only for trials that crossed
    idx = np.where(has_crossed)[0]
    crossing_times = t[first_crossing[idx]]


    # use nondecision times only for trials that crossed
    ndts_crossed = ndts[idx]
    rts[idx] = (crossing_times + ndts_crossed) / 1000  # convert to seconds

    # determine response type
    resp_hit = X[idx, first_crossing[idx]]
    resps[idx] = (resp_hit >= b).astype(int)

    return np.c_[rts, resps]


# Run DDMC experiment (4 conditions: c-c, c-i, i-c, i-i)
def ddmc_experiment(num_obs, muc, A1, A2, tau1, tau2, b, mu_r, sd_r): 
    max_time = 1500
    dt = 1

    # precompute time vector and noise
    t = np.arange(start=dt, stop=max_time+dt, step=dt)
    noise = np.random.normal(size = (num_obs, len(t)))
#    ndts = np.random.normal(size = num_obs, loc=mu_r, scale=sd_r)
    ndts = truncnorm.rvs((0 - mu_r) / sd_r, (np.inf - mu_r) / sd_r, mu_r, sd_r, size=num_obs)

    out = np.zeros((num_obs, 2))        # to store rt and resp

    # congruency conditions (equal split)
    quarter = int(np.ceil(num_obs / 4))
    conditions = np.repeat(np.arange(4), quarter)[:num_obs]

    # simulate CONG-CONG trials (A1, A2)
    out[:quarter] = ddmc_trial(
        muc=muc, A1=A1, A2=A2, tau1=tau1, tau2=tau2, b=b, t=t, 
        ndts=ndts[:quarter], noise=noise[:quarter]
    )
    # simulate CONG-INCONG trials (A1, -A2)
    out[quarter:quarter*2] = ddmc_trial(
        muc=muc, A1=A1, A2=-A2, tau1=tau1, tau2=tau2, b=b, t=t, 
        ndts=ndts[quarter:quarter*2], noise=noise[quarter:quarter*2]
    )
    # simulate INCONG-CONG trials (-A1, A2)
    out[quarter*2:quarter*3] = ddmc_trial(
        muc=muc, A1=-A1, A2=A2, tau1=tau1, tau2=tau2, b=b, t=t, 
        ndts=ndts[quarter*2:quarter*3], noise=noise[quarter*2:quarter*3]
    )
    # simulate INCONG-INCONG trials (-A1, -A2)
    out[quarter*3:] = ddmc_trial(
        muc=muc, A1=-A1, A2=-A2, tau1=tau1, tau2=tau2, b=b, t=t, 
        ndts=ndts[quarter*3:], noise=noise[quarter*3:]
    )
    
    # contaminate
    cont_rt = np.abs(np.random.standard_t(df = 1, size = num_obs))
    cont_resp = np.random.binomial(n = 1, p = 0.5, size = num_obs)

    cont_prob = 0.1
    replace = np.random.binomial(n = 1, p = cont_prob, size = num_obs)

    out[:, 0] = (1-replace)*out[:, 0] + replace*cont_rt
    out[:, 1] = (1-replace)*out[:, 1] + replace*cont_resp


    return dict(rt = out[:, 0], resp = out[:, 1], conditions = conditions, 
                num_obs = num_obs)

#ddmc_experiment(num_obs = 50, muc = 0.5, A1 = 60, A2 = 60, tau1 = 10, 
#                tau2 = 170, b = 50, mu_r = 300, sd_r = 36)


def meta(batch_size, num_obs = None): 
    if num_obs == None:
        num_obs = np.random.randint(50, 1000)
    return dict(num_obs = num_obs)

# define simulator
simulator = bf.simulators.make_simulator([ddmc_prior, ddmc_experiment], meta_fn = meta)

# define adapter
adapter = (bf.Adapter()
    .convert_dtype("float64", "float32")
    .sqrt("num_obs")
    .as_set(["rt", "resp", "conditions"])
    .concatenate(
        ["muc", "A1", "A2", "tau1", "tau2", "mu_r", "sd_r", "b"], 
        into = "inference_variables"
    )
    .concatenate(["rt", "resp", "conditions"], into = "summary_variables")
    .standardize(
        include="inference_variables", 
        mean=[0.5, 20, 20, 100, 100, 400, 30, 80],
        std=[0.22, 8, 8, 46, 46, 30, 10, 20]
    )
    .rename("num_obs", "inference_conditions")
)

# define workflow
workflow = bf.BasicWorkflow(
    simulator=simulator,
    adapter=adapter,
    initial_learning_rate=5e-4,
    inference_network=bf.networks.FlowMatching(
        dropout=0.01070354852467715
    ),
    summary_network=bf.networks.SetTransformer(
        dropout=0.01070354852467715, 
        num_seeds=7, 
        summary_dim=24,
        embed_dims=(128, 128)
    ),
    inference_variables = ["muc", "A1", "A2", "tau1", "tau2", "mu_r", "sd_r", "b"],
    inference_conditions = ["num_obs"], 
    summary_variables = ["rt", "resp", "conditions"],
    checkpoint_filepath = Path(os.getcwd()).resolve(),  # save in cwd
    checkpoint_name = "DDMC_1500RU"                 # file name
)


## Fit offline
#train_data = simulator.sample(5_00)
#validation_data = simulator.sample(1_00)
#
#history = workflow.fit_offline(
#    data=train_data,
#    epochs=5,
#    batch_size=32,
#    #num_batches_per_epoch=100,
#    validation_data=validation_data
#)

# Fit online
val_data = simulator.sample(200)
history = workflow.fit_online(
    epochs=250,
    num_batches_per_epoch=250,
    batch_size=64, 
    val_data=val_data
)


# Plot loss
f = bf.diagnostics.plots.loss(history = history)
plt.savefig("loss.pdf")

# Plot default diagnostics
figures = workflow.plot_default_diagnostics(test_data=val_data)
for k,i in figures.items():
    figures[k].savefig(k + '_posttraining.pdf')
