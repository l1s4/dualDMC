# note: this script is based on https://github.com/simschaefer/amortized-dmc
#Schaefer, S. B., Radev, S. T., Göttmann, J., & Schubert, A. -L. (2026). Amortized bayesian
#workflow for modeling congruency effects using the diffusion model for conflict
#tasks. https://doi.org/https://doi.org/10.31234/osf.io/dypcw_v2

import os
import math
import numpy as np  
import numba as nb
import seaborn as sns
import pandas as pd
#import keras
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import truncnorm
if "KERAS_BACKEND" not in os.environ:                # ensure the backend is set
    os.environ["KERAS_BACKEND"] = "torch"

import bayesflow as bf


# Priors
def mdmc_prior(): 
    lower = 0
    upper = np.inf

    means = np.array([360, 35, 60, 0.6, 20, 20, 120, 120])
    sds = np.array([25, 8, 13, 0.15, 8, 8, 40, 40])

    a = (lower - means) / sds
    b = (upper - means) / sds

    priors = truncnorm.rvs(a, b, loc=means, scale=sds)

    return dict(mu_r=priors[0], sd_r=priors[1], b=priors[2], muc=priors[3], 
                A1=priors[4], A2=priors[5], tau1=priors[6], tau2=priors[7])

#mdmc_prior()


# Simulate MDMC trials
def mdmc_trial(muc, A1, A2, tau1, tau2, b, ndts, noise, t, sigma = 4.0, dt = 1): 

    # Prepare output
    num_trials, _ = noise.shape
    rts = np.full(num_trials, -1.0)
    resps = np.full(num_trials, -1)

    # initial position for all trials
#    X0 = 0.0      # no starting point variability
    X0 = np.random.beta(3, 3, size=num_trials) * (2 * b) - b

    # with a = 2
    mu_t = (A1 / tau1 * np.exp(1 - t / tau1) * (1 - t / tau1)) +  (A2 / tau2 * np.exp(1 - t / tau2) * (1 - t / tau2)) + muc
    dX = mu_t[None, :] * dt + sigma * np.sqrt(dt) * noise

    X = np.cumsum(dX, axis=1) + X0[:, None]

    crossed_upper = X >= b
    crossed_lower = X <= -b
    crossed_any = crossed_upper | crossed_lower

    # First crossing index for each trial
    first_crossing = np.argmax(crossed_any, axis=1)
    has_crossed = np.any(crossed_any, axis=1)

    # Fill only for trials that crossed
    idx = np.where(has_crossed)[0]
    crossing_times = t[first_crossing[idx]]


    # use nondecision times only for trials that crossed
    ndts_crossed = ndts[idx]
#    rts[idx] = np.log((crossing_times + ndts_crossed) / 1000)  # convert to seconds and log transform for faster convergence
    rts[idx] = (crossing_times + ndts_crossed) / 1000  # convert to seconds

    # Determine response type
    resp_hit = X[idx, first_crossing[idx]]
    resps[idx] = (resp_hit >= b).astype(int)

    return np.c_[rts, resps]#, X, t

num_obs = 2
max_time = 1500
dt = 1

t = np.linspace(start=dt, stop=max_time, num=int(max_time/dt))  # time [ms]
noise = np.random.normal(size = (num_obs, len(t)))
ndts = np.random.normal(size = num_obs, loc=300, scale=0)

mdmc_trial(muc=0.5, A1=20, A2=20, tau1=30, tau2=30, b=50, t=t, 
           ndts=ndts, noise = noise)
           
#reacts, traj, t = mdmc_trial(muc=0.5, A1=20, A2=20, tau1=30, tau2=30, b=75, t=t, 
#            ndts=ndts, noise = noise)
#plt.plot(t / 1000, traj[0])        # time in seconds


# Run MDMC experiment (4 conditions: c-c, c-i, i-c, i-i)
def mdmc_experiment(num_obs, muc, A1, A2, tau1, tau2, b, mu_r, sd_r): 
    max_time = 1500
    dt = 1

    # precompute time vector and noise
    t = np.linspace(start=dt, stop=max_time, num=int(max_time/dt))
    noise = np.random.normal(size = (num_obs, len(t)))
    ndts = np.random.normal(size = num_obs, loc=mu_r, scale=sd_r)

    out = np.zeros((num_obs, 2))        # to store rt and resp

    # congruency conditions (equal split)
    quarter = int(np.ceil(num_obs / 4)) + 1
    conditions = np.repeat(np.arange(4), quarter)[:num_obs]

    # simulate CONG-CONG trials (A1, A2)
    out[:quarter] = mdmc_trial(
        muc=muc, A1=A1, A2=A2, tau1=tau1, tau2=tau2, b=b, t=t, ndts=ndts[:quarter], noise=noise[:quarter]
    )
    # simulate CONG-INCONG trials (A1, -A2)
    out[quarter:quarter*2] = mdmc_trial(
        muc=muc, A1=A1, A2=-A2, tau1=tau1, tau2=tau2, b=b, t=t, ndts=ndts[quarter:quarter*2], noise=noise[quarter:quarter*2]
    )
    # simulate INCONG-CONG trials (-A1, A2)
    out[quarter*2:quarter*3] = mdmc_trial(
        muc=muc, A1=-A1, A2=A2, tau1=tau1, tau2=tau2, b=b, t=t, ndts=ndts[quarter*2:quarter*3], noise=noise[quarter*2:quarter*3]
    )
    # simulate INCONG-INCONG trials (-A1, -A2)
    out[quarter*3:] = mdmc_trial(
        muc=muc, A1=-A1, A2=-A2, tau1=tau1, tau2=tau2, b=b, t=t, ndts=ndts[quarter*3:], noise=noise[quarter*3:]
    )
    

    return dict(rt = out[:, 0], resp = out[:, 1], conditions = conditions, num_obs = num_obs)

#mdmc_experiment(num_obs = 50, muc = 0.5, A1 = 60, A2 = 60, tau1 = 10, tau2 = 170, b = 50, mu_r = 300, sd_r = 36)


def meta(batch_size, num_obs = None): 
    if num_obs == None:
        num_obs = np.random.randint(50, 1000)
    return dict(num_obs = num_obs)

# define simulator
simulator = bf.simulators.make_simulator([mdmc_prior, mdmc_experiment], meta_fn = meta)

# define adapter
adapter = (bf.Adapter()
    .convert_dtype("float64", "float32")
    .sqrt("num_obs")
    .as_set(["rt", "resp", "conditions"])
    .concatenate(["muc", "A1", "A2", "tau1", "tau2", "mu_r", "sd_r", "b"], into = "inference_variables")
    .concatenate(["rt", "resp", "conditions"], into = "summary_variables")
    .standardize(include="inference_variables", mean=[0.6, 20, 20, 120, 120, 360, 35, 60], std=[0.15, 8, 8, 40, 40, 360, 35, 60])
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
        summary_dim=22,
        embed_dims=(128, 128)
    ),
    inference_variables = ["muc", "A1", "A2", "tau1", "tau2", "mu_r", "sd_r", "b"],
    inference_conditions = ["num_obs"], 
    summary_variables = ["rt", "resp", "conditions"],
    checkpoint_filepath = Path(os.getcwd()).resolve(),  # save in cwd
    checkpoint_name = "MDMC_2603_fitoff"                 # file name
)


# To fit offline
train_data = simulator.sample(5_00)
validation_data = simulator.sample(1_00)

history = workflow.fit_offline(
    data=train_data,
    epochs=5,
    batch_size=32,
    #num_batches_per_epoch=100,
    validation_data=validation_data
)

## To fit online
##TODO: check paper for specs
#val_data = simulator.sample(200)
#history = workflow.fit_online(
#    epochs=200,
#    num_batches_per_epoch=250,
#    batch_size=64,
#    validation_data=val_data
#)




# Plot loss
f = bf.diagnostics.plots.loss(history = history)
plt.savefig("loss.pdf")


val_data = simulator.sample(200)
# Plot default diagnostics
figures = workflow.plot_default_diagnostics(test_data=val_data)
for k,i in figures.items():
    figures[k].savefig(k + '_posttraining.pdf')
