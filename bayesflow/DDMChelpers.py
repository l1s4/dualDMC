import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def ddmc_trial(muc, A1, A2, tau1, tau2, b, ndts, noise, t, sigma = 4.0, dt = 1): 
    # Simulate multiple DDMC trials in parallel

    # Prepare output
    num_trials, _ = noise.shape
    rts = np.full(num_trials, -1.0)
    resps = np.full(num_trials, -1)

    X0 = np.random.beta(3, 3, size=num_trials) * (2 * b) - b    # initial pos

    # with a1 = a2 = 2
    mu_a1 = A1 / tau1 * np.exp(1 - t / tau1) * (1 - t / tau1)
    mu_a2 = A2 / tau2 * np.exp(1 - t / tau2) * (1 - t / tau2)
    mu_t = mu_a1 + mu_a2 + muc

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


    ndts_crossed = ndts[idx]    # nondecision times only for trials that crossed
    rts[idx] = (crossing_times + ndts_crossed) / 1000  # convert to seconds

    # get response type
    resp_hit = X[idx, first_crossing[idx]]
    resps[idx] = (resp_hit >= b).astype(int)

    return np.c_[rts, resps]


def ddmc_experiment(num_obs, muc, A1, A2, tau1, tau2, b, mu_r, sd_r): 
    max_time = 2000
    dt = 1

    # precompute time vector and noise
    t = np.linspace(start=dt, stop=max_time, num=int(max_time/dt))
    noise = np.random.normal(size = (num_obs, len(t)))
    ndts = np.random.normal(size = num_obs, loc=mu_r, scale=sd_r)

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

    return dict(rt = out[:, 0], resp = out[:, 1], conditions = conditions, 
    num_obs = num_obs)


def read_ReyMermetGade(directory, filelist, colnames, first, second):
    """
    Read in combi data from: Rey-Mermet, A., & Gade, M. (2016). Contextual 
    within-trial adaptation of cognitive control: Evidence from the combination 
    of conflict tasks. Journal of Experimental Psychology: Human Perception and 
    Performance, 42(10), 1505–1532. https://doi.org/10.1037/xhp0000229
    https://osf.io/ptg4n/wiki?wiki=r8bfj

    directory: path to combi files
    filelist: os.listdir(directory)
    colnames: colnames for data, see osf
    first: first conflict task, must match colname
    second: second conflict task, must match colname

    Returns: dataframe including congruency condition (CC / CI / IC / II)
    """
    df = pd.concat(
        [pd.read_csv(os.path.join(directory, f), sep=r"\s+", header=None)
        for f in filelist], 
        ignore_index=True
    )

    df.columns = colnames

    # drop unused columns
    #df = df.drop(df[df.catch == 1].index)
    #df.drop(['block', 'catch', 're', 'xr'], inplace=True, axis=1)

    df["rt"] = df["rt"] / 1000.0      # convert rt to s

    df["CI1"] = np.where(df[first] == 1, "congruent", "incongruent")
    df["CI2"] = np.where(df[second] == 1, "congruent", "incongruent")
    
    # condition column
    df["CI"] = np.where(
        (df["CI1"] == "congruent") & (df["CI2"] == "congruent"), "CC",
        np.where(
            (df["CI1"] == "incongruent") & (df["CI2"] == "congruent"), "IC",
            np.where(
                (df["CI1"] == "congruent") & (df["CI2"] == "incongruent"), "CI", 
                "II"
            )
        )
    )
    
    df["conditions"] = df["CI"].astype("category").cat.codes

    return df

def pred_samples_participants(data, approximator, n_resims=200, num_obs=200):
    """
    Iterate over participants, sample using approximator and predict
    Returns: tuple of combined df with samples and combined df with predictions
    """
    participants = data["subject"].unique()
    
    all_resims_list = []
    part_samples_list = []
    
    # iterate participants
    for i, subj in enumerate(participants):
        subj_data = data[data["subject"] == subj]

        subj_inference_dict = {
            key: np.array([subj_data[key].values.reshape(len(subj_data), 1)
        ])
        for key in ['rt', 'resp', 'conditions']
        }
        subj_inference_dict["num_obs"] = np.sum(subj_inference_dict["rt"], axis=1)

        samples = approximator.sample(
            conditions=subj_inference_dict, num_samples=200
        )

        # save samples
        samples_flat = {k: v.flatten() for k, v in samples.items()}
        subj_data_samples = pd.DataFrame(samples_flat)
        subj_data_samples["subject"] = subj
        part_samples_list.append(subj_data_samples)


        # shuffle
        for k, val in samples_flat.items():
            np.random.shuffle(val)
        
        # resimulate data
        resims_lst = []
        for j in range(n_resims): 
            
            iter_dict = {key: val[j] for key, val in samples_flat.items() if val[j] > 0}
            resim = ddmc_experiment(**iter_dict | {'num_obs': num_obs})

            resim_df = pd.DataFrame(resim)
            resim_df["n_resim"] = j
            resim_df["subject"] = subj

            resims_lst.append(pd.DataFrame(resim_df))

        resim_full = pd.concat(resims_lst)
        all_resims_list.append(resim_full)

    all_resims = pd.concat(all_resims_list)
    all_resims = all_resims[all_resims["rt"] != -1.0]   # exclude non-resps
    part_samples = pd.concat(part_samples_list)

    return (part_samples, all_resims)

def pred_samples_conditions(data, approximator, n_resims = 200, num_obs = 200): 
    """
    Iterate over conditions, sample using approximator and predict
    Returns: tuple of combined df with samples and combined df with predictions
    """
    conditions = data["CI"].unique()
    
    all_resims_list = []
    cond_samples_list = []
    
    # iterate participants
    for i, cond in enumerate(conditions):
        cond_data = data[data["CI"] == cond]

        cond_inference_dict = {
            key: np.array([cond_data[key].values.reshape(len(cond_data), 1)
        ])
        for key in ['rt', 'resp', 'conditions']
        }
        cond_inference_dict["num_obs"] = np.sum(cond_inference_dict["rt"], axis=1)

        samples = approximator.sample(
            conditions=cond_inference_dict, num_samples=200
        )

        # save samples
        samples_flat = {k: v.flatten() for k, v in samples.items()}
        cond_data_samples = pd.DataFrame(samples_flat)
        cond_data_samples["CI"] = cond
        cond_samples_list.append(cond_data_samples)

        # shuffle
        for k, val in samples_flat.items():
            np.random.shuffle(val)
        
        # resimulate data
        resims_lst = []
        for j in range(n_resims): 
            iter_dict = {key: val[j] for key, val in samples_flat.items()}
            resim = ddmc_experiment(**iter_dict | {'num_obs': num_obs})

            resim_df = pd.DataFrame(resim)
            resim_df["n_resim"] = j
            resim_df["CI"] = cond
        
            resims_lst.append(pd.DataFrame(resim_df))

        resim_full = pd.concat(resims_lst)
        all_resims_list.append(resim_full)

    all_resims = pd.concat(all_resims_list)
    all_resims = all_resims[all_resims["rt"] != -1.0]   # exclude non-resps
    cond_samples = pd.concat(cond_samples_list)

    return (cond_samples, all_resims)


def plt_subjs_emp (emp_data, pred_data):
    """
    Plot predicted rts vs empirical rts for each subject
    """

    participants = emp_data["subject"].unique()

    f, axarr = plt.subplots(3, 8, figsize=(12, 4), sharex=True, sharey=True)
    for idx, subj in enumerate(participants): 
        ax = axarr.flat[idx]
        sns.kdeplot(
            pred_data[pred_data["subject"] == subj]["rt"], 
            ax=ax, 
            color="blue"
        )
        sns.kdeplot(
            emp_data[emp_data["subject"] == subj]["rt"], 
            ax=ax, 
            color="orange"
        )
        sns.despine(ax=ax)
        ax.set_title(f"Subject {subj}")
        ax.set_xlabel("")
        ax.set_ylabel("")
    f.tight_layout(rect=[0, 0, 1, 1])


def plt_conds_emp (emp_data, pred_data):
    """
    Plot predicted rts vs empirical rts for each condition
    """

    conditions = emp_data["CI"].unique()

    f, axarr = plt.subplots(2, 2, figsize=(12, 4), sharex=True, sharey=True)
    for idx, cond in enumerate(conditions): 
        ax = axarr.flat[idx]
        sns.kdeplot(
            pred_data[pred_data["CI"] == cond]["rt"], 
            ax=ax, 
            color="blue"
        )
        sns.kdeplot(
            emp_data[emp_data["CI"] == cond]["rt"], 
            ax=ax, 
            color="orange"
        )
        sns.despine(ax=ax)
        ax.set_title(f"Condition {cond}")
        ax.set_xlabel("")
        ax.set_ylabel("")
    f.tight_layout(rect=[0, 0, 1, 1])
