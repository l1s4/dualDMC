import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Tuple, Optional, Mapping, Sequence, Union, Dict, List, Any, Iterable, Hashable, Literal


def read_ReyMermetGade(
    directory : str, 
    colnames : list, 
    first : str, 
    second : str, 
    exclude_catch : bool = True, 
    exclude_follow_catch : bool = False, 
    rt_min : int = 0, 
    rt_max : int = 2 
) -> pd.DataFrame: 
    """
    Read in combi data from: Rey-Mermet, A., & Gade, M. (2016). Contextual 
    within-trial adaptation of cognitive control: Evidence from the combination 
    of conflict tasks. Journal of Experimental Psychology: Human Perception and 
    Performance, 42(10), 1505-1532. https://doi.org/10.1037/xhp0000229
    https://osf.io/ptg4n/wiki?wiki=r8bfj

    Parameters: 
    -----------
    directory : str 
        Path to combi files

    colnames : list 
        Colnames for data as given on osf 

    first : str 
        Name of the first conflict task, must match one colname

    second : str
        Name of the second conflict task, must match one colname

    exclude_catch: boolean
        Whether or not to exclude catch trials. Defaults to True 

    exclude_follow_catch: boolean
        Whether or not to exclude trials following catch trials. 
        Defaults to False 

    rt_min : int
        Lower RT bound for outlier exclusion in seconds. Defaults to 0
    
    rt_max : int
        Upper RT bound for outlier exclusion in seconds. Defaults to 2.5

    Returns: 
    --------
    return_data : pandas.DataFrame
        DataFrame containing data from the respective experiment
        - ``id``        : int, subject identifier
        - ``rt``        : float, Reaction Time
        - ``corr``      : float, accuracy
        - ``condition`` : float, condition
        - ``CI``        : str, condition
    """

    filelist = os.listdir(directory)

    df = pd.concat(
        [pd.read_csv(os.path.join(directory, f), sep=r"\s+", header=None)
        for f in filelist], 
        ignore_index=True
    )
    df.columns = colnames
    
    df["rt"] = df["rt"] / 1000.0        # convert rt to s
    
    df = df[df["resp"] != 0]            # exclude timeout responses
    df = df[(df["rt"] < rt_max) & (df["rt"] > rt_min)] # exclude rt outliers
    if exclude_catch and "catch" in df: 
        df = df[df["catch"] != 1]       # exclude catch
    #TODO: exclude trials following catch if exclude_follow_catch

    # create condition column
    ci_map = {
        (1, 1) : "congruent-congruent", 
        (1, 0) : "congruent-incongruent", 
        (0, 1) : "incongruent-congruent", 
        (0, 0) : "incongruent-incongruent"
    }
    df["CI"] = pd.Series(zip(df[first], df[second]), index=df.index).map(ci_map)
    
    df["condition"] = df["CI"].astype("category").cat.codes
    df["id"] = pd.factorize(df["subject"])[0]

    return_df = df[["id", "corr", "rt", "CI", "condition"]]

    return return_df


def format_sim_data(
    sim_data: dict, 
    sample_parallel: bool = False, 
    congruency_coding: dict = {
        0 : "congruent-congruent", 
        1 : "congruent-incongruent", 
        2 : "incongruent-congruent", 
        3 : "incongruent-incongruent"
    }, 
    only_convergents: bool = True
) -> pd.DataFrame:
    """
    Format simulated behavioral data into long-format pandas DataFrame

    sim_data : Dict[str, np.ndarray]
        Dictionary containing simulation outputs (resp, rt, condition) with 
        key = "sim_data"

    sample_parallel : bool, optional
        Whether the simulator was called with .sample(n_samples) or 
        .sample_parallel(n_batches, n_subjects, n_trials)

    congruency_coding : dict, optional
        Dictionary mapping integer codes to conditions 

    only_convergents : bool, optional
        If True, remove trials with reaction time equal to -1 (timeout). 
        Default is True 

    Returns 
    -------
    pd.DataFrame
        Long-format DataFrame with one row per trial and columns: 
        - 'rt': Reaction time
        - 'response': Accuracy value 
        - 'condition': Condition code 
        - 'batch': Batch nr
        - 'subject_id': Subject number
        - 'trial': trial nr
    """
    
    if not sample_parallel: 
        df_complete = pd.DataFrame(
            sim_data["sim_data"], 
            columns=["resp", "rt", "condition"]
        )
    
    if sample_parallel: 
        n_batches = sim_data["sim_data"].shape[0]
        n_subjects = sim_data["sim_data"].shape[1]
        n_trials = sim_data["sim_data"].shape[2]

        df_complete = pd.DataFrame(
            sim_data["sim_data"].reshape(-1, 3), 
            columns=["response", "rt", "condition"]
        )
        df_complete["batch"] = np.repeat(
            np.arange(n_batches), n_subjects * n_trials
        )
        df_complete["subject_id"] = np.tile(
            np.repeat(np.arange(n_subjects), n_trials), n_batches
        )
        df_complete["trial"] = np.tile(
            np.arange(n_trials), n_batches * n_subjects
        )

    if only_convergents:
        df_complete = df_complete[df_complete["rt"] != -1]
        
    df_complete['congruency'] = [
        congruency_coding[x] for x in df_complete['condition']
    ]

    return df_complete
        

def plot_subj_params(
    sim_data: dict, 
    group_level: bool = True, 
    params_names: list = [
        "muc", "b", "A1", "A2", "tau1", "tau2", "t0"
    ], 
    params_labels: list = [
        r"$\mu_c$", r"$b$", r"$A_1$", r"$A_2$", r"$\tau_1$", r"$\tau_2$", 
        r"$t_0$"
    ], 
    plot_type: str = "hist"
): 

    # get min and max of A1, A2
    x_min_A = min([min(sim_data["A1"].flatten()), min(sim_data["A2"].flatten())])
    x_max_A = max([max(sim_data["A1"].flatten()), max(sim_data["A2"].flatten())])
    # get min and max of tau1, tau2
    x_min_tau = min([min(sim_data["tau1"].flatten()), min(sim_data["tau2"].flatten())])
    x_max_tau = max([max(sim_data["tau1"].flatten()), max(sim_data["tau2"].flatten())])

    
    plt_func = sns.histplot if plot_type == "hist" else sns.kdeplot
   
    f, axarr = plt.subplots(2, 4, figsize=(12, 4))
    f.suptitle("Subject-Parameter Distributions")
    for i in range(len(params_names)): 
        ax = axarr.flat[i]
        plt_func(sim_data[params_names[i]].flatten(), alpha=0.75, ax=ax)
        sns.despine(ax=ax)
        ax.set_xlabel(params_labels[i])
        ax.set_yticks([])
        if (params_names[i] == "A1") or (params_names[i] == "A2"): 
            ax.set_xlim(x_min_A, x_max_A)
        if (params_names[i] == "tau1") or (params_names[i] == "tau2"): 
            ax.set_xlim(x_min_tau, x_max_tau)

    for ax in axarr.flat[len(params_names):]: 
        ax.set_visible(False)

    f.tight_layout()

    return f, axarr


def plot_group_parameters(
    samples: dict, 
    params_names: list = [
        "mu_muc", "mu_A1", "mu_A2", "mu_tau1", "mu_tau2", "mu_t0", "mu_b"
    ], 
    params_labels: list = [
        r"$\mu_c$", r"$A_1$", r"$A_2$", r"$\tau_1$", r"$\tau_2$", 
        r"$t_0$", r"$b$"
    ], 
    plot_type: str = "kde", 
    mean_line: bool = True, 
    suptitle: str = ""
):
    """
    Plot distributions of group level parameters

    samples: dictionary containing params as keys
    params: List of parameter names. Optional
    params_labels: to write pretty greek letters. Optional
    type: "kde" or "hist"
    mean_line: adds line at mean
    """

    plt_func = sns.histplot if plot_type == "hist" else sns.kdeplot
    samples_flat = {k: v.flatten() for k, v in samples.items()}

    f, axarr = plt.subplots(
        2, int(np.ceil(len(params_names)/2)), figsize=(12, 4), tight_layout=True, 
    )
    f.suptitle(suptitle)
    for i in range(len(params_names)):
        ax = axarr.flat[i]
        plt_func(samples_flat[params_names[i]], ax=ax)
        ax.set_title(f"{params_labels[i]}")
        ax.set_xlabel(params_labels[i])
        if mean_line: 
            ax.axvline(x=np.mean(samples_flat[params_names[i]]), color='red')

    for ax in axarr.flat[len(params_names):]: 
        ax.set_visible(False)

    return f, axarr


def plot_rts(
    data: pd.DataFrame, 
    data_2: Optional[pd.DataFrame] = None, 
    congruency: str = "congruency", 
    congruency_2: str = "congruency", 
    accuracy: str = "accuracy", 
    accuracy_2: str = "accuracy", 
    id_col: str = None, 
    id_col_2: str = None, 
    exclude_errors: bool = True, 
    x_min: float = 0, 
    x_max: float = 2, 
    label: str = None, 
    label_2: str = None, 
    panel_header: str = None, 
    n_rows: int = 2,  
    plot_type_kde: bool = True
): 
    """
    Plot (aggregated) RT distributions

    data : DataFrame. The rt data to be plotted
        Must contain the following columns: 
        - congruency: Column indicating the condition/congruency
        - 'rt': Column containing reaction time
        - accuracy: Column containing accuracy values: 1 = correct, 0 = error 
        - id_col: Column containing id to group by. If None: RTs will be plotted 
              across every subject, condition, batch, ...

    data_2 : DataFrame, Optional. Additional rt data to be added to the plot
        Must contain the following columns: 
        - congruency: Column indicating the condition/congruency
        - 'rt': Column containing reaction time
        - accuracy: Column containing accuracy values: 1 = correct, 0 = error 
        - id_col: Column containing id to group by. If None: RTs will be plotted 
              across every subject, condition, batch, ...

    congruency : string, Column name of congruency/condition column in data
    congruency_2 : string, Column name of congruency/condition column in data_2

    accuracy : string, Column name of accuracy column in data
    accuracy_2 : string, Column name of accuracy column in data_2

    id_col : string, Column name of id column in data
        RT in 'data' will be plotted by grouping on unique values in this column
    id_col_2 : string, Column name of id column in data_2
        RT in 'data_2' will be plotted by grouping on unique values in this 
        column
    
    exclude_errors : boolean
        Whether error trials should be excluded or not. Defaults to True. 

    x_min : int, lower x-axis limit. Default is 0. 
    x_max : int, upper x-axis limit. Default is 2. 

    label : str, label of plotted rt data in 'data'
    label_2 : str, label of plotted rt data in 'data_2'

    panel_header : str, common part of the header of each panel in subplot. 
        The unique id is added to this, separated by " ", 
        e.g. panel_header = 'subject' results in 'subject 1', 'subject 2', ...
    
    n_rows : int, number of rows in case of multiple subplots

    plot_type_kde : boolean
        Whether sns.kdeplot should be used. Defaults to True. If False, 
        sns.histplot is used. 

    Returns: 
    Plot, possibly including subplots
    """

    if exclude_errors: 
        data = data[data[accuracy] == 1]
        if data_2 is not None: 
            data_2 = data_2[data_2[accuracy_2] == 1]

    plt_func = sns.kdeplot if plot_type_kde else sns.histplot

    if id_col == None:       # plot overall RT distribution(s)

        f = plt_func(data["rt"], color="blue", label=label)
        if data_2 is not None: 
            plt_func(data_2["rt"], color="orange", label=label_2)

        f.set_xlabel("RT")
        f.set_xlim(x_min, x_max) 
        f.legend()
        f.set_title(panel_header)
        return f

    # plot by id in id_col
    unique_ids = np.unique(data[id_col])
    panel_header = panel_header if panel_header is not None else id_col

    n_cols = int(np.ceil(len(unique_ids)/n_rows))
    f, axarr = plt.subplots(
        n_rows, n_cols, figsize=(2*n_cols, 2*n_rows), 
        constrained_layout=True
    )

    for idx, i in enumerate(unique_ids):
        ax = axarr.flat[idx]
        plt_func(
            data[data[id_col] == i]["rt"], ax=ax, 
            label=label, color="blue"
        )
        if data_2 is not None: 
            plt_func(
                data_2[data_2[id_col_2] == i]["rt"], ax=ax, 
                label=label_2, color="orange"
            )

        sns.despine(ax=ax)
        ax.set_xlim(x_min, x_max)
        ax.set_title(f"{panel_header} {i}")
        ax.legend()

    for ax in axarr.flat[len(unique_ids):]: 
        ax.set_visible(False)

    return f, axarr


def ddmc_stats(
    data: pd.DataFrame, 
    id_name: str = "id", 
    n_bins: int = 5, 
    rt: str = "rt", 
    accuracy: str = "accuracy", 
    congruency: str = "congruency", 
    quantiles: Union[np.ndarray, Sequence[float]] = np.arange(0.1, 1.0, 0.1)
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:

    """
    Compute distributional summary statistics for RT data, producing inputs 
    suitable for plot_ddmc_stats function
    """

    data[rt] = pd.to_numeric(data[rt], errors="coerce")

    delta_data = (
        data[data[accuracy] == 1]
        .groupby([id_name, congruency])[rt]
        .quantile(quantiles)
        .reset_index()
        .rename(columns={"level_2": "quantile"})
        .pivot(index=[id_name, "quantile"], columns=[congruency], values=rt)
        .reset_index()
        .assign(delta_1=lambda df: df["incongruent-incongruent"] - df["congruent-incongruent"])
        .assign(delta_2=lambda df: df["incongruent-congruent"] - df["congruent-congruent"])
        .assign(mean_qu_1=lambda df: (df["incongruent-incongruent"] + df["congruent-incongruent"]) / 2)
        .assign(mean_qu_2=lambda df: (df["incongruent-congruent"] + df["congruent-congruent"]) / 2)
    )

    df = data.copy()

    df["rt_bin"] = pd.qcut(df[rt], q=n_bins, labels=False)

    caf_data = (
        df.groupby([id_name, congruency, "rt_bin"])[accuracy]
        .mean()
        .reset_index()
        .rename(columns={"level_2": "quantile"})
        .reset_index()
    )

    cdf_data = pd.melt(
        delta_data,
        id_vars=[id_name, "quantile"],
        value_vars=[
            "congruent-congruent", "congruent-incongruent", 
            "incongruent-congruent", "incongruent-incongruent"
        ],
        var_name=congruency,
        value_name=rt,
    )
    
    return delta_data, cdf_data, caf_data


def plot_stats( 
    caf_data: pd.DataFrame,
    cdf_data: pd.DataFrame,
    delta_data: pd.DataFrame,
    caf_data_2: Optional[pd.DataFrame] = None,
    cdf_data_2: Optional[pd.DataFrame] = None,
    delta_data_2: Optional[pd.DataFrame] = None,
    alpha: float = 0.05,
    id_name: str = "id",
    congruency: str = "congruency",
    congruency_2: str = "congruency",
    accuracy: str = "accuracy", 
    accuracy_2: str = "accuracy", 
    rt : str = 'rt',
    n_delta_bins: int = 10,
    fontsize: int = 24,
    fontsize_axes: int = 15,
    delta_ylim: Tuple[float, float] = None, 
    delta_xlim: Tuple[float, float] = None,
    caf_ylim: Tuple[float, float] = (0, 1), 
    cdf_ylim: Tuple[float, float] = (0, 1)
): 
    """
    Plot three standard distributional diagnostics for reaction-time (RT) data 
    for up to two datasets:
    (1) conditional accuracy function (CAF), 
    (2) cumulative distribution function (CDF),
    (3) a delta-function summary of condition differences across the RT 
        distribution, conditioned on task-two congruency.

    The function creates a single figure with the three subplots arranged 
    horizontally. Each subplots contain one or two datasets, e.g. of predicted 
    and/or empirical data. 

    Parameters
    ----------
    delta_data : pandas.DataFrame
        Long-format data required for the Δ-function panel. Must contain at least:

        - ``'quantile'``: Quantile index/label (used for aggregation).
        - ``'mean_qu_1'``: Mean RT associated with each quantile for task-2 
            incongruent trials (x-axis of Δ-function).
        - ``'mean_qu_2'``: Mean RT associated with each quantile for task-2 
            congruent trials (x-axis of Δ-function).
        - ``'delta_1'``: Difference metric to plot (y-axis of Δ-function) for 
            task 2 - incongruent trials.
        - ``'delta_2'``: Difference metric to plot (y-axis of Δ-function) for 
            task 2 - congruent trials.
        - A column named by ``id_name``: Identifier for individual trajectories.

        Notes
        -----
        The function will add a temporary column ``'mean_qu_bins'`` via ``pd.cut``.
        (It is overwritten if already present.)

    caf_data : pandas.DataFrame
        Data for the CAF panel. Must contain at least:

        - ``'rt_bin'``: RT bin index/label (x-axis of CAF).
        - ``'accuracy'``: Accuracy per bin (y-axis of CAF).
        - A column named by ``congruency``: Grouping variable for CAF lines.

    cdf_data : pandas.DataFrame
        Long-format data for the CDF panel. Must contain at least:

        - ``'rt'``: Reaction times in seconds (x-axis of CDF).
        - ``'quantile'``: CDF quantiles (y-axis of CDF).
        - ``'condition'``: Condition label for grouping/colouring.
        - A column named by ``id_name``: Identifier for individual trajectories.

    
    alpha : float, default=0.05
        Opacity for individual CDF trajectories (panel 2). The mean CDF is 
        plotted with opacity 1.0.

    id_name : str, default='id'
        Column name used as an identifier for individual trajectories in the 
        CDF and Δ-function panels.

    congruency : str, default='congruency'
        Column name used to stratify the CAF panel.

    rt: str, default='rt'
        Columns name of the reaction time variable.

    n_delta_bins : int, default=10
        Number of bins used when discretizing ``delta_data['mean_qu']`` into
        ``'mean_qu_bins'``. (The function currently computes a binned summary, 
        but then replaces it with a mean-by-quantile aggregation for plotting.)

    fontsize : int, default=24
        Font size for subplot titles.

    fontsize_axes : int, default=20
        Font size for axis labels.

    delta_ylim : tuple[float, float] | None, default=None
        If provided (truthy), apply a fixed y-axis range to the Δ-function panel.

    delta_xlim : tuple[float, float] | None, default=None
        If provided (truthy), apply a fixed x-axis range to the Δ-function panel.

    caf_ylim : tuple[float, float] | None, default=None
        If provided (truthy), apply a fixed y-axis range to the caf-function panel.
    
    cdf_ylim : tuple[float, float] | None, default=None
        If provided (truthy), apply a fixed y-axis range to the cdf-function panel.
    

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created matplotlib figure.

    axes : numpy.ndarray of matplotlib.axes.Axes
        Array of axes in the order ``[CAF, CDF, Δ-function]``.

    Examples
    --------
    >>> fig, axes = plot_stats(
            caf_data, cdf_data, delta_data, 
            caf_data_2, cdf_data_2, delta_data_2, 
            id_name="subject")
    """


    fig, axes = plt.subplots(1, 3, figsize=(16, 4))

    # CAF
    sns.lineplot(caf_data, x="rt_bin", y=accuracy, hue=congruency, ax=axes[0])
    if caf_data_2 is not None: 
        sns.lineplot(
            caf_data_2, x="rt_bin", y=accuracy_2, hue=congruency_2, ax=axes[0], 
            marker="o", linestyle="--"
        )

    axes[0].set_title("CAF", fontsize=fontsize)
    axes[0].set_ylabel("CAF", fontsize=fontsize_axes)
    axes[0].set_xlabel("Bins", fontsize=fontsize_axes)
    axes[0].set(ylim=caf_ylim)
    axes[0].legend(title="", loc="lower right")

    # single CDF
    sns.lineplot(
        cdf_data, x=rt, y="quantile", ax=axes[1], hue=congruency, 
        style=id_name, legend=False, alpha=alpha
    )
    
    # mean CDF
    mean_data = cdf_data.groupby(
        ["quantile", congruency]
    )[rt].mean().reset_index()
    sns.lineplot(
        mean_data, x=rt, y="quantile", hue=congruency, alpha=1, ax=axes[1]
    )

    if cdf_data_2 is not None: 
        mean_data_2 = cdf_data_2.groupby(
            ["quantile", congruency_2]
        )[rt].mean().reset_index()

        sns.lineplot(
            mean_data_2, x=rt, y="quantile", ax=axes[1], hue=congruency_2, 
            alpha=1, marker="o", linestyle="--"
        )

    axes[1].set_title("CDF", fontsize=fontsize)
    axes[1].set_xlabel("RT[s]", fontsize=fontsize_axes)
    axes[1].set_ylabel('Cumulative Density', fontsize=fontsize_axes)
    axes[1].set(ylim=cdf_ylim)
    axes[1].get_legend().remove()



    # single Deltas
    plot_data = pd.concat(
    [
        delta_data[[id_name, "mean_qu_1", "delta_1"]]
        .rename(columns={"mean_qu_1": "mean_qu", "delta_1": "delta"})
        .assign(series=1),
        delta_data[[id_name, "mean_qu_2", "delta_2"]]
        .rename(columns={"mean_qu_2": "mean_qu", "delta_2": "delta"})
        .assign(series=2),
    ],
    ignore_index=True
    )

    sns.lineplot(
        data=plot_data, x="mean_qu", y="delta", hue=id_name, ax=axes[2], 
        units="series", estimator=None, legend=False, linewidth=0.5, 
        linestyle="--", marker="o", alpha=alpha
    )


    def aggregate_delta(data):
        data["mean_qu_bins_1"] = pd.cut(data["mean_qu_1"], bins=n_delta_bins)
        data["mean_qu_bins_2"] = pd.cut(data["mean_qu_2"], bins=n_delta_bins)

        delta_bins_1 = data.groupby(
            "mean_qu_bins_1", observed=False
        )["delta_1"].mean().reset_index()
        delta_bins_2 = data.groupby(
            "mean_qu_bins_2", observed=False
        )["delta_2"].mean().reset_index()

        delta_bins_1["bin_mid_1"] = delta_bins_1["mean_qu_bins_1"].apply(lambda x: x.mid)
        delta_bins_2["bin_mid_2"] = delta_bins_2["mean_qu_bins_2"].apply(lambda x: x.mid)

        delta_bins_1 = (
            data.groupby("quantile")[["mean_qu_1", "delta_1"]]
            .mean().reset_index().sort_values("mean_qu_1")
        )
        delta_bins_2 = (
            data.groupby("quantile")[["mean_qu_2", "delta_2"]]
            .mean().reset_index().sort_values("mean_qu_2")
        )

        return delta_bins_1, delta_bins_2
    


    # aggregated Deltas
    delta_bins_1, delta_bins_2 = aggregate_delta(delta_data)
    sns.lineplot(
        data=delta_bins_1, x="mean_qu_1", y="delta_1", ax=axes[2],
        linewidth=1, legend=False, color="orange"
    )

    sns.lineplot(
        data=delta_bins_2, x="mean_qu_2", y="delta_2", ax=axes[2],
        linewidth=1, legend=False, color="blue"
    )

    if delta_data_2 is not None: 
        delta_bins_1_2, delta_bins_2_2 = aggregate_delta(delta_data_2)
        sns.lineplot(
            delta_bins_1_2, x="mean_qu_1", y="delta_1", ax=axes[2],
            linewidth=0.5, linestyle="--", marker="o", 
            legend=False, color="orange"
        )

        sns.lineplot(
            data=delta_bins_2_2, x="mean_qu_2", y="delta_2", ax=axes[2],
            linewidth=0.5, linestyle="--", marker="o",
            legend=False, color="blue",
        )

    axes[2].set_ylabel("$\\Delta$", fontsize=fontsize_axes)
    axes[2].set_xlabel("RT[s]", fontsize=fontsize_axes)
    axes[2].set_title("$\\Delta$-Function", fontsize=fontsize)

    if delta_ylim is not None:
        axes[2].set(ylim=delta_ylim)
    if delta_xlim is not None:
        axes[2].set(xlim=delta_xlim)

    fig.tight_layout()

    return fig, axes


def plot_subj_stats(
    data_1: pd.DataFrame, 
    data_2: Optional[pd.DataFrame] = None, 
    id_col_1: str = "id", 
    id_col_2: str = None, 
    accuracy_1: str = "accuracy", 
    accuracy_2: str = "accuracy", 
    congruency_1: str = "congruency", 
    congruency_2: str = "congruency", 
    delta_ylim: tuple = None, 
    delta_xlim: tuple = None, 
    panel_header: str = "Subj", 
    n_rows: int = 3, 
    n_cols: int = 8, 
    n_delta_bins: int = 5, 
    plot_type: str = ["caf", "cdf", "delta"], 
    x_label: str = "", 
    y_label: str = ""
): 
    """
    """

    unique_ids = np.unique(data_1[id_col_1])
    panel_header = panel_header if panel_header is not None else id_col_1

    n_cols = int(np.ceil(len(unique_ids)/n_rows))
    f, axarr = plt.subplots(
        n_rows, n_cols, figsize=(2*n_cols, 2*n_rows), 
        constrained_layout=True
    )

    for idx, i in enumerate(unique_ids):
        ax = axarr.flat[idx]
        plot_data_1 = data_1[data_1[id_col_1] == i]

        match (plot_type):
            case ("caf"):
                sns.lineplot(
                    plot_data_1, x="rt_bin", y=accuracy_1, hue=congruency_1, 
                    ax=ax, legend=False
                )

            case ("cdf"):
                sns.lineplot(
                    plot_data_1, ax=ax, x="rt", y="quantile", 
                    hue=congruency_1, legend=False
                )
        
            case ("delta"): 
                sns.lineplot(
                    plot_data_1, ax=ax, 
                    x="mean_qu_1", y="delta_1",
                    linewidth=1, legend=False,
                    color="orange", label="TODO"
                )
                sns.lineplot(
                    plot_data_1, ax=ax, 
                    x="mean_qu_2", y="delta_2",
                    linewidth=1, legend=False,
                    color="blue", label="TODO"
                )

            case (_): 
                return None

        if data_2 is not None: 
            plot_data_2 = data_2[data_2[id_col_2] == i]
            match (plot_type):
                case ("caf"):
                    sns.lineplot(
                        plot_data_2, x="rt_bin", ax=ax, y=accuracy_2, 
                        hue=congruency_2, marker="o", linestyle="--", 
                        legend=False
                    )
                case ("cdf"):
                    sns.lineplot(
                        plot_data_2, ax=ax, x="rt", y="quantile", 
                        hue=congruency_2, legend=False, 
                        marker="o", linestyle="--"
                    )

                case ("delta"):
                    sns.lineplot(
                        plot_data_2, ax=ax, linestyle="--", marker="o",
                        x="mean_qu_1", y="delta_1", legend=False, 
                        color="orange", linewidth=0.5
                    )
                    sns.lineplot(
                        plot_data_2, ax=ax, linestyle="--", marker="o",
                        x="mean_qu_2", y="delta_2", legend=False, 
                        color="blue", linewidth=0.5
                    )
                case (_):
                    return None


        # axes limits for delta plot        
        if delta_ylim is not None:
            ax.set(ylim=delta_ylim)
        if delta_xlim is not None:
            ax.set(xlim=delta_xlim)

        sns.despine(ax=ax)
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)       
        ax.set_title(f"{panel_header} {i}")

        if plot_type == "delta": 
            ax.legend()


    for ax in axarr.flat[len(unique_ids):]: 
        ax.set_visible(False)

    return f, axarr



def format_empirical_data(
    data: pd.DataFrame,
    rt: str = None,
    accuracy: str = None,
    congruency: str = None,
    format_type: str = "group",
    subj_id_col: str = "subject_id"
) -> Dict[str, np.ndarray]:
    """
    Format empirical data
    """

    if format_type == "group": 
        # extract relevant variables
        obs_data = data[[accuracy, rt, congruency]].values
        obs_data = obs_data[None, :, :]     # add dimension

    else: 
        subjects = []
        for subj_id, subj_df in data.groupby(subj_id_col): 
            obs_data_subj = subj_df[[accuracy, rt, congruency]].values
            subjects.append(obs_data_subj)
            obs_data = np.array(subjects, dtype=object)

    return obs_data


def score_log_norm(x, m, s):
    """ 
    Score function for gaussian priors
    """
    return -(x-m) / s**2

def score_log_mixture(x, m1, s1, m2, s2):
    """
    Score function for normal mixture priors
    """
    n1 = np.exp(-0.5 * ((x - m1) / s1) ** 2) / (np.sqrt(2*np.pi) * s1)
    n2 = np.exp(-0.5 * ((x - m2) / s2) ** 2) / (np.sqrt(2*np.pi) * s2)

    p = 0.5 * n1 + 0.5 * n2

    return (
        0.5 * n1 * (-(x - m1) / s1**2)
        + 0.5 * n2 * (-(x - m2) / s2**2)
    ) / p

def prior_global_score(
    x: dict[str, np.ndarray], 
    global_prior: dict
) -> dict[str, np.ndarray]:
    """
    Prior global score function
    
    Returns: dictionary, containing prior global scores based on priors
    """
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


def fit_empirical_data(
    data: pd.DataFrame, 
    approximator_global: Any, 
    approximator_local: Optional[Any], 
    priors: dict, 
    num_samples: int = 100, 
    id_name: str = "subject_id", 
    rt: str = "rt", 
    accuracy: str = "response", 
    congruency: str = "congruency", 
    subject_level: bool = True, 
    batch_size: int = 64, 
    n_local_subjects: int = 5, 
) -> pd.DataFrame:
    """
    Samples posteriors for empirical data for each unique subject or group

    Iterates over unique identifiers in the input DataFrame (e.g. subjects), 
    formats their data, performs posterior sampling using the specified 
    approximators and aggregates the results into a combined DataFrame

    data : pandas DataFrame

    approximator_global : global bf Approximator

    approximator_local : local bf Approximator
        Not needed for group-level sampling

    priors : dictionary
        Must contain the global priors specified for training the model(s)

    num_samples : int
        How many samples should be drawn per id

    id_name : str
        Indicates a column in data with unique ids. Sampling is done per id in 
        this column
    
    rt : str
        Column in data that contains RT data
    
    accuracy : str 
        Column in data that contains accuracy data

    congruency : str 
        Column in data that contains congruency/condition data

    subject_level : boolean
        Whether resampling should be done per subject id or not. Defaults to 
        True

    batch_size : int
        Batch size used during training
    
    n_local_subjects : int
        Number of local subjects used to train the model. 
    """

    if not subject_level: 
        group_data = format_empirical_data(
            data=data, rt=rt, accuracy=accuracy, congruency=congruency
        )
        samples = approximator.sample(
            num_samples=num_samples, 
            conditions={"sim_data": group_data}, 
            batch_size=batch_size*10
        )

        # reformat 
        samples_flat={k: v.flatten() for k, v in samples.items()}
        data_samples_complete=pd.DataFrame(samples_flat)

    else: 
        subjects = []
        min_trials = None
        for sid, sdf in data.groupby(id_name): 
            arr = np.column_stack([sdf[accuracy], sdf[rt], sdf[congruency]])
            
            subjects.append(arr)

            if min_trials is None or arr.shape[0] < min_trials: 
                min_trials = arr.shape[0]

        # truncate trials
        subjects = [arr[:min_trials] for arr in subjects]

        obs_data = np.stack(subjects)
        obs_data = obs_data[None, :, :]
        obs_data = {"sim_data": obs_data}

        n_subj = len(np.unique(data[id_name]))

        if n_subj == 23: 
            dummy_subject = obs_data["sim_data"][:, -1:, :, :]
            obs_data["sim_data"] = np.concatenate(
                [obs_data["sim_data"], dummy_subject], axis=1
            )

        obs_data["sim_data"] = obs_data["sim_data"].reshape(
            (1, 24 // n_local_subjects, min_trials * 6, 3)
        )
        summaries = None

        n_resims = num_samples
        global_posterior = approximator_global.compositional_sample(
            num_samples=n_resims,
            conditions={'sim_data': obs_data['sim_data']}, # get rid of all unnecessary keys
            compute_prior_score=lambda x: prior_global_score(x, global_prior=priors),
            compositional_bridge_d1=0.1,      # hyperparameter, control error accumulation
            compositional_bridge_d0=0.1,      # hyperparameter, control error accumulation
            method='two_step_adaptive',       # hyperparameter, try "euler_maruyama" for non-adaptive method
            steps='adaptive',                 # hyperparameter, speed vs accuracy/stability
            mini_batch_size=10,               # hyperparameter, memory/speed vs accuracy
            batch_size=batch_size*10,
            return_summaries=True,
            summaries=summaries,
        )
        summaries = global_posterior.pop("_summaries")

        obs_data["sim_data"] = obs_data["sim_data"].reshape((1, 24, min_trials, 3))

        samples = approximator_local.ancestral_sample(
            conditions={'sim_data': obs_data['sim_data']},
            ancestral_conditions=global_posterior,
            batch_size=batch_size*num_samples,
        )

        samples_flat = {k: v.flatten() for k, v in samples.items()}
        data_samples_complete = pd.DataFrame(samples_flat)
        data_samples_complete["id"] = np.repeat(np.arange(1, 25), n_resims)

    return data_samples_complete


def resim_data_id(
    param_samples: pd.DataFrame, 
    simulator: Any, 
    params_names: list, 
    num_resims: int, 
    n_resim_trials: int, 
    congruency_coding: dict = {
        0 : "congruent-congruent", 
        1 : "congruent-incongruent", 
        2 : "incongruent-congruent", 
        3 : "incongruent-incongruent"
    }, 
): 
    """
    Resimulate data based on different parameter sets. 

    param_samples : pandas df 
        must at least contain columns named identical to the parameter names. 
    
    simulator : bf Simulator
        Must return dictionary with the key 'sim_data'. This must have three 
        columns (accuracy, rt, condition)
    
    params_names : list
        must contain parameter names matching the arguments passed to the 
        simulator

    num_resims : int 
        number of resimuation-batches

    n_resim_trials : int 
        number of trials to be resimulated

    Returns: 
    -------
    pandas DataFrame, with the following columns: 
        - 'corr': float
        - 'rt': float
        - 'condition': float 
        - 'congruency': str
    """

    id_resims = []
    for i in range(num_resims): 
        iteration_dict = {key: values.iloc[i] for key, values in param_samples.items() if key in params_names}
        simulated = simulator.simulate_ddmc(
            **iteration_dict | {'n_trials': n_resim_trials}#, 'contamination_probability':0}
        )
        resim_batch = pd.DataFrame(simulated["sim_data"])
        resim_batch.columns = ["corr", "rt", "condition"]
        resim_batch["resim_id"] = i
        
        id_resims.append(resim_batch)

    all_resim_id = pd.concat(id_resims)
    all_resim_id["congruency"] = [
        congruency_coding[x] for x in all_resim_id['condition']
    ]
    return all_resim_id


def resim_data(
    empirical_data: pd.DataFrame, 
    post_samples: pd.DataFrame, 
    simulator: Any, 
    id_name: str = "subect_id", 
    num_resims: int = 100, 
    n_resim_trials: int = 4, 
    resim_per_id: bool = True, 
    params_names: list = ["muc", "A1", "A2", "tau1", "tau2", "t0", "b"], 
    exclude_nonconvergents: bool = True
):
    """
    Resimulate DDMC data based on posterior samples of parameters

    empirical_data : pandas DataFrame 
        Must contain 

    post_samples : pandas DataFrame
        Must contain columns with identical names as those in params_names

    simuator : bf Simulator

    num_resims : int
        How many resimulations should be run per participant. 

    n_resim_trials : int
        How many trials should be resampled per participant. 

    resim_per_id : boolean
        Whether to resample per id (e.g. subject) or not. Defaults to True. 

    params_names : list 
        Parameter names indicating the columns with the posterior samples to be 
        used for resimulating

    exclude_nonconvergents : boolean 
        Whether timeout respnses (indicated by -1) should be excluded. Defaults 
        to True. 

    id_name : string
        Column name indicating the (e.g. subject) id column in both 
        empirical_data and post_samples

    Returns: 
    --------
    pandas DataFrame containing resimulated DDMC data 

    """

    post_samples_params = post_samples[params_names]

    list_resim_dfs = []

    if resim_per_id: 
        ids = empirical_data[id_name].unique()

        for i in range(len(ids)): 
            id = ids[i]
            n_obs = empirical_data[(empirical_data[id_name] == id)].shape[0]
            subj_samples = post_samples[post_samples[id_name] == id+1]

            # resimulate
            resim_subj_df = resim_data_id(
                subj_samples, simulator=simulator, num_resims=num_resims, 
                params_names=params_names, n_resim_trials=n_resim_trials
            )

            resim_subj_df["subj_id"] = id
            list_resim_dfs.append(resim_subj_df)

    else: 
        resim_samples = dict(post_samples)

        for i in range(num_resims): 
            iteration_dict = {key: values[i] for key, values in resim_samples.items() if key in params_names}
            resim =  simulator.simulate_ddmc(
                **iteration_dict | {'n_trials': num_resims}
            )

            resim_df = pd.DataFrame(resim)
            resim_df["num_resim"] = i
            
            list_resim_dfs.append(resim_df)


    data_resim = pd.concat(list_resim_dfs)

    if exclude_nonconvergents: 
        data_resim = data_resim[data_resim["rt"] != -1]

    return data_resim