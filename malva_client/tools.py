import tempfile
import subprocess
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from gprofiler import GProfiler

def mask_sequence(sequence: str) -> str:
    """
    Apply repeat masking using dustmasker to the provided sequence.
    
    The function writes the input sequence to a temporary FASTA file, then
    calls dustmasker with:
      - Input format: FASTA (-infmt fasta)
      - Parse sequence IDs (-parse_seqids)
      - Output format: FASTA (-outfmt fasta) so that the resulting 
        masked sequence is easy to parse.
      
    If the dustmasker call fails, the original sequence is returned.
    """
    try:
        with tempfile.NamedTemporaryFile("w+", delete=False) as tmp_in:
            tmp_in.write(">seq\n" + sequence)
            tmp_in_path = tmp_in.name
        
        tmp_out_path = tmp_in_path + ".masked"
        
        # Build the dustmasker command.
        # (If ASN.1 binary output is preferred, change "-outfmt fasta" 
        # to "-outfmt maskinfo_asn1_bin" and then parse accordingly.)
        cmd = [
            "dustmasker",
            "-in", tmp_in_path,
            "-infmt", "fasta",
            "-parse_seqids",
            "-outfmt", "fasta",
            "-out", tmp_out_path
        ]
        
        subprocess.run(cmd, check=True)

        masked_seq = ""
        with open(tmp_out_path, "r") as f:
            lines = f.readlines()
            if lines and lines[0].startswith(">"):
                masked_seq = "".join(line.strip() for line in lines[1:])
            else:
                masked_seq = "".join(line.strip() for line in lines)
        
        os.remove(tmp_in_path)
        os.remove(tmp_out_path)
        
        return masked_seq
    except Exception as e:
        print(f"Warning: dustmasker failed with error {e}. Returning unmasked sequence.")
        return sequence

def run_go_enrichment(
    gene_list: list[str],
    *,
    organism: str = "hsapiens",
    background: list[str] | None = None,
    sources: list[str] | None = ("GO:BP", "GO:MF", "GO:CC"),
    significance: float = 0.05,
    top_n: int = 10,
    plot: bool = True,
    figsize: tuple[int, int] = (6, 4),
    verbose: bool = True,
):
    """
    Simple GO enrichment for a list of genes.

    Parameters
    ----------
    gene_list : list[str]
        Symbols of genes of interest (e.g. markers from Scanpy).
    organism : str, default "hsapiens"
        Species code used by g:Profiler (mmusculus, drerio, …).
    background : list[str] | None
        Custom universe. If None, g:Profiler uses its own default.
    sources : list[str] | None
        Which annotation namespaces to query.  Default = all three GO branches.
    significance : float, default 0.05
        Adjusted-p-value cut-off (g_SCS multiple-testing).
    top_n : int, default 10
        How many significant terms to keep (after filtering).
    plot : bool, default True
        If True, draw a horizontal bar-plot (-log₁₀ adj-p).
    figsize : tuple[int, int]
        Size of the matplotlib figure (only if `plot=True`).
    verbose : bool
        Print the chosen terms to stdout.

    Returns
    -------
    go_df : pandas.DataFrame
        Enrichment table sorted by adjusted p-value.
    """
    gp = GProfiler(return_dataframe=True)
    res = gp.profile(
        organism=organism,
        query=gene_list,
        user_threshold=significance,
        significance_threshold_method="g_SCS",
        background=background,
        sources=sources,
    )

    if res.empty:
        if verbose:
            print("No GO term reached the chosen significance threshold.")
        return res

    # sort and take top_n
    res = res.sort_values("p_value").head(top_n).reset_index(drop=True)
    res["minus_log10_padj"] = -res["p_value"].apply(lambda p: np.log10(p))

    if verbose:
        display_cols = ["native", "name", "p_value", "intersection_size"]
        print(res[display_cols].to_string(index=False))

    if plot:
        plt.figure(figsize=figsize)
        ylabels = res["native"] + "  "  # small padding
        plt.barh(ylabels, res["minus_log10_padj"])
        plt.xlabel(r"$-\log_{10}(\mathrm{adj}\,p)$")
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.title("GO enrichment (top terms)")
        plt.tight_layout()
        plt.show()

    return res

def score_correlated_features(
    adata,
    feature_key: str,
    *,
    method: str = "gmm",
    gmm_components: int = 3,
    positive_label: str = "True",
    negative_label: str = "False",
    n_markers: int = 5,
    significance: float = 0.05,
    rank_kwargs: dict | None = None,
    umap_kwargs: dict | None = None,
    random_state: int | None = 0,
    show: bool = True,
    smooth_neighbors: bool = True,
    smooth_alpha: float = 0.5,
    smooth_iterations: int = 1,
):
    """
    Detects a data-driven threshold for a continuous per-cell feature, bins cells
    into positive/negative categories, and ranks/visualises genes correlated with
    the positive group. Optionally smooths the feature across transcriptomic neighbors
    before thresholding.

    Parameters
    ----------
    adata : AnnData
        Object containing the feature in `adata.obs[feature_key]` and an
        existing UMAP in `adata.obsm['X_umap']`. Should have a neighbors graph
        computed if smooth_neighbors=True.
    feature_key : str
        Column in `adata.obs` with the per-cell score to threshold.
    method : {"gmm", "otsu"}, default "gmm"
        Thresholding strategy:
        * "gmm" - "gmm_components"-component Gaussian mixture (scikit-learn)
    gmm_components : int
        Number of components to fit the GMM model
    positive_label / negative_label : str
        Category names for the new boolean flag.
    n_markers : int, default 5
        How many significant markers (per group) to visualise.
    significance : float, default 0.05
        Benjamini-Hochberg FDR cut-off for marker selection.
    rank_kwargs : dict | None
        Extra keyword arguments passed to `sc.tl.rank_genes_groups`.
    umap_kwargs : dict | None
        Extra keyword arguments passed to `sc.pl.umap`.
    random_state : int | None
        Ensures reproducibility for the GMM.
    show : bool
        Whether to immediately display the UMAP plots.
    smooth_neighbors : bool, default True
        Whether to smooth the feature across transcriptomic neighbors before thresholding.
    smooth_alpha : float, default 0.5
        Smoothing factor: 0 = no smoothing, 1 = complete neighbor averaging.
        Final value = (1-alpha) * original + alpha * neighbor_mean
    smooth_iterations : int, default 1
        Number of smoothing iterations to apply.

    Returns
    -------
    threshold : float
        Numeric threshold chosen by the algorithm.
    pos_markers, neg_markers : list[str], list[str]
        Lists of the marker gene names visualised.
    """
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp
    from sklearn.mixture import GaussianMixture

    if feature_key not in adata.obs:
        raise KeyError(f"'{feature_key}' not found in adata.obs")

    # Get the original feature values
    feature_values = np.nan_to_num(adata.obs[feature_key].values.astype(float))
    
    # Apply neighbor smoothing if requested
    if smooth_neighbors:
        # Check if neighbors graph exists
        if 'neighbors' not in adata.uns:
            raise ValueError(
                "Neighbors graph not found in adata.uns['neighbors']. "
                "Please compute neighbors first using:\n"
                "  sc.pp.neighbors(adata)\n"
                "or with custom parameters:\n"
                "  sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)"
            )
        
        # Get the connectivity matrix
        if 'connectivities' not in adata.obsp:
            raise ValueError(
                "Connectivity matrix not found in adata.obsp['connectivities']. "
                "Please recompute neighbors using sc.pp.neighbors(adata)"
            )
        
        connectivities = adata.obsp['connectivities']
        
        # Ensure it's a sparse matrix
        if not sp.issparse(connectivities):
            connectivities = sp.csr_matrix(connectivities)
        
        # Normalize connectivity matrix (row-wise) to get averaging weights
        # Each row sums to 1, representing the weight distribution among neighbors
        row_sums = np.array(connectivities.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        connectivities_norm = connectivities.multiply(1 / row_sums[:, np.newaxis])
        
        # Apply smoothing iterations
        smoothed_values = feature_values.copy()
        
        for iteration in range(smooth_iterations):
            # Compute neighbor-averaged values
            # Matrix multiplication: each cell gets weighted average of its neighbors
            neighbor_means = connectivities_norm.dot(smoothed_values)
            
            # Apply smoothing: blend original with neighbor average
            smoothed_values = (1 - smooth_alpha) * smoothed_values + smooth_alpha * neighbor_means
        
        # Store smoothed values for potential inspection
        smoothed_key = f"{feature_key}_smoothed"
        adata.obs[smoothed_key] = smoothed_values
        
        # Use smoothed values for thresholding
        x = smoothed_values.reshape(-1, 1)
        
        if show:
            print(f"Applied neighbor smoothing with alpha={smooth_alpha}, iterations={smooth_iterations}")
            print(f"Smoothed values stored in adata.obs['{smoothed_key}']")
    else:
        x = feature_values.reshape(-1, 1)

    # Apply thresholding method
    if method == "gmm":
        gmm = GaussianMixture(n_components=gmm_components, random_state=random_state).fit(x)
        means = gmm.means_.flatten()
        labels = gmm.predict(x)
        threshold = means.mean()          # midpoint between the component means
        positive_mask = labels == means.argmax()  # component with higher mean
    else:
        raise ValueError("method must be 'gmm'. Others are not supported for now")

    # Create categorical flag
    flag_key = f"{feature_key}_positive"
    adata.obs[flag_key] = pd.Categorical(
        np.where(positive_mask, positive_label, negative_label)
    )

    if umap_kwargs is None:
        umap_kwargs = dict(frameon=False)

    # Plot the categorical result and optionally the smoothed feature
    if show:
        plot_features = [flag_key]
        if smooth_neighbors:
            plot_features.append(smoothed_key)
        
        sc.pl.umap(adata, color=plot_features, **umap_kwargs, show=True)

    # Differential expression analysis
    if rank_kwargs is None:
        rank_kwargs = dict()

    sc.tl.rank_genes_groups(adata, groupby=flag_key, **rank_kwargs)

    pos_df = sc.get.rank_genes_groups_df(adata, group=positive_label)
    neg_df = sc.get.rank_genes_groups_df(adata, group=negative_label)

    pos_markers = (
        pos_df.loc[pos_df["pvals_adj"] < significance, "names"]
        .head(n_markers)
        .tolist()
    )
    neg_markers = (
        neg_df.loc[neg_df["pvals_adj"] < significance, "names"]
        .head(n_markers)
        .tolist()
    )

    if show:
        sc.pl.umap(
            adata,
            color=pos_markers + neg_markers,
            ncols=5,
            frameon=False,
            show=True,
        )

    return threshold, pos_markers, neg_markers