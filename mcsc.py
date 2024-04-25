VERSION = '0.2'





######################
### Internal stuff ###
######################

def _check_version():
    # check if there is a new version of the file

    # imports
    try:
        import urllib.request
    except ImportError as e:
        return

    # hosted file to compare to
    url = 'https://raw.githubusercontent.com/michaelcotner/mcsc/main/mcsc.py'

    try:
        with urllib.request.urlopen(url) as response:
            # get version from the first line
            file_head = response.readline().decode('utf-8')
            if file_head.startswith('VERSION = '):
                ref_version = file_head[len('VERSION = '):].strip('\'\n')
                print('Your mcsc version is', VERSION, 'but version', ref_version, 'exists on github.')
                print('You can update your mcsc by running:')
                print('\t', 'wget -N -O mcsc.py https://raw.githubusercontent.com/michaelcotner/mcsc/main/mcsc.py')
                return
    except Exception as e:
        return
    
    return

_check_version()





#################
### Functions ###
#################

def rgg_logfoldchange(adata:'anndata.AnnData', *, counts_layer='counts', uns_key_added='logfoldchanges', key='rank_genes_groups', copy=False):
    """
    Updates a specified `scanpy.tl.rank_genes_groups` `uns` layer, calculated from a counts layer.

    Since scanpy.get.rank_genes_groups expects logarithmized data, the log fold change calculation they give is meaningless. This function calculates the log fold change based on a counts layer and updates the `rank_genes_groups` `uns`.

    Parameters
    ----------
    adata
        The annotated data matrix, on which scanpy.tl.rank_genes_groups has already been run with a uns key specified in the `key` parameter
    counts_layer
        The layer (in `adata.layers`) containing raw counts or normalized counts which are *not* logarithmized, on which to do the log fold change calculation. `counts` by default.
    key
        The uns key that contains `rank_genes_group` data; specified as the `key_added` parameter in `scanpy.tl.rank_genes_groups`. The default value is `rank_genes_groups`, which is the default value in `sc.tl.rank_genes_groups` and is what you should use if you didn't specify a key.
    uns_key_added
        The name of the key within the `rank_genes_group` uns. If kept as its default value, `logfoldchanges`, this data will overwrite the default log fold change values from `scanpy.tl.rank_genes_groups`. If you do not want the old log fold change values to be overwritten, you should change this key.
    
    Returns
    -------
    `None` if `copy` = `False`. Updates the anndata object passed. If `copy` = `True`, returns an updated copy of the anndata object passed.

    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.rank_genes_groups(adata, 'bulk_labels', method='wilcoxon')
    >>> sc.get.rank_genes_groups_dotplot(adata, groups='Dendritic') # logfoldchange column is full of NaN
    >>> rgg_logfoldchange(adata)  # this dataset doesn't actually have a counts layer so this won't run. pretend it does.
    >>> sc.get.rank_genes_groups_dotplot(adata, groups='Dendritic') # logfoldchange contains meaningful values based on raw counts
    """

    # imports
    try:
        import anndata
        import numpy as np
        import pandas as pd
    except ImportError as e:
        print('One or more packages required for rgg_logfoldchange are not installed.\n', e)
        return None

    # check parameter values
    if key not in adata.uns_keys():
        if key == 'rank_genes_groups':
            raise ValueError('specified key', key, 'not found in uns. Ensure that you have run rank_genes_groups already and are using the same key value, if applicable.')
        else:
            raise ValueError('specified key', key, 'not found in uns. Ensure this value matches the value passed to \'key_added\' when running \'rank_genes_groups\'')
    
    if counts_layer not in adata.layers:
        raise ValueError('specified layer', counts_layer, 'not found in layers. Ensure you are specifying your desired counts layer if it is not named \`counts\`.')

    if uns_key_added != 'logfoldchanges':
        adata.uns[key][uns_key_added] = adata.uns['rank_genes_groups']['logfoldchanges'].copy()

    if copy:
        adata = adata.copy()

    groups = adata.uns['rank_genes_groups']['scores'].dtype.names  # list of groups in obs
    obs_name = adata.uns['rank_genes_groups']['params']['groupby']  # obs name that the groups come from


    # build a dataframe for log fold changes for each group
    lfcs = pd.DataFrame(index=adata.var_names, columns=groups)

    if adata.uns[key]['params']['reference'] == 'rest':
        for group in groups:
            lfcs[group] = np.log2(
                np.average(adata.layers[counts_layer][adata.obs[obs_name] == group].toarray(), axis=0)/
                np.average(adata.layers[counts_layer][adata.obs[obs_name] != group].toarray(), axis=0)
            )
    else:
        for group in groups:
            lfcs[group] = np.log2(
                np.average(adata.layers[counts_layer][adata.obs[obs_name] == group].toarray(), axis=0)/
                np.average(adata.layers[counts_layer][adata.obs[obs_name] == adata.uns[key]['params']['reference']].toarray(), axis=0)
            )

    # change the logfoldchange values in the uns
    for i in range(adata.uns[key][uns_key_added].shape[0]):
        for j, group in enumerate(groups):
            adata.uns[key][uns_key_added][i][j] = lfcs.loc[adata.uns[key]['names'][i][j], group]

    if copy:
        return adata
    else:
        return None

    '''
    
    this old looping method had a very long runtime
    the uns that scanpy adds makes this difficult - the genes are in their ranked order, which means a different order for each group
    we'll have to do the log fold change for each gene first, and then re-order them.

    # loop through ranked genes
    for i in range(adata.uns[key][uns_key_added].shape[0]):
        # loop through groups
        for j, group in enumerate(groups):
            # gene = adata.uns[key]['names'][i][j]  # gene name
            gene_idx = np.where(adata.var_names == adata.uns[key]['names'][i][j])[0][0] # faster to index off of gene index rather than gene name
            
            group_idx = np.where(adata.obs['sample'] == group)[0]  # indices where a cell is in our group
            if adata.uns[key]['params']['reference'] == 'rest':
                ref_idx = np.where(adata.obs['sample'] != group)[0]  # indices where a cell isn't in our group
            else:
                ref_idx = np.where(adata.obs['sample'] == adata.uns[key]['params']['reference'])[0] # indices where a cell is in the reference group

            exp_avg = np.mean(adata.layers[counts_layer][group_idx, gene_idx])  # average expression in group
            ref_exp_avg = np.mean(adata.layers[counts_layer][ref_idx, gene_idx])  # average expression in reference
            # print(exp_avg, ref_exp_avg)
            adata.uns[key][uns_key_added][i][j] = np.log2(exp_avg/ref_exp_avg)

    return adata
    '''
        
def volcano_plot():
    pass




