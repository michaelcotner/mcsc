# mcsc

A collection of functions that I've made over time for single cell analysis. Mostly things that I made that I have used or likely will use more than once. Use alongside `scanpy`. This is mostly for my own use and reference, but you're welcome to be here.

## use

Navigate to any directory accessible by your python script and run:

```
wget https://raw.githubusercontent.com/michaelcotner/mcsc/main/mcsc.py
```

In your python script, import the file with:

```py
from path_to_file import mcsc
```

### updating

If I've added new functions to `mcsc` that you want to use, you can update the file with:

```
wget -N -O mcsc.py https://raw.githubusercontent.com/michaelcotner/mcsc/main/mcsc.py
```

`mcsc` will warn you if there is a new version available.

# documentation

## rgg_logfoldchange
Updates a specified `scanpy.tl.rank_genes_groups` `uns` layer, calculated from a counts layer.

Since scanpy.get.rank_genes_groups expects logarithmized data, the log fold change calculation they give is meaningless. This function calculates the log fold change based on a counts layer and updates the `rank_genes_groups` `uns`.

### Parameters
 - adata
     - The annotated data matrix, on which scanpy.tl.rank_genes_groups has already been run with a uns key specified in the `key` parameter
 - counts_layer
     - default value: `'counts'`
     - The layer (in `adata.layers`) containing raw counts or normalized counts which are *not* logarithmized, on which to do the log fold change calculation.
 - key
     - default value: `'rank_genes_groups'`
     - The uns key that contains `rank_genes_group` data; specified as the `key_added` parameter in `scanpy.tl.rank_genes_groups`. The default value is `rank_genes_groups`, which is the default value in `sc.tl.rank_genes_groups` and is what you should use if you didn't specify a key.
 - uns_key_added
     - default value: `'logfoldchanges'`
     - The name of the key within the `rank_genes_group` uns. If kept as its default value, `logfoldchanges`, this data will overwrite the default log fold change values from `scanpy.tl.rank_genes_groups`. If you do not want the old log fold change values to be overwritten, you should change this key.
 - copy
     - default value: `False`
    Returns an updated copy of the matrix if `True`. Otherwise the passed anndata object is updated and `None` is returned.
        
### Returns
 - `None` if `copy` = `False`. Updates the anndata object passed. If `copy` = `True`, returns an updated copy of the anndata object passed.

### Examples
```py
>>> import mcsc
>>> import scanpy as sc
>>> adata
AnnData object with n_obs × n_vars = 12050 × 13818
    obs: 'sample'
    layers: 'counts'
>>> sc.tl.rank_genes_groups(adata, 'sample', method='wilcoxon')
>>> sc.get.rank_genes_groups_dotplot(adata, groups='pretreat') # logfoldchange column is full of NaN
>>> rgg_logfoldchange(adata)
>>> sc.get.rank_genes_groups_dotplot(adata, groups='pretreat') # logfoldchange contains meaningful values based on raw counts
```

## volcano_plot
Updates a specified `scanpy.tl.rank_genes_groups` `uns` layer, calculated from a counts layer.

Since scanpy.get.rank_genes_groups expects logarithmized data, the log fold change calculation they give is meaningless. This function calculates the log fold change based on a counts layer and updates the `rank_genes_groups` `uns`.

### Parameters

 - df
     - The output of `sc.get.rank_genes_groups_df()`, after running `mcsc.rgg_logfoldchange()` to fix log fold change values.
     - Realistically, this can be any dataframe with columns containing gene names, log fold change values and p-values, whose names are passed to `genes_key`, `lfc_key` and `p_key`, respectively.
 - lfc_key
     - default value `'logfoldchanges'`
     - The name of the column in `df` that contains the log2 fold change values for each gene. The default value is the column name from `sc.tl.rank_genes_groups()`
 - p_key
     - default value `'pvals_adj'`
     - The name of the column in `df` that contains the p-values for each gene. The default value is the column name from `sc.tl.rank_genes_groups()`
 - genes_key
     - default value `'names'`
     - The name of the column in `df` that contains the names of each gene. The default value is the column name from `sc.tl.rank_genes_groups()`
 - lfc_cutoff
     - default value `0.5`
     - The minimum log2 fold change in gene expression to highlight as upregulated or downregulated genes. Downregulated genes must be lower than the negative of this value.
 - p_cutoff
     - default value `'0.05'`
     - The maximum p-value to highlight as statistically significant upregulated or downregulated genes. 
 - cutoff_style
     - default value `'and'`
     - The style to color upregulated and downregulated genes. The styles are:
         - `'and'`: Genes must have a log2 fold change greater than the cutoff *and* a p-value smaller than the cutoff.
         - `'or'`: Genes must have a log2 fold change greater than the cutoff *or* a p-value smaller than the cutoff.
         - `'radius'`: Genes must meet ((log2foldchange-mean(log2foldchange))/`lfc_cutoff`)^2 + ((pvalue-mean(pvalue))/`p_cutoff`)^2 > 1
 - title
     - default value `None`
     - Adds this title to the axis if not `None`.
 - dot_size
     - default value `1`
     - The dot size passed as `matplotlib.pyplot.scatter(s=dot_size)`
 - upreg_color
     - default value `'#aa0000'`
     - The color to highlight upregulated genes as. Follows `matplotlib`'s [color specification](https://matplotlib.org/stable/users/explain/colors/colors.html).
 - downreg_color
     - default value `'#0000aa'`
     - The color to highlight upregulated genes as. Follows `matplotlib`'s [color specification](https://matplotlib.org/stable/users/explain/colors/colors.html).
 - grid
     - default value `False`
     - Draws a grid on the axis when `True`.
 - annotate_genes
     - default value `False`
     - When `True`, uses the `gene_annotation_cutoff_style` to annotate gene names that meet the cutoff.
     - When a list of genes, annotates those genes on the volcano plot.
 - gene_annotation_lfc_cutoff
     - default value `1.0`
     - The minimum log2 fold change in gene expression a gene must have to be annotated.
 - gene_annotation_p_cutoff
     - default value `0.001`
     - The minimum p-value a gene must have to be annotated.
 - gene_annotation_cutoff_style
     - default value `and`
     - The style to annotate genes with. See `cutoff_style` for a list of styles and their descriptions.
 - annotate_fontsize
     - default value `small`.
     - The font size passed as `matplotlib.pyplot.annotate(fontsize=annotate_fontsize)`.
     - Can be `xx-small`, `x-small`, `small`, `medium`, `large`, `x-large`, or `xx-large`.
 - plot
     - default value `True`
     - Shows the plot if `True`. Returns the `matplotlib.pyplot.ax` object if `False`.
 - return_genes
     - default value `False`
     - Returns the dataframe passed to `df=` with only the genes that met the log2 fold change and p-value cutoffs.

### Returns
 - df
     - If `return_genes` is `True`, returns an updated dataframe including only the genes that meet the passed log2 fold change and p-value cutoffs.
 - ax
     - If `plot` is `False`, returns a matplotlib.pyplot axis object containing the volcano plot.

### Examples
```py
>>> import mcsc
>>> import scanpy as sc
>>> adata
AnnData object with n_obs × n_vars = 12050 × 13818
    obs: 'sample'
    layers: 'counts'
>>> sc.tl.rank_genes_groups(adata, 'sample', method='wilcoxon')
>>> mcsc.rgg_logfoldchange(adata)
>>> mcsc.volcano_plot(
>>>     sc.get.rank_genes_groups_dotplot(adata, groups='pretreat')
>>> ) # displays the volcano plot
>>> mcsc.volcano_plot(
>>>     sc.get.rank_genes_groups_df(adata, group='PT'),
>>>     lfc_cutoff=1, 
>>>     p_cutoff=1e-10,
>>>     annotate_genes=True,
>>>     annotate_fontsize='x-small',
>>>     gene_annotation_cutoff_style='radius',
>>>     gene_annotation_lfc_cutoff=4,
>>>     gene_annotation_p_cutoff=1e-200,
>>>     plot=False
>>> )   # example with more options - returns an axis with annotated genes
```