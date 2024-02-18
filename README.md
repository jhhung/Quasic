# Quasic

Quasic is a single-cell RNA-seq quantification framework that considers the potential cell subpopulation structures within the data when estimating gene expression profiles. Quasic enables the analysis of cell subpopulation information, even in datasets with high dropout events. This empowers researchers to conduct more reliable cell subpopulation-related studies.

## Information

Workflow of Quasic:

![image](https://hackmd.io/_uploads/Bycojqk3a.png)

Three main steps are included in Quasic:

1. **Preprocessing**: An initial gene expression profile is estimated by Alevin (https://salmon.readthedocs.io/en/latest/alevin.html) to search as the first gauss of the prior
2. **Cenctroid Estimation**: Cells are clustered by Louvain clustering. The centroid for each cluster and the cluster assigning scores are then estimated
3. **EM-algorithm with constraints**:  The subpopulation-related constraints integrated into the EM algorithm, penalizing GEP estimations that deviate from the GEP associated with each centroid. 

Note that the code implementation of Qusic is built based on Salmon (https://github.com/COMBINE-lab/salmon) and tximport (https://github.com/thelovelab/tximport)

---
 
## Requirement

*  C++11 compiler (current test on gcc 7.5.0)
*  R (version >= 4.0)
*  [Cmake](https://cmake.org/)
*  [Boost](https://www.boost.org/) (current test on 1.71.0)
*  [Seurat](https://satijalab.org/seurat/) (current test on version 4.1.1)

---

## Installtion
```shell=
git clone https://github.com/jhhung/Quasic.git
cd Quasic/salmon
mkdir build && cd build
cmake ..
make
make install
```
---

## Execution

### Build salmon index

```shell
salmon index -t [transcripts_fasta_path] -i [transcripts_index_name]
```
### Run Quasic
```shell
./salmon alevin \
	--chromiumV3 \
	-l ISR -1 [inputR1] \
	-2 [inputR2] \
	-i [transcripts_index_name] \
  -o [output_dir] \
	--tgMap [transcript_to_gene_map] \
	-p [thread_num] \
	--regularization_parameter [regularization_intensity] \
	--rare_gene_cut_off [rare_gene_threshold] \
	--cluster_resolution [cluster_resolution] \
	--geneBlacklist [gene_blacklist]
```

#### Parameter explanation:

* Required
	* **-1 [inputR1]** : Input scRNA-seq R1 fastq reads。
	* **-2 [inputR2]**: Input scRNA-seq R2 fastq reads。
	* **-i [transcripts_index_name]** : Salmon index directory built by salmon index。See `data/salmon_index_filtered`.
	* **-o [output_dir]** : Output directory path
	* **--tgMap [transcript_to_gene_map]** : Mapping file for transcript and gene. See `data/txp2gene.tsv`.
* Option
	* **-p [thread_num] (default = 1)**: Number of using thread.
	* **--regularization_parameter [regularization_intensity] (default = 0.01)** : Intensity of subpopulation constraints. The value should be set between 0.0 and 1.0. Higher intensity represent the gene expression profile for each cell may close to its corresponding cell type signature.
	* **--rare_gene_cut_off [rare_gene_threshold] (default = 0.2)** : A gene may considered as a 'rare gene' for a cluster if the proportion of cells in the cluster that has expression this gene is lower than the threshold. During subpopulation-aware quantification, the expression of rare genes would not be regularized.
	* **--cluster_resolution (default = 0.2)[cluster_resolution]** : Resolution of Lovain clustering. Using higher resolution if you think there are more subpopulation in your scRNA-seq data.  
	* **--geneBlacklist (default = none) [gene_blacklist]** : List of transcripts that should not include in subpopulation-aware quantification. See `data/cell_cycle_related_transcript.txt`
