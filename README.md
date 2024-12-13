# DeepSME

*De Novo* Non-Canonical Nanopore Basecalling Unlocks Private Communication using Heavily-modified DNA Data at Single-Molecule Level

---

In short, this basecaller model can basecall the fully 5hmC modified DNA sequence with high accuracy.

## Install

### Clone the repository

``` bash
git clone https://github.com/sparkcyf/DeepSME.git
```

### Prepare the conda environment

> [!NOTE]
> You need a CUDA compatible GPU to train and infer the model. We recommend using a GPU with at least 12GB of memory.

The following conda environment may be required to install:

1. **bonito-py38** environment for train and infer the Enhanced Basecaller and Reinforced Basecaller. plase refer to [Document of Bonito](https://github.com/nanoporetech/bonito) to install the environment. We have tested the model on `bonito 0.81` with `python 3.8`, `cuda 11.8`.
2. **tombo-py37** environment for generate the k-mer table for 5hmC modified DNA. plase refer to [Document of Tombo](https://github.com/nanoporetech/tombo) to install the environment.

## Model weight and 5hmC-modified DNA Storage Datasets

You can view and download the config, weight of the model for 5hmC DeepSME, sequencing pod5 file (raw sequencing current) and basecalled fastq files from this link (https://mirrors.sustech.edu.cn/site/datasets-share/deepsme/), then use `tar -xvf reinforced_basecaller_model_5hmc.tar.gz` to extract the model. You can also check the raw sequence data and decode scripts of 5hmC-modified DNA Storage Datasets at https://github.com/sparkcyf/DeepSME_DNA_Storage_Decode_scripts .


## Usage of 5hmC DeepSME Basecaller

### Basecalling

> [!TIP]
> The model architecture and weight of DeepSME is compatible with Bonito. You may use bonito or bonito-compatible basecaller to basecall the sequence current.

#### Aligned (Output BAM)
``` python3
bonito basecaller \
--chunksize 3600 \
--device cuda:0 \
--reference reference.fasta \
reinforced_basecaller_model/ \
fast5_or_pod5_reads/ > basecalling.bam
```

#### Unaligned (Output fastq)
``` python
bonito basecaller \
--chunksize 3600 \
--device cuda:0 \
reinforced_basecaller_model/ \
fast5_or_pod5_reads/ > basecalling.fastq
```

#### (Optional) K-mer model extraction use uncalled4 based on 5hmC-modified 6-mer model

``` python
uncalled4 train reference.fasta \
gDNA_fast5/ \
--bam-in basecalling.bam \
-m dna_r9.4.1_450bps_6mer_5hmc.npz \
-p 16 -k 6 --kmer-shift 3 --out-dir kmer_5hmC.k6 \
--train-iterations 5 --init-mode moves
```

## Citation

If you find DeepSME useful in your research, please cite:

``` bibtex
@article {Fan2024.08.15.606762,
	author = {Fan, Qingyuan and Zhao, Xuyang and Li, Junyao and Liu, Ronghui and Liu, Ming and Long, Yanping and Fu, Yang and Feng, Qishun and Zhai, Jixian and Pan, Qing and Li, Yi},
	title = {De Novo Non-Canonical Nanopore Basecalling Unlocks Private Communication using Heavily-modified DNA Data at Single-Molecule Level},
	elocation-id = {2024.08.15.606762},
	year = {2024},
	doi = {10.1101/2024.08.15.606762},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/08/17/2024.08.15.606762},
	eprint = {https://www.biorxiv.org/content/early/2024/08/17/2024.08.15.606762.full.pdf},
	journal = {bioRxiv}
}
```

