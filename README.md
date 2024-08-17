# DeepSME

*De Novo Nanopore Basecalling of Motif Insensitive DNA Modifications And Alignment-free Digital Information Decryptions at Single-Molecule Level*

---
We have proposed and verified a scheme for high-informative secure communication combining DNA methylation induced information encryption and nanopore sequencing based decryption, via de novo assembly of a private DeepSME basecaller.

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
2. **uncalled4-py310** environment for generate the k-mer table of 5hmC modified DNA. plase refer to [Document of Uncalled4](https://github.com/skovaka/uncalled4) to install the environment.

## Model weight

You can view and download the config and the weight of the model for 5hmC from [this link](https://assets.sparktour.me/doc/publication/reinforced_basecaller_model_5hmc.tar.gz), then use `tar -xvf reinforced_basecaller_model_5hmc.tar.gz` to extract the model.

## Basecalling

> [!TIP]
> The model architecture and weight of DeepSME is compatible with Bonito. You may use bonito or bonito-compatible basecaller to basecall the sequence current.

### Aligned (Output BAM)
``` python3
bonito basecaller \
--chunksize 3600 \
--device cuda:0 \
--reference reference.fasta \
reinforced_basecaller_model/ \
fast5_or_pod5_reads/ > basecalling.bam
```

#### K-mer model extraction use uncalled4
``` python
uncalled4 train reference.fasta \
gDNA_fast5/ \
--bam-in basecalling.bam \
-m dna_r9.4.1_400bps_tombo.DNA.model.20240419_5hmC_central3mer.npz \
-p 16 -k 6 --kmer-shift 3 --out-dir kmer_5hmC.k6 \
--train-iterations 5 --init-mode moves
```

### Unaligned (Output fastq)
``` python
bonito basecaller \
--chunksize 3600 \
--device cuda:0 \
reinforced_basecaller_model/ \
fast5_or_pod5_reads/ > basecalling.fastq
```

## Citation

If you use DeepSME in your research, please cite:

``` bibtex
@article {Fan2024.08.15.606762,
	author = {Fan, Qingyuan and Zhao, Xuyang and Li, Junyao and Liu, Ronghui and Liu, Ming and Long, Yanping and Fu, Yang and Feng, Qishun and Zhai, Jixian and Pan, Qing and Li, Yi},
	title = {DeepSME: De Novo Nanopore Basecalling of Motif-insensitive DNA Methylation and Alignment-free Digital Information Decryptions at Single-Molecule Level},
	elocation-id = {2024.08.15.606762},
	year = {2024},
	doi = {10.1101/2024.08.15.606762},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/08/17/2024.08.15.606762},
	eprint = {https://www.biorxiv.org/content/early/2024/08/17/2024.08.15.606762.full.pdf},
	journal = {bioRxiv}
}
```

