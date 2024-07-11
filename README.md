# DeepSME

*De Novo Nanopore Basecalling of Motif Insensitive DNA Modifications And Alignment-free Digital Information Decryptions at Single-Molecule Level*

---
We have proposed and verified a scheme for high-informative secure communication combining DNA methylation induced information encryption and nanopore sequencing based decryption, via de novo assembly of a private DeepSME basecaller.

## TOC and folder structure

## Install

Four conda environment may be required to install:

1. **tombo-py37** environment for align the reference sequence to current. plase refer to [Document of Tombo](https://github.com/nanoporetech/tombo). You may need to lock the python version to 3.7 to successfully install it.
2. **preliminary_basecaller-py311** environment for train and infer the Preliminary Basecaller
3. **bonito-py38** environment for train and infer the Enhanced Basecaller and Reinforced Basecaller. plase refer to [Document of Bonito](https://github.com/nanoporetech/bonito)
4. **uncalled4-py310** environment for generate the k-mer table of 5hmC modified DNA. plase refer to [Document of Uncalled4](https://github.com/skovaka/uncalled4)

## Datasets and model weight

You can download the model weight and dataset chunks from [https://doi.org/10.5281/zenodo.12704171](https://doi.org/10.5281/zenodo.12704171).

Raw fast5 data are available on reasonable requests.

If you just want to test the basecaller, you can download the model weight and jump to the [Basecalling](#Basecalling) section.

## Train

### Preliminary Basecaller

Train the Preliminary Basecaller:

``` python
python ./scripts/train_original.py \
--data-dir preliminary_datasets_dir \
--output-dir model_dir \
--model bonito \
--window-size 2000 \
--batch-size 128 \
--use-scaler \
--num-epochs 50 \
--starting-lr 0.0005 \
--overwrite
```
Use Preliminary Basecaller to basecall the 6-mer QC Sequence get the CTC data:

``` python
python ./scripts/basecall_original.py \
--fast5-dir fast5_dir  \
--checkpoint model_dir/<best_checkpoint.pt> \
--model bonito \
--chunk-size 2000 \
--window-overlap 400 \
--batch-size 32 \
--output-file fastq_folder/basecalling_preliminary.fastq
```

#### Backtrace the CTC data

See `train/preliminary/kmer_extraction/extract_kmer_hpc.py` and `train/preliminary/kmer_extraction/extract_kmer_hpc_merge_pkl.py` . You may need to run it with a `mpi4py` and `dtaidistance` (we recommend you install these package in `tombo-py37` environment). 

After the backtrace and alignment, you shall get the kmers of the 5hmC modified DNA. You can use this kmer to replace the original kmer in the `legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model` from [https://github.com/nanoporetech/kmer_models/](https://github.com/nanoporetech/kmer_models/) . You can also use this kmer to generate the simulated current and train the Enhanced Basecaller.

### Enhanced Basecaller

#### Generate simulated current data

You can use squigulator or other current simulator to generate the simulated current data. Then you need to use tombo to generate chunks from the simulated current data.
``` bash
BIN_PATH="squigulator-v0.3.0/squigulator"
FASTA_PATH="reference.fasta"
PARAM="-x dna-r9-min --kmer-model kmer_model.csv -o blow5/simulated_current.blow5 -q blow5/simulated_current.fasta -f 2 --seed 1711350749  -t 32"
# RUN
$BIN_PATH $FASTA_PATH $PARAM
```
#### Train Enhanced Basecaller
``` python
bonito train \
--epochs 1 \
--lr 5e-6 \
--pretrained dna_r9.4.1_e8_sup@v3.3 \
--directory squigulator_chunks \
enhanced_basecaller_model
```

### Reinforced Basecaller

#### Build Datasets

``` python
bonito basecaller \
enhanced_basecaller_model \
--device cuda:0 \
--save-ctc \
--reference gDNA_reference.fasta \
--min-accuracy-save-ctc 0.7 \
--chunksize 3600 \
gDNA_fast5 > reinforced_basecaller_chunks/basecalls.sam
```

#### Train Reinforced Basecaller

``` python
bonito train \
--epochs 5 \
--lr 5e-4 \
--batch 128 \
--device cuda:0 \
--pretrained dna_r9.4.1_e8_sup@v3.3 \
--directory reinforced_basecaller_chunks/ \
reinforced_basecaller_model/
```

## Basecalling
The model architecture and weight of DeepSME is compatible with Bonito. You may use bonito or bonito-compatible basecaller to basecall the sequence current.

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