## Train the Preliminary Basecaller

1. clone https://github.com/sparkcyf/deepsme_preliminary_basecaller.git

```bash
git clone https://github.com/sparkcyf/deepsme_preliminary_basecaller.git
```

2. create the conda environment according to https://github.com/sparkcyf/deepsme_preliminary_basecaller/blob/7170375754d8c3b31a4b89ee3001bbf3315c5dac/environment.yml

```bash
conda env create -f environment.yml
```

3. Training Preliminary Basecallers

Once gathered enough current and reference chunks (in format of [marcpaga/basecalling_architectures](https://github.com/marcpaga/basecalling_architectures/?tab=readme-ov-file#chunk-the-raw-signal-and-save-it-numpy-arrays)) by using “constant velocity assumption” to segment the current and reference sequence, you can use the command and the code on the GitHub to train the preliminary basecaller and extract CTC data from them as follows:

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

4. Basecalling with Preliminary Basecallers

If you want to save the basecalling results to a fastq file, you can use the following command:

> [!TIP]  
> If you want to save the CTC data, you need to specify `CTC_TRACE_SAVE_FOLDER` environment variable to the folder where you want to save the CTC data for extracting the k-mer model.

``` python
export CTC_TRACE_SAVE_FOLDER=your_ctc_trace_folder/

python ./scripts/basecall_original.py \
--fast5-dir 6mer_QC_sequence_fast5_dir/  \
--checkpoint model_dir/<best_checkpoint.pt> \
--model bonito \
--chunk-size 2000 \
--window-overlap 400 \
--batch-size 32 \
--output-file fastq_folder/basecalling_preliminary.fastq
```

## Extract the k-mer table

You may need to run the script on HPC cluster with `train/preliminary/kmer_extraction/extract_kmer_hpc.lsf` to extract the k-mer table.

If you just want to see the extracted k-mer table, you can view the `train/preliminary/kmer_extraction/kmer_models/` folder.

- the `train/preliminary/kmer_extraction/kmer_models/dna_r9.4.1_400bps_6mer_5hmc.csv` can be used for [squigulator](https://github.com/hasindu2008/squigulator).
- the `train/preliminary/kmer_extraction/kmer_models/dna_r9.4.1_400bps_6mer_5hmc.npz` can be used for [Uncalled4](https://github.com/skovaka/uncalled4/).
- the `train/preliminary/kmer_extraction/kmer_models/tombo.DNA.5hmC.DeepSME.model` can be used for [Tombo](https://github.com/nanoporetech/tombo).
