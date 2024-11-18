"""
To run this script, you need to set the following environment variables:
- BASECALLING_SUMMARY_FN: the path to the basecalling summary file
- Q_SCORE_THRESHOLD: the threshold of q score
- JOB_MAGICK_NAME: the name of the job
- FAST5_FOLDER: the folder containing fast5 files
- CTC_TRACE_FOLDER: the folder containing ctc trace files (npz)
- REFERENCE_FASTA_FN: the path to the reference fasta file

to run the script (take caution to the tqdm and pkl output files, you may want to change the path):
mpirun -np 160 python3 extract_kmer_hpc.py

"""
import random
import numpy as np
from scipy import signal
from scipy import stats
from lib import fast5_processing
import matplotlib.pyplot as plt
import pandas as pd
from tombo import tombo_helper, tombo_stats, resquiggle
import pandas as pd
from tqdm import trange, tqdm
import os
import h5py, mappy
import time
import ctc_segmentation
from mpi4py import MPI
import sys
import pickle

# 初始化 MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
# print(rank, size)

# load input fn
input_basecalling_summary_fn = os.getenv('BASECALLING_SUMMARY_FN')
q_score_threshold = float(os.getenv('Q_SCORE_THRESHOLD'))
job_magick_name = os.getenv('JOB_MAGICK_NAME')
fast5_folder = os.getenv('FAST5_FOLDER')
ctc_trace_folder = os.getenv('CTC_TRACE_FOLDER')
reference_fasta_fn = os.getenv('REFERENCE_FASTA_FN')
# 读取数据（只在主进程）
if rank == 0:
    print("process size is: ", size)
    basecalling_summay = pd.read_csv(input_basecalling_summary_fn, sep='\t')

    # filter mean_q_score < 17.5, reset index
    basecalling_summay = basecalling_summay[basecalling_summay['mean_q_score'] > q_score_threshold].reset_index(drop=True)
    # filter sequence length < 1000, reset index
    basecalling_summay = basecalling_summay[basecalling_summay['sequence_length_template'] > 900].reset_index(drop=True) # you can change the length here
    basecalling_summay = basecalling_summay[basecalling_summay['sequence_length_template'] < 1400].reset_index(drop=True)

    print("loaded basecalling_summay, shape: ", basecalling_summay.shape)

    # 将数据分割并分发到各个进程
    # 确保将数据平均分配到所有进程
    # 计算每个chunk的基本大小和多出来的行数
    total_rows = len(basecalling_summay)
    rows_per_chunk = total_rows // size
    extra_rows = total_rows % size

    # 创建chunks
    chunks = []
    start = 0
    for i in range(size):
        # 为前extra_rows个chunk分配额外一行
        end = start + rows_per_chunk + (1 if i < extra_rows else 0)
        chunks.append(basecalling_summay[start:end])
        start = end

else:
    basecalling_summay = None
    chunks = None

# 将数据分割并分发到各个进程
chunk = comm.scatter(chunks, root=0)

char_list = [ "", "A", "C", "G", "T"]
encoding_dict = {'':  0, 'A': 1 , 'C':  2 , 'G':  3 , 'T':  4 }
# base_ctc_data_dir = "/data/SUSTech03/Cate-5KDNA-nomod/single_fast5_ctc_trace/"
base_ctc_data_dir = ctc_trace_folder
base_fast5_dir = fast5_folder
config = ctc_segmentation.CtcSegmentationParameters(char_list=char_list, min_window_size = 8000, max_window_size = 8000)
config.index_duration = 1/4000

# define mappy aligner
aligner = mappy.Aligner(reference_fasta_fn, preset=str('map-ont'), best_n=1)

def align_seq(aligner,seq_seq_str):
    try:
        ref_seq = aligner.seq("Cate_NAN")
        alignment = next(aligner.map(seq_seq_str))
        align_start = alignment.r_st
        align_end = alignment.r_en
        align_start_seq = alignment.q_st
        align_end_seq = alignment.q_en
        align_strand = alignment.strand
        if align_strand == -1:
            ref_seq = mappy.revcomp(ref_seq)
            align_start = len(ref_seq) - align_end
            align_end = len(ref_seq) - align_start
        aligned_ref_seq = ref_seq[align_start:align_end]
        return aligned_ref_seq, align_start, align_end, align_start_seq, align_end_seq, alignment.strand
    except StopIteration:
        return None, None, None, None, None, None

def extract_kmer_data(df, current):
    results = []
    # for i in range(2, len(df) - 3):
    for i in range(3, len(df) - 2):
        row = df.iloc[i]
        if row['conf'] > 0.8:
            kmer = ''.join(df['text'][i-3:i+3])
            start, end = int(row['start']), int(row['end'])
            current_segment = current[start:end]
            mean, std = np.mean(current_segment), np.std(current_segment)
            if std <= 0.4:
                results.append((kmer, mean))
    return results

def ctc_align_and_extract_current(row):
    read_id = read_id = row['read_id']
    # print("processing ", read_id)
    current = fast5_processing.extract_raw_signal_from_fast5(base_fast5_dir + read_id + ".fast5")
    probs_stack_len, probs_stack, trace_path = np.load(base_ctc_data_dir + read_id + ".npz").values()
    trace_path_start, trace_path_end = trace_path[0], trace_path[-1]
    stub_first = len(current) - probs_stack_len
    # call align_seq
    aligned_ref_seq, align_start, align_end, align_start_seq, align_end_seq, align_strand = align_seq(aligner, row['sequence'])
    # if return None, return empty list, end this function
    if aligned_ref_seq is None or len(aligned_ref_seq) < 900:
        return []
    

    # clip probs_stack according to trace_path_end
    probs_stack_start = np.max([0,trace_path[align_start_seq]-20])
    probs_stack_end = np.min([trace_path[align_end_seq-1] + 20, probs_stack.shape[0]])
    probs_stack = probs_stack[probs_stack_start:probs_stack_end]
    tokens = [np.array([encoding_dict[i]]) for i in aligned_ref_seq]
    bases_list = [i for i in aligned_ref_seq]

    # ctc seg
    # ground_truth_mat, utt_begin_indices = ctc_segmentation.prepare_token_list(config, tokens)
    # timings, char_probs, state_list = ctc_segmentation.ctc_segmentation(config, probs_stack, ground_truth_mat)
    # segments = ctc_segmentation.determine_utterance_segments(config, utt_begin_indices, char_probs, timings, bases_list)
    try:
        # ctc seg
        ground_truth_mat, utt_begin_indices = ctc_segmentation.prepare_token_list(config, tokens)
        timings, char_probs, state_list = ctc_segmentation.ctc_segmentation(config, probs_stack, ground_truth_mat)
        
        # 如果 ctc_segmentation 抛出异常，以下代码将不会执行
        segments = ctc_segmentation.determine_utterance_segments(config, utt_begin_indices, char_probs, timings, bases_list)
    except Exception as e:
        # 在这里处理异常，当 ctc_segmentation 函数抛出异常时返回空数组
        print("An error occurred: ", str(e))
        return []

    # convert seg
    expdict = [{"text" : t, "start" : p[0], "end" : p[1], "conf" : p[2]} for t,p in zip(bases_list, segments)]
    expdict_df = pd.DataFrame(expdict)
    expdict_df['start'] = expdict_df['start']*4000 + stub_first + probs_stack_start
    expdict_df['end'] = expdict_df['end']*4000 + stub_first + probs_stack_start

    # norm current
    norm_current_start, norm_current_end = int(expdict_df['start'][2]), int(expdict_df['end'][-3:-2])
    norm_current = tombo_stats.normalize_raw_signal(current[norm_current_start:norm_current_end])[0]
    current[norm_current_start:norm_current_end] = norm_current
    return extract_kmer_data(expdict_df, current)


# 每个进程处理自己的部分
# TODO add logging
# results = [ctc_align_and_extract_current(row) for index, row in chunks.iterrows()]


# 每个进程处理自己的部分
# 每个进程处理自己的部分，并将结果保存为pickle文件
def process_and_save(chunk, rank):
    results = []
    progress_filename = f"tqdm_log/{job_magick_name}_process_{rank}_progress.txt"
    with open(progress_filename, "w") as tqdm_file:
        for index, row in tqdm(chunk.iterrows(), total=len(chunk), file=tqdm_file, desc=f"Process {rank}"):
            result = ctc_align_and_extract_current(row)
            results.append(result)
    
    # 将结果保存为pickle文件
    pickle_filename = f"intermediate_pkl/{job_magick_name}_process_{rank}_results.pkl"
    print("save pkl: ", pickle_filename)
    with open(pickle_filename, "wb") as f:
        pickle.dump(results, f)

# 每个进程调用process_and_save函数
process_and_save(chunk, rank)

# 主进程合并结果并保存
if rank == 0:
    print("save pkl done.")


# 结束MPI
MPI.Finalize()