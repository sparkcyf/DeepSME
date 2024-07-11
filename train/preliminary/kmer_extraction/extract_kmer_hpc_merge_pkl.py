import pandas as pd
import os
import pickle

all_results = []
size = 160
for i in range(size):  # 假设size是进程数
    pickle_filename = f"intermediate_pkl/process_{i}_results.pkl"
    print("load pkl: ", pickle_filename)
    if os.path.exists(pickle_filename):
        with open(pickle_filename, "rb") as f:
            results = pickle.load(f)
            all_results.extend(results)
# 合并结果
print("-------merging result---------")
# 展平嵌套列表， flatten twice
# combined_results = [item for sublist in all_results for sublist2 in sublist for item in sublist2]
combined_results = [item for sublist in all_results for item in sublist]
# print(combined_results)
print("len of combined result: ", len(combined_results))
result_df = pd.DataFrame(combined_results, columns=['kmer', 'mean'])

# save to csv
result_df.to_csv('kmer_data/kmer_data_raw.csv', index=False)

def aggregate_kmer_means(df):
    # 按照 kmer 分组，并计算每个组的 mean 平均值
    aggregated_df = df.groupby('kmer')['mean'].mean().reset_index()
    return aggregated_df

# 调用函数并获取合并后的 DataFrame
aggregated_result_df = aggregate_kmer_means(result_df)

# save to csv
aggregated_result_df.to_csv('kmer_data/kmer_data_aggregated.csv', index=False)

print("save csv procedure done.")