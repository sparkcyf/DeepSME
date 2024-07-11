import numpy as np
import h5py
from tombo import tombo_helper, tombo_stats, resquiggle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from ont_fast5_api.fast5_interface import get_fast5_file
import h5py

def extract_raw_signal_from_fast5(fast5_path):
    with h5py.File(fast5_path, 'r') as f:
        # Get the list of all reads in the file
        read_names = list(f['Raw/Reads'].keys())
        
        # Assuming you're interested in the first read for simplicity
        read_name = read_names[0]
        
        # Extract the raw signal
        raw_signal = f[f'Raw/Reads/{read_name}/Signal'][:]
        
        # Extract the conversion factors (if they exist)
        try:
            digitisation = f[f'UniqueGlobalKey/channel_id'].attrs['digitisation']
            offset = f[f'UniqueGlobalKey/channel_id'].attrs['offset']
            range_ = f[f'UniqueGlobalKey/channel_id'].attrs['range']
            sampling_rate = f[f'UniqueGlobalKey/channel_id'].attrs['sampling_rate']
            
            # Convert the raw signal to pA
            normalized_signal = (raw_signal + offset) * range_ / digitisation
        except:
            normalized_signal = raw_signal  # If conversion factors don't exist, return the raw signal as is

    return normalized_signal

def extract_tombo_data_from_fast5(filename: str):
    with h5py.File(filename, 'r') as fast5_file:
        # Extract current array
        current_array = extract_raw_signal_from_fast5(filename)
        
        # Extract aligned data
        events_data = fast5_file['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'][:]
        read_start_rel_to_raw = fast5_file['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'].attrs['read_start_rel_to_raw']
        
        # Using NumPy vectorized operations to improve performance
        starts = events_data['start'] + read_start_rel_to_raw
        bases = events_data['base'].astype(str)
        lengths = events_data['length']
        norm_means = events_data['norm_mean']

        # Creating a DataFrame in one go to save time
        aligned_data = pd.DataFrame({'start': starts, 'base': bases, 'length': lengths, 'norm_mean': norm_means})

        # Adding the mapped index column
        alignment = fast5_file['Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment']
        mapped_start = alignment.attrs['mapped_start']
        mapped_end = alignment.attrs['mapped_end']
        mapped_index = np.arange(mapped_start, mapped_end)
        aligned_data['mapped_index'] = mapped_index[:len(aligned_data)]
        
    return current_array, aligned_data

def extract_tombo_data_only_from_fast5(filename: str):
    with h5py.File(filename, 'r') as fast5_file:

        # Extract aligned data
        events_data = fast5_file['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'][:]
        read_start_rel_to_raw = fast5_file['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'].attrs['read_start_rel_to_raw']
        
        # Using NumPy vectorized operations to improve performance
        starts = events_data['start'] + read_start_rel_to_raw
        bases = events_data['base'].astype(str)
        lengths = events_data['length']
        norm_means = events_data['norm_mean']

        # Creating a DataFrame in one go to save time
        aligned_data = pd.DataFrame({'start': starts, 'base': bases, 'length': lengths, 'norm_mean': norm_means})

        # Adding the mapped index column
        alignment = fast5_file['Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment']
        mapped_start = alignment.attrs['mapped_start']
        mapped_end = alignment.attrs['mapped_end']
        mapped_index = np.arange(mapped_start, mapped_end)
        aligned_data['mapped_index'] = mapped_index[:len(aligned_data)]
        
    return aligned_data

def plot_current_with_aligned_bases(current_array, aligned_data, x_start: int, x_end: int):
    """
    Plots the current signal with aligned bases from a Fast5 file with Tombo-aligned data.
    
    Parameters:
    - x_start (int): Starting point of the x-axis range.
    - x_end (int): Ending point of the x-axis range.
    """
    
    # 1. Extract Data
    #current_array, aligned_data = extract_tombo_data_from_fast5(filename)
    
    # 2. Filter Aligned Data
    filtered_aligned_data = aligned_data[(aligned_data['start'] >= x_start) & (aligned_data['start'] <= x_end)]
    
    # 3. Adjust Indices for Current Array
    adjusted_indices = filtered_aligned_data['start']
    
    # Ensure the adjusted range is within the bounds of current_array
    adjusted_range_start = max(0, x_start)
    adjusted_range_end = min(len(current_array), x_end)
    
    # Adjusted x-values for plotting
    x_values = np.arange(x_start, x_start + (adjusted_range_end - adjusted_range_start))
    
    # 4. Plot Current Signal
    plt.figure(figsize=(15, 5))
    plt.plot(x_values, current_array[adjusted_range_start:adjusted_range_end], label="Current Signal")
    
    # 5. Annotate with Aligned Bases
    for start, base, adj_index in zip(filtered_aligned_data['start'], filtered_aligned_data['base'], adjusted_indices):
        if adjusted_range_start <= adj_index < adjusted_range_end:  # Ensure the annotation is within the plot range
            #plt.annotate(base, (start, current_array[adj_index]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=9, color='red')
            plt.annotate(base, (start+1, 2), textcoords="offset points", xytext=(0,10), ha='center', fontsize=9, color='red')
            plt.axvline(start, color='grey', linestyle='--', linewidth=0.5)
    
    # 6. Customize and Display Plot
    plt.xlim(x_start, x_end)
    plt.xlabel("Position in Raw Signal")
    plt.ylabel("Normalized Current")
    plt.title("Current Signal with Aligned Bases")
    plt.legend()
    plt.tight_layout()
    plt.show()