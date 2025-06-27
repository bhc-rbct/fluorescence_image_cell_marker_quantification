import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import gaussian_kde
from scipy.signal import savgol_filter, find_peaks
from skimage.filters import threshold_otsu, threshold_triangle, threshold_minimum
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
np.random.seed(42) 


PANCK_THRESHOLD = 0.0716
KRT56_THRESHOLD = 0.1259 


EPSILON = 1e-6  # Smallest detectable signal

def process_all_files(results_folder):
    """
    Process KRT56 positive and negative CSV files from quantifications subfolders.
    
    Args:
        results_folder (str): Path to the root folder containing tissue folders
        
    Returns:
        tuple: (krt56_pos_df, krt56_neg_df) - DataFrames containing combined data
    """
    # Initialize empty lists to store DataFrames before concatenation
    all_cells_dfs = []
    
    for folder_name in os.listdir(results_folder):
        folder_path = os.path.join(results_folder, folder_name, "quantifications")
        
        if os.path.isdir(folder_path):
            print(f"Processing folder: {folder_path}")
            
            # Get all files in the quantifications folder
            files = os.listdir(folder_path)

            all_cells_files = [f for f in files 
                             if f.startswith('Tissue_') and f.endswith('_allCells_allMarkers.csv')]
            
            # Process positive files
            for f in all_cells_files:
                file_path = os.path.join(folder_path, f)
                df = pd.read_csv(file_path)
                
                # Add metadata columns
                df['source_tissue'] = folder_name
                df['source_file'] = f
                
                all_cells_dfs.append(df)
            
    
    # Combine all DataFrames
    panCK_df = pd.concat(all_cells_dfs, ignore_index=True) if all_cells_dfs else pd.DataFrame()
    
    return panCK_df

def calculate_thresholds(data):
    """Calculate multiple thresholds for a given dataset."""
    thresholds = {}
    
    # Otsu
    thresholds['Otsu'] = threshold_otsu(data)
    
    # Triangle
    thresholds['Triangle'] = threshold_triangle(data)

    # Minimal
    try:
        thresholds['Minimum'] = threshold_minimum(data)
    except:
        pass

    # K-means (K=2)
    kmeans = KMeans(n_clusters=2, n_init=50, max_iter=500).fit(data.reshape(-1, 1))
    thresholds['K-means'] = np.mean(kmeans.cluster_centers_)
    
    # GMM (Gaussian Mixture Model)
    gmm = GaussianMixture(n_components=2, n_init=50, max_iter=500, init_params='kmeans').fit(data.reshape(-1, 1))
    thresholds['GMM'] = np.mean(gmm.means_)
    
    return thresholds
def create_mutiple_dist_histogram(data_list, label_list, color_list, plot_name, column, tick_space = 0.5, transform_method = "log10"):
    plt.figure(figsize=(10, 6))
    ax = plt.gca()
    plot_title = f'{column} Distribution'

    for i, (data, label, color) in enumerate(zip(data_list, label_list, color_list)):

        if transform_method: 
            # Data processing
            transformed_data = safe_transform(data, transform_method)
            if i == 0:
                plot_name = f'{transform_method}_{plot_name}'
                ax.set_xlabel(f'{transform_method}(Intensity)')
                if transform_method == "log10":
                    ax.set_xlabel(f'log10(Intensity + 1e-6)')
        # Plot multiple histograms with different colors, transparency, and labels
        sns.histplot(transformed_data, color=color, alpha=0.5, bins=75, kde=True, ax=ax, label=label)
    
    ax.set_title(plot_title)

    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f'{plot_name}.png'), dpi=300)
    plt.close()

def safe_transform(data, method="arcsinh", scale=10):
    if isinstance(data, dict):
        return {k: safe_transform(v, method, scale) for k, v in data.items()}
    if method == "arcsinh":
        return np.arcsinh(scale * data)
    elif method == "log10":
        return np.log10(data + EPSILON)  # Pseudocount
    elif method == "sqrt":
        return np.sqrt(data)
    else:
        raise ValueError("Choose: arcsinh/log10/sqrt")

def create_single_intensity_histogram(data_raw, label, plot_name, column, get_sattelpunkt=True, tick_space = 0.5, transform_method = "log10", show_thresholds = True):
    print(column)
    try:
        title = f'{plot_name}.png'
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), gridspec_kw={'height_ratios': [2, 1]})
        plot_title = f'{column} Distribution'
        
        if transform_method: 
            # Data processing
            
            transformed_data = safe_transform(data_raw, transform_method)

            title = f'{transform_method}_{title}'
            ax2.set_xlabel(f'{transform_method}(Intensity)')
            if transform_method == "log10":
                ax2.set_xlabel(f'log10(Intensity + 1e-6)')
        else:
            transformed_data = data_raw

            
        if show_thresholds:

            thresholds = calculate_thresholds(transformed_data.to_numpy())

            print(thresholds)
            for method, value in thresholds.items():
                ax1.axvline(value, color='black', linestyle='--', linewidth=1, alpha=0.7)
                ax1.text(value, ax1.get_ylim()[1]*0.9, f'{method}: {value:.3f}', 
                        rotation=90, va='top', ha='right', fontsize=8, backgroundcolor='white')
        
        if len(transformed_data) == 0:
            print(f"Skipping {plot_name} - no valid data")
            return
            
        # --- KDE Calculation ---
        kde = gaussian_kde(transformed_data)
        x = np.linspace(transformed_data.min(), transformed_data.max(), 1000)
        y = kde(x) * len(transformed_data) * (x[1]-x[0])
        
        # --- Derivative Analysis ---
        dy = np.gradient(y, x)
        dy_smooth = savgol_filter(dy, 51, 3)

        d2y = np.gradient(dy_smooth, x)
        
        # Find zero-crossings where derivative changes from negative to positive
        zero_crossings = np.where((np.diff(np.sign(dy_smooth)) > 0))[0]
       
       # 2. Finde Nullstellen der zweiten Ableitung (Hochpunkte in erster Ableitung)
        hochpunkte = np.where(np.diff(np.sign(d2y)))[0]
        
        # If we found any zero crossings, select the most prominent one (highest y value)
        if len(zero_crossings) > 0 and get_sattelpunkt:
            main_zero_crossing = zero_crossings[np.argmax(y[zero_crossings])]
            marked_x = x[main_zero_crossing]
            zero_crossing_y = y[main_zero_crossing]
        else:
            marked_x = None
            
        # --- Top Plot ---
        sns.histplot(transformed_data, color='red', alpha=0.5, bins=75, kde=True, ax=ax1, label=label)

        # Mark the zero crossing point if found
        if marked_x is not None:
            ax1.scatter(marked_x, zero_crossing_y, color='blue', s=100, 
                        label=f'Local Min (x={marked_x:.2f})')
        ax1.set_xlabel('')
        ax1.set_xticks(np.arange(np.floor(x.min()), np.ceil(x.max()) + tick_space, tick_space))
        ax1.set_xticks(np.arange(np.floor(x.min()), np.ceil(x.max()) + 0.01, 0.01), minor=True)
        ax1.set_title(plot_title)
        ax1.legend()
        
        # --- Bottom Plot ---
        ax2.plot(x, dy_smooth, color='purple', label='First Derivative')
        ax2.axhline(0, color='gray', linestyle=':')
        
        # Mark the zero crossing point if found
        if marked_x is not None:
            ax2.axvline(marked_x, color='blue', linestyle='--')
            ax2.scatter(marked_x, 0, color='blue', s=100)
            ax2.text(marked_x, ax2.get_ylim()[0]*0.9, f'{marked_x:.2f}', 
                    color='blue', ha='center', va='center')
        
            
        ax2.set_xticks(np.arange(np.floor(x.min()), np.ceil(x.max()) + tick_space, tick_space))
        ax2.set_xticks(np.arange(np.floor(x.min()), np.ceil(x.max()) + 0.01, 0.01), minor=True)
        ax2.set_ylabel('Slope')
        ax2.legend()

        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, title), dpi=300)
        plt.close()

        
    except Exception as e:
        print(f"Error in {column}: {str(e)}")
        plt.close()

def create_all_histograms(df, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    
    create_single_intensity_histogram(df["panCK cMean"],"All detected cells", "panCK", "panCK", True, 0.5)
    #create_single_intensity_histogram(df["panCK cMean"], None, "All detected cells", "panCK", "panCK", True, 0.5, None)

    
    pancK_pos_cells_df = df[df["panCK cMean"] >  PANCK_THRESHOLD]
    pancK_neg_cells_df = df[df["panCK cMean"] <= PANCK_THRESHOLD]

    create_single_intensity_histogram(pancK_pos_cells_df["KRT56 cMean"], f"panCK+ cells (> {PANCK_THRESHOLD})", "KRT56_in_panCKpos", "KRT56 in panCK+ cells", 0.25)
    create_mutiple_dist_histogram([df["KRT56 cMean"], pancK_pos_cells_df["KRT56 cMean"], pancK_neg_cells_df["KRT56 cMean"]], ["All detected cells", f"panCK+ cells (> {PANCK_THRESHOLD})", f"panCK- cells (<= {PANCK_THRESHOLD})"], ["yellow", "blue", "red"], "KRT56_in_all_panCK_pos_neg", "KRT56 Mean Intesities per cell")
    
    #panCK_pos_KRT56_pos = pancK_pos_cells_df[pancK_pos_cells_df["KRT56 cMean"] > KRT56_THRESHOLD ]
    #panCK_pos_KRT56_neg = pancK_pos_cells_df[pancK_pos_cells_df["KRT56 cMean"] <= KRT56_THRESHOLD ]
    
    create_single_intensity_histogram(pancK_pos_cells_df["HNF4A nMean"], f"panCK+ cells (> {PANCK_THRESHOLD})", "HNF4A_in_panCKpos", "HNF4A in panCK+ cells", True, 0.25)
    
    create_single_intensity_histogram(pancK_pos_cells_df["KLF5 nMean"], f"panCK+ cells (> {PANCK_THRESHOLD})", "KLF5_in_panCKpos", "KLF5 in panCK+ cells", False, 0.25)

    create_single_intensity_histogram(pancK_pos_cells_df["Bcl-xL nMean"], f"panCK+ cells (> {PANCK_THRESHOLD})", "Bcl-xL_in_panCKpos", "Bcl-xL in panCK+ cells", True, 0.25)


if __name__ == "__main__":
    results_folder = r""
    output_folder = r""
    
    # Get the combined DataFrames
    all_cells_df = process_all_files(results_folder)
    
    # Create and save histograms
    create_all_histograms(all_cells_df, output_folder)
    
    # Print summary information
    print("\nDataFrame:")
    print(f"Total rows: {len(all_cells_df)}")
    print(f"Columns: {list(all_cells_df.columns)}")
    
