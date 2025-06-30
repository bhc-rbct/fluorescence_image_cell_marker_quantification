


from datetime import date
import os

import pandas as pd

PANCK_THRESHOLD = 0.026829
KRT56_THRESHOLD = 0.0770
KLF5_THRESHOLD = 0.076599

HNF4A_THRESHOLD = 0.0786
BCLXL_THRESHOLD = 0.0794318


def get_all_data_per_image(folder_path):
    whole_image_df = pd.DataFrame()
    if os.path.isdir(folder_path):
        all_tiles_dfs = []
        # Get all files in the quantifications folder
        files = os.listdir(folder_path)

        all_cells_files = [f for f in files 
                            if f.startswith('Tissue_') and f.endswith('_allCells_allMarkers.csv')]
        
        # Process positive files
        for f in all_cells_files:
            file_path = os.path.join(folder_path, f)
            df = pd.read_csv(file_path)

            all_tiles_dfs.append(df)

        whole_image_df = pd.concat(all_tiles_dfs, ignore_index=True) if all_tiles_dfs else pd.DataFrame()
    return whole_image_df

def create_summary_csv(folder_path, output_csv_path):
    # Load all data
    whole_image_df = get_all_data_per_image(folder_path)
    
    # Create a copy of panCK positive cells to avoid SettingWithCopyWarning
    panCK_pos = whole_image_df[whole_image_df["panCK cMean"] > PANCK_THRESHOLD].copy()
    
    # Add new columns using .loc to avoid warnings
    panCK_pos.loc[:, 'KRT56_pos'] = panCK_pos['KRT56 cMean'] > KRT56_THRESHOLD
    panCK_pos.loc[:, 'HNF4A_pos'] = panCK_pos['HNF4A nMean'] > HNF4A_THRESHOLD
    panCK_pos.loc[:, 'KLF5_pos'] = panCK_pos['KLF5 nMean'] > KLF5_THRESHOLD
    panCK_pos.loc[:, 'Bcl_xL_pos'] = panCK_pos['Bcl-xL nMean'] > BCLXL_THRESHOLD
    
    # Calculate all required counts
    summary_data = {
        'Total_cells': [len(whole_image_df)],
        'Total_panCK_pos': [len(panCK_pos)],
        'KLF5_Bcl_xL_double_pos': [len(panCK_pos[panCK_pos['KLF5_pos'] & panCK_pos['Bcl_xL_pos']])],
        'KLF5_single_pos': [len(panCK_pos[panCK_pos['KLF5_pos'] & ~panCK_pos['Bcl_xL_pos']])],
        'Bcl_xL_single_pos': [len(panCK_pos[~panCK_pos['KLF5_pos'] & panCK_pos['Bcl_xL_pos']])],
        
        'KRT56_pos_HNF4A_neg': [len(panCK_pos[panCK_pos['KRT56_pos'] & ~panCK_pos['HNF4A_pos']])],
        'KRT56_neg_HNF4A_pos': [len(panCK_pos[~panCK_pos['KRT56_pos'] & panCK_pos['HNF4A_pos']])],
        
        # Subgroup counts
        'KRT56_pos_HNF4A_neg_KLF5_Bcl_xL_double_pos': [
            len(panCK_pos[panCK_pos['KRT56_pos'] & ~panCK_pos['HNF4A_pos'] & 
                         panCK_pos['KLF5_pos'] & panCK_pos['Bcl_xL_pos']])
        ],
        'KRT56_pos_HNF4A_neg_KLF5_single_pos': [
            len(panCK_pos[panCK_pos['KRT56_pos'] & ~panCK_pos['HNF4A_pos'] & 
                         panCK_pos['KLF5_pos'] & ~panCK_pos['Bcl_xL_pos']])
        ],
        'KRT56_pos_HNF4A_neg_Bcl_xL_single_pos': [
            len(panCK_pos[panCK_pos['KRT56_pos'] & ~panCK_pos['HNF4A_pos'] & 
                         ~panCK_pos['KLF5_pos'] & panCK_pos['Bcl_xL_pos']])
        ],
        
        'KRT56_neg_HNF4A_pos_KLF5_Bcl_xL_double_pos': [
            len(panCK_pos[~panCK_pos['KRT56_pos'] & panCK_pos['HNF4A_pos'] & 
                         panCK_pos['KLF5_pos'] & panCK_pos['Bcl_xL_pos']])
        ],
        'KRT56_neg_HNF4A_pos_KLF5_single_pos': [
            len(panCK_pos[~panCK_pos['KRT56_pos'] & panCK_pos['HNF4A_pos'] & 
                         panCK_pos['KLF5_pos'] & ~panCK_pos['Bcl_xL_pos']])
        ],
        'KRT56_neg_HNF4A_pos_Bcl_xL_single_pos': [
            len(panCK_pos[~panCK_pos['KRT56_pos'] & panCK_pos['HNF4A_pos'] & 
                         ~panCK_pos['KLF5_pos'] & panCK_pos['Bcl_xL_pos']])
        ]
    }
    
    # Create DataFrame and save to CSV
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_csv_path, index=False)
    return summary_df

def process_all_files(results_folder):
    for experiment_folder in os.listdir(results_folder):
        quant_folder = os.path.join(results_folder, experiment_folder, "quantifications")
        
        if os.path.isdir(quant_folder):
            print(f"Processing folder: {quant_folder}")
            # Create output path for this experiment's summary
            current_date = date.today().strftime('%Y%m%d')  # Changed to YYYYMMDD format
            output_csv = os.path.join(results_folder, experiment_folder, f"{current_date}_{experiment_folder}_summary.csv")
            create_summary_csv(quant_folder, output_csv)
        else: 
            print(f"{quant_folder} is not a folder")


process_all_files(r'')
