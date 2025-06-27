# Image Analysis Pipeline for Marker Quantification

This pipeline processes multiplex immunofluorescence images to quantify marker expression patterns in epithelial cells. The workflow consists of three main steps:

## 1. Cell Detection and Initial Quantification (QuPath)

**Script:** `quantify_classify.groovy`

### Requirements:
- QuPath (v0.5.0 or later)
- StarDist extension for QuPath
- Pre-trained StarDist model (`dsb2018_heavy_augment.pb`)

### What this script does:
1. Performs nucleus detection using StarDist on DAPI channel
2. Expands nuclei to approximate cell boundaries
3. Measures mean intensity for each marker in each cell
4. Performs initial classification based on panCK and KRT56 expression
5. Exports measurements for all cells to CSV files (organized by tile region)

### Key outputs:
- `Tissue_X_allCells_allMarkers.csv` files containing:
  - Normalized intensity measurements for all markers (between 0 and 1)

## 2. Threshold Determination and Visualization (Python)

**Script:** `plots_threshold_detection.py`

### Requirements:
- Python 3.8+
- Packages: numpy, pandas, matplotlib, seaborn, scipy, scikit-image

### What this script does:
1. Aggregates data from all tile regions
2. Calculates multiple potential thresholds for each marker using:
   - Otsu's method
   - Triangle method
   - K-means clustering
   - Gaussian Mixture Models
3. Generates diagnostic plots showing:
   - Intensity distributions
   - Automated threshold suggestions
   - Derivative analysis to identify inflection points

### Key outputs:
- Histogram plots for each marker (`sqrt_panCK.png`, `sqrt_KRT56_in_panCKpos.png`, etc.)
- Visualizations of different thresholding methods
- Suggested thresholds based on distribution analysis

## 3. Final Classification and Summary Statistics (Python)

**Script:** `summary_file_creation.py`

### Requirements:
- Same as step 2

### What this script does:
1. Applies final thresholds to classify cells as positive/negative for each marker
2. Calculates counts for all biologically relevant combinations:
   - KLF5/Bcl-xL co-expression patterns
   - KRT56/HNF4A mutual exclusion
   - Subpopulation frequencies
3. Generates a comprehensive summary CSV for each image

### Key outputs:
- `YYYYMMDD_[image_name]_summary.csv` containing:
  - Total cell counts per category
  - Subpopulation frequencies
  - Co-expression statistics

## Workflow Order:
1. First run `quantify_classify.groovy` in QuPath
2. Then run `plots_threshold_detection.py` to determine optimal thresholds
3. Finally run `summary_file_creation.py` with your chosen thresholds

## Customization:
- Update paths in the Groovy script for:
  - Output directory (`generalOutputDir`)
  - StarDist model path (`starDistPathModel`)
- Adjust thresholds in both Python scripts as needed based on diagnostic plots
- Modify marker combinations in `summary_file_creation.py` for different biological questions

## Notes:
- The pipeline assumes a specific channel order: DAPI, p40, KLF5, KRT56, HNF4A, panCK, Bcl-xL
- All scripts include extensive comments for parameter adjustment

## License

MIT License

## Contact

For questions or issues, please open an issue on the repository.

*This tool is intended for research use and should be validated for your specific imaging setup and markers.*
