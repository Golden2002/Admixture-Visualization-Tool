# Admixture-Visualization-Tool
Multi-K ADMIXTURE visualization pipeline with hierarchical clustering and cross-K ancestry mapping, supporting partial topological sorting and ancestry component alignment.

---

# Multi-K ADMIXTURE Plot Script - Sanitized and Fully Documented

## Description

This R script generates admixture plots for multiple K values, mapping ancestry components across different K values to maintain consistency. It uses hierarchical clustering to order populations, a partial topological sorting algorithm to handle component swaps with potential cycles, and a combination of direct renaming and pairwise exchanges to ensure ancestry consistency across K.

## Features

* Reads ADMIXTURE Q files and sample information.
* Supports multiple K values.
* Maps ancestry components from lower K to the reference K (maximum K).
* Handles cycles in component swaps using partial topological sorting.
* Orders populations using hierarchical clustering.
* Generates aligned bar plots for each K and side plots for target populations.
* Fully customizable color palette.

## Input Files

* `.fam` file containing sample identifiers.
* `.info` file containing sample metadata (including population labels).
* ADMIXTURE `.Q` files for each K in the specified range.

## Output

* A multi-K ADMIXTURE plot in PNG format.
* Plots include:

  * Left: full sample ancestry proportion bars.
  * Right: average ancestry proportions for specified target populations.

## Usage

```R
# Set directories and file paths
run_dir <- "YOUR_RUN_DIRECTORY"
fam_file <- file.path(run_dir, "samples.fam")
info_file <- file.path(run_dir, "sample_info.txt")
K_range <- 3:7
target_pops <- c("Pop1", "Pop2", "Pop3")
output_png <- "Admixture_Plot_K3-7.png"

# Read and process data
admix_data <- read_admixture_data(run_dir, fam_file, info_file, K_range)

# Perform population clustering
admix_data <- cluster_populations(admix_data, K_range)

# Map ancestry components across K
admix_mapped <- map_ancestries(admix_data, K_range)

# Create plots
main_plots <- create_main_plots(admix_mapped, color_dict, K_range)
pop_axis <- create_population_axis(sample_order_df)
target_plots <- create_target_population_plots(admix_mapped, target_pops, color_dict, K_range)

# Save final plot
final_plot <- combine_plots(main_plots, pop_axis, target_plots)
ggsave(output_png, final_plot, width=plot_width, height=plot_height, units="in")
```

## Notes

* Hierarchical clustering is performed using `ward.D2` method.
* Color palette can be customized with `RColorBrewer` schemes.
* Partial topological sort ensures safe execution of component swaps even when cycles exist.
* Unmapped components are renamed with prefix `Extra_` to avoid conflicts.
* Compatible with R >= 4.0.

## Example Color Palette

```R
c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F")
```

---

# Detailed Methods
———— Multi-K ADMIXTURE Plotter

A sophisticated R script for creating publication-quality ADMIXTURE plots across multiple K values with consistent ancestry coloring and intelligent component mapping.

## Features

### 🎨 **Visualization Features**
- **Multi-K Comparison**: Simultaneous visualization of ancestry proportions across multiple K values (e.g., K=3-7)
- **Consistent Color Scheme**: Automatic color consistency for the same ancestry components across different K values
- **Dual-Panel Layout**: 
  - Left panel: Complete population structure with hierarchical clustering
  - Right panel: Focused view of target populations
- **Publication-Ready**: High-resolution output with customizable dimensions and fonts

### 🔬 **Advanced Algorithms**
- **Hierarchical Clustering**: Ward's method clustering for optimal population ordering
- **Cross-K Ancestry Mapping**: Intelligent matching of ancestry components across different K values
- **Opposite Logic Exchange**: Sophisticated component swapping with opposite pairing detection
- **Topological Sorting**: Dependency-aware exchange order to handle complex mapping scenarios

### 📊 **Data Processing**
- **Automatic Sample Ordering**: Logical ordering based on population clustering
- **Component Reconciliation**: Handles missing, extra, and unmapped ancestry components
- **Robust Error Handling**: Graceful handling of missing files and edge cases

## Installation

### Prerequisites
- R (≥ 4.0.0)
- Required R packages:
  ```r
  install.packages(c("dplyr", "tidyr", "ggplot2", "patchwork", 
                     "tibble", "RColorBrewer", "igraph"))
  ```

### Quick Start
1. Clone or download this repository
2. Prepare your ADMIXTURE results and sample metadata
3. Modify the paths in the script to point to your data
4. Run the script:
   ```bash
   Rscript admixture_plotter.R
   ```

## Input Files

### Required Files
```
project_directory/
├── admixture_results/
│   ├── input.fam                  # PLINK .fam file with sample IDs
│   ├── sample_info.txt           # Sample metadata (tab-delimited)
│   ├── input.3.Q                 # ADMIXTURE Q file for K=3
│   ├── input.4.Q                 # ADMIXTURE Q file for K=4
│   └── ...                       # Additional Q files for other K values
└── admixture_plotter.R           # This script
```

### File Formats

#### 1. **sample_info.txt** (tab-delimited)
```
sample  pop         region
S001    Han_N       East_Asia
S002    Japanese    East_Asia
S003    Tibetan     Central_Asia
...
```
- **Required columns**: `sample` (matching .fam file), `pop` (population label)
- **Optional columns**: Any additional metadata

#### 2. **ADMIXTURE Q Files**
- Standard ADMIXTURE output format (space-delimited ancestry proportions)
- File naming convention: `input.[K].Q` (adjustable in script)

## Configuration

### Main Parameters (edit in script)
```r
# File paths
run_dir <- "/path/to/your/admixture_results"
fam_file <- file.path(run_dir, "input.fam")
info_file <- file.path(run_dir, "sample_info.txt")

# K-value range
K_range <- 3:7

# Target populations (highlighted in right panel)
target_pops <- c("Han_N", "Japanese", "Tibetan", "Korean", "Mongolian")

# Output settings
output_png <- "Admixture_Plot_K3-7.png"
plot_width <- 16
plot_height_base <- 8
color_palette <- "Paired"  # RColorBrewer palette
```

### Customization Options
- **Color Palette**: Any RColorBrewer palette (Set3, Set2, Accent, Dark2, etc.)
- **Plot Dimensions**: Adjust `plot_width` and `plot_height_base` for different aspect ratios
- **Font Sizes**: Modify `size` parameters in `geom_text()` calls
- **Margins**: Adjust `plot.margin` values in theme settings

## Algorithm Details

### 1. **Population Ordering**
- Uses hierarchical clustering (Ward.D2 method) on population-level ancestry means at the highest K
- Creates biologically meaningful grouping of populations

### 2. **Cross-K Ancestry Mapping**
The core innovation of this script is maintaining color consistency across K values:

#### Step 1: Reference Establishment
- Uses the highest K value as reference
- Identifies representative populations for each ancestry component

#### Step 2: Direct Renaming
When a target ancestry name doesn't exist in the current K:
- Directly rename the component to match the reference

#### Step 3: Component Exchange
When both components exist in the current K:
- **Regular Exchange**: Swap labels between two components
- **Opposite Logic**: If component A has already exchanged with B, and now C wants to exchange with A, exchange C with A's opposite (B) instead
- **Dependency-Aware**: Uses topological sorting to determine safe exchange order

#### Step 4: Unmapped Components
- Label as "Extra_[original_name]" to distinguish from mapped components

### 3. **Safe Exchange Order Generation**
- Builds dependency graph between components
- Uses topological sorting to find acyclic ordering
- Handles cycles by separating cyclic nodes
- Ensures stable and reproducible mapping

## Output

### Primary Output
- **High-resolution PNG**: Publication-ready multi-panel plot
- **File name**: Configurable (default: `Admixture_Plot_K3-7.png`)

### Supplementary Outputs (optional)
- `ancestry_mapping_summary.csv`: Mapping between original and plotted ancestry labels
- `color_mapping.csv`: Color assignments for each ancestry component

### Plot Structure
```
┌─────────────────────────────────────────┬─────────┐
│                                         │         │
│  K=3 ████████████████████████████████  │ Target  │
│                                         │ Pop 1   │
│  K=4 ████████████████████████████████  │ ██████  │
│                                         │         │
│  K=5 ████████████████████████████████  │ Target  │
│                                         │ Pop 2   │
│  K=6 ████████████████████████████████  │ ██████  │
│                                         │         │
│  K=7 ████████████████████████████████  │ Target  │
│                                         │ Pop 3   │
│                                         │ ██████  │
│  Population  Population  Population     │         │
│  Labels      Labels      Labels         │ Labels  │
└─────────────────────────────────────────┴─────────┘
```

## Usage Examples

### Basic Usage
```r
# Edit parameters in script, then:
source("admixture_plotter.R")
```

### Different K Ranges
```r
K_range <- 2:10  # Plot K=2 through K=10
target_pops <- c("European", "African", "Asian", "Native_American")
```

### Different Color Scheme
```r
color_palette <- "Set3"  # 12-color qualitative palette
```

## Troubleshooting

### Common Issues

1. **"File not found" warnings**
   - Check `run_dir` path and Q file naming convention
   - Ensure all K values in `K_range` have corresponding .Q files

2. **Missing target populations**
   - Verify population names in `target_pops` match those in `sample_info.txt`
   - Check for typos or case sensitivity

3. **Color consistency issues**
   - The script uses the highest K as reference; ensure K_max has meaningful ancestry structure
   - Check the ancestry mapping logs printed during execution

4. **Plot labels overlapping**
   - Adjust `plot_width` and `plot_height_base` for better spacing
   - Modify label positioning parameters in `create_population_axis()` and `create_target_population_axis()`

### Debug Mode
The script includes built-in debugging output:
- Component mapping decisions at each step
- Exchange operations performed
- Final mapping relationships for each K

## Citation

If you use this script in your research, please cite:

```
ADMIXTURE Plotter: A tool for consistent multi-K ancestry visualization.
GitHub: https://github.com/yourusername/admixture-plotter
```

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with detailed description

## Acknowledgments

- Built with R and ggplot2
- Inspired by the need for consistent ancestry visualization across K values
- Thanks to the ADMIXTURE developers for the underlying analysis tool

## Contact

For questions, issues, or feature requests:
- Open an issue on GitHub
- Email: your.email@example.com
```

## 3. 其他支持文件

### `LICENSE` (可选)
```text
MIT License

Copyright (c) 2025 ADMIXTURE Plotter Developers

Permission is hereby granted...
```

### `example_config.R` (配置示例)
```r
# Example configuration file
# Copy this to a new file and customize for your project

run_dir <- "/path/to/your/admixture_results"
fam_file <- file.path(run_dir, "your_data.fam")
info_file <- file.path(run_dir, "sample_metadata.txt")

# ADMIXTURE Q file pattern (adjust based on your naming)
Q_file_pattern <- "your_data.%d.Q"  # Will become your_data.3.Q, etc.

K_range <- 3:7
target_pops <- c("Population1", "Population2", "Population3")

output_png <- "My_Admixture_Plot.png"
plot_width <- 16
plot_height_base <- 8
color_palette <- "Set3"
```

### `requirements.txt` (R包依赖)
```r
# Required R packages
# Install with: install.packages(c("dplyr", "tidyr", ...))

dplyr
tidyr
ggplot2
patchwork
tibble
RColorBrewer
igraph
```

## 4. 使用说明总结

1. **下载所有文件**到本地目录
2. **修改配置文件**：
   - 更新`admixture_plotter.R`中的文件路径
   - 设置适当的K范围和目标群体
3. **准备输入文件**：
   - ADMIXTURE的.Q文件
   - 样本元数据文件
   - PLINK .fam文件
4. **运行脚本**：
   ```bash
   Rscript admixture_plotter.R
   ```
5. **检查输出**：
   - 查看生成的PNG图像
   - 检查映射表（如果启用了）

这个完整的GitHub仓库包含：
- 去敏感化的专业脚本
- 详尽的英文注释
- 完整的用户文档
- 示例配置文件
- 依赖说明

这样您的代码就适合分享到GitHub，其他人可以轻松理解和使用您的工具。
