## Visualizing spatial population structure with estimated effective migration surfaces(EEMS)

**## a visual representation of population structure that can highlight potential regions of higher-than-average and lower-than-average historic gene flow.**

##Reference:(Petkova et al. 2016) https://www.nature.com/articles/ng.3464

##Authors: Bhuwan Singh Bist
##Date:10/14/2024

# Step 1: Load required packages in R
# Install, and load 'vcfR' for reading the VCF file
```
if (!require("vcfR")) install.packages("vcfR", dependencies=TRUE)
library(vcfR)
```
# Step 2: Read in the VCF file (replace with your actual VCF file path)
```
vcf_data <- read.vcfR("Armadillos_filtered.DP5g50maf05.recode.vcf")
```

# Step 3: Extract genotype matrix and sample coordinates for EEMS input
# This generates the genetic differences matrix that will be used in EEMS
```
library(adegenet)  # Use this to process genotypes if needed
genind_obj <- vcfR2genind(vcf_data)
dist_matrix <- dist(tab(genind_obj))
```
## Save distance matrix in the format required by EEMS
write.table(as.matrix(dist_matrix), "output_path/output.diffs", col.names=FALSE, row.names=FALSE)

# Step 4: Prepare coordinates file for EEMS
# Assuming you have a file with the sample coordinates in this format: longitude latitude (one sample per line)
# Replace 'coords_file.txt' with your actual file
```
coords <- read.table("coords_file.txt", header=FALSE)
write.table(coords, "output_path/output.coord", col.names=FALSE, row.names=FALSE)
```
# Step 5: Prepare habitat outer boundary (polygon file)
# EEMS requires the coordinates outlining the habitat as a closed polygon
# Replace 'polygon_file.txt' with your actual habitat boundary file
```
polygon <- read.table("polygon_file.txt", header=FALSE)
write.table(polygon, "output_path/output.outer", col.names=FALSE, row.names=FALSE)
```
# Step 6: Update params.ini file for EEMS
# Ensure your params.ini file is properly configured with:
# - nIndiv = 88 (for 88 samples)
# - nDemes = 50 (for moderate grid resolution)
# - datapath = 'output_path/output'
# - mcmcpath = 'output_path/output-results'
# Example content of params.ini:
# datapath = output_path/output
# mcmcpath = output_path/output-results
# nIndiv = 88
# nSites = (number of SNPs in the VCF)
# nDemes = 50
# diploid = true
# numMCMCIter = 2000000
# numBurnIter = 1000000
# numThinIter = 9999

# Step 7: Run EEMS in the Linux cluster (from the command line)
# After preparing your files, run EEMS with the following command:
```
./runeems_snps --params params.ini --seed 123
```
# Explanation:
# ./runeems_snps        : EEMS executable for SNP data
# --params params.ini   : Path to your EEMS configuration file (params.ini)
# --seed 123            : Optional seed for reproducibility

# Step 8: Install and load rEEMSplots for visualizing EEMS output
```
if (!require("rEEMSplots")) install.packages("rEEMSplots", repos = NULL, type = "source")
```
# Step 9: Generate visualizations for migration surface
```
library(rEEMSplots)
mcmcpath <- "output_path/output-results"
plotpath <- "output_path/plots"
eems.plots(mcmcpath, plotpath, longlat = TRUE)
```
# Step 10: Reproduce example plots using ggplot2
# Load data saved by rEEMSplots and create custom scatter plots
```
load(paste0(plotpath, "-rdist.RData"))
```

# Reproduce plot showing observed vs fitted dissimilarities between demes
```
library(ggplot2)
ggplot(B.component %>% filter(size > 1), aes(fitted, obsrvd)) +
  geom_point() +
  labs(title = "Observed vs Fitted Dissimilarities Between Demes")
```
# Reproduce plot showing observed vs fitted dissimilarities within demes
```
ggplot(W.component %>% filter(size > 1), aes(fitted, obsrvd)) +
  geom_point() +
  labs(title = "Observed vs Fitted Dissimilarities Within Demes")
```
# Step 11: Optional - Interactive Shiny app
# Visualize pairwise dissimilarities and locations on the map using a Shiny app.


