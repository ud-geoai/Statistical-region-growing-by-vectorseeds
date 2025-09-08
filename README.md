# Statistical region growing by vectorseeds

This R script provides functions to perform seed-based statistical region growing on raster data. It includes utilities for extracting neighboring cells, growing regions from seed points, and converting raster regions to polygons. The script is designed for geospatial analysis using the terra package. The script is designed to help collect reference pixels from rasters around point objects based on the statistical similarity of the neighboring pixel values. The function calculates the relative standard deviation and stops at a user-defined threshold, enabling automated selection of spectrally homogeneous areas for analysis. Accordingly, the function incorporates new pixels step-by-step and calculates the mean ± z × SD values where the "z" defines the threshold: when z=1, the pixel values have to be within 68.2% range around the mean; when z=0.5, the range is 38.3%. This iterative approach ensures consistent and reliable pixel sampling for remote sensing applications and spatial analysis.

## Contents

- **get_neighbors** – Returns neighboring cell indices for a raster cell.
- **seed_based_region_growing** – Performs statistical region growing using seed points.
- **regions_to_polygons** – Converts raster regions into polygon features, optionally saving them to a file.

## Requirements

- **R (≥ 4.0 recommended)**
- **terra package**

---

### Function: get_neighbors

**Description**\
*Determines the neighboring cell indices for a given cell within a raster grid, based on either 4-connectivity or 8-connectivity.*

**Usage**
> get_neighbors(cell, nrows, ncols, connectivity = 8)

**Arguments**

> - cell (integer): The index of the target cell.
> - nrows (integer): Total number of rows in the raster grid.
> - ncols (integer): Total number of columns in the raster grid.
> - connectivity (integer): Neighborhood type; 4 for orthogonal neighbors only, 8 to include diagonal neighbors. Defaults to 8.

**Value**\
*Returns an integer vector containing the cell indices of valid neighbors.*

### Function: seed_based_region_growing

**Description**\
*Performs statistical region growing on a raster using vector seed points. Regions are expanded iteratively based on pixel similarity to a local reference (seed neighborhood), measured in standard deviations.*

**Usage**
> seed_based_region_growing(raster, seed_points, threshold_value = 1.0, connectivity = 8, max_iterations = 1000)

**Arguments**

> - raster (SpatRaster): Input raster to segment.
> - seed_points (SpatVector): Vector points defining initial seeds for region growing.
> - threshold_value (numeric): Maximum allowed z-score for including neighboring pixels. Defaults to 1.0.
> - connectivity (integer): 4 or 8 connectivity for neighbor evaluation. Defaults to 8.
> - max_iterations (integer): Maximum number of iterations per region. Defaults to 1000.

**Value**\
*A raster with the same dimensions as the input. Each cell that belongs to a grown region is labeled with the index of the seed point that generated it. Cells not assigned to any region are 0.*

### Function: region_to_polygons

**Description**\
*Converts labeled raster regions into polygon features. Optionally merges adjacent cells of the same region (dissolve) and saves the result to a shapefile.*

**Usage**
> regions_to_polygons(region_raster, dissolve = TRUE, output_path = NULL)

**Arguments**

> - region_raster (SpatRaster): Labeled raster where each region has a unique integer value.
> - dissolve (logical): If TRUE, merges contiguous cells of the same region into single polygons. Defaults to TRUE.
> - output_path (character): Optional path to save the polygons as a shapefile. Defaults to NULL (no file saved).

**Value**\
*A SpatVector containing polygon features with an attribute _region_id_ corresponding to the raster region. Returns NULL if no regions are found.*

---

# Complete workflow example

```
# Load the required package
source("statistical_region_growing_by_vectorseeds.R")
library(terra)

# Step 1: Load raster and seed points
raster_path <- "example_raster.tif"
seed_points_path <- "seed_points.shp"

r <- rast(raster_path)
seeds <- vect(seed_points_path)

# Step 2: Perform seed-based statistical region growing
result_raster <- seed_based_region_growing(
	raster = r,
	seed_points = seeds,
	threshold_value = 1.0,
	connectivity = 8,
	max_iterations = 1000
)

# Step 3: Convert the raster regions into polygons
regions_polygons <- regions_to_polygons(
	region_raster = result_raster,
	dissolve = TRUE,
	output_path = "regions.shp"
)

# Step 4: Plot the polygons
plot(regions_polygons)

```
