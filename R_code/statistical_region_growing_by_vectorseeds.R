get_neighbors <- function(cell, nrows, ncols, connectivity = 8) {
  row <- ((cell - 1) %/% ncols) + 1
  col <- ((cell - 1) %% ncols) + 1
  
  offsets <- if (connectivity == 4) {
    list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
  } else {
    list(c(-1, -1), c(-1, 0), c(-1, 1),
         c(0, -1),             c(0, 1),
         c(1, -1),  c(1, 0),  c(1, 1))
  }
  
  neighbors <- c()
  for (offset in offsets) {
    new_row <- row + offset[1]
    new_col <- col + offset[2]
    if (new_row >= 1 && new_row <= nrows && new_col >= 1 && new_col <= ncols) {
      neighbors <- c(neighbors, (new_row - 1) * ncols + new_col)
    }
  }
  neighbors
}

regions_to_polygons <- function(region_raster, dissolve = TRUE, output_path = NULL) {
  unique_regions <- unique(values(region_raster))
  unique_regions <- unique_regions[unique_regions > 0 & !is.na(unique_regions)]
  polygons_list <- list()
  
  for (region_id in unique_regions) {
    binary_raster <- region_raster
    values(binary_raster) <- ifelse(values(region_raster) == region_id, 1, NA)
    region_polygons <- as.polygons(binary_raster, dissolve = dissolve)
    if (nrow(region_polygons) > 0) {
      region_polygons$region_id <- region_id
      polygons_list[[length(polygons_list) + 1]] <- region_polygons
    }
  }
  
  if (length(polygons_list) > 0) {
    merged_polygons <- do.call(rbind, polygons_list)
    if (!is.null(output_path)) {
      writeVector(merged_polygons, output_path)
    }
    return(merged_polygons)
  } else {
    return(NULL)
  }
}

seed_based_region_growing <- function(raster, seed_points, 
                                      threshold_value = 1.0,
                                      connectivity = 8,
                                      max_iterations = 1000) {
  
  if (!same.crs(seed_points, raster)) {
    message("CRS mismatch - reprojecting seed points to match raster")
    seed_points <- project(seed_points, raster)
  }
  
  raster_values <- values(raster)
  nrows <- nrow(raster)
  ncols <- ncol(raster)
  
  result_raster <- raster[[1]]
  values(result_raster) <- 0
  seed_cells <- cells(raster, seed_points)
  
  for (seed_idx in 1:length(seed_cells)) {
    seed_cell <- seed_cells[seed_idx]
    if (is.na(seed_cell)) next
    
    seed_values <- if (nlyr(raster) > 1) raster_values[seed_cell, ] else raster_values[seed_cell]
    if (any(is.na(seed_values))) next
    
    seed_neighbors <- get_neighbors(seed_cell, nrows, ncols, connectivity)
    if (nlyr(raster) > 1) {
      neighbor_vals <- raster_values[seed_neighbors, ]
      neighbor_vals <- neighbor_vals[!apply(is.na(neighbor_vals), 1, any), ]
      if (nrow(neighbor_vals) > 1) {
        reference_mean <- colMeans(neighbor_vals, na.rm = TRUE)
        reference_std <- apply(neighbor_vals, 2, sd, na.rm = TRUE)
        reference_std[is.na(reference_std) | reference_std == 0] <- 0.1
      } else {
        reference_mean <- seed_values
        reference_std <- rep(0.1, length(seed_values))
      }
    } else {
      neighbor_vals <- raster_values[seed_neighbors]
      neighbor_vals <- neighbor_vals[!is.na(neighbor_vals)]
      if (length(neighbor_vals) > 1) {
        reference_mean <- mean(neighbor_vals, na.rm = TRUE)
        reference_std <- sd(neighbor_vals, na.rm = TRUE)
        if (is.na(reference_std) || reference_std == 0) reference_std <- 0.1
      } else {
        reference_mean <- seed_values
        reference_std <- 0.1
      }
    }
    
    region <- c(seed_cell)
    checked <- rep(FALSE, ncell(raster))
    checked[seed_cell] <- TRUE
    queue <- c(seed_cell)
    iteration <- 0
    
    while (length(queue) > 0 && iteration < max_iterations) {
      iteration <- iteration + 1
      current_cell <- queue[1]
      queue <- queue[-1]
      neighbors <- get_neighbors(current_cell, nrows, ncols, connectivity)
      
      for (neighbor in neighbors) {
        if (!checked[neighbor] && neighbor > 0 && neighbor <= ncell(raster)) {
          neighbor_values <- if (nlyr(raster) > 1) raster_values[neighbor, ] else raster_values[neighbor]
          if (any(is.na(neighbor_values))) {
            checked[neighbor] <- TRUE
            next
          }
          
          similarity <- if (nlyr(raster) > 1) {
            max(abs((neighbor_values - reference_mean) / reference_std), na.rm = TRUE)
          } else {
            abs((neighbor_values - reference_mean) / reference_std)
          }
          
          if (!is.na(similarity) && similarity <= threshold_value) {
            region <- c(region, neighbor)
            queue <- c(queue, neighbor)
          }
          
          checked[neighbor] <- TRUE
        }
      }
    }
    
    result_values <- values(result_raster)
    result_values[region] <- seed_idx
    values(result_raster) <- result_values
  }
  
  result_raster
}