# ============================================
# Multi-K ADMIXTURE Plotting Script
# Version: 1.2.0
# Last Updated: 2025-12-30
# ============================================

# =========================
# 0. Load Required Libraries
# =========================
suppressPackageStartupMessages({
  library(dplyr)      # Data manipulation
  library(tidyr)      # Data reshaping
  library(ggplot2)    # Plotting
  library(patchwork)  # Plot composition
  library(tibble)     # Modern data frames
  library(RColorBrewer) # Color schemes
  library(igraph)     # Graph theory for topological sorting
})

# =========================
# 1. Global Parameter Settings
# =========================

# NOTE: Replace these paths with your own data directory structure
# Example structure:
# project_dir/
#   ├── admixture_results/
#   │   ├── input.fam
#   │   ├── sample_info.txt
#   │   └── *.Q files (K3.Q, K4.Q, etc.)
#   └── admixture_plotter.R

run_dir <- "/path/to/your/admixture_results"  # Directory containing ADMIXTURE results
fam_file <- file.path(run_dir, "input.fam")   # PLINK .fam file with sample IDs
info_file <- file.path(run_dir, "sample_info.txt")  # Sample metadata (sample, pop, etc.)

# K-value range to plot (e.g., K=3 to K=7)
K_range <- 3:7

# Target populations to highlight in the right panel
# These should match population names in your sample_info.txt
target_pops <- c("Population1", "Population2", "Population3", "Population4", "Population5")

# Output settings
output_png <- "Admixture_Plot_K3-7.png"
plot_width <- 16                    # Plot width in inches
plot_height_base <- 8               # Base height for 3 K values
height_per_K <- 1.2                 # Additional height per K value
color_palette <- "Paired"           # RColorBrewer palette name

# =========================
# 2. Data Reading Function
# =========================
read_admixture_data <- function(run_dir, fam_file, info_file, K_range) {
  # Read PLINK .fam file to get sample IDs
  fam <- read.table(fam_file, stringsAsFactors = FALSE)
  # Read sample metadata (must contain 'sample' and 'pop' columns)
  sample_info <- read.table(info_file, header = TRUE, stringsAsFactors = FALSE)
  samples <- fam$V1  # Sample IDs from the first column of .fam file
  
  all_data <- list()
  
  for (k in K_range) {
    # Construct ADMIXTURE Q file name (adjust pattern as needed)
    Q_file <- file.path(run_dir, sprintf("input.%d.Q", k))
    
    if (!file.exists(Q_file)) {
      warning(sprintf("File not found: %s", Q_file))
      next
    }
    
    # Read ancestry proportions
    Q <- read.table(Q_file)
    colnames(Q) <- sprintf("A%02d", 1:k)  # Format as A01, A02, etc.
    Q$sample <- samples  # Add sample IDs
    
    # Merge with sample metadata
    Q_merged <- Q %>%
      left_join(sample_info, by = "sample")
    
    # Reshape from wide to long format
    Q_long <- Q_merged %>%
      pivot_longer(
        cols = starts_with("A"),
        names_to = "Ancestry",
        values_to = "Proportion"
      ) %>%
      mutate(K = factor(sprintf("K=%d", k), 
                        levels = sprintf("K=%d", K_range)))
    
    all_data[[as.character(k)]] <- Q_long
  }
  
  # Combine all K values into single data frame
  bind_rows(all_data)
}

# Read the ADMIXTURE data
admix_data <- read_admixture_data(run_dir, fam_file, info_file, K_range)

# =========================
# 3. Population Clustering for Ordering
# =========================

# Use the highest K value as reference for population ordering
K_ref <- max(K_range)
ref_K_label <- sprintf("K=%d", K_ref)
ref_data <- admix_data %>% filter(K == ref_K_label)

# Create population-by-ancestry matrix for hierarchical clustering
ref_matrix <- ref_data %>%
  group_by(pop, Ancestry) %>%
  summarise(mean_prop = mean(Proportion), .groups = "drop") %>%
  pivot_wider(
    names_from = Ancestry,
    values_from = mean_prop,
    values_fill = 0
  ) %>%
  column_to_rownames("pop")

# Perform hierarchical clustering using Ward's method
hc <- hclust(dist(ref_matrix), method = "ward.D2")
pop_order <- hc$labels[hc$order]
cat("Population hierarchical clustering order:\n")
print(pop_order)

# Apply the population order factor
admix_data$pop <- factor(admix_data$pop, levels = pop_order)

# =========================
# 4. Sample Ordering and Coordinate Assignment
# =========================

# Create a consistent sample order based on population clustering
sample_order_df <- admix_data %>%
  distinct(sample, pop) %>%
  arrange(pop) %>%
  mutate(sample_order = row_number())

# Add sample order to main data
admix_data <- admix_data %>%
  left_join(sample_order_df, by = c("sample", "pop"))

# =========================
# 5. Cross-K Ancestry Component Mapping
# =========================

# ------------------------------------------------------------
# 5.1 Build Reference Ancestry Entities (based on K_max)
# ------------------------------------------------------------
K_max <- max(K_range)
ref_K_label <- sprintf("K=%d", K_max)

# Use the highest K as reference for ancestry labeling
ref_data <- admix_data %>% filter(K == ref_K_label)

# Identify representative population for each ancestry component
# (the population with the highest proportion of that ancestry)
ref_top <- ref_data %>%
  group_by(pop, Ancestry) %>%
  summarise(mean_prop = mean(Proportion), .groups = "drop") %>%
  group_by(pop) %>%
  slice_max(mean_prop, n = 1, with_ties = FALSE)

ancestry_representative <- ref_top %>%
  group_by(Ancestry) %>%
  slice_max(mean_prop, n = 1, with_ties = FALSE) %>%
  select(Ancestry, rep_pop = pop, rep_prop = mean_prop)

cat("Reference K (K=", K_max, ") ancestry component representatives:\n")
print(ancestry_representative)

# Save original data for mapping
admix_original <- admix_data

# Create mapped data frame with initial identical ancestry labels
admix_mapped <- admix_data %>%
  mutate(Ancestry_plot = Ancestry)

# ------------------------------------------------------------
# 5.2 Process Each K Value Independently (K < K_max)
# ------------------------------------------------------------
for (k_value in K_range[K_range < K_max]) {
  cat("\n", strrep("=", 60), "\n")
  cat("Processing K =", k_value, "\n")
  cat(strrep("=", 60), "\n")
  
  # Extract data for current K value
  dat <- admix_original %>% filter(K == paste0("K=", k_value))
  orig_components <- unique(dat$Ancestry)
  
  # Initialize mapping: each component maps to itself initially
  final_map <- setNames(orig_components, orig_components)
  
  # ------------------------------------------------------------
  # 5.3 Build component_mapping for current K
  # ------------------------------------------------------------
  component_mapping <- data.frame()
  for (i in seq_len(nrow(ancestry_representative))) {
    rep_pop <- ancestry_representative$rep_pop[i]
    ref_anc <- ancestry_representative$Ancestry[i]
    
    # Find the dominant ancestry component in the representative population
    top_in_pop <- dat %>%
      filter(pop == rep_pop) %>%
      group_by(Ancestry) %>%
      summarise(mean_prop = mean(Proportion), .groups = "drop") %>%
      slice_max(mean_prop, n = 1, with_ties = FALSE)
    
    if (nrow(top_in_pop) > 0) {
      component_mapping <- bind_rows(component_mapping, data.frame(
        current_anc = top_in_pop$Ancestry[1],
        ref_anc = ref_anc,
        rep_pop = rep_pop,
        mean_prop = top_in_pop$mean_prop[1]
      ))
    }
  }
  
  # For each current_anc, keep only the "best" mapping (highest mean_prop)
  component_mapping <- component_mapping %>%
    group_by(current_anc) %>%
    slice_max(mean_prop, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # ------------------------------------------------------------
  # 5.4 Step 1: Direct Renaming (when target doesn't exist)
  # ------------------------------------------------------------
  if (nrow(component_mapping) > 0) {
    cat("\n[Step 1] Direct renaming (when target name doesn't exist):\n")
    for (i in seq_len(nrow(component_mapping))) {
      current_anc <- component_mapping$current_anc[i]
      ref_anc <- component_mapping$ref_anc[i]
      
      if (!(ref_anc %in% orig_components)) {
        final_map[current_anc] <- ref_anc
        cat(sprintf("  ✓ Original %s -> Final %s (%.2f%% in %s)\n",
                    current_anc, ref_anc,
                    component_mapping$mean_prop[i] * 100,
                    component_mapping$rep_pop[i]))
      }
    }
  }
  
  # ------------------------------------------------------------
  # 5.5 Build Safe Exchange Order Using Topological Sorting
  # ------------------------------------------------------------
  
  # Components already renamed in Step 1
  already_renamed <- names(final_map)[final_map != names(final_map)]
  
  # Components that still need swapping
  swap_needs <- component_mapping %>%
    filter(current_anc != ref_anc & 
             ref_anc %in% orig_components & 
             !(current_anc %in% already_renamed))
  
  # Generate safe exchange order using topological sorting
  if (nrow(swap_needs) == 0) {
    safe_order <- character(0)
  } else {
    # Build dependency graph
    edges <- data.frame(from=character(), to=character(), stringsAsFactors = FALSE)
    for (i in seq_len(nrow(swap_needs))) {
      a <- swap_needs$current_anc[i]
      b <- swap_needs$ref_anc[i]
      if (a %in% swap_needs$ref_anc[-i]) {
        # b depends on a, process a before b
        edges <- rbind(edges, data.frame(from=a, to=b, stringsAsFactors = FALSE))
      } else {
        edges <- rbind(edges, data.frame(from=b, to=a, stringsAsFactors = FALSE))
      }
    }
    
    # Create directed graph
    g <- graph_from_data_frame(edges, directed = TRUE)
    
    # Find minimum feedback arc set and remove cycles
    fb_edges <- feedback_arc_set(g)
    g_acyclic <- delete_edges(g, fb_edges)
    
    # Topological sort for acyclic part
    topo_order <- names(topo_sort(g_acyclic, mode = "out"))
    
    # Nodes in cycles = all nodes - topologically sortable nodes
    all_nodes <- swap_needs$current_anc
    cycle_nodes <- setdiff(all_nodes, topo_order)
    
    # Order cycle nodes by original order
    cycle_order <- all_nodes[all_nodes %in% cycle_nodes]
    
    # Final safe order = topological order + cycle nodes
    safe_order <- c(topo_order, cycle_order)
  }
  
  cat("Generated safe exchange order:\n")
  print(safe_order)
  
  # ------------------------------------------------------------
  # 5.6 Step 2: Component Exchange (with opposite logic)
  # ------------------------------------------------------------
  if (nrow(swap_needs) > 0) {
    cat("\n[Step 2] Processing component exchanges (with opposite logic):\n")
    
    pair_map <- list()      # Track pairing relationships a <-> b
    processed_pairs <- list()  # Track processed exchange pairs
    
    # Process in safe order to avoid dependency issues
    for (a_orig in safe_order) {
      
      # Find the row index for this component
      i <- which(swap_needs$current_anc == a_orig)
      if (length(i) == 0) next
      
      b_ref <- swap_needs$ref_anc[i]
      a <- final_map[a_orig]
      b <- b_ref
      
      # Check if components are already paired
      a_paired <- !is.null(pair_map[[a_orig]])
      b_paired <- !is.null(pair_map[[b]])
      same_pair <- FALSE
      
      if (a_paired && b_paired) {
        same_pair <- pair_map[[a_orig]] == b || pair_map[[b]] == a_orig
      }
      
      # Case 1: Both already paired and are the same pair → already exchanged
      if (a_paired && b_paired && same_pair) {
        cat("  ⓘ Already exchanged, skipping\n")
        next
      }
      
      cat(sprintf("\nExchange need %d: Original %s (%s) -> Reference %s\n",
                  i, a_orig, a, b))
      
      # Find which original component currently has label b
      target_owner <- NULL
      for (orig_comp in orig_components) {
        if (final_map[orig_comp] == b) {
          target_owner <- orig_comp
          break
        }
      }
      
      if (is.null(target_owner)) {
        # Current K doesn't have target component → unilateral renaming
        cat(sprintf("  ⓘ Target component not in current K, renaming: %s -> %s\n", a_orig, b))
        final_map[a_orig] <- b
        pair_map[[b]] <- a_orig
        next
      }
      
      # Check if b has already been exchanged with another component
      b_paired <- pair_map[[b]]
      if (!is.null(b_paired)) {
        cat(sprintf("  ⓘ %s already exchanged with %s, using opposite\n", b, b_paired))
        b_owner <- NULL
        for (orig_comp in orig_components) {
          if (final_map[orig_comp] == b_paired) {
            b_owner <- orig_comp
            break
          }
        }
        if (is.null(b_owner)) {
          cat("  ✗ Cannot find opposite of b, skipping\n")
          next
        }
        b <- final_map[b_owner]
      } else {
        b_owner <- target_owner
      }
      
      # Check if a_orig has already participated in an exchange
      if (!is.null(pair_map[[a_orig]])) {
        # a_orig already exchanged with another component
        a_opposite <- pair_map[[a_orig]]
        cat(sprintf("  ⓘ %s already exchanged with %s, letting %s exchange with %s's opposite\n",
                    a_orig, a_opposite, b, a_orig))
        
        # Find component currently labeled as a_opposite
        a_opposite_owner <- NULL
        for (orig_comp in orig_components) {
          if (final_map[orig_comp] == a_opposite) {
            a_opposite_owner <- orig_comp
            break
          }
        }
        
        if (is.null(a_opposite_owner)) {
          cat("  ✗ Cannot find opposite of a_orig, skipping\n")
          next
        }
        
        # Check if b is already a_opposite
        if (b == a_opposite) {
          cat("  ⓘ b is already opposite of a_orig, no action needed\n")
          next
        }
        
        # Find component currently labeled as b
        b_owner <- NULL
        for (orig_comp in orig_components) {
          if (final_map[orig_comp] == b) {
            b_owner <- orig_comp
            break
          }
        }
        
        if (is.null(b_owner)) {
          # Unilateral rename a_opposite_owner to b
          cat(sprintf("  ⓘ Target component not in current K, renaming: %s -> %s\n", 
                      a_opposite_owner, b))
          final_map[a_opposite_owner] <- b
          next
        }
        
        # Exchange a_opposite_owner with b_owner
        tmp <- final_map[a_opposite_owner]
        final_map[a_opposite_owner] <- final_map[b_owner]
        final_map[b_owner] <- tmp
        cat(sprintf("  ✓ Executed exchange: %s <-> %s\n", 
                    final_map[a_opposite_owner], final_map[b_owner]))
        
        # Record pairing
        pair_map[[final_map[a_opposite_owner]]] <- final_map[b_owner]
        pair_map[[final_map[b_owner]]] <- final_map[a_opposite_owner]
        next
      }
      
      # Regular exchange between a_orig and b_owner
      tmp <- final_map[a_orig]
      final_map[a_orig] <- final_map[b_owner]
      final_map[b_owner] <- tmp
      cat(sprintf("  ✓ Executed exchange: %s <-> %s\n", final_map[a_orig], final_map[b_owner]))
      
      # Record pairing
      pair_map[[final_map[a_orig]]] <- final_map[b_owner]
      pair_map[[final_map[b_owner]]] <- final_map[a_orig]
    }
  } else {
    cat("\n[Step 2] No component exchanges needed\n")
  }
  
  # ------------------------------------------------------------
  # 5.7 Step 3: Handle Unmapped Components
  # ------------------------------------------------------------
  mapped_components <- component_mapping$current_anc
  unmapped_components <- setdiff(orig_components, mapped_components)
  
  if (length(unmapped_components) > 0) {
    cat("\n[Step 3] Processing unmapped components:\n")
    for (comp in unmapped_components) {
      new_name <- paste0("Extra_", comp)
      final_map[comp] <- new_name
      cat(sprintf("  ⓘ Renaming unmapped component: %s -> %s\n", comp, new_name))
    }
  }
  
  # ------------------------------------------------------------
  # 5.8 Apply Final Mapping to Current K Value
  # ------------------------------------------------------------
  k_label <- paste0("K=", k_value)
  
  # Apply mapping using the Ancestry column
  current_indices <- which(admix_mapped$K == k_label)
  current_ancestries <- admix_mapped$Ancestry[current_indices]
  
  admix_mapped$Ancestry_plot[current_indices] <- final_map[current_ancestries]
  
  cat("\n[Final mapping] K=", k_value, ":\n")
  print(final_map)
}

# ------------------------------------------------------------
# 5.9 Ensure K_max (maximum K) Ancestry_plot Column is Correct
# ------------------------------------------------------------
max_k_indices <- which(admix_mapped$K == ref_K_label)
admix_mapped$Ancestry_plot[max_k_indices] <- admix_mapped$Ancestry[max_k_indices]

# ------------------------------------------------------------
# 5.10 Build Complete Ancestry Component List for Color Mapping
# ------------------------------------------------------------
all_ancestries <- unique(admix_mapped$Ancestry_plot)
ref_ancestries <- admix_mapped %>%
  filter(K == ref_K_label) %>%
  distinct(Ancestry_plot) %>%
  arrange(Ancestry_plot) %>%
  pull(Ancestry_plot)

# Ensure all components are included in color mapping
ancestries_for_color <- unique(c(ref_ancestries, all_ancestries))
n_ancestries <- length(ancestries_for_color)

# Generate colors
if (n_ancestries > 12) {
  color_vector <- colorRampPalette(brewer.pal(12, color_palette))(n_ancestries)
} else {
  color_vector <- brewer.pal(max(3, n_ancestries), color_palette)[1:n_ancestries]
}

color_dict <- setNames(color_vector, ancestries_for_color)

cat("\n=== Cross-K ancestry component mapping completed ===\n")
cat(sprintf("Total mapped ancestry components: %d\n", length(all_ancestries)))
cat(sprintf("Reference K (K=%d) has %d components\n", K_max, length(ref_ancestries)))

# =========================
# 6. Create Main Plots (Left Panel)
# =========================
create_main_plots <- function(admix_mapped, color_dict, K_range) {
  plot_list <- list()
  
  for (k_val in K_range) {
    k_label <- sprintf("K=%d", k_val)
    
    plot_data <- admix_mapped %>%
      filter(K == k_label) %>%
      arrange(sample_order) %>%
      mutate(
        Ancestry_plot = factor(Ancestry_plot, levels = names(color_dict))
      )
    
    # Calculate x-axis range
    total_samples <- max(plot_data$sample_order)
    
    p <- ggplot(plot_data, aes(x = sample_order, y = Proportion, fill = Ancestry_plot)) +
      geom_bar(stat = "identity", width = 1, linewidth = 0) +
      scale_fill_manual(values = color_dict) +
      scale_x_continuous(
        expand = expansion(mult = c(0, 0.005)),
        limits = c(0.5, total_samples + 0.5)
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0))) +
      labs(title = NULL, x = NULL, y = NULL) +
      # Add K-value label on left side
      annotate("text", 
               x = -total_samples * 0.045,
               y = 0.5,
               label = k_label,
               hjust = 1,
               vjust = 0.5,
               size = 4.0,
               fontface = "bold") +
      coord_cartesian(xlim = c(0.5, total_samples + 0.5), clip = "off") +
      theme_void() +
      theme(
        plot.margin = margin(1, 2, 1, 5, "pt"),
        legend.position = "none"
      )
    
    plot_list[[k_label]] <- p
  }
  
  plot_list
}

main_plots <- create_main_plots(admix_mapped, color_dict, K_range)

# =========================
# 7. Create Population Axis (Bottom of Left Panel)
# =========================
create_population_axis <- function(sample_order_df) {
  # Calculate population boundaries
  pop_positions <- sample_order_df %>%
    group_by(pop) %>%
    summarise(
      start = min(sample_order), 
      end = max(sample_order), 
      .groups = "drop"
    ) %>%
    mutate(mid = (start + end) / 2)
  
  # Calculate x-axis range matching bar plot
  x_min <- 0.5
  x_max <- max(sample_order_df$sample_order) + 0.5
  
  axis_plot <- ggplot(pop_positions) +
    # Baseline covering the exact bar plot range
    annotate("segment", 
             x = x_min, 
             xend = x_max, 
             y = 0, 
             yend = 0, 
             linewidth = 0.3) +
    # Small tick marks at population boundaries
    geom_segment(aes(x = start + 0.5, xend = start + 0.5, 
                     y = 0, yend = -0.03), 
                 linewidth = 0.2) +
    # Last boundary tick mark
    annotate("segment",
             x = x_max,
             xend = x_max,
             y = 0,
             yend = -0.03,
             linewidth = 0.2) +
    # Population labels (rotated 90 degrees)
    geom_text(aes(x = mid, y = -0.08, label = pop), 
              angle = 90, hjust = 1, size = 3,
              check_overlap = FALSE) +
    scale_x_continuous(
      limits = c(x_min - 0.1, x_max + 0.1),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(
      xlim = c(x_min, x_max),
      ylim = c(-0.15, 0.05), 
      clip = "off"
    ) +
    theme_void() +
    theme(
      plot.margin = margin(5, 5, 5, 5, "pt")
    )
  
  return(axis_plot)
}

population_axis <- create_population_axis(sample_order_df)

# =========================
# 8. Create Target Population Plots (Right Panel)
# =========================
create_target_population_plots <- function(admix_mapped, target_pops, color_dict, K_range) {
  
  # Calculate average ancestry proportions for target populations
  target_avg <- admix_mapped %>%
    filter(pop %in% target_pops) %>%
    group_by(K, pop, Ancestry_plot) %>%
    summarise(mean_prop = mean(Proportion), .groups = "drop") %>%
    mutate(
      pop = factor(pop, levels = target_pops),
      Ancestry_plot = factor(Ancestry_plot, levels = names(color_dict))
    )
  
  # Create one ggplot per K value
  right_plots <- list()
  
  for (k_val in K_range) {
    k_label <- paste0("K=", k_val)
    
    plot_data <- target_avg %>% filter(K == k_label)
    
    p_k <- ggplot(
      plot_data,
      aes(x = pop, y = mean_prop, fill = Ancestry_plot)
    ) +
      geom_bar(
        stat = "identity",
        width = 0.8,
        linewidth = 0,
        show.legend = FALSE
      ) +
      scale_fill_manual(values = color_dict, drop = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0, 0))) +
      # Add K-value label on right side
      annotate(
        "text",
        x = length(target_pops) + 0.5,
        y = 1.05,
        label = k_label,
        hjust = 0,
        vjust = 4,
        size = 4,
        fontface = "bold"
      ) +
      coord_cartesian(ylim = c(0, 1.1), clip = "off") +
      theme_void() +
      theme(
        plot.margin = margin(1, 5, 1, 5, "pt"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank()
      )
    
    right_plots[[k_label]] <- p_k
  }
  
  return(right_plots)
}

right_plots <- create_target_population_plots(
  admix_mapped,
  target_pops,
  color_dict,
  K_range
)

# =========================
# 9. Create Target Population Axis
# =========================
create_target_population_axis <- function(target_pops) {
  
  axis_df <- data.frame(
    pop = factor(target_pops, levels = target_pops),
    x = seq_along(target_pops)
  )
  
  ggplot(axis_df, aes(x = x, y = 0, label = pop)) +
    geom_text(
      angle = 90, 
      hjust = 0.5,
      vjust = -0.1,
      size = 4, 
      fontface = "bold"
    ) +
    scale_x_continuous(
      limits = c(0.5, length(target_pops) + 0.5),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(ylim = c(-0.2, 0.2), clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(1, 1, 1, 1, "pt")
    )
}

target_axis <- create_target_population_axis(target_pops)

# =========================
# 10. Assemble Final Plot
# =========================

# Combine right panel plots with appropriate heights
n_K <- length(K_range)
left_heights <- c(rep(10, n_K), 12)  # Extra height for population axis

right_panel <- wrap_plots(
  c(right_plots, list(target_axis)),
  ncol = 1,
  heights = left_heights
)

# Add margin to right panel for label spacing
right_panel <- right_panel +
  theme(plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt"))

# Combine left and right panels
left_panel <- wrap_plots(
  c(main_plots, list(population_axis)),
  ncol = 1,
  heights = left_heights
)

final_plot <- left_panel | right_panel

# Set width ratio (left:right = 5:1.2)
final_plot <- final_plot + 
  plot_layout(widths = c(5, 1.2))

# Add plot annotations
final_plot <- final_plot +
  plot_annotation(
    title = "Admixture Analysis of Target and Reference Populations",
    subtitle = sprintf("K values: %s", paste(K_range, collapse = ", ")),
    caption = "Note: Colors are consistent across K values for the same ancestry component",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5,
                                margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0.5,
                                   margin = margin(b = 5)),
      plot.caption = element_text(size = 9, hjust = 0.5, face = "italic",
                                  margin = margin(t = 15, b = 5)),
      plot.margin = margin(10, 10, 25, 10, "pt")
    )
  )

# Display the plot
print(final_plot)

# =========================
# 11. Save the Plot
# =========================
plot_height <- plot_height_base + (length(K_range) - 3) * height_per_K

ggsave(
  filename = output_png,
  plot = final_plot,
  width = plot_width,
  height = plot_height,
  dpi = 300,
  bg = "white"
)

cat(sprintf("Plot saved to: %s\n", output_png))
cat(sprintf("Plot dimensions: %.1f x %.1f inches\n", plot_width, plot_height))

# =========================
# 12. Optional: Save Ancestry Mapping Tables
# =========================

# Save ancestry mapping summary
ancestry_summary <- admix_mapped %>%
  distinct(K, Ancestry, Ancestry_plot) %>%
  arrange(K, Ancestry)

write.csv(ancestry_summary, "ancestry_mapping_summary.csv", row.names = FALSE)
cat("Ancestry mapping summary saved to ancestry_mapping_summary.csv\n")

# Save color mapping
color_df <- data.frame(
  Ancestry = names(color_dict),
  Color = color_dict,
  stringsAsFactors = FALSE
)
write.csv(color_df, "color_mapping.csv", row.names = FALSE)
cat("Color mapping saved to color_mapping.csv\n")

cat("\n=== Plotting completed successfully ===\n")
cat(sprintf("K-value range: %d-%d\n", min(K_range), max(K_range)))
cat(sprintf("Total samples: %d\n", length(unique(admix_mapped$sample))))
cat(sprintf("Total populations: %d\n", length(unique(admix_mapped$pop))))
cat(sprintf("Target populations: %s\n", paste(target_pops, collapse = ", ")))