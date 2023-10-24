# Define a function to create Fst plots
createFstPlots <- function(input_file, output_dir, chrom_to_remove, start_column, end_column, plotting_params) {
  # Read the input file
  data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Create a new column 'win_mid' as the average of 'Start' and 'End'
  data$win_mid <- (data[, start_column] + data[, end_column]) / 2
  
  # Remove lines with a specific chromosome value
  data <- data[!data$Chromosome %in% chrom_to_remove, ]
  
  # Create individual Fst plots for each column containing Fst values
  for (col_name in colnames(data)) {
    if (col_name != "Chromosome" && col_name != start_column && col_name != end_column && col_name != "features" && col_name != "win_mid") {
      plot_data <- data[, c("win_mid", col_name)]
      plot <- ggplot(plot_data, aes(x = win_mid, y = .data[[col_name]])) +
        geom_line(size = plotting_params$line_size) +
        geom_hline(yintercept = plotting_params$hline_y) +
        geom_vline(xintercept = plotting_params$vline_x) +
        labs(title = plotting_params$plot_title, x = plotting_params$x_label, y = plotting_params$y_label) +
        scale_y_continuous(limits = c(plotting_params$ymin, plotting_params$ymax)) +
        theme_bw() +
        theme(
          text = element_text(
            family = plotting_params$font_family,
            size = plotting_params$font_size,
            face = plotting_params$font_face
          ),
          strip.text = element_text(size = plotting_params$strip_text_size),
          axis.title.x = element_text(size = plotting_params$axis_title_size),
          axis.title.y = element_text(size = plotting_params$axis_title_size),
          axis.text.x = element_text(size = plotting_params$axis_text_size),
          axis.text.y = element_text(size = plotting_params$axis_text_size),
          legend.text = element_text(size = plotting_params$legend_text_size),
          legend.title = element_text(size = plotting_params$legend_title_size)
        ) +
        facet_grid(facet_rows = plotting_params$facet_rows, facet_cols = plotting_params$facet_cols, scales = plotting_params$facet_scales, space = plotting_params$facet_space)
      
      # Save the plot in the specified output directory
      output_file <- file.path(output_dir, paste0("Fst_Plot_", col_name, ".png"))
      ggsave(output_file, plot, width = plotting_params$width, height = plotting_params$height)
   

# Define command line options
option_list <- list(
  make_option("--input", type = "character", help = "Input file with Fst values"),
  make_option("--output-dir", type = "character", help = "Output directory for saving plots"),
  make_option("--chrom-to-remove", type = "character", default = "RCWQ01006588", help = "Chromosome to remove"),
  make_option("--start-column", type = "character", default = "Start", help = "Name of the Start column"),
  make_option("--end-column", type = "character", default = "End", help = "Name of the End column"),
  make_option("--title-size", type = "integer", default = 12, help = "Title font size"),
  make_option("--label-size", type = "integer", default = 10, help = "Label font size"),
  make_option("--width", type = "numeric", default = 8, help = "Plot width"),
  make_option("--height", type = "numeric", default = 6, help = "Plot height"),
  make_option("--facet-rows", type = "character", default = NULL, help = "Facet grid rows"),
  make_option("--facet-cols", type = "character", default = NULL, help = "Facet grid columns"),
  make_option("--facet-scales", type = "character", default = NULL, help = "Facet grid scales"),
  make_option("--facet-space", type = "character", default = NULL, help = "Facet grid space"),
  make_option("--plot-title", type = "character", default = NULL, help = "Main plot title"),
  make_option("--x-label", type = "character", default = NULL, help = "X-axis label"),
  make_option("--y-label", type = "character", default = NULL, help = "Y-axis label"),
  make_option("--line-size", type = "numeric", default = 1, help = "Line plot width"),
  make_option("--hline-y", type = "numeric", default = NULL, help = "Horizontal line position"),
  make_option("--vline-x", type = "numeric", default = NULL, help = "Vertical line position"),
  make_option("--ymin", type = "numeric", default = NULL, help = "Minimum value for the y-axis"),
  make_option("--ymax", type = "numeric", default = NULL, help = "Maximum value for the y-axis"),
  make_option("--font-family", type = "character", default = NULL, help = "Font family"),
  make_option("--font-size", type = "integer", default = 12, help = "Font size"),
  make_option("--font-face", type = "character", default = "plain", help = "Font face"),
  make_option("--strip-text-size", type = "integer", default = 10, help = "Strip text size"),
  make_option("--axis-title-size", type = "integer", default = 12, help = "Axis title size"),
  make_option("--axis-text-size", type = "integer", default = 10, help = "Axis text size"),
  make_option("--legend-text-size", type = "integer", default = 10, help = "Legend text size"),
  make_option("--legend-title-size", type = "integer", default = 12, help = "Legend title size")
)

# Parse command line options
opt <- parse_args(OptionParser(option_list = option_list))

# Call the function to create Fst plots
createFstPlots(
  input_file = opt$input,
  output_dir = opt$output_dir,
  chrom_to_remove = opt$chrom_to_remove,
  start_column = opt$start_column,
  end_column = opt$end_column,
  plotting_params = list(
    line_size = opt$line_size,
    hline_y = opt$hline_y,
    vline_x = opt$vline_x,
    plot_title = opt$plot_title,
    x_label = opt$x_label,
    y_label = opt$y_label,
    ymin = opt$ymin,
    ymax = opt$ymax,
    font_family = opt$font_family,
    font_size = opt$font_size,
    font_face = opt$font_face,
    strip_text_size = opt$strip_text_size,
    axis_title_size = opt$axis_title_size,
    axis_text_size = opt$axis_text_size,
    legend_text_size = opt$legend_text_size,
    legend_title_size = opt$legend_title_size,
    facet_rows = opt$facet_rows,
    facet_cols = opt$facet_cols,
    facet_scales = opt$facet_scales,
    facet_space = opt$facet_space,
    width = opt$width,
    height = opt$height
  ))

    }
  }
}



