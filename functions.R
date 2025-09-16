#' Calculate Shannon entropy
#'
#' @param frequencies A vector of nucleotide frequencies
#'
#' @return A numerical Shannon entropy value
#' @export
calculate_shannon <- function(frequencies) {
  # If frequencies don't sum to 1, calculate reference frequency
  if (sum(frequencies) < 1) {
    ref <- 1 - sum(frequencies)
    
    # Add reference frequency to list of frequencies
    frequencies %<>% c(ref)
  }
  
  # Calculate shannon entropy
  se <- -sum(frequencies * log2(frequencies))
  
  return(se)
}

#' Assign position to feature in SARS-CoV-2 genome
#'
#' @param POS a numerical genome position
#' @param anno_db a data frame containing at least the following columns: start, end, gene, feature
#'
#' @return
#' @export
assign_feature <- function(pos_val, anno_db) {
  x <- anno_db[start <= pos_val & end >= pos_val, feature]
  if (length(x) == 0) x <- NA
  return(x)
}


#' Generate a linear mixed effect model or lm if lme4 triggers errors
#'
#' @param df a data frame containing all variables specified in other arguments
#' @param y character string specifying outcome variable
#' @param predictors a character vector specifying predictor variables
#' @param random character string specifying random effect variables
#' @param add_RDV logical, whether to add RDV as predictor
#'
#' @return an lme4 or lm class object
#' @export
generate_model <- function(df, 
                           y,
                           predictors,
                           interaction=NULL,
                           random = NULL,
                           add_RDV = T) {
  
  # If data frame contains any variability in RDV status, add as predictor
  if (length(unique(df$RDV)) > 1 & add_RDV) {predictors %<>% c("RDV")}
  
  # If a predictor var is a character type...
  cols <- predictors[lapply(df[, ..predictors], is.character) %>% unlist()]
  
    # ...convert into factor
    df %<>% .[, (cols) := lapply(.SD, as.factor),
       .SDcols = cols]
  
  # Generate formula
  formula <- create_formula(y = y,
                            predictors = c(predictors,interaction),
                            random = random)
  
  # Define safe functions for error handling
  possibly_lme4 <- possibly(\(x) lmerTest::lmer(formula, data = x))
  
  possibly_lm <- possibly(\(x) lm(create_formula(y = y,
                                                 predictors = predictors), 
                                  data = x))
  # Generate model
  model <- df |> possibly_lme4()
  
  # If no output from lme4, try lm
  if (is.null(model)) {
    model <- df |> possibly_lm()
  }
  
  return(model)
}

#' Helper function to generate QQ plot for model
#'
#' @param model A lm or lmer class object
#' @param title Character string for title to print on plot
#'
#' @return a base R QQ plot
#' @export
assess_QQ_plot <- function(model,
                           title) {
  par(mfrow = c(1, 1))
  qqnorm(resid(model), main = title)
  qqline(resid(model), col = "dodgerblue", lwd = 2)
}

#' Save a ggplot
#'
#' @param plot ggplot object
#' @param filename path for saving plot 
#'
#' @return saves plot at filename
#' @export
save_plot <- function(plot, filename) {
  # Determine dimensions for saving
  out <- ggplot_build(plot)
  
  # Scale plot dimensions according to rows and columns of plot
  ggsave(filename, 
         plot,
         width = 800 * max(out$layout$layout$COL),
         height = 600 * max(out$layout$layout$ROW),
         scale = c(3:1)[max(out$layout$layout$ROW)],
         units = "px")
}

#' Generate scatter plot with SARS-CoV-2 genome annotation
#'
#' @param df a data frame containing all variables listed in x, y, color, facet_var
#' @param x character string specifying X axis variable
#' @param y character string specifying Y axis variable
#' @param color character string specifying variable to use to color points
#' @param facet_var character string specifying variable to use to facet plots
#' @param xlab character string specifying X axis title
#' @param ylab character string specifying Y axis title
#' @param color_lab character string specifying label for color legend
#' @param title character string specifying plot title
#' @param anno a data frame containing at least the following columns: start, end, gene
#'
#' @return a ggplot object
#' @export
scatter_plot_with_genome <- function(df,
                                     x,
                                     y,
                                     color,
                                     facet_var,
                                     xlab,
                                     ylab,
                                     color_lab = color,
                                     title = element_blank(),
                                     anno = NULL) {
  ggplot(df, 
         mapping = aes(x = .data[[x]], 
                       y = .data[[y]],
                       color = .data[[color]])) +
    geom_point() +
    theme_bw() +
    facet_wrap(vars(.data[[facet_var]])) +
    geom_gene_arrow(data = anno,
                    aes(xmin = start, 
                        xmax = end, 
                        y = -0.05,
                        fill = gene),
                    arrowhead_height = unit(1, "char"), 
                    arrow_body_height = unit(1, "char"),
                    arrowhead_width = unit(0, "char"),
                    inherit.aes = F) +
    scale_fill_paletteer_d(palette) +
    geom_gene_label(data = anno,
                    aes(xmin = start, 
                        xmax = end, 
                        y = -0.05 , 
                        label = gene),
                    inherit.aes = F,
                    align = "centre",
                    grow = T) +
    labs(x = xlab,
         y = ylab,
         title = title,
         color = color_lab) -> plot
  
  return(plot)
}

#' Generate and save multiple plots generated with scatter_plot_with_genome()
#'
#' @param df a data frame containing all variables listed in x, y, color, facet_var
#' @param plot_var a character string specifying variable to split df by for separate plots
#' @param save_path path to save generated figures
#' @param x character string specifying X axis variable
#' @param y character string specifying Y axis variable
#' @param color character string specifying variable to use to color points
#' @param facet_var character string specifying variable to use to facet plots
#' @param xlab character string specifying X axis title
#' @param ylab character string specifying Y axis title
#' @param color_lab character string specifying label for color legend
#' @param title character string specifying plot title
#' @param anno a data frame containing at least the following columns: start, end, gene
#'
#' @return saved plots at save_path
#' @export
make_and_save_scatter_plots <- function(df, 
                                        plot_var,
                                        save_path = "plots/",
                                        x,
                                        y,
                                        color,
                                        facet_var,
                                        xlab,
                                        ylab,
                                        color_lab,
                                        anno = anno_db_gene) {
  
  # Split data frames by variable
  df %<>% split(f = .[[plot_var]])
  
  # Get title names
  title <- names(df)
  
  # Generate plots
  list(df, title) |> 
    purrr::pmap(\(d, t) scatter_plot_with_genome(d, 
                                                 x = x,
                                                 y = y,
                                                 color = color,
                                                 facet_var = facet_var,
                                                 xlab = xlab,
                                                 ylab = ylab, 
                                                 color_lab = color_lab,
                                                 title = t,
                                                 anno = anno)) -> plots
  
  # Generate save file names
  filenames <- paste(save_path, title, ".png", sep = "")
  
  # Save plots
  list(plots, filenames) |> 
    purrr::pmap(\(plot, filename) save_plot(plot, filename))
  
}

#' Calculate days since last Remdesivir treatment
#'
#' @param df a data frame containing at least the following columns: MRN, Days
#' @param treatment a data frame containing at least the following columns: MRN, start, end, Medication
#'
#' @return numerical value indicating days since end of last RDV treatment
#' @export
days_since_RDV <- function(df, treatment) {
  # Take only entries for relevant patient started before sample data
  treatment %<>% 
    .[MRN == unique(df[["MRN"]]) & 
        Medication == "Remdesivir" &
        start <= df[["Days"]],]
  
  days_since <- as.numeric(df[["Days"]]) - as.numeric(treatment[["end"]])
  
  # If any remdesivir treatments administered
  if (length(days_since) > 0) {
    # return days since most recent administration
    days_since %<>% min()
  } else {
    # If no RDV ever administered, return NA
    days_since <- NA
  }
  
  # If during administration (i.e. negative end date), return 0
  if (!is.na(days_since) & days_since < 0) {days_since <- 0}
  
  return(days_since)
}

#' Create formula for lmem or glm
#'
#' @param y character string specifying outcome variable
#' @param predictors a character vector specifying predictor variables
#' @param random character string specifying random effect variables
#'
#' @return a class formula object
#' @export
create_formula <- function(y, 
                           predictors, 
                           random = NULL) {
  # Collate all predictors 
  predictors %<>% paste(., collapse = "+")
  formula <- paste(y, "~", predictors)
  
  # If random variable provided
  if (!is.null(random)) {
    # Format random variable
    random <- paste("(1|", random, ")", sep = "")
    
    # Include in formula
    formula %<>% paste(random, sep = " + ")
  }
  
  # Convert string to formula class object
  formula %<>% as.formula()
  
  return(formula)
}