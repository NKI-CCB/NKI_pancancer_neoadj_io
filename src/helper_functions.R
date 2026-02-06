# title: "General utility functions"
# author: "Pedro Batista Tan"
# date: "09/10/2023"

# Function to require packages, and install missing ones
using_packages <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs, require,character.only = TRUE))
  not_installed <- libs[req == FALSE]
  if(length(not_installed) > 0){ 
    install.packages(not_installed)
    lapply(not_installed, require, character.only = TRUE)
  }
}

using_bioconductor_packages <- function(...) {
  if (!requireNamespace("BiocManager")) {
    install.packages("BiocManager")
  }
  libs <- unlist(list(...))
  req <- unlist(lapply(libs, require,character.only = TRUE))
  not_installed <- libs[req == FALSE]
  if(length(not_installed) > 0){ 
    BiocManager::install(not_installed)
    lapply(not_installed, require, character.only = TRUE)
  }
}

save_forestplot_pdf <- function(forest_plot_data,
                                filename, width = 12, height = 7,
                                title_text = NULL, 
                                title_x = 0.5,
                                title_y = 0.775,
                                ...) {
  stopifnot(!missing(forest_plot_data))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  forest(forest_plot_data, studlab = TRUE, common = FALSE, random = TRUE, ...)
  
  # Add title if provided
  if (!is.null(title_text)) {
    grid::grid.text(title_text,
                    x = title_x, y = title_y,
                    gp = grid::gpar(fontsize = 14, fontface = "bold"))
  }
  # fp1 <-  grid.grabExpr(print(forest(forest_plot_data, studlab = TRUE)))
  # gridExtra::grid.arrange(fp1, fp1)
  dev.off()
}

# Function to require packages, and install missing ones
using_packages("ggplot2", "gridExtra", "reshape2", "pheatmap")

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}


factor_summary <- function(object, maxsum = 20) {
  object %>% mutate_if(is.character, as.factor) %>% summary(maxsum = maxsum)
}

floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)

ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)


FormatDecimals <- function(values) {
  result <- numeric(length(values))  # Initialize an empty numeric vector
  
  for (i in seq_along(values)) {
    value <- values[i]
    
    if (is.na(value)) {
      result[i] <- value
    } else if (abs(value) >= 10) {
      result[i] <- round(value, digits = 1)
    } else if (abs(value) >= 0.1) {
      result[i] <- round(value, digits = 2)
    } else if (abs(value) >= 0.01) {
      result[i] <- round(value, digits = 3)
    } else if (abs(value) >= 0.001) {
      result[i] <- round(value, digits = 4)
    } else {
      result[i] <- format(value, scientific = TRUE, digits = 2)
    }
  }
  
  return(result)
}

select_top_extremes <- function(all_scores, top_n = 100) {
  # function that returns TRUE/FALSE for rankings within the top_n values
  increasing <- sort(all_scores, decreasing = FALSE)[1:top_n]
  decreasing <- sort(all_scores, decreasing = TRUE)[1:top_n]
  return(all_scores %in% increasing | all_scores %in% decreasing)
}

select_pval <- function(all_pvals, threshold = 0.01) {
  # function that returns TRUE/FALSE for pvalues below threshold
  return(all_pvals < threshold)
}

# Function to extract correlation coefficient and p-values
corr_func <- function(var1, var2, data) {
  res = cor.test(data[,var1], data[,var2])
  tidy(res) %>% mutate(var1 = var1, var2 = var2)
}

save_pdf <- function(plot, filename, width = 12, height = 7) {
  stopifnot(!missing(plot))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid.arrange(plot)
  # fp1 <-  grid.grabExpr(print(forest(forest_plot_data, studlab = TRUE)))
  # gridExtra::grid.arrange(fp1, fp1)
  dev.off()
}

# TODO Test Complex heatmap? Different behaviour than pheatmap?
# save_pheatmap_pdf <- function(x, filename, width = 7, height = 7) {
#   stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   pdf(filename, width = width, height = height)
#   grid::grid.newpage()
#   draw(x)
#   dev.off()
# }

save_pheatmap_pdf <- function(x, filename, width = 7, height = 7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_complex_heatmap <- function(x, filename, width = 7, height = 7, res = 300) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width = width, height = height, units = "in", res = res)
  draw(x)
  dev.off()
}

counts_CI_for_frequency_tables <- function(data, column_name, value) {
  # Check how many times the value appears on the corresponding column
  n_value <- sum(data[[column_name]] == value, na.rm = TRUE)
  # Calculate the confidence intervals for proportions
  question_ci <- freqCI(data[[column_name]], level = 0.95)
  # Get the index of the category that corresponds to the value queried
  index <- which(question_ci$cat_names == value)
  
  question_ci_count <- c(Number_of_values = n_value,
                         Total = nrow(data),
                         percentage_of_total = round(
                           question_ci$rel_freq[index] * 100, 1),
                         CI_lower = round(question_ci$CIs_low[index] * 100, 1),
                         CI_upper = round(question_ci$CIs_high[index] * 100, 1)
  )
  if (n_value == 0) {
    question_ci_count <- c(Number_of_values = 0,
                           Total = nrow(data),
                           percentage_of_total = 0,
                           CI_lower = 0,
                           CI_upper = 0
    )
  }
  
  return(list(question_ci_count))
}


run_enrichr <- function(deg_df, dbs, output_filename) {
  
  enrichr_res <- enrichr(deg_df$gene, dbs)
  saveRDS(enrichr_res, here(output_dir, output_filename))
  
  for (db in dbs){
    print(db)
    enrichr_res_db <- as.data.frame(
      enrichr_res[[db]][enrichr_res[[db]]$Adjusted.P.value < 0.1, ])
    
    if (nrow(enrichr_res_db) != 0) {
      plot <- try(plotEnrich(enrichr_res_db,
                             numChar = 40, showTerms = min(25, nrow(enrichr_res_db)),
                             y = "Count", orderBy = "Adjusted.P.value", title = db)) +
        theme_Publication() +
        theme(legend.position = "right",
              legend.direction = "vertical") +
        guides(fill = guide_colourbar(direction = "vertical", position = "right",
                                      barwidth = 0.5, barheight = 5))
      
      print(plot)
    } else {
      print("No results below FDRq threshold")
    }
  }
}

filter_enrichr <- function(enrichr_data, dbs, fdrq_threshold = 0.1) {
  db_data <- data.frame(matrix(ncol = 10, nrow = 0)) %>%
    setNames(c("Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value",
               "Old.Adjusted.P.value", "Odds.Ratio", "Combined.Score",
               "Genes", "kb_db"))
  
  for (current_db in dbs) {
    db_data_f <- enrichr_data[[current_db]] %>%
      filter(Adjusted.P.value < fdrq_threshold) %>%
      mutate(kb_db = current_db)
    
    db_data <- db_data %>% rbind(db_data_f)
  }
  # db_data <- db_data %>% select(Term, Overlap, Adjusted.P.value,
  # Odds.Ratio, Genes, kb_db)
  
  db_data <- db_data %>%
    separate(Overlap, c("Term_gene_count", "Term_total"),
             sep = "/", remove = FALSE, convert = TRUE) %>%
    mutate(Overlap_n = Term_gene_count / Term_total)
  return(db_data)
}

plot_enrichr_gene_count_fraction <- function(enrichr_filtered_res) {
  fig <- enrichr_filtered_res  %>%
    ggplot(aes(x = reorder(Term, Adjusted.P.value), y = Term_gene_count,
               fill = Overlap_n)) +
    geom_col()  +
    facet_wrap(~kb_db, scales = "free") +
    theme_Publication() +
    theme(axis.text.x = element_text(
      size = 9, angle = 90, vjust = 0.5, hjust = 1))
  print(fig)
  
  fig <- enrichr_filtered_res  %>%
    ggplot(aes(x = reorder(Term, Adjusted.P.value), y = Overlap_n,
               fill = log10(Term_gene_count))) +
    geom_col()  +
    facet_wrap(~kb_db, scales = "free") +
    theme_Publication() +
    theme(axis.text.x = element_text(
      size = 9, angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = "", y = "Representend Fraction\nof Ontology Term", fill = "log10(Count)")
  print(fig)
  
}

plot_enrichr_term_overlap <- function(enrichr_filtered_res){
  db_terms <- enrichr_filtered_res %>%
    select(Term, Genes) %>%
    deframe() %>%
    strsplit(";") %>%
    enframe() 
  db_terms <- unnest(db_terms) %>%
    setNames(c("Term", "Gene"))
  
  db_terms_l <- split(
    db_terms$Gene, db_terms$Term)
  db_terms_l <- db_terms_l[
    db_terms$Term %>% unique()]
  
  term_overlap <- computeGeneSetsOverlap(db_terms_l,
                                         db_terms$Gene,
                                         min.sz = 1, max.sz = Inf)
  
  term_overlap_l <- melt(term_overlap)
  term_overlap_l <- term_overlap_l[term_overlap_l$value != 0, ]
  
  term_sizes <- enrichr_filtered_res %>%
    select(Term, Term_gene_count, Term_total)
  
  term_overlap_l <- term_overlap_l %>%
    left_join(term_sizes %>%
                select(Term, Term_gene_count) %>%
                unique() %>%
                mutate(Var1 = Term,
                       geneset_1 = sprintf("%s (%s)", Term, Term_gene_count)),
              by = "Var1") %>%
    left_join(term_sizes %>%
                select(Term, Term_gene_count) %>%
                unique() %>%
                mutate(Var2 = Term,
                       geneset_2 = sprintf(
                         "%s (%s)", Term, Term_gene_count)) %>%
                select(Var2, geneset_2),
              by = "Var2")
  
  plot <- ggplot(term_overlap_l, aes(x = geneset_2, y = geneset_1)) +
    geom_raster(aes(fill = value)) +
    scale_fill_gradient(low = "grey90", high = "red") +
    labs(x = "", y = "", title = "Gene set overlap") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 11))
  print(plot)
  
  term_overlap_l
}

check_enrichr_term_genes <- function(enrichr_filtered_res, terms) {
  db_terms <- enrichr_filtered_res %>%
    select(Term, Genes) %>%
    deframe() %>%
    strsplit(";") %>%
    enframe() 
  db_terms <- unnest(db_terms) %>%
    setNames(c("Term", "Gene")) %>%
    left_join(enrichr_filtered_res %>%
                select(Term, kb_db)) %>%
    left_join(gene_descriptions %>% rename(Gene = gene)) %>%
    filter(Term %in% terms)
  
  db_terms
  
}


# title: "ggplot theme for publication ready Plots"
# author: "Koundinya Desiraju"
# date: "04/07/2015"
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}


calc_stat_groups <- function(data, treatment, column, group1, group2,
                             contrast, current_gene_signature) {
  
  timepoint_group1 <- data[data$Treatment == treatment &
                             data[, column] == group1,
                           current_gene_signature, drop = TRUE]
  
  mean_group1 <- mean(timepoint_group1)
  sd_group1 <- sd(timepoint_group1)
  median_group1 <- median(timepoint_group1)
  n_group1 <- length(timepoint_group1)
  
  timepoint_group2 <- data[data$Treatment == treatment &
                             data[, column] == group2,
                           current_gene_signature, drop = TRUE]
  
  mean_group2 <- mean(timepoint_group2)
  sd_group2 <- sd(timepoint_group2)
  median_group2 <- median(timepoint_group2)
  n_group2 <- length(timepoint_group2)
  
  wilcoxon_test <- wilcox.test(timepoint_group1, timepoint_group2,
                               alternative = "two.sided")
  n_total <- n_group1 + n_group2
  
  stat_df <- wilcoxon_test %>%
    broom::tidy() %>%
    mutate(
      "contrast" = contrast,
      "group1" = group1,
      "group2" = group2,
      "signature" = current_gene_signature,
      "n" = n_total,
      "mean_group1" = mean_group1,
      "sd_group1" = sd_group1,
      "median_group1" = median_group1,
      "n_group1" = n_group1,
      "mean_group2" = mean_group2,
      "sd_group2" = sd_group2,
      "median_group2" = median_group2,
      "n_group2" = n_group2)
  
  return(stat_df)
}

print_stat_groups <- function(data){
  print(data[1, ] %>%
          mutate(p.value = FormatDecimals(p.value),
                 "n [group1/2]" = sprintf("%s [%s/%s]",
                                          n, n_group1, n_group2),
                 "Mean (SD) group1" = sprintf(
                   "%s (%s)", FormatDecimals(mean_group1),
                   FormatDecimals(sd_group1)),
                 "Median group1" = FormatDecimals(median_group1),
                 "Mean (SD) group2" = sprintf(
                   "%s (%s)", FormatDecimals(mean_group2),
                   FormatDecimals(sd_group2)),
                 "Median group2" = FormatDecimals(
                   median_group2)) %>%
          select(c("statistic", "p.value", "contrast", "n [group1/2]",
                   "Mean (SD) group1", "Mean (SD) group2",
                   "Median group1", "Median group2")) %>%
          rbind(
            data[2, ] %>%
              mutate(p.value = FormatDecimals(p.value),
                     "n [group1/2]" = sprintf("%s [%s/%s]",
                                              n, n_group1, n_group2),
                     "Mean (SD) group1" = sprintf(
                       "%s (%s)", FormatDecimals(mean_group1),
                       FormatDecimals(sd_group1)),
                     "Median group1" = FormatDecimals(median_group1),
                     "Mean (SD) group2" = sprintf(
                       "%s (%s)", FormatDecimals(mean_group2),
                       FormatDecimals(sd_group2)),
                     "Median group2" = FormatDecimals(
                       median_group2)) %>%
              select(c("statistic", "p.value", "contrast", "n [group1/2]",
                       "Mean (SD) group1", "Mean (SD) group2",
                       "Median group1", "Median group2"))
          ) %>%
          knitr::kable())
}

# Scatter <- ggplot(mtcars, aes(mpg,disp,color=factor(carb))) + geom_point(size=3) + labs(title="Scatter Plot")
# grid.arrange(Scatter,(Scatter +scale_colour_Publication()+ theme_Publication()),nrow=1)
# 
# Bar <- ggplot(mtcars, aes(factor(carb),fill=factor(carb))) + geom_bar(alpha=0.7) + labs(title="Bar Plot")
# 
# grid.arrange(Bar,(Bar + scale_fill_Publication() +theme_Publication()),nrow=1)
# 
# Bubble <- ggplot(mtcars, aes(mpg,disp,color=factor(carb),size=hp)) + geom_point(alpha=0.7) + labs(title="Bubble Plot") + scale_size_continuous(range = c(3,10))
# 
# grid.arrange(Bubble,(Bubble +scale_colour_Publication()+ theme_Publication()),nrow=1)
# 
# mtcars$Index <- 1:nrow(mtcars)
# dat <- melt(mtcars,id.vars = c("Index"),measure.vars = c("drat","wt"))
# Line <- ggplot(dat,aes(Index,value,colour=variable))+geom_line(size=1.3) + labs(title="Line Plot") 
# grid.arrange(Line,(Line +scale_colour_Publication()+ theme_Publication()),nrow=1)
# 
# P <- ggplot(data = mpg,aes(cty, hwy,color=class))+geom_point(size=3) + facet_wrap(~ manufacturer,scales="free")+
#   labs(title="Plot With Facets")
# P
# 
# P + scale_colour_Publication() + theme_Publication()

