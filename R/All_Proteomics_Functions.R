#' Expression Level Phenotype Analysis
#'
#' This function will take expression data directly from ProteomeDiscoverer (v. 2.3) and generate a volcano plot and hit list
#' based on user-defined SD and p-values
#'
#' @param Expression_data Raw ProteomeDiscoverer (v. 2.3) data
#' @param SD_cutoff A positive integer (default = 2)
#' @param p_cutoff A positive integer between 0 and 1 (default = 0.05)
#' @param TMTplex 8, 10, 16 TMTplex (default = 10)
#' @param plot Display volcano plot (default = TRUE)
#' @param xlab X-axis Label (default = "Z Score")
#' @param ylab Y-axis Label (default = "- log10 (p value)")
#' @param labelcolor Color of Not significant and Significant data points (default = c("grey", "red"))
#' @param alpha Transparency of plot (default = 0.5)
#' @param plottextsize Text Size of plot (default = 12)
#' @param correctedpvalue Bonferroni corrected p-value (default = FALSE)
#' @param interactiveplot Plotly interactive Volcano plot (default = FALSE)
#' @return Volcano plot, Hit list, Total proteins assayed, and Unique protein hits
#' @export

ExpressionLevel_Phenotype.Fun <- function(Expression_data, SD_cutoff = 2, p_cutoff = 0.05, TMTplex = 10, plot = TRUE, xlab = "Z Score", ylab = "- log10 (p value)", labelcolor = c("grey", "red"), alpha = 0.5, plottextsize = 12, correctedpvalue = FALSE, interactiveplot = FALSE){

  if (TMTplex !=10)
    stop("This TMTplex is not yet supported in this version.\n")

  if (SD_cutoff < 0)
    stop("SD_cutoff must be non-negative.\n")

  if (p_cutoff < 0)
    stop("p_value cutoff must be non-negative.\n")


  High_Confidence_Exp_Level <-
    dplyr::filter(Expression_data, Expression_data[,grepl("Confidence", names(Expression_data))] == "High")

  # This dataframe takes all useful columns from the raw data and will be outputed at the end for reference

  Exp_Level_Concise <- magrittr::"%>%" (High_Confidence_Exp_Level,
                                        dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Description", "Accession", (grep(("^Abundances\\s\\(Grouped):"), names(High_Confidence_Exp_Level), value = TRUE))))))

  # Remove all blank rows in Dataframe

  NA_Removed <- tidyr::drop_na(Exp_Level_Concise)

  # Step 1

  Exp_Level_1 <- magrittr::"%>%" (NA_Removed,
                                  dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Description", "Accession"))))
  Exp_Level_2 <- magrittr::"%>%" (NA_Removed,
                                  dplyr::select(tidyselect::any_of(grep(("^Abundances\\s\\(Grouped):"), names(High_Confidence_Exp_Level), value = TRUE))))
  step1 <- cbind(Exp_Level_1, Exp_Level_2)

  if ("Master Protein Accessions" %in% names(step1)){
    step1 <- dplyr::rename(step1, "Accession" = "Master Protein Accessions")
  }

  # Step 2

  if (("Accession" %in% names(step1) && "Description" %in% names(step1))) {

    step1_real <-
      dplyr::filter(
        step1,
        step1[,3] > 0,
        step1[,4] > 0,
        step1[,5] > 0,
        step1[,6] > 0,
        step1[,7] > 0,
        step1[,8] > 0,
        step1[,9] > 0,
        step1[,10] > 0,
        step1[,11] > 0,
        step1[,12] > 0)

  Normalized_Exp_1 <- (as.numeric(unlist(step1_real[,3])) * mean(as.numeric(unlist(step1_real[,8])))) / ((as.numeric(unlist(step1_real[,8]))) * mean(as.numeric(unlist(step1_real[,3]))))
  Normalized_Exp_2 <- (as.numeric(unlist(step1_real[,4])) * mean(as.numeric(unlist(step1_real[,9])))) / ((as.numeric(unlist(step1_real[,9]))) * mean(as.numeric(unlist(step1_real[,4]))))
  Normalized_Exp_3 <- (as.numeric(unlist(step1_real[,5])) * mean(as.numeric(unlist(step1_real[,10])))) / ((as.numeric(unlist(step1_real[,10]))) * mean(as.numeric(unlist(step1_real[,5]))))
  Normalized_Exp_4 <- (as.numeric(unlist(step1_real[,6])) * mean(as.numeric(unlist(step1_real[,11])))) / ((as.numeric(unlist(step1_real[,11]))) * mean(as.numeric(unlist(step1_real[,6]))))
  Normalized_Exp_5 <- (as.numeric(unlist(step1_real[,7])) * mean(as.numeric(unlist(step1_real[,12])))) / ((as.numeric(unlist(step1_real[,12]))) * mean(as.numeric(unlist(step1_real[,7]))))

  }

  if (!"Description" %in% names(step1)){

    step1_real <-
      dplyr::filter(
        step1,
        step1[,2] > 0,
        step1[,3] > 0,
        step1[,4] > 0,
        step1[,5] > 0,
        step1[,6] > 0,
        step1[,7] > 0,
        step1[,8] > 0,
        step1[,9] > 0,
        step1[,10] > 0,
        step1[,11] > 0)

    Normalized_Exp_1 <- (as.numeric(unlist(step1_real[,2])) * mean(as.numeric(unlist(step1_real[,7])))) / ((as.numeric(unlist(step1_real[,7]))) * mean(as.numeric(unlist(step1_real[,2]))))
    Normalized_Exp_2 <- (as.numeric(unlist(step1_real[,3])) * mean(as.numeric(unlist(step1_real[,8])))) / ((as.numeric(unlist(step1_real[,8]))) * mean(as.numeric(unlist(step1_real[,3]))))
    Normalized_Exp_3 <- (as.numeric(unlist(step1_real[,4])) * mean(as.numeric(unlist(step1_real[,9])))) / ((as.numeric(unlist(step1_real[,9]))) * mean(as.numeric(unlist(step1_real[,4]))))
    Normalized_Exp_4 <- (as.numeric(unlist(step1_real[,5])) * mean(as.numeric(unlist(step1_real[,10])))) / ((as.numeric(unlist(step1_real[,10]))) * mean(as.numeric(unlist(step1_real[,5]))))
    Normalized_Exp_5 <- (as.numeric(unlist(step1_real[,6])) * mean(as.numeric(unlist(step1_real[,11])))) / ((as.numeric(unlist(step1_real[,11]))) * mean(as.numeric(unlist(step1_real[,6]))))

  }

  step2 <-
    cbind(step1_real,
          Normalized_Exp_1,
          Normalized_Exp_2,
          Normalized_Exp_3,
          Normalized_Exp_4,
          Normalized_Exp_5)

  # Step 3

  Log_2_Exp_1 <- log(Normalized_Exp_1, 2)
  Log_2_Exp_2 <- log(Normalized_Exp_2, 2)
  Log_2_Exp_3 <- log(Normalized_Exp_3, 2)
  Log_2_Exp_4 <- log(Normalized_Exp_4, 2)
  Log_2_Exp_5 <- log(Normalized_Exp_5, 2)

  step3 <- cbind(step2, Log_2_Exp_1, Log_2_Exp_2, Log_2_Exp_3, Log_2_Exp_4, Log_2_Exp_5)

  # Step 4

  Exp_Log_Avg <- rowMeans(step3[, c("Log_2_Exp_1", "Log_2_Exp_2", "Log_2_Exp_3", "Log_2_Exp_4", "Log_2_Exp_5")], na.rm = T)
  step4 <- cbind(step3, Exp_Log_Avg)

  # Step 5

  c1 <- step4[, c("Log_2_Exp_1", "Log_2_Exp_2", "Log_2_Exp_3", "Log_2_Exp_4", "Log_2_Exp_5")]
  SD_Exp <- apply(c1, 1, sd)
  step5 <- cbind(step4, SD_Exp)

  # Step 6

  f1 <- step5$Exp_Log_Avg
  mean_log_avg_Exp <- mean(f1)

  sd_log_avg_Exp <- sd(f1)
  z_score_Exp <- ((f1 - mean_log_avg_Exp)) / (sd_log_avg_Exp)

  step6 <- cbind(step5, z_score_Exp)

  # Step 7

  g1 <- abs(Exp_Log_Avg)
  g2 <- step6$SD_Exp

  T_value_Exp <- (g1 * sqrt(5) / g2)
  step7 <- cbind(step6, T_value_Exp)

  # Step 8

  h1 <- step7$Exp_Log_Avg
  h2 <- mean(h1)
  h3 <- step7$T_value_Exp
  p_value_Exp <- 2 * pt(h3, 4, lower.tail = F)
  step8 <- cbind(step7, p_value_Exp)

  # Step 9

  j1 <- step8$p_value_Exp
  log_10_Exp <- -log10(j1)
  step9 <- cbind(step8, log_10_Exp)

  # Step 10

  Hit_ident1_Exp <- magrittr::"%>%" (step9,
                                     dplyr::select(tidyselect::any_of(c("Description", "Accession"))))
  Log_Avg_Exp <- step9$Exp_Log_Avg
  Log_10_Exp <- step9$log_10_Exp
  Hit_ident_Exp <- cbind(Hit_ident1_Exp, Log_Avg_Exp, Log_10_Exp, Log_2_Exp_1, Log_2_Exp_2, Log_2_Exp_3, Log_2_Exp_4, Log_2_Exp_5, z_score_Exp)

  Expunique <- dplyr::n_distinct(step9$Accession)

  if (correctedpvalue == FALSE){
  n1 <-
    (Hit_ident_Exp$z_score_Exp > SD_cutoff &
       Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff)) |
    (Hit_ident_Exp$z_score_Exp < -SD_cutoff & Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff))
  n2 <- as.data.frame(n1)
  }

  if (correctedpvalue == TRUE){
    n1 <-
      (Hit_ident_Exp$z_score_Exp > SD_cutoff &
         Hit_ident_Exp$Log_10_Exp > -log10((p_cutoff/Expunique))) |
      (Hit_ident_Exp$z_score_Exp < -SD_cutoff & Hit_ident_Exp$Log_10_Exp > -log10((p_cutoff/Expunique)))
    n2 <- as.data.frame(n1)
  }

  Sig_Check_Exp <-
    dplyr::if_else(n2$n1 == "TRUE", "Significant", "Not Significant")
  step10 <- cbind(Hit_ident_Exp, p_value_Exp, Sig_Check_Exp)

  ddkd <- step10$Sig_Check_Exp

  as.data.frame(ddkd)
  abkk <- (ddkd == "Significant")
  Exp_Number_Of_Hits <- length(which(abkk))


  Hit_List <- step10[step10$Sig_Check_Exp == "Significant", ]
  ExpExport <- unique(Hit_List$Accession)
  ExpSigExport <- magrittr::"%>%" (step10,
                                   dplyr::select(tidyselect::any_of(c("Accession", "Description", "z_score_Exp", "p_value_Exp", "Sig_Check_Exp"))))
  Hit_ListExport <- magrittr::"%>%" (Hit_List,
                                     dplyr::select(tidyselect::any_of(c("Accession", "Description", "z_score_Exp", "p_value_Exp", "Sig_Check_Exp"))))

  print(ExpSigExport)
  print(paste("There are", Exp_Number_Of_Hits, "Hits"))
  print(paste(Expunique, "Unique Proteins were assayed"))
  print(Hit_ListExport)


  utils::write.table(ExpExport, file = "ExpressionLevelPhenotypeUniqueHits.csv", row.names = FALSE, col.names = FALSE)
  utils::write.table(ExpSigExport, file = "ExpressionLevelPhenotypeOut.csv", row.names = FALSE)
  utils::write.table(Hit_ListExport, file = "ExpressionLevelPhenotypeHitsOut.csv", row.names = FALSE)


  if (plot == TRUE){

  Volcano_Plot <- ggplot2::ggplot(step10,
                                  ggplot2::aes (
                                    x = z_score_Exp ,
                                    y = Log_10_Exp ,
                                    label = Accession,
                                    colour = Sig_Check_Exp
                                  )) + ggplot2::geom_point(size = 1.5, alpha = alpha) +
    ggplot2::labs(x= xlab,y= ylab) + ggplot2::scale_colour_manual(breaks = c("Not Significant", "Significant") , values = labelcolor) +
    ggplot2::expand_limits(x=0, y=0) +
    ggplot2::scale_x_continuous(breaks = -12:12) +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = plottextsize)) + ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  print(Volcano_Plot)
  }

ggplot2::ggsave("Volcano_Plot_Expression_Phenotype.png", plot = Volcano_Plot)

if (interactiveplot == TRUE){
  plotly::ggplotly(Volcano_Plot)
}

}

#' OnePotTPP Phenotype Analysis
#'
#' This function will take expression data and OnePotTPP (TMT10plex) data directly from ProteomeDiscoverer (v. 2.3) and generate a volcano plot and hit list
#' based on user-defined SD and p-values
#'
#' @param Expression_data Raw ProteomeDiscoverer (v. 2.3) data
#' @param TPP_data Raw OnePotTPP (TMT10plex) ProteomeDiscoverer (v. 2.3) data
#' @param SD_cutoff A positive integer (default = 2)
#' @param p_cutoff A positive integer between 0 and 1 (default = 0.05)
#' @param TMTplex 8, 10, 16 TMTplex (default = 10)
#' @param plot Display volcano plot (default = TRUE)
#' @param xlab X-axis Label (default = "Z Score")
#' @param ylab Y-axis Label (default = "- log10 (p value)")
#' @param labelcolor Color of Not significant and Significant data points (default = c("grey", "red"))
#' @param alpha Transparency of plot (default = 0.5)
#' @param plottextsize Text Size of plot (default = 12)
#' @param correctedpvalue Bonferroni corrected p-value (default = FALSE)
#' @param interactiveplot Plotly interactive Volcano plot (default = FALSE)
#' @return Volcano plot, Hit list, Total proteins assayed, and Unique protein hits
#' @export

OnePotTPP_Phenotype.Fun <- function(Expression_data, TPP_Raw, SD_cutoff = 2, p_cutoff = 0.05, TMTplex = 10, plot = TRUE, xlab = "Z Score", ylab = "- log10 (p value)", labelcolor = c("grey", "red"), alpha = 0.5, plottextsize = 12, correctedpvalue = FALSE, interactiveplot = FALSE){

  if (TMTplex !=10)
    stop("This TMTplex is not yet supported in this version.\n")

  if (SD_cutoff < 0)
    stop("SD_cutoff must be non-negative.\n")

  if (p_cutoff < 0)
    stop("p_value cutoff must be non-negative.\n")

# Define SD and p-value cutoffs
# Load "Expression_data" and "TPP_Raw"

options(max.print = 25000)

# We can vary this to select for High | Medium | Low or both

High_Confidence_Exp_Level <-
  dplyr::filter(Expression_data, Expression_data[,grepl(("Confidence"), names(Expression_data))] == "High")

# This dataframe takes all useful columns from the raw data and will be outputed at the end for reference

Exp_Level_Concise <- magrittr::"%>%" (High_Confidence_Exp_Level,
                                      dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Description", "Accession", (grep(("^Abundances\\s\\(Grouped):"), names(High_Confidence_Exp_Level), value = TRUE))))))

# Remove all blank rows in Dataframe

NA_Removed <- tidyr::drop_na(Exp_Level_Concise)

# Step 1

Exp_Level_1 <- magrittr::"%>%" (NA_Removed,
                                dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Description", "Accession"))))
Exp_Level_2 <- magrittr::"%>%" (NA_Removed,
                                dplyr::select(tidyselect::any_of(grep(("^Abundances\\s\\(Grouped):"), names(High_Confidence_Exp_Level), value = TRUE))))
step1 <- cbind(Exp_Level_1, Exp_Level_2)

if ("Master Protein Accessions" %in% names(step1)){
  step1 <- dplyr::rename(step1, "Accession" = "Master Protein Accessions")
}


# Step 2

if (("Accession" %in% names(step1) && "Description" %in% names(step1))) {

  step1_real <-
    dplyr::filter(
      step1,
      step1[,3] > 0,
      step1[,4] > 0,
      step1[,5] > 0,
      step1[,6] > 0,
      step1[,7] > 0,
      step1[,8] > 0,
      step1[,9] > 0,
      step1[,10] > 0,
      step1[,11] > 0,
      step1[,12] > 0)

  Normalized_Exp_1 <- (as.numeric(unlist(step1_real[,3])) * mean(as.numeric(unlist(step1_real[,8])))) / ((as.numeric(unlist(step1_real[,8]))) * mean(as.numeric(unlist(step1_real[,3]))))
  Normalized_Exp_2 <- (as.numeric(unlist(step1_real[,4])) * mean(as.numeric(unlist(step1_real[,9])))) / ((as.numeric(unlist(step1_real[,9]))) * mean(as.numeric(unlist(step1_real[,4]))))
  Normalized_Exp_3 <- (as.numeric(unlist(step1_real[,5])) * mean(as.numeric(unlist(step1_real[,10])))) / ((as.numeric(unlist(step1_real[,10]))) * mean(as.numeric(unlist(step1_real[,5]))))
  Normalized_Exp_4 <- (as.numeric(unlist(step1_real[,6])) * mean(as.numeric(unlist(step1_real[,11])))) / ((as.numeric(unlist(step1_real[,11]))) * mean(as.numeric(unlist(step1_real[,6]))))
  Normalized_Exp_5 <- (as.numeric(unlist(step1_real[,7])) * mean(as.numeric(unlist(step1_real[,12])))) / ((as.numeric(unlist(step1_real[,12]))) * mean(as.numeric(unlist(step1_real[,7]))))

}

if (!"Description" %in% names(step1)){

  step1_real <-
    dplyr::filter(
      step1,
      step1[,2] > 0,
      step1[,3] > 0,
      step1[,4] > 0,
      step1[,5] > 0,
      step1[,6] > 0,
      step1[,7] > 0,
      step1[,8] > 0,
      step1[,9] > 0,
      step1[,10] > 0,
      step1[,11] > 0)

  Normalized_Exp_1 <- (as.numeric(unlist(step1_real[,2])) * mean(as.numeric(unlist(step1_real[,7])))) / ((as.numeric(unlist(step1_real[,7]))) * mean(as.numeric(unlist(step1_real[,2]))))
  Normalized_Exp_2 <- (as.numeric(unlist(step1_real[,3])) * mean(as.numeric(unlist(step1_real[,8])))) / ((as.numeric(unlist(step1_real[,8]))) * mean(as.numeric(unlist(step1_real[,3]))))
  Normalized_Exp_3 <- (as.numeric(unlist(step1_real[,4])) * mean(as.numeric(unlist(step1_real[,9])))) / ((as.numeric(unlist(step1_real[,9]))) * mean(as.numeric(unlist(step1_real[,4]))))
  Normalized_Exp_4 <- (as.numeric(unlist(step1_real[,5])) * mean(as.numeric(unlist(step1_real[,10])))) / ((as.numeric(unlist(step1_real[,10]))) * mean(as.numeric(unlist(step1_real[,5]))))
  Normalized_Exp_5 <- (as.numeric(unlist(step1_real[,6])) * mean(as.numeric(unlist(step1_real[,11])))) / ((as.numeric(unlist(step1_real[,11]))) * mean(as.numeric(unlist(step1_real[,6]))))

}

step2 <-
  cbind(step1_real,
        Normalized_Exp_1,
        Normalized_Exp_2,
        Normalized_Exp_3,
        Normalized_Exp_4,
        Normalized_Exp_5)

# Step 3

Log_2_Exp_1 <- log(Normalized_Exp_1, 2)
Log_2_Exp_2 <- log(Normalized_Exp_2, 2)
Log_2_Exp_3 <- log(Normalized_Exp_3, 2)
Log_2_Exp_4 <- log(Normalized_Exp_4, 2)
Log_2_Exp_5 <- log(Normalized_Exp_5, 2)

step3 <- cbind(step2, Log_2_Exp_1, Log_2_Exp_2, Log_2_Exp_3, Log_2_Exp_4, Log_2_Exp_5)

# Step 4

Exp_Log_Avg <- rowMeans(step3[, c("Log_2_Exp_1", "Log_2_Exp_2", "Log_2_Exp_3", "Log_2_Exp_4", "Log_2_Exp_5")], na.rm = T)
step4 <- cbind(step3, Exp_Log_Avg)

# Step 5

c1 <- step4[, c("Log_2_Exp_1", "Log_2_Exp_2", "Log_2_Exp_3", "Log_2_Exp_4", "Log_2_Exp_5")]
SD_Exp <- apply(c1, 1, sd)
step5 <- cbind(step4, SD_Exp)

# Step 6

f1 <- step5$Exp_Log_Avg
mean_log_avg_Exp <- mean(f1)
sd_log_avg_Exp <- sd(f1)

z_score_Exp <- ((f1 - mean_log_avg_Exp)) / (sd_log_avg_Exp)
step6 <- cbind(step5, z_score_Exp)

# Step 7

g1 <- abs(Exp_Log_Avg)
g2 <- step6$SD_Exp

T_value_Exp <- (g1 * sqrt(5) / g2)
step7 <- cbind(step6, T_value_Exp)

# Step 8

h1 <- step7$Exp_Log_Avg
h2 <- mean(h1)
h3 <- step7$T_value_Exp
p_value_Exp <- 2 * pt(h3, 4, lower.tail = F)
step8 <- cbind(step7, p_value_Exp)

# Step 9

j1 <- step8$p_value_Exp
log_10_Exp <- -log10(j1)
step9 <- cbind(step8, log_10_Exp)

# Step 10

Hit_ident1_Exp <- magrittr::"%>%" (step9,
                                   dplyr::select(tidyselect::any_of(c("Description", "Accession"))))
Log_Avg_Exp <- step9$Exp_Log_Avg
Log_10_Exp <- step9$log_10_Exp
Hit_ident_Exp <- cbind(Hit_ident1_Exp, Log_Avg_Exp, Log_10_Exp, Log_2_Exp_1, Log_2_Exp_2, Log_2_Exp_3, Log_2_Exp_4, Log_2_Exp_5, z_score_Exp)

n1 <-
  (Hit_ident_Exp$z_score_Exp > SD_cutoff &
     Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff)) |
  (Hit_ident_Exp$z_score_Exp < -SD_cutoff & Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff))
n2 <- as.data.frame(n1)

Sig_Check_Exp <-
  dplyr::if_else(n2$n1 == "TRUE", "Significant", "Not Significant")
step10 <- cbind(Hit_ident_Exp, p_value_Exp, Sig_Check_Exp)

# TPP_Data

TPP_Confidence <- dplyr::filter(TPP_Raw, TPP_Raw[,grepl(("Confidence"), names(TPP_Raw))] == "High")

Data2 <- magrittr::"%>%" (TPP_Confidence, dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Accession", "Description", (grep(("^Abundances\\s\\(Grouped):"), names(TPP_Confidence), value = TRUE))))))

# Remove all blank rows in Dataframe

NA_Removed_TPP <- tidyr::drop_na(Data2)

TPP_1 <- magrittr::"%>%" (NA_Removed_TPP,
                          dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Accessions", "Description"))))
TPP_2 <- magrittr::"%>%" (NA_Removed_TPP,
                          dplyr::select(tidyselect::any_of(grep(("^Abundances\\s\\(Grouped):"), names(NA_Removed_TPP), value = TRUE))))
step1_TPP <- cbind(TPP_1, TPP_2)

if ("Master Protein Accessions" %in% names(step1_TPP)){
  step1_TPP <- dplyr::rename(step1_TPP, "Accession" = "Master Protein Accessions")
}

TPP_merger <- merge(Hit_ident_Exp, step1_TPP)

TPP_merge <-
  dplyr::filter(
    TPP_merger,
    TPP_merger[,11] > 0,
    TPP_merger[,12] > 0,
    TPP_merger[,13] > 0,
    TPP_merger[,14] > 0,
    TPP_merger[,15] > 0,
    TPP_merger[,16] > 0,
    TPP_merger[,17] > 0,
    TPP_merger[,18] > 0,
    TPP_merger[,19] > 0,
    TPP_merger[,20] > 0)


# Step 11

Normalized_TPP_T1 <- (as.numeric(unlist(TPP_merge[,11])) * mean(as.numeric(unlist(TPP_merge[,16])))) / ((as.numeric(unlist(TPP_merge[,16]))) * mean(as.numeric(unlist(TPP_merge[,11]))))
Normalized_TPP_T2 <- (as.numeric(unlist(TPP_merge[,12])) * mean(as.numeric(unlist(TPP_merge[,17])))) / ((as.numeric(unlist(TPP_merge[,17]))) * mean(as.numeric(unlist(TPP_merge[,12]))))
Normalized_TPP_T3 <- (as.numeric(unlist(TPP_merge[,13])) * mean(as.numeric(unlist(TPP_merge[,18])))) / ((as.numeric(unlist(TPP_merge[,18]))) * mean(as.numeric(unlist(TPP_merge[,13]))))
Normalized_TPP_T4 <- (as.numeric(unlist(TPP_merge[,14])) * mean(as.numeric(unlist(TPP_merge[,19])))) / ((as.numeric(unlist(TPP_merge[,19]))) * mean(as.numeric(unlist(TPP_merge[,14]))))
Normalized_TPP_T5 <- (as.numeric(unlist(TPP_merge[,15])) * mean(as.numeric(unlist(TPP_merge[,20])))) / ((as.numeric(unlist(TPP_merge[,20]))) * mean(as.numeric(unlist(TPP_merge[,15]))))

step11 <-
  cbind(TPP_merge,
        Normalized_TPP_T1,
        Normalized_TPP_T2,
        Normalized_TPP_T3,
        Normalized_TPP_T4,
        Normalized_TPP_T5)

step111 <- tidyr::drop_na(step11)

# Step 12

Log_2_TPP_T1 <- log(Normalized_TPP_T1, 2)
Log_2_TPP_T2 <- log(Normalized_TPP_T2, 2)
Log_2_TPP_T3 <- log(Normalized_TPP_T3, 2)
Log_2_TPP_T4 <- log(Normalized_TPP_T4, 2)
Log_2_TPP_T5 <- log(Normalized_TPP_T5, 2)

step12 <- cbind(step111, Log_2_TPP_T1, Log_2_TPP_T2, Log_2_TPP_T3, Log_2_TPP_T4, Log_2_TPP_T5)

# Step 13

Log_Avg_TPP <- rowMeans(step12[, c("Log_2_TPP_T1", "Log_2_TPP_T2", "Log_2_TPP_T3", "Log_2_TPP_T4", "Log_2_TPP_T5")])
step13 <- cbind(step12, Log_Avg_TPP)

# Step  14

rep1_exp1 <- step13$Log_2_TPP_T1 - step13$Log_2_Exp_1
rep2_exp2 <- step13$Log_2_TPP_T2 - step13$Log_2_Exp_2
rep3_exp3 <- step13$Log_2_TPP_T3 - step13$Log_2_Exp_3
rep4_exp4 <- step13$Log_2_TPP_T4 - step13$Log_2_Exp_4
rep5_exp5 <- step13$Log_2_TPP_T5 - step13$Log_2_Exp_5
step14 <- cbind(step13, rep1_exp1, rep2_exp2, rep3_exp3, rep4_exp4, rep5_exp5)

# Step 15

TPP_Avg <- rowMeans(step14[, c("rep1_exp1", "rep2_exp2", "rep3_exp3", "rep4_exp4", "rep5_exp5")], na.rm = T)
step15 <- cbind(step14, TPP_Avg)

# Step 16

k1 <- step15[, c("rep1_exp1", "rep2_exp2", "rep3_exp3", "rep4_exp4", "rep5_exp5")]
SD_TPP <- apply(k1, 1, sd)
step16 <- cbind(step15, SD_TPP)

# Step 17

step17 <- tidyr::drop_na(step16)
TPPunique <- dplyr::n_distinct(step17$Accession)
# Step 18

y1 <- step17$TPP_Avg
TPP_mean_avg <- mean(y1)
TPP_sd_avg <- sd(y1)
z_score_TPP <- ((y1 - TPP_mean_avg)) / (TPP_sd_avg)
step18 <- cbind(step17, z_score_TPP)

# Step 19

hg1 <- abs(TPP_Avg)
hg2 <- step18$SD_TPP

TPP_T_value <- (hg1 * sqrt(5) / hg2)
step19 <- cbind(step18, TPP_T_value)

# Step 20

hk1 <- step19$TPP_Avg
hk2 <- mean(hg1)
hk3 <- step19$TPP_T_value
p_value_TPP <- 2 * pt(hk3, 4, lower.tail = F)
step20 <- cbind(step19, p_value_TPP)

# Step 21

jk1 <- step20$p_value_TPP
TPP_log_10 <- -log10(jk1)
step21 <- cbind(step20, TPP_log_10)

TPP_Hit_ident1 <- magrittr::"%>%" (step21,
                                   dplyr::select(tidyselect::any_of(c("Accession", "Description"))))
TPP_Log_Avg <- step21[, c("z_score_TPP", "p_value_TPP")]
TPP_Hit_ident <- cbind(TPP_Hit_ident1, TPP_Log_Avg, TPP_log_10)


# Step 22

if (correctedpvalue == FALSE){

nk1 <-
  (TPP_Hit_ident$z_score_TPP > SD_cutoff &
     TPP_Hit_ident$TPP_log_10 > -log10(p_cutoff)) |
  (TPP_Hit_ident$z_score_TPP < -SD_cutoff & TPP_Hit_ident$TPP_log_10 > -log10(p_cutoff))
nk2 <- as.data.frame(nk1)
}

if (correctedpvalue == TRUE){
  nk1 <-
    (TPP_Hit_ident$z_score_TPP > SD_cutoff &
       TPP_Hit_ident$TPP_log_10 > -log10((p_cutoff/TPPunique))) |
    (TPP_Hit_ident$z_score_TPP < -SD_cutoff & TPP_Hit_ident$TPP_log_10 > -log10((p_cutoff/TPPunique)))
  nk2 <- as.data.frame(nk1)

}

Sig_Check_TPP <-
  dplyr::if_else(nk2$nk1 == "TRUE", "Significant", "Not Significant")
step22 <- cbind(TPP_Hit_ident, Sig_Check_TPP)

ddk <- step22$Sig_Check_TPP

as.data.frame(ddk)
abk <- (ddk == "Significant")
TPP_Number_Of_Hits <- length(which(abk))

Hit_List <- step22[step22$Sig_Check_TPP == "Significant", ]
TPPuniquehits <- dplyr::n_distinct(Hit_List$Accession)
TPPExport <- unique(Hit_List$Accession)
TPPSigExport <- magrittr::"%>%" (step22,
                                 dplyr::select(tidyselect::any_of(c("Accession", "Description", "z_score_TPP", "p_value_TPP", "Sig_Check_TPP"))))
Hit_ListExportTPP <- magrittr::"%>%" (Hit_List,
                                      dplyr::select(tidyselect::any_of(c("Accession", "Description", "z_score_TPP", "p_value_TPP", "Sig_Check_TPP"))))


print(TPPSigExport)
print(paste("There are", TPP_Number_Of_Hits, "Hits"))
print(paste("There are", TPPuniquehits, "Unique Hits"))
print(paste(TPPunique, "Unique Proteins were assayed"))
print(Hit_ListExportTPP)

utils::write.table(TPPExport, file = "OnePotTPPPhenotypeUniqueHits.csv", row.names = FALSE, col.names = FALSE)
utils::write.table(TPPSigExport, file = "OnePotTPPPhenotypeOut.csv", row.names = FALSE)
utils::write.table(Hit_ListExportTPP, file = "OnePotTPPPhenotypeHitsOut.csv", row.names = FALSE)

# Visualizing Data

if (plot == TRUE){

Volcano_Plot <- ggplot2::ggplot(step22,
         ggplot2::aes (
           x = z_score_TPP ,
           y = TPP_log_10 ,
           label = Accession,
           colour = Sig_Check_TPP
         )) + ggplot2::geom_point(size = 1.5, alpha = alpha) +
  ggplot2::labs(x= xlab, y= ylab) + ggplot2::scale_colour_manual(breaks = c("Not Significant", "Significant") ,values = labelcolor) +
  ggplot2::expand_limits(x=0, y=0) +
  ggplot2::scale_x_continuous(breaks = -12:12) +
  ggplot2::theme_bw()+ ggplot2::theme(text = ggplot2::element_text(size = plottextsize)) + ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
  ggplot2::theme(legend.title = ggplot2::element_blank())

print(Volcano_Plot)
}

ggplot2::ggsave("Volcano_Plot_TPP_Phenotype.png", plot = Volcano_Plot)



if (interactiveplot == TRUE){
  plotly::ggplotly(Volcano_Plot)
}

}


#' OnePotTPP Ligand Analysis
#'
#' This function will take OnePotTPP (TMT10plex) data directly from ProteomeDiscoverer (v. 2.3) and generate a volcano plot and hit list
#' based on user-defined SD and p-values
#'
#' @param Expression_data Raw ProteomeDiscoverer (v. 2.3) data
#' @param TPP_data Raw OnePotTPP (TMT10plex) ProteomeDiscoverer (v. 2.3) data
#' @param SD_cutoff A positive integer (default = 2)
#' @param p_cutoff A positive integer between 0 and 1 (default = 0.05)
#' @param TMTplex 8, 10, 16 TMTplex (default = 10)
#' @param plot Display volcano plot (default = TRUE)
#' @param xlab X-axis Label (default = "Z Score")
#' @param ylab Y-axis Label (default = "- log10 (p value)")
#' @param labelcolor Color of Not significant and Significant data points (default = c("grey", "red"))
#' @param alpha Transparency of plot (default = 0.5)
#' @param plottextsize Text Size of plot (default = 12)
#' @param correctedpvalue Bonferroni corrected p-value (default = FALSE)
#' @param interactiveplot Plotly interactive Volcano plot (default = FALSE)
#' @return Volcano plot, Hit list, Total proteins assayed, and Unique protein hits
#' @export

OnePotTPP_Ligand.Fun <- function(TPP_Raw, SD_cutoff = 2, p_cutoff = 0.05, TMTplex = 10, plot = TRUE, xlab = "Z Score", ylab = "- log10 (p value)", labelcolor = c("grey", "red"), alpha = 0.5, plottextsize = 12, correctedpvalue = FALSE, interactiveplot = FALSE){

  if (TMTplex !=10)
    stop("This TMTplex is not yet supported in this version.\n")

  if (SD_cutoff < 0)
    stop("SD_cutoff must be non-negative.\n")

  if (p_cutoff < 0)
    stop("p_value cutoff must be non-negative.\n")

  # Define SD and p-value cutoffs
  # Load "TPP_Raw"

  options(max.print = 25000)

  # We can vary this to select for High | Medium | Low or both

  # TPP_Data

  TPP_Confidence <- dplyr::filter(TPP_Raw, TPP_Raw[,grepl("Confidence", names(TPP_Raw))] == "High")

  Data2 <- magrittr::"%>%" (TPP_Confidence,
                            dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Accession", "Description", (grep(("^Abundances\\s\\(Grouped):"), names(TPP_Confidence), value = TRUE))))))

  if ("Master Protein Accessions" %in% names(Data2)){
    Data2 <- dplyr::rename(Data2, "Accession" = "Master Protein Accessions")
  }

  # Remove all blank rows in Dataframe

  NA_Removed_TPP <- tidyr::drop_na(Data2)

  TPP_1 <- magrittr::"%>%" (NA_Removed_TPP,
                            dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Accession", "Description"))))
  TPP_2 <- magrittr::"%>%" (NA_Removed_TPP,
                            dplyr::select(tidyselect::any_of(grep("^Abundances\\s\\(Grouped):", names(NA_Removed_TPP), value = TRUE))))
  step1_TPP <- cbind(TPP_1, TPP_2)

  if (("Accession" %in% names(step1_TPP) && "Description" %in% names(step1_TPP))) {

    step1_real <-
      dplyr::filter(
        step1_TPP,
        step1_TPP[,3] > 0,
        step1_TPP[,4] > 0,
        step1_TPP[,5] > 0,
        step1_TPP[,6] > 0,
        step1_TPP[,7] > 0,
        step1_TPP[,8] > 0,
        step1_TPP[,9] > 0,
        step1_TPP[,10] > 0,
        step1_TPP[,11] > 0,
        step1_TPP[,12] > 0)

    Normalized_TPP_T1 <- (as.numeric(unlist(step1_real[,3])) * mean(as.numeric(unlist(step1_real[,8])))) / ((as.numeric(unlist(step1_real[,8]))) * mean(as.numeric(unlist(step1_real[,3]))))
    Normalized_TPP_T2 <- (as.numeric(unlist(step1_real[,4])) * mean(as.numeric(unlist(step1_real[,9])))) / ((as.numeric(unlist(step1_real[,9]))) * mean(as.numeric(unlist(step1_real[,4]))))
    Normalized_TPP_T3 <- (as.numeric(unlist(step1_real[,5])) * mean(as.numeric(unlist(step1_real[,10])))) / ((as.numeric(unlist(step1_real[,10]))) * mean(as.numeric(unlist(step1_real[,5]))))
    Normalized_TPP_T4 <- (as.numeric(unlist(step1_real[,6])) * mean(as.numeric(unlist(step1_real[,11])))) / ((as.numeric(unlist(step1_real[,11]))) * mean(as.numeric(unlist(step1_real[,6]))))
    Normalized_TPP_T5 <- (as.numeric(unlist(step1_real[,7])) * mean(as.numeric(unlist(step1_real[,12])))) / ((as.numeric(unlist(step1_real[,12]))) * mean(as.numeric(unlist(step1_real[,7]))))

  }

  if (!"Description" %in% names(step1_TPP)){

    step1_real <-
      dplyr::filter(
        step1_TPP,
        step1_TPP[,2] > 0,
        step1_TPP[,3] > 0,
        step1_TPP[,4] > 0,
        step1_TPP[,5] > 0,
        step1_TPP[,6] > 0,
        step1_TPP[,7] > 0,
        step1_TPP[,8] > 0,
        step1_TPP[,9] > 0,
        step1_TPP[,10] > 0,
        step1_TPP[,11] > 0)

    Normalized_TPP_T1 <- (as.numeric(unlist(step1_real[,2])) * mean(as.numeric(unlist(step1_real[,7])))) / ((as.numeric(unlist(step1_real[,7]))) * mean(as.numeric(unlist(step1_real[,2]))))
    Normalized_TPP_T2 <- (as.numeric(unlist(step1_real[,3])) * mean(as.numeric(unlist(step1_real[,8])))) / ((as.numeric(unlist(step1_real[,8]))) * mean(as.numeric(unlist(step1_real[,3]))))
    Normalized_TPP_T3 <- (as.numeric(unlist(step1_real[,4])) * mean(as.numeric(unlist(step1_real[,9])))) / ((as.numeric(unlist(step1_real[,9]))) * mean(as.numeric(unlist(step1_real[,4]))))
    Normalized_TPP_T4 <- (as.numeric(unlist(step1_real[,5])) * mean(as.numeric(unlist(step1_real[,10])))) / ((as.numeric(unlist(step1_real[,10]))) * mean(as.numeric(unlist(step1_real[,5]))))
    Normalized_TPP_T5 <- (as.numeric(unlist(step1_real[,6])) * mean(as.numeric(unlist(step1_real[,11])))) / ((as.numeric(unlist(step1_real[,11]))) * mean(as.numeric(unlist(step1_real[,6]))))

  }

  # Step 11

  step11 <-
    cbind(step1_real,
          Normalized_TPP_T1,
          Normalized_TPP_T2,
          Normalized_TPP_T3,
          Normalized_TPP_T4,
          Normalized_TPP_T5)

  step111 <- tidyr::drop_na(step11)

  # Step 12

  Log_2_TPP_T1 <- log(Normalized_TPP_T1, 2)
  Log_2_TPP_T2 <- log(Normalized_TPP_T2, 2)
  Log_2_TPP_T3 <- log(Normalized_TPP_T3, 2)
  Log_2_TPP_T4 <- log(Normalized_TPP_T4, 2)
  Log_2_TPP_T5 <- log(Normalized_TPP_T5, 2)

  step12 <- cbind(step111, Log_2_TPP_T1, Log_2_TPP_T2, Log_2_TPP_T3, Log_2_TPP_T4, Log_2_TPP_T5)

  # Step 13

  Log_Avg_TPP <- rowMeans(step12[, c("Log_2_TPP_T1", "Log_2_TPP_T2", "Log_2_TPP_T3", "Log_2_TPP_T4", "Log_2_TPP_T5")])
  step13 <- cbind(step12, Log_Avg_TPP)

  # Step  14

  step14 <- step13

  # Step 15

  step15 <- step14

  # Step 16

  k1 <- step15[, c("Log_2_TPP_T1", "Log_2_TPP_T2", "Log_2_TPP_T3", "Log_2_TPP_T4", "Log_2_TPP_T5")]
  SD_TPP <- apply(k1, 1, sd)
  step16 <- cbind(step15, SD_TPP)


  # Step 17

  step17 <- tidyr::drop_na(step16)
  TPPunique <- dplyr::n_distinct(step17$Accession)

  # Step 18

  y1 <- step17$Log_Avg_TPP
  TPP_mean_avg <- mean(y1)
  TPP_sd_avg <- sd(y1)
  z_score_TPP <- ((y1 - TPP_mean_avg)) / (TPP_sd_avg)
  step18 <- cbind(step17, z_score_TPP)

  # Step 19

  hg1 <- abs(step18$Log_Avg_TPP)
  hg2 <- step18$SD_TPP

  TPP_T_value <- (hg1 * sqrt(5) / hg2)
  step19 <- cbind(step18, TPP_T_value)

  # Step 20

  hk1 <- step19$Log_Avg_TPP
  hk2 <- mean(hg1)
  hk3 <- step19$TPP_T_value
  p_value_TPP <- 2 * pt(hk3, 4, lower.tail = F)
  step20 <- cbind(step19, p_value_TPP)

  # Step 21

  jk1 <- step20$p_value_TPP
  TPP_log_10 <- -log10(jk1)
  step21 <- cbind(step20, TPP_log_10)

  TPP_Hit_ident1 <- magrittr::"%>%" (step21,
                                     dplyr::select(tidyselect::any_of(c("Accession", "Description"))))
  TPP_Log_Avg <- step21[, c("z_score_TPP", "p_value_TPP")]
  TPP_Hit_ident <- cbind(TPP_Hit_ident1, TPP_Log_Avg, TPP_log_10)


  # Step 22

if (correctedpvalue == FALSE){

  nk1 <-
    (TPP_Hit_ident$z_score_TPP > SD_cutoff &
       TPP_Hit_ident$TPP_log_10 > -log10(p_cutoff)) |
    (TPP_Hit_ident$z_score_TPP < -SD_cutoff & TPP_Hit_ident$TPP_log_10 > -log10(p_cutoff))
  nk2 <- as.data.frame(nk1)

}

  if (correctedpvalue == TRUE){
    nk1 <-
      (TPP_Hit_ident$z_score_TPP > SD_cutoff &
         TPP_Hit_ident$TPP_log_10 > -log10((p_cutoff/TPPunique))) |
      (TPP_Hit_ident$z_score_TPP < -SD_cutoff & TPP_Hit_ident$TPP_log_10 > -log10((p_cutoff/TPPunique)))
    nk2 <- as.data.frame(nk1)
  }

    Sig_Check_TPP <-
    dplyr::if_else(nk2$nk1 == "TRUE", "Significant", "Not Significant")
  step22 <- cbind(TPP_Hit_ident, Sig_Check_TPP)

  ddk <- step22$Sig_Check_TPP

  as.data.frame(ddk)
  abk <- (ddk == "Significant")
  TPP_Number_Of_Hits <- length(which(abk))

  Hit_List <- step22[step22$Sig_Check_TPP == "Significant", ]
  TPPuniquehits <- dplyr::n_distinct(Hit_List$Accession)
  TPPExport <- unique(Hit_List$Accession)
  TPPSigExport <- magrittr::"%>%" (step22,
                                   dplyr::select(tidyselect::any_of(c("Accession", "Description", "z_score_TPP", "p_value_TPP", "Sig_Check_TPP"))))
  Hit_ListExportTPP <- magrittr::"%>%" (Hit_List,
                                        dplyr::select(tidyselect::any_of(c("Accession", "Description", "z_score_TPP", "p_value_TPP", "Sig_Check_TPP"))))


  print(TPPSigExport)
  print(paste("There are", TPP_Number_Of_Hits, "Hits"))
  print(paste("There are", TPPuniquehits, "Unique Hits"))
  print(paste(TPPunique, "Unique Proteins were assayed"))
  print(Hit_ListExportTPP)

  utils::write.table(TPPExport, file = "OnePotTPPLigandUniqueHits.csv", row.names = FALSE, col.names = FALSE)
  utils::write.table(TPPSigExport, file = "OnePotTPPLigandOut.csv", row.names = FALSE)
  utils::write.table(Hit_ListExportTPP, file = "OnePotTPPLigandHitsOut.csv", row.names = FALSE)

  # Visualizing Data

  if (plot == TRUE){

    Volcano_Plot <- ggplot2::ggplot(step22,
                                    ggplot2::aes (
                                      x = z_score_TPP ,
                                      y = TPP_log_10 ,
                                      label = Accession,
                                      colour = Sig_Check_TPP
                                    )) + ggplot2::geom_point(size = 1.5, alpha = alpha) +
      ggplot2::labs(x= xlab, y= ylab) + ggplot2::scale_colour_manual(breaks = c("Not Significant", "Significant") ,values = labelcolor) +
      ggplot2::expand_limits(x=0, y=0) +
      ggplot2::scale_x_continuous(breaks = -12:12) +
      ggplot2::theme_bw()+ ggplot2::theme(text = ggplot2::element_text(size = plottextsize)) + ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
      ggplot2::theme(legend.title = ggplot2::element_blank())

    print(Volcano_Plot)
  }

  ggplot2::ggsave("Volcano_Plot_TPP_Ligand.png", plot = Volcano_Plot)

  if (interactiveplot == TRUE){
    plotly::ggplotly(Volcano_Plot)
  }


}


#' OnePotSPROX Phenotype Analysis
#'
#' This function will take expression data and OnePotSPROX (TMT10plex) data directly from ProteomeDiscoverer (v. 2.3) and generate a volcano plot and hit list
#' based on user-defined SD and p-values
#'
#' @param Expression_data Raw ProteomeDiscoverer (v. 2.3) data
#' @param TPP_data Raw OnePotSPROX (TMT10plex) ProteomeDiscoverer (v. 2.3) data
#' @param SD_cutoff A positive integer (default = 2)
#' @param p_cutoff A positive integer between 0 and 1 (default = 0.05)
#' @param TMTplex 8, 10, 16 TMTplex (default = 10)
#' @param plot Display volcano plot (default = TRUE)
#' @param xlab X-axis Label (default = "Z Score")
#' @param ylab Y-axis Label (default = "- log10 (p value)")
#' @param labelcolor Color of Not significant and Significant data points (default = c("grey", "red"))
#' @param alpha Transparency of plot (default = 0.5)
#' @param plottextsize Text Size of plot (default = 12)
#' @param correctedpvalue Bonferroni corrected p-value (default = FALSE)
#' @param interactiveplot Plotly interactive Volcano plot (default = FALSE)
#' @return Volcano plot, Hit list, Total proteins assayed and Unique protein hits
#' @export

OnePotSPROX_Phenotype.Fun <- function(Expression_data, SPROX_Raw, SD_cutoff = 2, p_cutoff = 0.05, TMTplex = 10, plot = TRUE, xlab = "Z Score", ylab = "- log10 (p value)", labelcolor = c("grey", "red"), alpha = 0.5, plottextsize = 12, correctedpvalue = FALSE, interactiveplot = FALSE){

  if (TMTplex !=10)
    stop("This TMTplex is not yet supported in this version.\n")

  if (SD_cutoff < 0)
    stop("SD_cutoff must be non-negative.\n")

  if (p_cutoff < 0)
    stop("p_value cutoff must be non-negative.\n")

  High_Confidence_Exp_Level <-
    dplyr::filter(Expression_data, Expression_data[,grepl("Confidence", names(Expression_data))] == "High")

  # This dataframe takes all useful columns from the raw data and will be outputed at the end for reference

  Exp_Level_Concise <- magrittr::"%>%" (High_Confidence_Exp_Level,
                                        dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Description", "Accession", (grep(("^Abundances\\s\\(Grouped):"), names(High_Confidence_Exp_Level), value = TRUE))))))

  # Remove all blank rows in Dataframe

  NA_Removed <- tidyr::drop_na(Exp_Level_Concise)

  # Step 1

  Exp_Level_1 <- magrittr::"%>%" (NA_Removed,
                                  dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Accession", "Description"))))
  Exp_Level_2 <- magrittr::"%>%" (NA_Removed,
                                  dplyr::select(tidyselect::any_of(grep("^Abundances\\s\\(Grouped):", names(NA_Removed), value = TRUE))))
  step1 <- cbind(Exp_Level_1, Exp_Level_2)

  if ("Master Protein Accessions" %in% names(step1)){
    step1 <- dplyr::rename(step1, "Accession" = "Master Protein Accessions")
  }

  # Step 2

  if (("Accession" %in% names(step1) && "Description" %in% names(step1))) {

    step1_real <-
      dplyr::filter(
        step1,
        step1[,3] > 0,
        step1[,4] > 0,
        step1[,5] > 0,
        step1[,6] > 0,
        step1[,7] > 0,
        step1[,8] > 0,
        step1[,9] > 0,
        step1[,10] > 0,
        step1[,11] > 0,
        step1[,12] > 0)

    Normalized_Exp_1 <- (as.numeric(unlist(step1_real[,3])) * mean(as.numeric(unlist(step1_real[,8])))) / ((as.numeric(unlist(step1_real[,8]))) * mean(as.numeric(unlist(step1_real[,3]))))
    Normalized_Exp_2 <- (as.numeric(unlist(step1_real[,4])) * mean(as.numeric(unlist(step1_real[,9])))) / ((as.numeric(unlist(step1_real[,9]))) * mean(as.numeric(unlist(step1_real[,4]))))
    Normalized_Exp_3 <- (as.numeric(unlist(step1_real[,5])) * mean(as.numeric(unlist(step1_real[,10])))) / ((as.numeric(unlist(step1_real[,10]))) * mean(as.numeric(unlist(step1_real[,5]))))
    Normalized_Exp_4 <- (as.numeric(unlist(step1_real[,6])) * mean(as.numeric(unlist(step1_real[,11])))) / ((as.numeric(unlist(step1_real[,11]))) * mean(as.numeric(unlist(step1_real[,6]))))
    Normalized_Exp_5 <- (as.numeric(unlist(step1_real[,7])) * mean(as.numeric(unlist(step1_real[,12])))) / ((as.numeric(unlist(step1_real[,12]))) * mean(as.numeric(unlist(step1_real[,7]))))

  }

  if (!"Description" %in% names(step1)){

    step1_real <-
      dplyr::filter(
        step1,
        step1[,2] > 0,
        step1[,3] > 0,
        step1[,4] > 0,
        step1[,5] > 0,
        step1[,6] > 0,
        step1[,7] > 0,
        step1[,8] > 0,
        step1[,9] > 0,
        step1[,10] > 0,
        step1[,11] > 0)

    Normalized_Exp_1 <- (as.numeric(unlist(step1_real[,2])) * mean(as.numeric(unlist(step1_real[,7])))) / ((as.numeric(unlist(step1_real[,7]))) * mean(as.numeric(unlist(step1_real[,2]))))
    Normalized_Exp_2 <- (as.numeric(unlist(step1_real[,3])) * mean(as.numeric(unlist(step1_real[,8])))) / ((as.numeric(unlist(step1_real[,8]))) * mean(as.numeric(unlist(step1_real[,3]))))
    Normalized_Exp_3 <- (as.numeric(unlist(step1_real[,4])) * mean(as.numeric(unlist(step1_real[,9])))) / ((as.numeric(unlist(step1_real[,9]))) * mean(as.numeric(unlist(step1_real[,4]))))
    Normalized_Exp_4 <- (as.numeric(unlist(step1_real[,5])) * mean(as.numeric(unlist(step1_real[,10])))) / ((as.numeric(unlist(step1_real[,10]))) * mean(as.numeric(unlist(step1_real[,5]))))
    Normalized_Exp_5 <- (as.numeric(unlist(step1_real[,6])) * mean(as.numeric(unlist(step1_real[,11])))) / ((as.numeric(unlist(step1_real[,11]))) * mean(as.numeric(unlist(step1_real[,6]))))

  }


  step2 <-
    cbind(step1_real,
          Normalized_Exp_1,
          Normalized_Exp_2,
          Normalized_Exp_3,
          Normalized_Exp_4,
          Normalized_Exp_5)

  # Step 3

  Log_2_Exp_1 <- log(Normalized_Exp_1, 2)
  Log_2_Exp_2 <- log(Normalized_Exp_2, 2)
  Log_2_Exp_3 <- log(Normalized_Exp_3, 2)
  Log_2_Exp_4 <- log(Normalized_Exp_4, 2)
  Log_2_Exp_5 <- log(Normalized_Exp_5, 2)

  step3 <- cbind(step2, Log_2_Exp_1, Log_2_Exp_2, Log_2_Exp_3, Log_2_Exp_4, Log_2_Exp_5)

  # Step 4

  Exp_Log_Avg <- rowMeans(step3[, c("Log_2_Exp_1", "Log_2_Exp_2", "Log_2_Exp_3", "Log_2_Exp_4", "Log_2_Exp_5")], na.rm = T)
  step4 <- cbind(step3, Exp_Log_Avg)

  # Step 5

  c1 <- step4[, c("Log_2_Exp_1", "Log_2_Exp_2", "Log_2_Exp_3", "Log_2_Exp_4", "Log_2_Exp_5")]
  SD_Exp <- apply(c1, 1, sd)
  step5 <- cbind(step4, SD_Exp)

  # Step 6

  f1 <- step5$Exp_Log_Avg
  mean_log_avg_Exp <- mean(f1)
  sd_log_avg_Exp <- sd(f1)

  z_score_Exp <- ((f1 - mean_log_avg_Exp)) / (sd_log_avg_Exp)
  step6 <- cbind(step5, z_score_Exp)

  # Step 7

  g1 <- abs(Exp_Log_Avg)
  g2 <- step6$SD_Exp

  T_value_Exp <- (g1 * sqrt(5) / g2)
  step7 <- cbind(step6, T_value_Exp)

  # Step 8

  h1 <- step7$Exp_Log_Avg
  h2 <- mean(h1)
  h3 <- step7$T_value_Exp
  p_value_Exp <- 2 * pt(h3, 4, lower.tail = F)
  step8 <- cbind(step7, p_value_Exp)

  # Step 9

  j1 <- step8$p_value_Exp
  log_10_Exp <- -log10(j1)
  step9 <- cbind(step8, log_10_Exp)

  # Step 10

  Hit_ident1_Exp <- magrittr::"%>%" (step9,
                                     dplyr::select(tidyselect::any_of(c("Accession", "Description"))))
  Log_Avg_Exp <- step9$Exp_Log_Avg
  Log_10_Exp <- step9$log_10_Exp
  Hit_ident_Exp <- cbind(Hit_ident1_Exp, Log_Avg_Exp, Log_10_Exp, Log_2_Exp_1, Log_2_Exp_2, Log_2_Exp_3, Log_2_Exp_4, Log_2_Exp_5, z_score_Exp)

  n1 <-
    (Hit_ident_Exp$z_score_Exp > SD_cutoff &
       Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff)) |
    (Hit_ident_Exp$z_score_Exp < -SD_cutoff & Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff))
  n2 <- as.data.frame(n1)

  Sig_Check_Exp <-
    dplyr::if_else(n2$n1 == "TRUE", "Significant", "Not Significant")
  step10 <- cbind(Hit_ident_Exp, p_value_Exp, Sig_Check_Exp)

  # Remove [K] and [R] [+], [-], All bracketed Amino Acids
  # If contains oxidation modification remove it

  SPROX_Confidence <- dplyr::filter(SPROX_Raw, SPROX_Raw[,grepl("Confidence", names(SPROX_Raw))] == "High")

  SPROX_No_Oxid <- dplyr::filter(SPROX_Confidence, !grepl('Oxidation', Modifications))

  SPROX_removed_AA1 <- dplyr::mutate_all(SPROX_No_Oxid, ~gsub("\\[.*?].", "", .))
  SPROX_removed_AA2 <- dplyr::mutate_all(SPROX_removed_AA1, ~gsub("\\.\\[.*?]", "", .))


  # If contains M filter into new column

  SPROX_Just_M <- SPROX_removed_AA2[grep("M", SPROX_removed_AA2$`Annotated Sequence`), ]


  # Go through TPP analysis

  SPROX_1 <- magrittr::"%>%" (SPROX_Just_M,
                              dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Positions in Master Proteins"))))
  SPROX_2 <- magrittr::"%>%" (SPROX_Just_M,
                              dplyr::select(tidyselect::any_of(c("Annotated Sequence", "Modifications"))))
  SPROX_3 <- magrittr::"%>%" (SPROX_Just_M,
                              dplyr::select(tidyselect::any_of(grep("^Abundances\\s\\(Grouped):", names(SPROX_Just_M), value = TRUE))))
  step1_SPROX <- cbind(SPROX_1, SPROX_2, SPROX_3)

  if ("Master Protein Accessions" %in% names(step1_SPROX)){
    step1_SPROX <- dplyr::rename(step1_SPROX, "Accession" = "Master Protein Accessions")
  }

  NA_Removed_SPROX <- tidyr::drop_na(step1_SPROX)

  SPROX_merger <- merge(Hit_ident_Exp, NA_Removed_SPROX)


  if (("Accession" %in% names(NA_Removed_SPROX) && "Positions in Master Proteins" %in% names(NA_Removed_SPROX))) {

    step1_real <-
      dplyr::filter(
        SPROX_merger,
        SPROX_merger[,14] > 0,
        SPROX_merger[,15] > 0,
        SPROX_merger[,16] > 0,
        SPROX_merger[,17] > 0,
        SPROX_merger[,18] > 0,
        SPROX_merger[,19] > 0,
        SPROX_merger[,20] > 0,
        SPROX_merger[,21] > 0,
        SPROX_merger[,22] > 0,
        SPROX_merger[,23] > 0)

    Normalized_SPROX_T1 <- (as.numeric(unlist(step1_real[,14])) * mean(as.numeric(unlist(step1_real[,19])))) / ((as.numeric(unlist(step1_real[,19]))) * mean(as.numeric(unlist(step1_real[,14]))))
    Normalized_SPROX_T2 <- (as.numeric(unlist(step1_real[,15])) * mean(as.numeric(unlist(step1_real[,20])))) / ((as.numeric(unlist(step1_real[,20]))) * mean(as.numeric(unlist(step1_real[,15]))))
    Normalized_SPROX_T3 <- (as.numeric(unlist(step1_real[,16])) * mean(as.numeric(unlist(step1_real[,21])))) / ((as.numeric(unlist(step1_real[,21]))) * mean(as.numeric(unlist(step1_real[,16]))))
    Normalized_SPROX_T4 <- (as.numeric(unlist(step1_real[,17])) * mean(as.numeric(unlist(step1_real[,22])))) / ((as.numeric(unlist(step1_real[,22]))) * mean(as.numeric(unlist(step1_real[,17]))))
    Normalized_SPROX_T5 <- (as.numeric(unlist(step1_real[,18])) * mean(as.numeric(unlist(step1_real[,23])))) / ((as.numeric(unlist(step1_real[,23]))) * mean(as.numeric(unlist(step1_real[,18]))))

  }

  if (!"Positions in Master Proteins" %in% names(NA_Removed_SPROX)){

    step1_real <-
      dplyr::filter(
        SPROX_merger,
        SPROX_merger[,13] > 0,
        SPROX_merger[,14] > 0,
        SPROX_merger[,15] > 0,
        SPROX_merger[,16] > 0,
        SPROX_merger[,17] > 0,
        SPROX_merger[,18] > 0,
        SPROX_merger[,19] > 0,
        SPROX_merger[,20] > 0,
        SPROX_merger[,21] > 0,
        SPROX_merger[,22] > 0)

    Normalized_SPROX_T1 <- (as.numeric(unlist(step1_real[,13])) * mean(as.numeric(unlist(step1_real[,18])))) / ((as.numeric(unlist(step1_real[,18]))) * mean(as.numeric(unlist(step1_real[,13]))))
    Normalized_SPROX_T2 <- (as.numeric(unlist(step1_real[,14])) * mean(as.numeric(unlist(step1_real[,19])))) / ((as.numeric(unlist(step1_real[,19]))) * mean(as.numeric(unlist(step1_real[,14]))))
    Normalized_SPROX_T3 <- (as.numeric(unlist(step1_real[,15])) * mean(as.numeric(unlist(step1_real[,20])))) / ((as.numeric(unlist(step1_real[,20]))) * mean(as.numeric(unlist(step1_real[,15]))))
    Normalized_SPROX_T4 <- (as.numeric(unlist(step1_real[,16])) * mean(as.numeric(unlist(step1_real[,21])))) / ((as.numeric(unlist(step1_real[,21]))) * mean(as.numeric(unlist(step1_real[,16]))))
    Normalized_SPROX_T5 <- (as.numeric(unlist(step1_real[,17])) * mean(as.numeric(unlist(step1_real[,22])))) / ((as.numeric(unlist(step1_real[,22]))) * mean(as.numeric(unlist(step1_real[,17]))))

  }

  SPROX_merger <- merge(Hit_ident_Exp, step1_real)

  SPROX_merge <- SPROX_merger

  # Step 11

  step11 <-
    cbind(SPROX_merge,
          Normalized_SPROX_T1,
          Normalized_SPROX_T2,
          Normalized_SPROX_T3,
          Normalized_SPROX_T4,
          Normalized_SPROX_T5)

  step111 <- tidyr::drop_na(step11)

  # Step 12

  Log_2_SPROX_T1 <- log(Normalized_SPROX_T1, 2)
  Log_2_SPROX_T2 <- log(Normalized_SPROX_T2, 2)
  Log_2_SPROX_T3 <- log(Normalized_SPROX_T3, 2)
  Log_2_SPROX_T4 <- log(Normalized_SPROX_T4, 2)
  Log_2_SPROX_T5 <- log(Normalized_SPROX_T5, 2)

  step12 <- cbind(step111, Log_2_SPROX_T1, Log_2_SPROX_T2, Log_2_SPROX_T3, Log_2_SPROX_T4, Log_2_SPROX_T5)

  # Step 13

  Log_Avg_SPROX <- rowMeans(step12[, c("Log_2_SPROX_T1", "Log_2_SPROX_T2", "Log_2_SPROX_T3", "Log_2_SPROX_T4", "Log_2_SPROX_T5")])
  step13 <- cbind(step12, Log_Avg_SPROX)

  # Step  14

  rep1_exp1 <- step13$Log_2_SPROX_T1 - step13$Log_2_Exp_1
  rep2_exp2 <- step13$Log_2_SPROX_T2 - step13$Log_2_Exp_2
  rep3_exp3 <- step13$Log_2_SPROX_T3 - step13$Log_2_Exp_3
  rep4_exp4 <- step13$Log_2_SPROX_T4 - step13$Log_2_Exp_4
  rep5_exp5 <- step13$Log_2_SPROX_T5 - step13$Log_2_Exp_5
  step14 <- cbind(step13, rep1_exp1, rep2_exp2, rep3_exp3, rep4_exp4, rep5_exp5)

  # Step 15

  SPROX_Avg <- rowMeans(step14[, c("rep1_exp1", "rep2_exp2", "rep3_exp3", "rep4_exp4", "rep5_exp5")], na.rm = T)
  step15 <- cbind(step14, SPROX_Avg)

  # Step 16

  k1 <- step15[, c("rep1_exp1", "rep2_exp2", "rep3_exp3", "rep4_exp4", "rep5_exp5")]
  SD_SPROX <- apply(k1, 1, sd)
  step16 <- cbind(step15, SD_SPROX)

  # Step 17

  step17 <- tidyr::drop_na(step16)
  SPROXunique <- dplyr::n_distinct(step17$Accession)
  SPROXuniquePep <- nrow(step17)
  # Step 18

  y1 <- step17$SPROX_Avg
  SPROX_mean_avg <- mean(y1)
  SPROX_sd_avg <- sd(y1)
  z_score_SPROX <- ((y1 - SPROX_mean_avg)) / (SPROX_sd_avg)
  step18 <- cbind(step17, z_score_SPROX)

  # Step 19

  hg1 <- abs(SPROX_Avg)
  hg2 <- step18$SD_SPROX

  SPROX_T_value <- (hg1 * sqrt(5) / hg2)
  step19 <- cbind(step18, SPROX_T_value)

  # Step 20

  hk1 <- step19$SPROX_Avg
  hk2 <- mean(hg1)
  hk3 <- step19$SPROX_T_value
  p_value_SPROX <- 2 * pt(hk3, 4, lower.tail = F)
  step20 <- cbind(step19, p_value_SPROX)

  # Step 21

  jk1 <- step20$p_value_SPROX
  SPROX_log_10 <- -log10(jk1)
  step21 <- cbind(step20, SPROX_log_10)

  SPROX_Hit_ident1 <- magrittr::"%>%" (step21,
                                       dplyr::select(tidyselect::any_of(c("Accession", "Description"))))
  SPROX_Hit_ident2 <- magrittr::"%>%" (step21,
                                       dplyr::select(tidyselect::any_of(c("Annotated Sequence", "Modifications"))))
  SPROX_Log_Avg <- step21[, c("SPROX_Avg", "SD_SPROX", "z_score_SPROX", "SPROX_T_value", "p_value_SPROX", "SPROX_log_10")]
  SPROX_Hit_ident <- cbind(SPROX_Hit_ident1, SPROX_Hit_ident2, SPROX_Log_Avg)

  # Step 22

  if (correctedpvalue == FALSE){

  nk1 <-
    (SPROX_Hit_ident$z_score_SPROX > SD_cutoff &
       SPROX_Hit_ident$SPROX_log_10 > -log10(p_cutoff)) |
    (SPROX_Hit_ident$z_score_SPROX < -SD_cutoff & SPROX_Hit_ident$SPROX_log_10 > -log10(p_cutoff))
  nk2 <- as.data.frame(nk1)

  }

  if (correctedpvalue == TRUE){
    nk1 <-
      (SPROX_Hit_ident$z_score_SPROX > SD_cutoff &
         SPROX_Hit_ident$SPROX_log_10 > -log10((p_cutoff/SPROXuniquePep))) |
      (SPROX_Hit_ident$z_score_SPROX < -SD_cutoff & SPROX_Hit_ident$SPROX_log_10 > -log10((p_cutoff/SPROXuniquePep)))
    nk2 <- as.data.frame(nk1)

  }

  Sig_Check_SPROX <-
    dplyr::if_else(nk2$nk1 == "TRUE", "Significant", "Not Significant")
  step22 <- cbind(SPROX_Hit_ident, Sig_Check_SPROX)

  ddk <- step22$Sig_Check_SPROX

  as.data.frame(ddk)
  abk <- (ddk == "Significant")
  SPROX_Number_Of_Hits <- length(which(abk))

  Hit_List <- step22[ step22$Sig_Check_SPROX == "Significant", ]
  SPROXuniqueHits <- dplyr::n_distinct(Hit_List$Accession)
  SPROXExport <- unique(Hit_List$Accession)
  SPROXSigExport <- magrittr::"%>%" (step22,
                                     dplyr::select(tidyselect::any_of(c("Accession", "Description", "Annotated Sequence", "Modifications", "z_score_SPROX", "p_value_SPROX", "Sig_Check_SPROX"))))
  Hit_ListExportSPROX <- magrittr::"%>%" (Hit_List,
                                          dplyr::select(tidyselect::any_of(c("Accession", "Description", "Annotated Sequence", "Modifications", "z_score_SPROX", "p_value_SPROX", "Sig_Check_SPROX"))))


  print(SPROXSigExport)
  print(paste("There are", SPROX_Number_Of_Hits, "Hits"))
  print(paste("There are", SPROXuniqueHits, "Unique Hits"))
  print(paste(SPROXuniquePep, "Unique Peptides were assayed"))
  print(paste(SPROXunique, "Unique Proteins were assayed"))
  print(Hit_ListExportSPROX)

  utils::write.table(SPROXExport, file = "OnePotSPROXPhenotypeUniqueHits.csv", row.names = FALSE, col.names = FALSE)
  utils::write.table(SPROXSigExport, file = "OnePotSPROXPhenotypeOut.csv", row.names = FALSE)
  utils::write.table(Hit_ListExportSPROX, file = "OnePotSPROXPhenotypeHitsOut.csv", row.names = FALSE)

  # Visualizing Data

  if (plot == TRUE){

  Volcano_Plot <- ggplot2::ggplot(step22,
                         ggplot2::aes (
                           x = z_score_SPROX ,
                           y = SPROX_log_10 ,
                           label = Accession,
                           colour = Sig_Check_SPROX
                         )) + ggplot2::geom_point(size = 1.5, alpha = alpha) +
    ggplot2::labs(x= xlab, y= ylab) + ggplot2::scale_colour_manual(breaks = c("Not Significant", "Significant"), values = labelcolor) +
    ggplot2::expand_limits(x=0, y=0) +
    ggplot2::scale_x_continuous(breaks = -12:12) +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = plottextsize)) + ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  print(Volcano_Plot)
  }

ggplot2::ggsave("Volcano_Plot_SPROX_Phenotype.png", plot = Volcano_Plot)

if (interactiveplot == TRUE){
  plotly::ggplotly(Volcano_Plot)
}

}


#' OnePotSPROX Ligand Analysis
#'
#' This function will take OnePotSPROX (TMT10plex) data directly from ProteomeDiscoverer (v. 2.3) and generate a volcano plot and hit list
#' based on user-defined SD and p-values
#'
#' @param Expression_data Raw ProteomeDiscoverer (v. 2.3) data
#' @param TPP_data Raw OnePotSPROX (TMT10plex) ProteomeDiscoverer (v. 2.3) data
#' @param SD_cutoff A positive integer (default = 2)
#' @param p_cutoff A positive integer between 0 and 1 (default = 0.05)
#' @param TMTplex 8, 10, 16 TMTplex (default = 10)
#' @param plot Display volcano plot (default = TRUE)
#' @param xlab X-axis Label (default = "Z Score")
#' @param ylab Y-axis Label (default = "- log10 (p value)")
#' @param labelcolor Color of Not significant and Significant data points (default = c("grey", "red"))
#' @param alpha Transparency of plot (default = 0.5)
#' @param plottextsize Text Size of plot (default = 12)
#' @param correctedpvalue Bonferroni corrected p-value (default = FALSE)
#' @param interactiveplot Plotly interactive Volcano plot (default = FALSE)
#' @return Volcano plot, Hit list, Total proteins assayed and Unique protein hits
#' @export

OnePotSPROX_Ligand.Fun <- function(SPROX_Raw, SD_cutoff = 2, p_cutoff = 0.05, TMTplex = 10, plot = TRUE, xlab = "Z Score", ylab = "- log10 (p value)", labelcolor = c("grey", "red"), alpha = 0.5, plottextsize = 12, correctedpvalue = FALSE, interactiveplot = FALSE){

  if (TMTplex !=10)
    stop("This TMTplex is not yet supported in this version.\n")

  if (SD_cutoff < 0)
    stop("SD_cutoff must be non-negative.\n")

  if (p_cutoff < 0)
    stop("p_value cutoff must be non-negative.\n")

  # Remove [K] and [R] [+], [-], All bracketed Amino Acids
  # If contains oxidation modification remove it

  SPROX_Confidence <- dplyr::filter(SPROX_Raw, SPROX_Raw[,grepl("Confidence", names(SPROX_Raw))] == "High")

  SPROX_No_Oxid <- dplyr::filter(SPROX_Confidence, !grepl('Oxidation', Modifications))

  SPROX_removed_AA1 <- dplyr::mutate_all(SPROX_No_Oxid, ~gsub("\\[.*?].", "", .))
  SPROX_removed_AA2 <- dplyr::mutate_all(SPROX_removed_AA1, ~gsub("\\.\\[.*?]", "", .))


  # If contains M filter into new column

  SPROX_Just_M <- SPROX_removed_AA2[grep("M", SPROX_removed_AA2$`Annotated Sequence`), ]


  # Go through TPP analysis

  SPROX_1 <- magrittr::"%>%" (SPROX_Just_M,
                              dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Positions in Master Proteins"))))
  SPROX_2 <- magrittr::"%>%" (SPROX_Just_M,
                              dplyr::select(tidyselect::any_of(c("Annotated Sequence", "Modifications"))))
  SPROX_3 <- magrittr::"%>%" (SPROX_Just_M,
                              dplyr::select(tidyselect::any_of(grep("^Abundances\\s\\(Grouped):", names(SPROX_Just_M), value = TRUE))))
  step1_SPROX <- cbind(SPROX_1, SPROX_2, SPROX_3)

  if ("Master Protein Accessions" %in% names(step1_SPROX)){
    step1_SPROX <- dplyr::rename(step1_SPROX, "Accession" = "Master Protein Accessions")
  }

  # Filter blanks

  NA_Removed_SPROX <- tidyr::drop_na(step1_SPROX)

  if (("Accession" %in% names(NA_Removed_SPROX) && "Positions in Master Proteins" %in% names(NA_Removed_SPROX))) {

    step1_real <-
      dplyr::filter(
        NA_Removed_SPROX,
        NA_Removed_SPROX[,5] > 0,
        NA_Removed_SPROX[,6] > 0,
        NA_Removed_SPROX[,7] > 0,
        NA_Removed_SPROX[,8] > 0,
        NA_Removed_SPROX[,9] > 0,
        NA_Removed_SPROX[,10] > 0,
        NA_Removed_SPROX[,11] > 0,
        NA_Removed_SPROX[,12] > 0,
        NA_Removed_SPROX[,13] > 0,
        NA_Removed_SPROX[,14] > 0)

    Normalized_SPROX_T1 <- (as.numeric(unlist(step1_real[,5])) * mean(as.numeric(unlist(step1_real[,10])))) / ((as.numeric(unlist(step1_real[,10]))) * mean(as.numeric(unlist(step1_real[,5]))))
    Normalized_SPROX_T2 <- (as.numeric(unlist(step1_real[,6])) * mean(as.numeric(unlist(step1_real[,11])))) / ((as.numeric(unlist(step1_real[,11]))) * mean(as.numeric(unlist(step1_real[,6]))))
    Normalized_SPROX_T3 <- (as.numeric(unlist(step1_real[,7])) * mean(as.numeric(unlist(step1_real[,12])))) / ((as.numeric(unlist(step1_real[,12]))) * mean(as.numeric(unlist(step1_real[,7]))))
    Normalized_SPROX_T4 <- (as.numeric(unlist(step1_real[,8])) * mean(as.numeric(unlist(step1_real[,13])))) / ((as.numeric(unlist(step1_real[,13]))) * mean(as.numeric(unlist(step1_real[,8]))))
    Normalized_SPROX_T5 <- (as.numeric(unlist(step1_real[,9])) * mean(as.numeric(unlist(step1_real[,14])))) / ((as.numeric(unlist(step1_real[,14]))) * mean(as.numeric(unlist(step1_real[,9]))))

  }

  if (!"Positions in Master Proteins" %in% names(NA_Removed_SPROX)){

    step1_real <-
      dplyr::filter(
        NA_Removed_SPROX,
        NA_Removed_SPROX[,4] > 0,
        NA_Removed_SPROX[,5] > 0,
        NA_Removed_SPROX[,6] > 0,
        NA_Removed_SPROX[,7] > 0,
        NA_Removed_SPROX[,8] > 0,
        NA_Removed_SPROX[,9] > 0,
        NA_Removed_SPROX[,10] > 0,
        NA_Removed_SPROX[,11] > 0,
        NA_Removed_SPROX[,12] > 0,
        NA_Removed_SPROX[,13] > 0)

    Normalized_SPROX_T1 <- (as.numeric(unlist(step1_real[,4])) * mean(as.numeric(unlist(step1_real[,9])))) / ((as.numeric(unlist(step1_real[,9]))) * mean(as.numeric(unlist(step1_real[,4]))))
    Normalized_SPROX_T2 <- (as.numeric(unlist(step1_real[,5])) * mean(as.numeric(unlist(step1_real[,10])))) / ((as.numeric(unlist(step1_real[,10]))) * mean(as.numeric(unlist(step1_real[,5]))))
    Normalized_SPROX_T3 <- (as.numeric(unlist(step1_real[,6])) * mean(as.numeric(unlist(step1_real[,11])))) / ((as.numeric(unlist(step1_real[,11]))) * mean(as.numeric(unlist(step1_real[,6]))))
    Normalized_SPROX_T4 <- (as.numeric(unlist(step1_real[,7])) * mean(as.numeric(unlist(step1_real[,12])))) / ((as.numeric(unlist(step1_real[,12]))) * mean(as.numeric(unlist(step1_real[,7]))))
    Normalized_SPROX_T5 <- (as.numeric(unlist(step1_real[,8])) * mean(as.numeric(unlist(step1_real[,13])))) / ((as.numeric(unlist(step1_real[,13]))) * mean(as.numeric(unlist(step1_real[,8]))))

  }



  # Step 11

  step11 <-
    cbind(step1_real,
          Normalized_SPROX_T1,
          Normalized_SPROX_T2,
          Normalized_SPROX_T3,
          Normalized_SPROX_T4,
          Normalized_SPROX_T5)

  step111 <- tidyr::drop_na(step11)

  # Step 12

  Log_2_SPROX_T1 <- log(Normalized_SPROX_T1, 2)
  Log_2_SPROX_T2 <- log(Normalized_SPROX_T2, 2)
  Log_2_SPROX_T3 <- log(Normalized_SPROX_T3, 2)
  Log_2_SPROX_T4 <- log(Normalized_SPROX_T4, 2)
  Log_2_SPROX_T5 <- log(Normalized_SPROX_T5, 2)

  step12 <- cbind(step111, Log_2_SPROX_T1, Log_2_SPROX_T2, Log_2_SPROX_T3, Log_2_SPROX_T4, Log_2_SPROX_T5)

  # Step 13

  Log_Avg_SPROX <- rowMeans(step12[, c("Log_2_SPROX_T1", "Log_2_SPROX_T2", "Log_2_SPROX_T3", "Log_2_SPROX_T4", "Log_2_SPROX_T5")])
  step13 <- cbind(step12, Log_Avg_SPROX)

  # Step  14

  step14 <- step13

  # Step 15

  step15 <- step14

  # Step 16

  k1 <- step15[, c("Log_2_SPROX_T1", "Log_2_SPROX_T2", "Log_2_SPROX_T3", "Log_2_SPROX_T4", "Log_2_SPROX_T5")]
  SD_SPROX <- apply(k1, 1, sd)
  step16 <- cbind(step15, SD_SPROX)

  # Step 17

  step17 <- tidyr::drop_na(step16)
  SPROXunique <- dplyr::n_distinct(step17$Accession)
  SPROXuniquePep <- nrow(step17)
  # Step 18

  y1 <- step17$Log_Avg_SPROX
  SPROX_mean_avg <- mean(y1)
  SPROX_sd_avg <- sd(y1)
  z_score_SPROX <- ((y1 - SPROX_mean_avg)) / (SPROX_sd_avg)
  step18 <- cbind(step17, z_score_SPROX)

  # Step 19

  hg1 <- abs(Log_Avg_SPROX)
  hg2 <- step18$SD_SPROX

  SPROX_T_value <- (hg1 * sqrt(5) / hg2)
  step19 <- cbind(step18, SPROX_T_value)

  # Step 20

  hk1 <- step19$Log_Avg_SPROX
  hk2 <- mean(hg1)
  hk3 <- step19$SPROX_T_value
  p_value_SPROX <- 2 * pt(hk3, 4, lower.tail = F)
  step20 <- cbind(step19, p_value_SPROX)

  # Step 21

  jk1 <- step20$p_value_SPROX
  SPROX_log_10 <- -log10(jk1)
  step21 <- cbind(step20, SPROX_log_10)
  SPROX_Hit_ident1 <- magrittr::"%>%" (step21,
                                       dplyr::select(tidyselect::any_of(c("Accession", "Annotated Sequence", "Modifications"))))
  SPROX_Log_Avg <- magrittr::"%>%" (step21,
                                    dplyr::select(tidyselect::any_of(c("Log_Avg_SPROX", "SD_SPROX", "z_score_SPROX", "SPROX_T_value", "p_value_SPROX", "SPROX_log_10"))))
  SPROX_Hit_ident <- cbind(SPROX_Hit_ident1, SPROX_Log_Avg)

  # Step 22

  if (correctedpvalue == FALSE){

  nk1 <-
    (SPROX_Hit_ident$z_score_SPROX > SD_cutoff &
       SPROX_Hit_ident$SPROX_log_10 > -log10(p_cutoff)) |
    (SPROX_Hit_ident$z_score_SPROX < -SD_cutoff & SPROX_Hit_ident$SPROX_log_10 > -log10(p_cutoff))
  nk2 <- as.data.frame(nk1)

  }

  if (correctedpvalue == TRUE){
    nk1 <-
      (SPROX_Hit_ident$z_score_SPROX > SD_cutoff &
         SPROX_Hit_ident$SPROX_log_10 > -log10((p_cutoff/SPROXuniquePep))) |
  (SPROX_Hit_ident$z_score_SPROX < -SD_cutoff & SPROX_Hit_ident$SPROX_log_10 > -log10((p_cutoff/SPROXuniquePep)))
nk2 <- as.data.frame(nk1)

  }

  Sig_Check_SPROX <-
    dplyr::if_else(nk2$nk1 == "TRUE", "Significant", "Not Significant")
  step22 <- cbind(SPROX_Hit_ident, Sig_Check_SPROX)

  ddk <- step22$Sig_Check_SPROX

  as.data.frame(ddk)
  abk <- (ddk == "Significant")
  SPROX_Number_Of_Hits <- length(which(abk))

  Hit_List <- step22[ step22$Sig_Check_SPROX == "Significant", ]
  SPROXuniqueHits <- dplyr::n_distinct(Hit_List$Accession)
  SPROXExport <- unique(Hit_List$Accession)
  SPROXSigExport <- magrittr::"%>%" (step22,
                                     dplyr::select(tidyselect::any_of(c("Accession", "Annotated Sequence", "Modifications", "z_score_SPROX", "p_value_SPROX", "Sig_Check_SPROX"))))
  Hit_ListExportSPROX <- magrittr::"%>%" (Hit_List,
                                          dplyr::select(tidyselect::any_of(c("Accession", "Annotated Sequence", "Modifications", "z_score_SPROX", "p_value_SPROX", "Sig_Check_SPROX"))))


  print(SPROXSigExport)
  print(paste("There are", SPROX_Number_Of_Hits, "Hits"))
  print(paste("There are", SPROXuniqueHits, "Unique Hits"))
  print(paste(SPROXuniquePep, "Unique Peptides were assayed"))
  print(paste(SPROXunique, "Unique Proteins were assayed"))
  print(Hit_ListExportSPROX)

  utils::write.table(SPROXExport, file = "OnePotSPROXLigandUniqueHits.csv", row.names = FALSE, col.names = FALSE)
  utils::write.table(SPROXSigExport, file = "OnePotSPROXLigandOut.csv", row.names = FALSE)
  utils::write.table(Hit_ListExportSPROX, file = "OnePotSPROXLigandHitsOut.csv", row.names = FALSE)

  # Visualizing Data

  if (plot == TRUE){

    Volcano_Plot <- ggplot2::ggplot(step22,
                                    ggplot2::aes (
                                      x = z_score_SPROX ,
                                      y = SPROX_log_10 ,
                                      label = Accession,
                                      colour = Sig_Check_SPROX
                                    )) + ggplot2::geom_point(size = 1.5, alpha = alpha) +
      ggplot2::labs(x= xlab, y= ylab) + ggplot2::scale_colour_manual(breaks = c("Not Significant", "Significant"), values = labelcolor) +
      ggplot2::expand_limits(x=0, y=0) +
      ggplot2::scale_x_continuous(breaks = -12:12) +
      ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = plottextsize)) + ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
      ggplot2::theme(legend.title = ggplot2::element_blank())

    print(Volcano_Plot)
  }

  ggplot2::ggsave("Volcano_Plot_SPROX_Ligand.png", plot = Volcano_Plot)

  if (interactiveplot == TRUE){
    plotly::ggplotly(Volcano_Plot)
  }

}



#' STEPP Phenotype Analysis
#'
#' This function will take expression data and STEPP/LiP (TMT10plex) data directly from ProteomeDiscoverer (v. 2.3) and generate a volcano plot and hit list
#' based on user-defined SD and p-values
#'
#' @param Expression_data Raw ProteomeDiscoverer (v. 2.3) data
#' @param TPP_data Raw STEPP (TMT10plex) ProteomeDiscoverer (v. 2.3) data
#' @param SD_cutoff A positive integer (default = 2)
#' @param p_cutoff A positive integer between 0 and 1 (default = 0.05)
#' @param TMTplex 8, 10, 16 TMTplex (default = 10)
#' @param plot Display volcano plot (default = TRUE)
#' @param xlab X-axis Label (default = "Z Score")
#' @param ylab Y-axis Label (default = "- log10 (p value)")
#' @param labelcolor Color of Not significant and Significant data points (default = c("grey", "red"))
#' @param alpha Transparency of plot (default = 0.5)
#' @param plottextsize Text Size of plot (default = 12)
#' @param correctedpvalue Bonferroni corrected p-value (default = FALSE)
#' @param interactiveplot Plotly interactive Volcano plot (default = FALSE)
#' @return Volcano plot, Hit list, Total proteins assayed and Unique protein hits
#' @export

STEPP_Phenotype.Fun <- function(Expression_data, STEPP_Raw, SD_cutoff = 2, p_cutoff = 0.05, TMTplex = 10, plot = TRUE, xlab = "Z Score", ylab = "- log10 (p value)", labelcolor = c("grey", "red"), alpha = 0.5, plottextsize = 12, correctedpvalue = FALSE, interactiveplot = FALSE){

  if (TMTplex !=10)
    stop("This TMTplex is not yet supported in this version.\n")

  if (SD_cutoff < 0)
    stop("SD_cutoff must be non-negative.\n")

  if (p_cutoff < 0)
    stop("p_value cutoff must be non-negative.\n")

  High_Confidence_Exp_Level <-
    dplyr::filter(Expression_data, Expression_data[,grepl("Confidence", names(Expression_data))] == "High")

  # This dataframe takes all useful columns from the raw data and will be outputed at the end for reference

  Exp_Level_Concise <- magrittr::"%>%" (High_Confidence_Exp_Level,
                                        dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Description", "Accession", (grep(("^Abundances\\s\\(Grouped):"), names(High_Confidence_Exp_Level), value = TRUE))))))


  # Remove all blank rows in Dataframe

  NA_Removed <- tidyr::drop_na(Exp_Level_Concise)

  # Step 1

  Exp_Level_1 <- magrittr::"%>%" (NA_Removed,
                                  dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Accession", "Description"))))
  Exp_Level_2 <- magrittr::"%>%" (NA_Removed,
                                  dplyr::select(tidyselect::any_of(grep("^Abundances\\s\\(Grouped):", names(High_Confidence_Exp_Level), value = TRUE))))
  step1 <- cbind(Exp_Level_1, Exp_Level_2)

  if ("Master Protein Accessions" %in% names(step1)){
    step1 <- dplyr::rename(step1, "Accession" = "Master Protein Accessions")
  }

  # Step 2

  if (("Accession" %in% names(step1) && "Description" %in% names(step1))) {

    step1_real <-
      dplyr::filter(
        step1,
        step1[,3] > 0,
        step1[,4] > 0,
        step1[,5] > 0,
        step1[,6] > 0,
        step1[,7] > 0,
        step1[,8] > 0,
        step1[,9] > 0,
        step1[,10] > 0,
        step1[,11] > 0,
        step1[,12] > 0)

    Normalized_Exp_1 <- (as.numeric(unlist(step1_real[,3])) * mean(as.numeric(unlist(step1_real[,8])))) / ((as.numeric(unlist(step1_real[,8]))) * mean(as.numeric(unlist(step1_real[,3]))))
    Normalized_Exp_2 <- (as.numeric(unlist(step1_real[,4])) * mean(as.numeric(unlist(step1_real[,9])))) / ((as.numeric(unlist(step1_real[,9]))) * mean(as.numeric(unlist(step1_real[,4]))))
    Normalized_Exp_3 <- (as.numeric(unlist(step1_real[,5])) * mean(as.numeric(unlist(step1_real[,10])))) / ((as.numeric(unlist(step1_real[,10]))) * mean(as.numeric(unlist(step1_real[,5]))))
    Normalized_Exp_4 <- (as.numeric(unlist(step1_real[,6])) * mean(as.numeric(unlist(step1_real[,11])))) / ((as.numeric(unlist(step1_real[,11]))) * mean(as.numeric(unlist(step1_real[,6]))))
    Normalized_Exp_5 <- (as.numeric(unlist(step1_real[,7])) * mean(as.numeric(unlist(step1_real[,12])))) / ((as.numeric(unlist(step1_real[,12]))) * mean(as.numeric(unlist(step1_real[,7]))))

  }

  if (!"Description" %in% names(step1)){

    step1_real <-
      dplyr::filter(
        step1,
        step1[,2] > 0,
        step1[,3] > 0,
        step1[,4] > 0,
        step1[,5] > 0,
        step1[,6] > 0,
        step1[,7] > 0,
        step1[,8] > 0,
        step1[,9] > 0,
        step1[,10] > 0,
        step1[,11] > 0)

    Normalized_Exp_1 <- (as.numeric(unlist(step1_real[,2])) * mean(as.numeric(unlist(step1_real[,7])))) / ((as.numeric(unlist(step1_real[,7]))) * mean(as.numeric(unlist(step1_real[,2]))))
    Normalized_Exp_2 <- (as.numeric(unlist(step1_real[,3])) * mean(as.numeric(unlist(step1_real[,8])))) / ((as.numeric(unlist(step1_real[,8]))) * mean(as.numeric(unlist(step1_real[,3]))))
    Normalized_Exp_3 <- (as.numeric(unlist(step1_real[,4])) * mean(as.numeric(unlist(step1_real[,9])))) / ((as.numeric(unlist(step1_real[,9]))) * mean(as.numeric(unlist(step1_real[,4]))))
    Normalized_Exp_4 <- (as.numeric(unlist(step1_real[,5])) * mean(as.numeric(unlist(step1_real[,10])))) / ((as.numeric(unlist(step1_real[,10]))) * mean(as.numeric(unlist(step1_real[,5]))))
    Normalized_Exp_5 <- (as.numeric(unlist(step1_real[,6])) * mean(as.numeric(unlist(step1_real[,11])))) / ((as.numeric(unlist(step1_real[,11]))) * mean(as.numeric(unlist(step1_real[,6]))))

  }


  step2 <-
    cbind(step1_real,
          Normalized_Exp_1,
          Normalized_Exp_2,
          Normalized_Exp_3,
          Normalized_Exp_4,
          Normalized_Exp_5)

  # Step 3

  Log_2_Exp_1 <- log(Normalized_Exp_1, 2)
  Log_2_Exp_2 <- log(Normalized_Exp_2, 2)
  Log_2_Exp_3 <- log(Normalized_Exp_3, 2)
  Log_2_Exp_4 <- log(Normalized_Exp_4, 2)
  Log_2_Exp_5 <- log(Normalized_Exp_5, 2)

  step3 <- cbind(step2, Log_2_Exp_1, Log_2_Exp_2, Log_2_Exp_3, Log_2_Exp_4, Log_2_Exp_5)

  # Step 4

  Exp_Log_Avg <- rowMeans(step3[, c("Log_2_Exp_1", "Log_2_Exp_2", "Log_2_Exp_3", "Log_2_Exp_4", "Log_2_Exp_5")], na.rm = T)
  step4 <- cbind(step3, Exp_Log_Avg)

  # Step 5

  c1 <- step4[, c("Log_2_Exp_1", "Log_2_Exp_2", "Log_2_Exp_3", "Log_2_Exp_4", "Log_2_Exp_5")]
  SD_Exp <- apply(c1, 1, sd)
  step5 <- cbind(step4, SD_Exp)

  # Step 6

  f1 <- step5$Exp_Log_Avg
  mean_log_avg_Exp <- mean(f1)
  sd_log_avg_Exp <- sd(f1)

  z_score_Exp <- ((f1 - mean_log_avg_Exp)) / (sd_log_avg_Exp)
  step6 <- cbind(step5, z_score_Exp)

  # Step 7

  g1 <- abs(Exp_Log_Avg)
  g2 <- step6$SD_Exp

  T_value_Exp <- (g1 * sqrt(5) / g2)
  step7 <- cbind(step6, T_value_Exp)

  # Step 8

  h1 <- step7$Exp_Log_Avg
  h2 <- mean(h1)
  h3 <- step7$T_value_Exp
  p_value_Exp <- 2 * pt(h3, 4, lower.tail = F)
  step8 <- cbind(step7, p_value_Exp)

  # Step 9

  j1 <- step8$p_value_Exp
  log_10_Exp <- -log10(j1)
  step9 <- cbind(step8, log_10_Exp)

  # Step 10

  Hit_ident1_Exp <- magrittr::"%>%" (step9,
                                     dplyr::select(tidyselect::any_of(c("Accession", "Description"))))
  Log_Avg_Exp <- step9$Exp_Log_Avg
  Log_10_Exp <- step9$log_10_Exp
  Hit_ident_Exp <- cbind(Hit_ident1_Exp, Log_Avg_Exp, Log_10_Exp, Log_2_Exp_1, Log_2_Exp_2, Log_2_Exp_3, Log_2_Exp_4, Log_2_Exp_5, z_score_Exp)

  n1 <-
    (Hit_ident_Exp$z_score_Exp > SD_cutoff &
       Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff)) |
    (Hit_ident_Exp$z_score_Exp < -SD_cutoff & Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff))
  n2 <- as.data.frame(n1)

  Sig_Check_Exp <-
    dplyr::if_else(n2$n1 == "TRUE", "Significant", "Not Significant")
  step10 <- cbind(Hit_ident_Exp, p_value_Exp, Sig_Check_Exp)


  # Only get Semi-Tryptic Peptides


  STEPP_Confidence <- dplyr::filter(STEPP_Raw, STEPP_Raw[,grepl("Confidence", names(STEPP_Raw))] == "High")

  STEPP_Semi_Tryp <- dplyr::filter(STEPP_Confidence, grepl("1xTMT6plex [N-Term]", Modifications, fixed = TRUE))

  STEPP_removed_AA1 <- dplyr::mutate_all(STEPP_Semi_Tryp, ~gsub("\\[.*?].", "", .))
  STEPP_removed_AA2 <- dplyr::mutate_all(STEPP_removed_AA1, ~gsub("\\.\\[.*?]", "", .))


  # Go through TPP analysis

  STEPP_1 <- magrittr::"%>%" (STEPP_removed_AA2,
                              dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Accession"))))
  STEPP_2 <- STEPP_removed_AA2[, c("Annotated Sequence", "Modifications")]
  STEPP_3 <- STEPP_removed_AA2[, grep("^Abundances\\s\\(Grouped):", names(STEPP_removed_AA2), value = TRUE)]
  step1_STEPP <- cbind(STEPP_1, STEPP_2, STEPP_3)

  if ("Master Protein Accessions" %in% names(step1_STEPP)){
    step1_STEPP <- dplyr::rename(step1_STEPP, "Accession" = "Master Protein Accessions")
  }

  # Filter blanks

  NA_Removed_STEPP <- tidyr::drop_na(step1_STEPP)

    step1_real <-
      dplyr::filter(
        NA_Removed_STEPP,
        NA_Removed_STEPP[,4] > 0,
        NA_Removed_STEPP[,5] > 0,
        NA_Removed_STEPP[,6] > 0,
        NA_Removed_STEPP[,7] > 0,
        NA_Removed_STEPP[,8] > 0,
        NA_Removed_STEPP[,9] > 0,
        NA_Removed_STEPP[,10] > 0,
        NA_Removed_STEPP[,11] > 0,
        NA_Removed_STEPP[,12] > 0,
        NA_Removed_STEPP[,13] > 0)


  STEPP_merger <- merge(Hit_ident_Exp, step1_real)

  STEPP_merge <- STEPP_merger

  # Step 11

  Normalized_STEPP_T1 <- (as.numeric(unlist(STEPP_merge[,13])) * mean(as.numeric(unlist(STEPP_merge[,18])))) / ((as.numeric(unlist(STEPP_merge[,18]))) * mean(as.numeric(unlist(STEPP_merge[,13]))))
  Normalized_STEPP_T2 <- (as.numeric(unlist(STEPP_merge[,14])) * mean(as.numeric(unlist(STEPP_merge[,19])))) / ((as.numeric(unlist(STEPP_merge[,19]))) * mean(as.numeric(unlist(STEPP_merge[,14]))))
  Normalized_STEPP_T3 <- (as.numeric(unlist(STEPP_merge[,15])) * mean(as.numeric(unlist(STEPP_merge[,20])))) / ((as.numeric(unlist(STEPP_merge[,20]))) * mean(as.numeric(unlist(STEPP_merge[,15]))))
  Normalized_STEPP_T4 <- (as.numeric(unlist(STEPP_merge[,16])) * mean(as.numeric(unlist(STEPP_merge[,21])))) / ((as.numeric(unlist(STEPP_merge[,21]))) * mean(as.numeric(unlist(STEPP_merge[,16]))))
  Normalized_STEPP_T5 <- (as.numeric(unlist(STEPP_merge[,17])) * mean(as.numeric(unlist(STEPP_merge[,22])))) / ((as.numeric(unlist(STEPP_merge[,22]))) * mean(as.numeric(unlist(STEPP_merge[,17]))))

  step11 <-
    cbind(STEPP_merge,
          Normalized_STEPP_T1,
          Normalized_STEPP_T2,
          Normalized_STEPP_T3,
          Normalized_STEPP_T4,
          Normalized_STEPP_T5)

  step111 <- tidyr::drop_na(step11)

  # Step 12

  Log_2_STEPP_T1 <- log(Normalized_STEPP_T1, 2)
  Log_2_STEPP_T2 <- log(Normalized_STEPP_T2, 2)
  Log_2_STEPP_T3 <- log(Normalized_STEPP_T3, 2)
  Log_2_STEPP_T4 <- log(Normalized_STEPP_T4, 2)
  Log_2_STEPP_T5 <- log(Normalized_STEPP_T5, 2)

  step12 <- cbind(step111, Log_2_STEPP_T1, Log_2_STEPP_T2, Log_2_STEPP_T3, Log_2_STEPP_T4, Log_2_STEPP_T5)

  # Step 13

  Log_Avg_STEPP <- rowMeans(step12[, c("Log_2_STEPP_T1", "Log_2_STEPP_T2", "Log_2_STEPP_T3", "Log_2_STEPP_T4", "Log_2_STEPP_T5")])
  step13 <- cbind(step12, Log_Avg_STEPP)

  # Step  14

  rep1_exp1 <- step13$Log_2_STEPP_T1 - step13$Log_2_Exp_1
  rep2_exp2 <- step13$Log_2_STEPP_T2 - step13$Log_2_Exp_2
  rep3_exp3 <- step13$Log_2_STEPP_T3 - step13$Log_2_Exp_3
  rep4_exp4 <- step13$Log_2_STEPP_T4 - step13$Log_2_Exp_4
  rep5_exp5 <- step13$Log_2_STEPP_T5 - step13$Log_2_Exp_5
  step14 <- cbind(step13, rep1_exp1, rep2_exp2, rep3_exp3, rep4_exp4, rep5_exp5)

  # Step 15

  STEPP_Avg <- rowMeans(step14[, c("rep1_exp1", "rep2_exp2", "rep3_exp3", "rep4_exp4", "rep5_exp5")], na.rm = T)
  step15 <- cbind(step14, STEPP_Avg)

  # Step 16

  k1 <- step15[, c("rep1_exp1", "rep2_exp2", "rep3_exp3", "rep4_exp4", "rep5_exp5")]
  SD_STEPP <- apply(k1, 1, sd)
  step16 <- cbind(step15, SD_STEPP)

  # Step 17

  step17 <- tidyr::drop_na(step16)
  STEPPunique <- dplyr::n_distinct(step17$Accession)
  STEPPuniquePeptide <- nrow(step17[, c("Annotated Sequence", "Accession")])

  # Step 18

  y1 <- step17$STEPP_Avg
  STEPP_mean_avg <- mean(y1)
  STEPP_sd_avg <- sd(y1)
  z_score_STEPP <- ((y1 - STEPP_mean_avg)) / (STEPP_sd_avg)
  step18 <- cbind(step17, z_score_STEPP)

  # Step 19

  hg1 <- abs(STEPP_Avg)
  hg2 <- step18$SD_STEPP

  STEPP_T_value <- (hg1 * sqrt(5) / hg2)
  step19 <- cbind(step18, STEPP_T_value)

  # Step 20

  hk1 <- step19$STEPP_Avg
  hk2 <- mean(hg1)
  hk3 <- step19$STEPP_T_value
  p_value_STEPP <- 2 * pt(hk3, 4, lower.tail = F)
  step20 <- cbind(step19, p_value_STEPP)

  # Step 21

  jk1 <- step20$p_value_STEPP
  STEPP_log_10 <- -log10(jk1)
  step21 <- cbind(step20, STEPP_log_10)

  STEPP_Hit_ident1 <- magrittr::"%>%" (step21,
                                       dplyr::select(tidyselect::any_of(c("Accession", "Description"))))
  STEPP_Hit_ident2 <- magrittr::"%>%" (step21,
                                       dplyr::select(tidyselect::any_of(c("Annotated Sequence", "Modifications"))))
  STEPP_Log_Avg <- magrittr::"%>%" (step21,
                                    dplyr::select(tidyselect::any_of(c("STEPP_Avg", "SD_STEPP", "z_score_STEPP", "STEPP_T_value", "p_value_STEPP", "STEPP_log_10"))))
  STEPP_Hit_ident <- cbind(STEPP_Hit_ident1, STEPP_Hit_ident2, STEPP_Log_Avg)

  # Step 22

  if (correctedpvalue == FALSE){

  nk1 <-
    (STEPP_Hit_ident$z_score_STEPP > SD_cutoff &
       STEPP_Hit_ident$STEPP_log_10 > -log10(p_cutoff)) |
    (STEPP_Hit_ident$z_score_STEPP < -SD_cutoff & STEPP_Hit_ident$STEPP_log_10 > -log10(p_cutoff))
  nk2 <- as.data.frame(nk1)

  }

  if (correctedpvalue == TRUE){
    nk1 <-
      (STEPP_Hit_ident$z_score_STEPP > SD_cutoff &
         STEPP_Hit_ident$STEPP_log_10 > -log10((p_cutoff/STEPPuniquePeptide))) |
  (STEPP_Hit_ident$z_score_STEPP < -SD_cutoff & STEPP_Hit_ident$STEPP_log_10 > -log10((p_cutoff/STEPPuniquePeptide)))
nk2 <- as.data.frame(nk1)

  }


  Sig_Check_STEPP <-
    dplyr::if_else(nk2$nk1 == "TRUE", "Significant", "Not Significant")
  step22 <- cbind(STEPP_Hit_ident, Sig_Check_STEPP)

  ddk <- step22$Sig_Check_STEPP

  as.data.frame(ddk)
  abk <- (ddk == "Significant")
  STEPP_Number_Of_Hits <- length(which(abk))

  Hit_List <- step22[ step22$Sig_Check_STEPP == "Significant", ]
  STEPPuniquehits <- dplyr::n_distinct(Hit_List$Accession)
  STEPPExport <- unique(Hit_List$Accession)
  STEPPSigExport <- magrittr::"%>%" (step22,
                                     dplyr::select(tidyselect::any_of(c("Accession", "Description", "Annotated Sequence", "Modifications", "z_score_STEPP", "p_value_STEPP", "Sig_Check_STEPP"))))
  Hit_ListExportSTEPP <- magrittr::"%>%" (Hit_List,
                                          dplyr::select(tidyselect::any_of(c("Accession", "Description", "Annotated Sequence", "Modifications", "z_score_STEPP", "p_value_STEPP", "Sig_Check_STEPP"))))

  print(STEPPSigExport)
  print(paste("There are", STEPP_Number_Of_Hits, "Hits"))
  print(paste("There are", STEPPuniquehits, "Unique Hits"))
  print(paste(STEPPuniquePeptide, "Unique Peptides were assayed"))
  print(paste(STEPPunique, "Unique Proteins were assayed"))
  print(Hit_ListExportSTEPP)

  utils::write.table(STEPPExport, file = "STEPPPhenotypeUniqueHits.csv", row.names = FALSE, col.names = FALSE)
  utils::write.table(STEPPSigExport, file = "STEPPPhenotypeOut.csv", row.names = FALSE)
  utils::write.table(Hit_ListExportSTEPP, file = "STEPPPhenotypeHitsOut.csv", row.names = FALSE)

  # Visualizing Data

  if (plot == TRUE){

  Volcano_Plot <- ggplot2::ggplot(step22,
                         ggplot2::aes (
                           x = z_score_STEPP ,
                           y = STEPP_log_10 ,
                           label = Accession,
                           colour = Sig_Check_STEPP
                         )) + ggplot2::geom_point(size = 1.5, alpha = alpha) +
    ggplot2::labs(x= xlab, y= ylab) + ggplot2::scale_colour_manual(breaks = c("Not Significant", "Significant"), values = labelcolor) +
    ggplot2::expand_limits(x=0, y=0) +
    ggplot2::scale_x_continuous(breaks = -12:12) +
    ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = plottextsize)) + ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  print(Volcano_Plot)
  }

ggplot2::ggsave("Volcano_Plot_STEPP_Phenotype.png", plot = Volcano_Plot)

if (interactiveplot == TRUE){
  plotly::ggplotly(Volcano_Plot)
}

}


#' STEPP Ligand Analysis
#'
#' This function will take STEPP/LiP (TMT10plex) data directly from ProteomeDiscoverer (v. 2.3) and generate a volcano plot and hit list
#' based on user-defined SD and p-values
#'
#' @param Expression_data Raw ProteomeDiscoverer (v. 2.3) data
#' @param TPP_data Raw STEPP (TMT10plex) ProteomeDiscoverer (v. 2.3) data
#' @param SD_cutoff A positive integer (default = 2)
#' @param p_cutoff A positive integer between 0 and 1 (default = 0.05)
#' @param TMTplex 8, 10, 16 TMTplex (default = 10)
#' @param plot Display volcano plot (default = TRUE)
#' @param xlab X-axis Label (default = "Z Score")
#' @param ylab Y-axis Label (default = "- log10 (p value)")
#' @param labelcolor Color of Not significant and Significant data points (default = c("grey", "red"))
#' @param alpha Transparency of plot (default = 0.5)
#' @param plottextsize Text Size of plot (default = 12)
#' @param correctedpvalue Bonferroni corrected p-value (default = FALSE)
#' @param interactiveplot Plotly interactive Volcano plot (default = FALSE)
#' @return Volcano plot, Hit list, Total proteins assayed and Unique protein hits
#' @export

STEPP_Ligand.Fun <- function(STEPP_Raw, SD_cutoff = 2, p_cutoff = 0.05, TMTplex = 10, plot = TRUE, xlab = "Z Score", ylab = "- log10 (p value)", labelcolor = c("grey", "red"), alpha = 0.5, plottextsize = 12, correctedpvalue = FALSE, interactiveplot = FALSE){

  if (TMTplex !=10)
    stop("This TMTplex is not yet supported in this version.\n")

  if (SD_cutoff < 0)
    stop("SD_cutoff must be non-negative.\n")

  if (p_cutoff < 0)
    stop("p_value cutoff must be non-negative.\n")

  # Only get Semi-Tryptic Peptides


  STEPP_Confidence <- dplyr::filter(STEPP_Raw, STEPP_Raw[,grepl("Confidence", names(STEPP_Raw))] == "High")

  STEPP_Semi_Tryp <- dplyr::filter(STEPP_Confidence, grepl("1xTMT6plex [N-Term]", Modifications, fixed = TRUE))

  STEPP_removed_AA1 <- dplyr::mutate_all(STEPP_Semi_Tryp, ~gsub("\\[.*?].", "", .))
  STEPP_removed_AA2 <- dplyr::mutate_all(STEPP_removed_AA1, ~gsub("\\.\\[.*?]", "", .))


  # Go through TPP analysis

  STEPP_1 <- magrittr::"%>%" (STEPP_removed_AA2,
                              dplyr::select(tidyselect::any_of(c("Master Protein Accessions", "Accession"))))
  STEPP_2 <- STEPP_removed_AA2[, c("Annotated Sequence", "Modifications")]
  STEPP_3 <- STEPP_removed_AA2[, grep("^Abundances\\s\\(Grouped):", names(STEPP_removed_AA2), value = TRUE)]
  step1_STEPP <- cbind(STEPP_1, STEPP_2, STEPP_3)

  if ("Master Protein Accessions" %in% names(step1_STEPP)){
    step1_STEPP <- dplyr::rename(step1_STEPP, "Accession" = "Master Protein Accessions")
  }

  # Filter blanks

  NA_Removed_STEPP <- tidyr::drop_na(step1_STEPP)

    step1_real <-
      dplyr::filter(
        NA_Removed_STEPP,
        NA_Removed_STEPP[,4] > 0,
        NA_Removed_STEPP[,5] > 0,
        NA_Removed_STEPP[,6] > 0,
        NA_Removed_STEPP[,7] > 0,
        NA_Removed_STEPP[,8] > 0,
        NA_Removed_STEPP[,9] > 0,
        NA_Removed_STEPP[,10] > 0,
        NA_Removed_STEPP[,11] > 0,
        NA_Removed_STEPP[,12] > 0,
        NA_Removed_STEPP[,13] > 0)

  STEPP_merge <- step1_real

  # Step 11

  Normalized_STEPP_T1 <- (as.numeric(unlist(STEPP_merge[,4])) * mean(as.numeric(unlist(STEPP_merge[,9])))) / ((as.numeric(unlist(STEPP_merge[,9]))) * mean(as.numeric(unlist(STEPP_merge[,4]))))
  Normalized_STEPP_T2 <- (as.numeric(unlist(STEPP_merge[,5])) * mean(as.numeric(unlist(STEPP_merge[,10])))) / ((as.numeric(unlist(STEPP_merge[,10]))) * mean(as.numeric(unlist(STEPP_merge[,5]))))
  Normalized_STEPP_T3 <- (as.numeric(unlist(STEPP_merge[,6])) * mean(as.numeric(unlist(STEPP_merge[,11])))) / ((as.numeric(unlist(STEPP_merge[,11]))) * mean(as.numeric(unlist(STEPP_merge[,6]))))
  Normalized_STEPP_T4 <- (as.numeric(unlist(STEPP_merge[,7])) * mean(as.numeric(unlist(STEPP_merge[,12])))) / ((as.numeric(unlist(STEPP_merge[,12]))) * mean(as.numeric(unlist(STEPP_merge[,7]))))
  Normalized_STEPP_T5 <- (as.numeric(unlist(STEPP_merge[,8])) * mean(as.numeric(unlist(STEPP_merge[,13])))) / ((as.numeric(unlist(STEPP_merge[,13]))) * mean(as.numeric(unlist(STEPP_merge[,8]))))


  step11 <-
    cbind(STEPP_merge,
          Normalized_STEPP_T1,
          Normalized_STEPP_T2,
          Normalized_STEPP_T3,
          Normalized_STEPP_T4,
          Normalized_STEPP_T5)

  step111 <- tidyr::drop_na(step11)

  # Step 12

  Log_2_STEPP_T1 <- log(Normalized_STEPP_T1, 2)
  Log_2_STEPP_T2 <- log(Normalized_STEPP_T2, 2)
  Log_2_STEPP_T3 <- log(Normalized_STEPP_T3, 2)
  Log_2_STEPP_T4 <- log(Normalized_STEPP_T4, 2)
  Log_2_STEPP_T5 <- log(Normalized_STEPP_T5, 2)

  step12 <- cbind(step111, Log_2_STEPP_T1, Log_2_STEPP_T2, Log_2_STEPP_T3, Log_2_STEPP_T4, Log_2_STEPP_T5)

  # Step 13

  Log_Avg_STEPP <- rowMeans(step12[, c("Log_2_STEPP_T1", "Log_2_STEPP_T2", "Log_2_STEPP_T3", "Log_2_STEPP_T4", "Log_2_STEPP_T5")])
  step13 <- cbind(step12, Log_Avg_STEPP)

  # Step  14

  step14 <- step13

  # Step 15

  step15 <- step14

  # Step 16

  k1 <- step15[, c("Log_2_STEPP_T1", "Log_2_STEPP_T2", "Log_2_STEPP_T3", "Log_2_STEPP_T4", "Log_2_STEPP_T5")]
  SD_STEPP <- apply(k1, 1, sd)
  step16 <- cbind(step15, SD_STEPP)

  # Step 17

  step17 <- tidyr::drop_na(step16)
  STEPPunique <- dplyr::n_distinct(step17$Accession)
  STEPPuniquePeptide <- nrow(step17[, c("Annotated Sequence", "Accession")])

  # Step 18

  y1 <- step17$Log_Avg_STEPP
  STEPP_mean_avg <- mean(y1)
  STEPP_sd_avg <- sd(y1)
  z_score_STEPP <- ((y1 - STEPP_mean_avg)) / (STEPP_sd_avg)
  step18 <- cbind(step17, z_score_STEPP)

  # Step 19

  hg1 <- abs(Log_Avg_STEPP)
  hg2 <- step18$SD_STEPP

  STEPP_T_value <- (hg1 * sqrt(5) / hg2)
  step19 <- cbind(step18, STEPP_T_value)

  # Step 20

  hk1 <- step19$Log_Avg_STEPP
  hk2 <- mean(hg1)
  hk3 <- step19$STEPP_T_value
  p_value_STEPP <- 2 * pt(hk3, 4, lower.tail = F)
  step20 <- cbind(step19, p_value_STEPP)

  # Step 21

  jk1 <- step20$p_value_STEPP
  STEPP_log_10 <- -log10(jk1)
  step21 <- cbind(step20, STEPP_log_10)

  STEPP_Hit_ident1 <- magrittr::"%>%" (step21,
                                       dplyr::select(tidyselect::any_of(c("Accession", "Annotated Sequence", "Modifications"))))
  STEPP_Log_Avg <- magrittr::"%>%" (step21,
                                    dplyr::select(tidyselect::any_of(c("Log_Avg_STEPP", "SD_STEPP", "z_score_STEPP", "STEPP_T_value", "p_value_STEPP", "STEPP_log_10"))))
  STEPP_Hit_ident <- cbind(STEPP_Hit_ident1, STEPP_Log_Avg)

  # Step 22

  if (correctedpvalue == FALSE){

  nk1 <-
    (STEPP_Hit_ident$z_score_STEPP > SD_cutoff &
       STEPP_Hit_ident$STEPP_log_10 > -log10(p_cutoff)) |
    (STEPP_Hit_ident$z_score_STEPP < -SD_cutoff & STEPP_Hit_ident$STEPP_log_10 > -log10(p_cutoff))
  nk2 <- as.data.frame(nk1)

  }

  if (correctedpvalue == TRUE){
    nk1 <-
      (STEPP_Hit_ident$z_score_STEPP > SD_cutoff &
         STEPP_Hit_ident$STEPP_log_10 > -log10((p_cutoff/STEPPuniquePeptide))) |
      (STEPP_Hit_ident$z_score_STEPP < -SD_cutoff & STEPP_Hit_ident$STEPP_log_10 > -log10((p_cutoff/STEPPuniquePeptide)))
    nk2 <- as.data.frame(nk1)

  }

  Sig_Check_STEPP <-
    dplyr::if_else(nk2$nk1 == "TRUE", "Significant", "Not Significant")
  step22 <- cbind(STEPP_Hit_ident, Sig_Check_STEPP)

  ddk <- step22$Sig_Check_STEPP

  as.data.frame(ddk)
  abk <- (ddk == "Significant")
  STEPP_Number_Of_Hits <- length(which(abk))

  Hit_List <- step22[ step22$Sig_Check_STEPP == "Significant", ]
  STEPPuniquehits <- dplyr::n_distinct(Hit_List$Accession)
  STEPPExport <- unique(Hit_List$Accession)
  STEPPSigExport <- magrittr::"%>%" (step22,
                                     dplyr::select(tidyselect::any_of(c("Accession", "Annotated Sequence", "Modifications", "z_score_STEPP", "p_value_STEPP", "Sig_Check_STEPP"))))
  Hit_ListExportSTEPP <- magrittr::"%>%" (Hit_List,
                                          dplyr::select(tidyselect::any_of(c("Accession", "Annotated Sequence", "Modifications", "z_score_STEPP", "p_value_STEPP", "Sig_Check_STEPP"))))

  print(STEPPSigExport)
  print(paste("There are", STEPP_Number_Of_Hits, "Hits"))
  print(paste("There are", STEPPuniquehits, "Unique Hits"))
  print(paste(STEPPuniquePeptide, "Unique Peptides were assayed"))
  print(paste(STEPPunique, "Unique Proteins were assayed"))
  print(Hit_ListExportSTEPP)

  utils::write.table(STEPPExport, file = "STEPPLigandUniqueHits.csv", row.names = FALSE, col.names = FALSE)
  utils::write.table(STEPPSigExport, file = "STEPPLigandOut.csv", row.names = FALSE)
  utils::write.table(Hit_ListExportSTEPP, file = "STEPPLigandHitsOut.csv", row.names = FALSE)

  # Visualizing Data

  if (plot == TRUE){

    Volcano_Plot <- ggplot2::ggplot(step22,
                                    ggplot2::aes (
                                      x = z_score_STEPP ,
                                      y = STEPP_log_10 ,
                                      label = Accession,
                                      colour = Sig_Check_STEPP
                                    )) + ggplot2::geom_point(size = 1.5, alpha = alpha) +
      ggplot2::labs(x= xlab, y= ylab) + ggplot2::scale_colour_manual(breaks = c("Not Significant", "Significant"), values = labelcolor) +
      ggplot2::expand_limits(x=0, y=0) +
      ggplot2::scale_x_continuous(breaks = -12:12) +
      ggplot2::theme_bw() + ggplot2::theme(text = ggplot2::element_text(size = plottextsize)) + ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
      ggplot2::theme(legend.title = ggplot2::element_blank())

    print(Volcano_Plot)
  }

  ggplot2::ggsave("Volcano_Plot_STEPP_Ligand.png", plot = Volcano_Plot)

  if (interactiveplot == TRUE){
    plotly::ggplotly(Volcano_Plot)
  }

}


#' Two Method Comparison Venn Diagram
#'
#' This function will make a Venn Diagram showing unique hits for two different proteomics methods
#'
#'
#' @param data1 Unqiue Hit .csv Output method 1
#' @param data2 Unique Hit .csv Output method 2
#' @param name1 Method 1 Name
#' @param name2 Method 2 Name
#' @param filename Filename for Output (default = 'Two_Comparison_Venn_Diagram')
#' @param imagetype Image Output Extension (default = png)
#' @param textcolor Color of Circle's Circumference (default = black)
#' @param outlinetype Circle Outline Dash pattern (default = dotted)
#' @param outlinewidth Circle Outline Thickness (default = 1)
#' @param areafill Color of Circles (default = c("#6495ed", "#bf3eff"))
#' @param alpha Transparency of Circles (default = 0.3)
#' @param arealabelsize Size of Area Labels (default = 0.65)
#' @param categorynamesize Size of Category Names (default = 0.5)
#' @param categorycolor Color of Category Names (default = black)
#' @param scaled If TRUE Circle area will scale based on value (default = FALSE)
#' @return Venn Diagram
#' @export


Venn2.Fun <- function(data1, data2, name1, name2, filename = 'Two_Comparison_Venn_Diagram.PNG', imagetype = "png", textcolor = "black", outlinetype = "dotted", outlinewidth =1, areafill = c("#6495ed", "#bf3eff"),
                      alpha = c(0.3,0.3), arealabelsize = 0.65, categorynamesize = 0.5, categorycolor = "black", scaled = FALSE){

  z <- utils::read.csv2(file = data1, header = FALSE)
  z1 <- utils::read.csv2(file = data2, header = FALSE)

  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

  VennDiagram::venn.diagram(x = c(z, z1),
               category.names = c(name1, name2),
               filename = filename,
               output = TRUE ,
               imagetype= imagetype ,
               height = 480 ,
               width = 480 ,
               resolution = 300,
               compression = "lzw",
               lwd = outlinewidth,
               col= categorycolor,
               lty = outlinetype,
               fill = areafill,
               alpha = alpha,
               cex = arealabelsize,
               fontfamily = "serif",
               cat.cex = categorynamesize,
               cat.default.pos = "outer",
               fontface = "bold",
               cat.pos = c(-1, 1),
               cat.dist = c(0.055, 0.055),
               cat.fontfamily = "serif",
               cat.col = textcolor,
               scaled = scaled

  )




}

#' Three Method Comparison Venn Diagram
#'
#' This function will make a Venn Diagram showing unique hits for three different proteomics methods
#'
#'
#' @param data1 Unqiue Hit .csv Output method 1
#' @param data2 Unique Hit .csv Output method 2
#' @param data3 Unique Hit .csv Output method 3
#' @param name1 Method 1 Name
#' @param name2 Method 2 Name
#' @param name3 Method 3 Name
#' @param filename Filename for Output (default = 'Three_Comparison_Venn_Diagram')
#' @param imagetype Image Output Extension (default = png)
#' @param textcolor Color of Circle's Circumference (default = black)
#' @param outlinetype Circle Outline Dash pattern (default = dotted)
#' @param outlinewidth Circle Outline Thickness (default = 1)
#' @param areafill Color of Circles (default = c("#6495ed", "#bf3eff", "#ffffba"))
#' @param alpha Transparency of Circles (default = 0.3)
#' @param arealabelsize Size of Area Labels (default = 0.5)
#' @param categorynamesize Size of Category Names (default = 0.5)
#' @param categorycolor Color of Category Names (default = black)
#' @param scaled If TRUE Circle area will scale based on value (default = FALSE)
#' @return Venn Diagram
#' @export


Venn3.Fun <- function(data1, data2, data3, name1, name2, name3, filename = 'Three_Comparison_Venn_Diagram.PNG', imagetype = "png", textcolor = "black", outlinetype = "dotted", outlinewidth = 1, areafill = c("#6495ed", "#bf3eff", "#ffffba"),
                      alpha = c(0.3,0.3, 0.3), arealabelsize = 0.5, categorynamesize = 0.5, categorycolor = "black", scaled = FALSE){

  z <- utils::read.csv2(file = data1, header = FALSE)
  z1 <- utils::read.csv2(file = data2, header = FALSE)
  z2 <- utils::read.csv2(file = data3, header = FALSE)

  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

  VennDiagram::venn.diagram(x = c(z, z1, z2),
                            category.names = c(name1, name2, name3),
                            filename = filename,
                            output = TRUE ,
                            imagetype= imagetype ,
                            height = 480 ,
                            width = 480 ,
                            resolution = 300,
                            compression = "lzw",
                            lwd = outlinewidth,
                            col= categorycolor,
                            lty = outlinetype,
                            fill = areafill,
                            alpha = alpha,
                            cex = arealabelsize,
                            fontfamily = "serif",
                            cat.cex = categorynamesize,
                            cat.default.pos = "outer",
                            fontface = "bold",
                            cat.pos = c(-27, 27, 135),
                            cat.dist = c(0.055, 0.055, 0.085),
                            cat.fontfamily = "serif",
                            cat.col = textcolor,
                            scaled = scaled

  )
}

#' Four Method Comparison Venn Diagram
#'
#' This function will make a Venn Diagram showing unique hits for three different proteomics methods
#'
#'
#' @param data1 Unqiue Hit .csv Output method 1
#' @param data2 Unique Hit .csv Output method 2
#' @param data3 Unique Hit .csv Output method 3
#' @param data4 Unique Hit .csv Output method 4
#' @param name1 Method 1 Name
#' @param name2 Method 2 Name
#' @param name3 Method 3 Name
#' @param name4 Method 4 Name
#' @param filename Filename for Output (default = 'Four_Comparison_Venn_Diagram')
#' @param imagetype Image Output Extension (default = png)
#' @param textcolor Color of Circle's Circumference (default = black)
#' @param outlinetype Circle Outline Dash pattern (default = dotted)
#' @param outlinewidth Circle Outline Thickness (default = 1)
#' @param areafill Color of Circles (default = c("#6495ed", "#bf3eff", "#ffffba", "#90ee90"))
#' @param alpha Transparency of Circles (default = 0.3)
#' @param arealabelsize Size of Area Labels (default = 0.5)
#' @param categorynamesize Size of Category Names (default = 0.5)
#' @param categorycolor Color of Category Names (default = black)
#' @param scaled If TRUE Circle area will scale based on value (default = FALSE)
#' @return Venn Diagram
#' @export


Venn4.Fun <- function(data1, data2, data3, data4, name1, name2, name3, name4, filename = 'Four_Comparison_Venn_Diagram.PNG', imagetype = "png", textcolor = "black", outlinetype = "dotted", outlinewidth = 1, areafill = c("#6495ed", "#bf3eff", "#ffffba", "#90ee90"),
                      alpha = c(0.3,0.3, 0.3, 0.3), arealabelsize = 0.5, categorynamesize = 0.5, categorycolor = "black", scaled = FALSE){

  z <- utils::read.csv2(file = data1, header = FALSE)
  z1 <- utils::read.csv2(file = data2, header = FALSE)
  z2 <- utils::read.csv2(file = data3, header = FALSE)
  z3 <- utils::read.csv2(file = data4, header = FALSE)

  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

  VennDiagram::venn.diagram(x = c(z, z1, z2, z3),
                            category.names = c(name1, name2, name3, name4),
                            filename = filename,
                            output = TRUE ,
                            imagetype= imagetype ,
                            height = 480 ,
                            width = 480 ,
                            resolution = 300,
                            compression = "lzw",
                            lwd = outlinewidth,
                            col= categorycolor,
                            lty = outlinetype,
                            fill = areafill,
                            alpha = alpha,
                            cex = arealabelsize,
                            fontfamily = "serif",
                            cat.cex = categorynamesize,
                            cat.default.pos = "outer",
                            fontface = "bold",
                            cat.fontfamily = "serif",
                            cat.col = textcolor,
                            scaled = scaled

  )
}






















