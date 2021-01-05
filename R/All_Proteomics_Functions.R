#' Expression Level Analysis
#'
#' This function will take expression data directly from ProteomeDiscoverer (v. 2.3) and generate a volcano plot and hit list
#' based on user-defined SD and p-values
#'
#' @param Expression_data Raw ProteomeDiscoverer (v. 2.3) data
#' @param SD_cutoff A number
#' @param p_cutoff A number
#' @return Volcano plot, Hit list, Total proteins assayed, and Unique protein hits
#' @export

ExpressionLevel.Fun <- function(Expression_data, SD_cutoff, p_cutoff){

  High_Confidence_Exp_Level <-
    Expression_data[Expression_data$`Protein FDR Confidence: Combined` == "High",]

  # This dataframe takes all useful columns from the raw data and will be outputed at the end for reference

  Exp_Level_Concise <- High_Confidence_Exp_Level[, 2:27]

  # Remove all blank rows in Dataframe

  NA_Removed <- tidyr::drop_na(Exp_Level_Concise)

  # Step 1

  Exp_Level_1 <- NA_Removed[, 3:4]
  Exp_Level_2 <- NA_Removed[, 17:26]
  step1 <- cbind(Exp_Level_1, Exp_Level_2)

  # Step 2

  Normalized_Exp_1 <- (step1$`Abundances (Grouped): 126` * mean(step1$`Abundances (Grouped): 129N`)) / (step1$`Abundances (Grouped): 129N` * mean(step1$`Abundances (Grouped): 126`))
  Normalized_Exp_2 <- (step1$`Abundances (Grouped): 127N` * mean(step1$`Abundances (Grouped): 129C`)) / (step1$`Abundances (Grouped): 129C` * mean(step1$`Abundances (Grouped): 127N`))
  Normalized_Exp_3 <- (step1$`Abundances (Grouped): 127C` * mean(step1$`Abundances (Grouped): 130N`)) / (step1$`Abundances (Grouped): 130N` * mean(step1$`Abundances (Grouped): 127C`))
  Normalized_Exp_4 <- (step1$`Abundances (Grouped): 128N` * mean(step1$`Abundances (Grouped): 130C`)) / (step1$`Abundances (Grouped): 130C` * mean(step1$`Abundances (Grouped): 128N`))
  Normalized_Exp_5 <- (step1$`Abundances (Grouped): 128C` * mean(step1$`Abundances (Grouped): 131`)) / (step1$`Abundances (Grouped): 131` * mean(step1$`Abundances (Grouped): 128C`))

  step2 <-
    cbind(step1,
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

  Exp_Log_Avg <- rowMeans(step3[, 18:22], na.rm = T)
  step4 <- cbind(step3, Exp_Log_Avg)

  # Step 5

  c1 <- step4[, 18:22]
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

  Hit_ident1_Exp <- step9[, 1:2]
  Log_Avg_Exp <- step9$Exp_Log_Avg
  Log_10_Exp <- step9$log_10_Exp
  Hit_ident_Exp <- cbind(Hit_ident1_Exp, Log_Avg_Exp, Log_10_Exp, Log_2_Exp_1, Log_2_Exp_2, Log_2_Exp_3, Log_2_Exp_4, Log_2_Exp_5, z_score_Exp)

  n1 <-
    (Hit_ident_Exp$z_score_Exp > SD_cutoff &
       Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff)) |
    (Hit_ident_Exp$Log_Avg < -SD_cutoff & Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff))
  n2 <- as.data.frame(n1)

  Sig_Check_Exp <-
    dplyr::if_else(n2$n1 == "TRUE", "Significant", "Not Significant")
  step10 <- cbind(Hit_ident_Exp, Sig_Check_Exp)

  ddkd <- step10$Sig_Check_Exp

  as.data.frame(ddkd)
  abkk <- (ddkd == "Significant")
  Exp_Number_Of_Hits <- length(which(abkk))

  Hit_List <- step10[ step10$Sig_Check_Exp == "Significant", ]
  Expunique <- dplyr::n_distinct(step10$Accession)
  ExpExport <- unique(Hit_List$Accession)

  print(step10)
  print(paste("There are", Exp_Number_Of_Hits, "Hits"))
  print(paste(Expunique, "Proteins were assayed"))
  print(Hit_List)

  utils::write.table(ExpExport, file = "ExpressionLevelUniqueHits.csv", row.names = FALSE, col.names = FALSE)

  Volcano_Plot <- ggplot2::ggplot(step10,
                                  ggplot2::aes (
                                    x = z_score_Exp ,
                                    y = Log_10_Exp ,
                                    label = Accession,
                                    colour = Sig_Check_Exp
                                  )) + ggplot2::geom_point(size = 1.5, alpha = 0.5) +
    ggplot2::labs(x="Z Score",y="- log10 (p value)") + ggplot2::scale_colour_manual(breaks = c("Not Significant","Significant"),values = c("grey","red")) +
    ggplot2::scale_x_continuous(breaks = seq(-8,8,by = 1), limits=c(-7, 7)) +
    ggplot2::scale_y_continuous(breaks = seq(0,8,by=1), limits=c(0, 7.5)) +
    ggplot2::expand_limits(x=0, y=0) +
    ggplot2::theme_bw()+ ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  print(Volcano_Plot)

}

#' OnePotTPP Analysis
#'
#' This function will take expression data and OnePotTPP (TMT10plex) data directly from ProteomeDiscoverer (v. 2.3) and generate a volcano plot and hit list
#' based on user-defined SD and p-values
#'
#' @param Expression_data Raw ProteomeDiscoverer (v. 2.3) data
#' @param TPP_data Raw OnePotTPP (TMT10plex) ProteomeDiscoverer (v. 2.3) data
#' @param SD_cutoff A number
#' @param p_cutoff A number
#' @return Volcano plot, Hit list, Total proteins assayed, and Unique protein hits
#' @return Hit list, volcano plot
#' @export

OnePotTPP.Fun <- function(Expression_data, TPP_Raw, SD_cutoff, p_cutoff){

# Define SD and p-value cutoffs
# Load "Expression_data" and "TPP_Raw"

options(max.print = 25000)

# We can vary this to select for High | Medium | Low or both

High_Confidence_Exp_Level <-
  Expression_data[Expression_data$`Protein FDR Confidence: Combined` == "High",]

# This dataframe takes all useful columns from the raw data and will be outputed at the end for reference

Exp_Level_Concise <- High_Confidence_Exp_Level[, 2:27]

# Remove all blank rows in Dataframe

NA_Removed <- tidyr::drop_na(Exp_Level_Concise)

# Step 1

Exp_Level_1 <- NA_Removed[, 3:4]
Exp_Level_2 <- NA_Removed[, 17:26]
step1 <- cbind(Exp_Level_1, Exp_Level_2)

# Step 2

Normalized_Exp_1 <- (step1$`Abundances (Grouped): 126` * mean(step1$`Abundances (Grouped): 129N`)) / (step1$`Abundances (Grouped): 129N` * mean(step1$`Abundances (Grouped): 126`))
Normalized_Exp_2 <- (step1$`Abundances (Grouped): 127N` * mean(step1$`Abundances (Grouped): 129C`)) / (step1$`Abundances (Grouped): 129C` * mean(step1$`Abundances (Grouped): 127N`))
Normalized_Exp_3 <- (step1$`Abundances (Grouped): 127C` * mean(step1$`Abundances (Grouped): 130N`)) / (step1$`Abundances (Grouped): 130N` * mean(step1$`Abundances (Grouped): 127C`))
Normalized_Exp_4 <- (step1$`Abundances (Grouped): 128N` * mean(step1$`Abundances (Grouped): 130C`)) / (step1$`Abundances (Grouped): 130C` * mean(step1$`Abundances (Grouped): 128N`))
Normalized_Exp_5 <- (step1$`Abundances (Grouped): 128C` * mean(step1$`Abundances (Grouped): 131`)) / (step1$`Abundances (Grouped): 131` * mean(step1$`Abundances (Grouped): 128C`))

step2 <-
  cbind(step1,
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

Exp_Log_Avg <- rowMeans(step3[, 18:22], na.rm = T)
step4 <- cbind(step3, Exp_Log_Avg)

# Step 5

c1 <- step4[, 18:22]
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

Hit_ident1_Exp <- step9[, 1:2]
Log_Avg_Exp <- step9$Exp_Log_Avg
Log_10_Exp <- step9$log_10_Exp
Hit_ident_Exp <- cbind(Hit_ident1_Exp, Log_Avg_Exp, Log_10_Exp, Log_2_Exp_1, Log_2_Exp_2, Log_2_Exp_3, Log_2_Exp_4, Log_2_Exp_5, z_score_Exp)

n1 <-
  (Hit_ident_Exp$z_score_Exp > SD_cutoff &
     Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff)) |
  (Hit_ident_Exp$Log_Avg < -SD_cutoff & Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff))
n2 <- as.data.frame(n1)

Sig_Check_Exp <-
  dplyr::if_else(n2$n1 == "TRUE", "Significant", "Not Significant")
step10 <- cbind(Hit_ident_Exp, Sig_Check_Exp)

# TPP_Data

TPP_Confidence <- TPP_Raw[TPP_Raw$`Protein FDR Confidence: Combined` == "High",]

Data2 <- TPP_Confidence[, 2:27]

# Remove all blank rows in Dataframe

NA_Removed_TPP <- tidyr::drop_na(Data2)

TPP_1 <- NA_Removed_TPP[, 3:4]
TPP_2 <- NA_Removed_TPP[, 17:26]
step1_TPP <- cbind(TPP_1, TPP_2)

TPP_merger <- merge(Hit_ident_Exp, step1_TPP)
TPP_merge <-
  dplyr::filter(
    TPP_merger,
    `Abundances (Grouped): 126` > 0,
    `Abundances (Grouped): 127N` > 0,
    `Abundances (Grouped): 127C` > 0,
    `Abundances (Grouped): 128N` > 0,
    `Abundances (Grouped): 128C` > 0,
    `Abundances (Grouped): 129N` > 0,
    `Abundances (Grouped): 129C` > 0,
    `Abundances (Grouped): 130N` > 0,
    `Abundances (Grouped): 130C` > 0,
    `Abundances (Grouped): 131` > 0
  )

# Step 11

Normalized_TPP_T1 <- (TPP_merge$`Abundances (Grouped): 126` * mean(TPP_merge$`Abundances (Grouped): 129N`)) / (TPP_merge$`Abundances (Grouped): 129N` * mean(TPP_merge$`Abundances (Grouped): 126`))
Normalized_TPP_T2 <- (TPP_merge$`Abundances (Grouped): 127N` * mean(TPP_merge$`Abundances (Grouped): 129C`)) / (TPP_merge$`Abundances (Grouped): 129C` * mean(TPP_merge$`Abundances (Grouped): 127N`))
Normalized_TPP_T3 <- (TPP_merge$`Abundances (Grouped): 127C` * mean(TPP_merge$`Abundances (Grouped): 130N`)) / (TPP_merge$`Abundances (Grouped): 130N` * mean(TPP_merge$`Abundances (Grouped): 127C`))
Normalized_TPP_T4 <- (TPP_merge$`Abundances (Grouped): 128N` * mean(TPP_merge$`Abundances (Grouped): 130C`)) / (TPP_merge$`Abundances (Grouped): 130C` * mean(TPP_merge$`Abundances (Grouped): 128N`))
Normalized_TPP_T5 <- (TPP_merge$`Abundances (Grouped): 128C` * mean(TPP_merge$`Abundances (Grouped): 131`)) / (TPP_merge$`Abundances (Grouped): 131` * mean(TPP_merge$`Abundances (Grouped): 128C`))

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

Log_Avg_TPP <- rowMeans(step12[, 25:29])
step13 <- cbind(step12, Log_Avg_TPP)

# Step  14

rep1_exp1 <- step13$Log_2_TPP_T1 - step13$Log_2_Exp_1
rep2_exp2 <- step13$Log_2_TPP_T2 - step13$Log_2_Exp_2
rep3_exp3 <- step13$Log_2_TPP_T3 - step13$Log_2_Exp_3
rep4_exp4 <- step13$Log_2_TPP_T4 - step13$Log_2_Exp_4
rep5_exp5 <- step13$Log_2_TPP_T5 - step13$Log_2_Exp_5
step14 <- cbind(step13, rep1_exp1, rep2_exp2, rep3_exp3, rep4_exp4, rep5_exp5)

# Step 15

TPP_Avg <- rowMeans(step14[, 31:35], na.rm = T)
step15 <- cbind(step14, TPP_Avg)

# Step 16

k1 <- step15[, 31:35]
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
TPP_p_value <- 2 * pt(hk3, 4, lower.tail = F)
step20 <- cbind(step19, TPP_p_value)

# Step 21

jk1 <- step20$TPP_p_value
TPP_log_10 <- -log10(jk1)
step21 <- cbind(step20, TPP_log_10)

TPP_Hit_ident1 <- step21[, 1:2]
TPP_Log_Avg <- step21[, 36:41]
TPP_Hit_ident <- cbind(TPP_Hit_ident1, TPP_Log_Avg, TPP_log_10)

# Step 22

nk1 <-
  (TPP_Hit_ident$z_score_TPP > SD_cutoff &
     TPP_Hit_ident$TPP_log_10 > -log10(p_cutoff)) |
  (TPP_Hit_ident$z_score_TPP < -SD_cutoff & TPP_Hit_ident$TPP_log_10 > -log10(p_cutoff))
nk2 <- as.data.frame(nk1)

TPP_Sig_Check <-
  dplyr::if_else(nk2$nk1 == "TRUE", "Significant", "Not Significant")
step22 <- cbind(TPP_Hit_ident, TPP_Sig_Check)

ddk <- step22$TPP_Sig_Check

as.data.frame(ddk)
abk <- (ddk == "Significant")
TPP_Number_Of_Hits <- length(which(abk))

Hit_List <- step22[ step22$TPP_Sig_Check == "Significant", ]
TPPuniquehits <- dplyr::n_distinct(Hit_List$Accession)
TPPExport <- unique(Hit_List$Accession)
print(step22)
print(paste("There are", TPP_Number_Of_Hits, "Hits"))
print(paste("There are", TPPuniquehits, "Unique Hits"))
print(paste(TPPunique, "Proteins were assayed"))
print(Hit_List)

utils::write.table(TPPExport, file = "OnePotTPPUniqueHits.csv", row.names = FALSE, col.names = FALSE)

# Visualizing Data

Volcano_Plot <- ggplot2::ggplot(step22,
         ggplot2::aes (
           x = z_score_TPP ,
           y = TPP_log_10 ,
           label = Accession,
           colour = TPP_Sig_Check
         )) + ggplot2::geom_point(size = 1.5, alpha = 0.5) +
  ggplot2::labs(x="Z Score",y="- log10 (p value)") + ggplot2::scale_colour_manual(breaks = c("Not Significant","Significant"),values = c("grey","red")) +
  ggplot2::scale_x_continuous(breaks = seq(-8,8,by = 1), limits=c(-7, 7)) +
  ggplot2::scale_y_continuous(breaks = seq(0,8,by=1), limits=c(0, 7.5)) +
  ggplot2::expand_limits(x=0, y=0) +
  ggplot2::theme_bw()+ ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
  ggplot2::theme(legend.title = ggplot2::element_blank())

print(Volcano_Plot)

}

#' OnePotSPROX Analysis
#'
#' This function will take expression data and OnePotSPROX (TMT10plex) data directly from ProteomeDiscoverer (v. 2.3) and generate a volcano plot and hit list
#' based on user-defined SD and p-values
#'
#' @param Expression_data Raw ProteomeDiscoverer (v. 2.3) data
#' @param TPP_data Raw OnePotSPROX (TMT10plex) ProteomeDiscoverer (v. 2.3) data
#' @param SD_cutoff A number
#' @param p_cutoff A number
#' @return Volcano plot, Hit list, Total proteins assayed and Unique protein hits
#' @return Hit list, volcano plot
#' @export

OnePotSPROX.Fun <- function(Expression_data, SPROX_Raw, SD_cutoff, p_cutoff){

  High_Confidence_Exp_Level <-
    Expression_data[Expression_data$`Protein FDR Confidence: Combined` == "High",]

  # This dataframe takes all useful columns from the raw data and will be outputed at the end for reference

  Exp_Level_Concise <- High_Confidence_Exp_Level[, 2:27]

  # Remove all blank rows in Dataframe

  NA_Removed <- tidyr::drop_na(Exp_Level_Concise)

  # Step 1

  Exp_Level_1 <- NA_Removed[, 3:4]
  Exp_Level_2 <- NA_Removed[, 17:26]
  step1 <- cbind(Exp_Level_1, Exp_Level_2)

  # Step 2

  Normalized_Exp_1 <- (step1$`Abundances (Grouped): 126` * mean(step1$`Abundances (Grouped): 129N`)) / (step1$`Abundances (Grouped): 129N` * mean(step1$`Abundances (Grouped): 126`))
  Normalized_Exp_2 <- (step1$`Abundances (Grouped): 127N` * mean(step1$`Abundances (Grouped): 129C`)) / (step1$`Abundances (Grouped): 129C` * mean(step1$`Abundances (Grouped): 127N`))
  Normalized_Exp_3 <- (step1$`Abundances (Grouped): 127C` * mean(step1$`Abundances (Grouped): 130N`)) / (step1$`Abundances (Grouped): 130N` * mean(step1$`Abundances (Grouped): 127C`))
  Normalized_Exp_4 <- (step1$`Abundances (Grouped): 128N` * mean(step1$`Abundances (Grouped): 130C`)) / (step1$`Abundances (Grouped): 130C` * mean(step1$`Abundances (Grouped): 128N`))
  Normalized_Exp_5 <- (step1$`Abundances (Grouped): 128C` * mean(step1$`Abundances (Grouped): 131`)) / (step1$`Abundances (Grouped): 131` * mean(step1$`Abundances (Grouped): 128C`))

  step2 <-
    cbind(step1,
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

  Exp_Log_Avg <- rowMeans(step3[, 18:22], na.rm = T)
  step4 <- cbind(step3, Exp_Log_Avg)

  # Step 5

  c1 <- step4[, 18:22]
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

  Hit_ident1_Exp <- step9[, 1:2]
  Log_Avg_Exp <- step9$Exp_Log_Avg
  Log_10_Exp <- step9$log_10_Exp
  Hit_ident_Exp <- cbind(Hit_ident1_Exp, Log_Avg_Exp, Log_10_Exp, Log_2_Exp_1, Log_2_Exp_2, Log_2_Exp_3, Log_2_Exp_4, Log_2_Exp_5, z_score_Exp)

  n1 <-
    (Hit_ident_Exp$z_score_Exp > SD_cutoff &
       Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff)) |
    (Hit_ident_Exp$Log_Avg < -z_score_Exp & Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff))
  n2 <- as.data.frame(n1)

  Sig_Check_Exp <-
    dplyr::if_else(n2$n1 == "TRUE", "Significant", "Not Significant")
  step10 <- cbind(Hit_ident_Exp, Sig_Check_Exp)


  # Remove [K] and [R] [+], [-], All bracketed Amino Acids
  # If contains oxidation modification remove it

  SPROX_Confidence <- SPROX_Raw[SPROX_Raw$Confidence == "High",]

  SPROX_No_Oxid <- dplyr::filter(SPROX_Confidence, !grepl('Oxidation', Modifications))

  SPROX_removed_AA1 <- dplyr::mutate_all(SPROX_No_Oxid, ~gsub("\\[[A-Z]].", "", .))
  SPROX_removed_AA2 <- dplyr::mutate_all(SPROX_removed_AA1, ~gsub("\\.\\[[A-Z]]", "", .))
  SPROX_removed_AA3 <- dplyr::mutate_all(SPROX_removed_AA2, ~gsub("\\[[+]].", "", .))
  SPROX_removed_AA4 <- dplyr::mutate_all(SPROX_removed_AA3, ~gsub("\\[[-]].", "", .))
  SPROX_removed_AA5 <- dplyr::mutate_all(SPROX_removed_AA4, ~gsub("\\.\\[[+]]", "", .))
  SPROX_removed_AA6 <- dplyr::mutate_all(SPROX_removed_AA5, ~gsub("\\.\\[[-]]", "", .))

  # If contains M filter into new column

  SPROX_Just_M <- SPROX_removed_AA6[grep("M", SPROX_removed_AA6$`Annotated Sequence`), ]

  # Go through TPP analysis

  SPROX_1 <- SPROX_Just_M[, 10:11]
  SPROX_2 <- SPROX_Just_M[, 3:4]
  SPROX_3 <- SPROX_Just_M[, 15:24]
  step1_SPROX <- cbind(SPROX_1, SPROX_2, SPROX_3)

  # Filter blanks

  NA_Removed_SPROX <- tidyr::drop_na(step1_SPROX)
  NA_Removed_SPROX_Renamed <- magrittr::"%>%" (NA_Removed_SPROX,
    dplyr::rename("Accession" = "Master Protein Accessions"))

  SPROX_merger <- merge(Hit_ident_Exp, NA_Removed_SPROX_Renamed)
  SPROX_merge <-
    dplyr::filter(
      SPROX_merger,
      `Abundances (Grouped): 126` > 0,
      `Abundances (Grouped): 127N` > 0,
      `Abundances (Grouped): 127C` > 0,
      `Abundances (Grouped): 128N` > 0,
      `Abundances (Grouped): 128C` > 0,
      `Abundances (Grouped): 129N` > 0,
      `Abundances (Grouped): 129C` > 0,
      `Abundances (Grouped): 130N` > 0,
      `Abundances (Grouped): 130C` > 0,
      `Abundances (Grouped): 131` > 0
    )

  # Step 11

  Normalized_SPROX_T1 <- (as.numeric(SPROX_merge$`Abundances (Grouped): 126`) * mean(as.numeric(SPROX_merge$`Abundances (Grouped): 129N`))) / (as.numeric(SPROX_merge$`Abundances (Grouped): 129N`) * mean(as.numeric(SPROX_merge$`Abundances (Grouped): 126`)))
  Normalized_SPROX_T2 <- (as.numeric(SPROX_merge$`Abundances (Grouped): 127N`) * mean(as.numeric(SPROX_merge$`Abundances (Grouped): 129C`))) / (as.numeric(SPROX_merge$`Abundances (Grouped): 129C`) * mean(as.numeric(SPROX_merge$`Abundances (Grouped): 127N`)))
  Normalized_SPROX_T3 <- (as.numeric(SPROX_merge$`Abundances (Grouped): 127C`) * mean(as.numeric(SPROX_merge$`Abundances (Grouped): 130N`))) / (as.numeric(SPROX_merge$`Abundances (Grouped): 130N`) * mean(as.numeric(SPROX_merge$`Abundances (Grouped): 127C`)))
  Normalized_SPROX_T4 <- (as.numeric(SPROX_merge$`Abundances (Grouped): 128N`) * mean(as.numeric(SPROX_merge$`Abundances (Grouped): 130C`))) / (as.numeric(SPROX_merge$`Abundances (Grouped): 130C`) * mean(as.numeric(SPROX_merge$`Abundances (Grouped): 128N`)))
  Normalized_SPROX_T5 <- (as.numeric(SPROX_merge$`Abundances (Grouped): 128C`) * mean(as.numeric(SPROX_merge$`Abundances (Grouped): 131`))) / (as.numeric(SPROX_merge$`Abundances (Grouped): 131`) * mean(as.numeric(SPROX_merge$`Abundances (Grouped): 128C`)))

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

  Log_Avg_SPROX <- rowMeans(step12[, 29:33])
  step13 <- cbind(step12, Log_Avg_SPROX)

  # Step  14

  rep1_exp1 <- step13$Log_2_SPROX_T1 - step13$Log_2_Exp_1
  rep2_exp2 <- step13$Log_2_SPROX_T2 - step13$Log_2_Exp_2
  rep3_exp3 <- step13$Log_2_SPROX_T3 - step13$Log_2_Exp_3
  rep4_exp4 <- step13$Log_2_SPROX_T4 - step13$Log_2_Exp_4
  rep5_exp5 <- step13$Log_2_SPROX_T5 - step13$Log_2_Exp_5
  step14 <- cbind(step13, rep1_exp1, rep2_exp2, rep3_exp3, rep4_exp4, rep5_exp5)

  # Step 15

  SPROX_Avg <- rowMeans(step14[, 34:38], na.rm = T)
  step15 <- cbind(step14, SPROX_Avg)

  # Step 16

  k1 <- step15[, 34:38]
  SD_SPROX <- apply(k1, 1, sd)
  step16 <- cbind(step15, SD_SPROX)

  # Step 17

  step17 <- tidyr::drop_na(step16)
  SPROXunique <- dplyr::n_distinct(step17$Accession)

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
  SPROX_p_value <- 2 * pt(hk3, 4, lower.tail = F)
  step20 <- cbind(step19, SPROX_p_value)

  # Step 21

  jk1 <- step20$SPROX_p_value
  SPROX_log_10 <- -log10(jk1)
  step21 <- cbind(step20, SPROX_log_10)

  SPROX_Hit_ident1 <- step21[, 1:2]
  SPROX_Hit_ident2 <- step21[, 12:13]
  SPROX_Log_Avg <- step21[, 40:45]
  SPROX_Hit_ident <- cbind(SPROX_Hit_ident1, SPROX_Hit_ident2, SPROX_Log_Avg)

  # Step 22

  nk1 <-
    (SPROX_Hit_ident$z_score_SPROX > SD_cutoff &
       SPROX_Hit_ident$SPROX_log_10 > -log10(p_cutoff)) |
    (SPROX_Hit_ident$z_score_SPROX < -SD_cutoff & SPROX_Hit_ident$SPROX_log_10 > -log10(p_cutoff))
  nk2 <- as.data.frame(nk1)

  SPROX_Sig_Check <-
    dplyr::if_else(nk2$nk1 == "TRUE", "Significant", "Not Significant")
  step22 <- cbind(SPROX_Hit_ident, SPROX_Sig_Check)

  ddk <- step22$SPROX_Sig_Check

  as.data.frame(ddk)
  abk <- (ddk == "Significant")
  SPROX_Number_Of_Hits <- length(which(abk))

  Hit_List <- step22[ step22$SPROX_Sig_Check == "Significant", ]
  SPROXuniqueHits <- dplyr::n_distinct(Hit_List$Accession)
  SPROXExport <- unique(Hit_List$Accession)
  print(step22)
  print(paste("There are", SPROX_Number_Of_Hits, "Hits"))
  print(paste("There are", SPROXuniqueHits, "Unique Hits"))
  print(paste(SPROXunique, "Proteins were assayed"))
  print(Hit_List)

  utils::write.table(SPROXExport, file = "OnePotSPROXUniqueHits.csv", row.names = FALSE, col.names = FALSE)

  # Visualizing Data

  Volcano_Plot <- ggplot2::ggplot(step22,
                         ggplot2::aes (
                           x = z_score_SPROX ,
                           y = SPROX_log_10 ,
                           label = Accession,
                           colour = SPROX_Sig_Check
                         )) + ggplot2::geom_point(size = 1.5, alpha = 0.5) +
    ggplot2::labs(x="Z Score",y="- log10 (p value)") + ggplot2::scale_colour_manual(breaks = c("Not Significant","Significant"),values = c("grey","red")) +
    ggplot2::scale_x_continuous(breaks = seq(-8,8,by = 1), limits=c(-12, 12)) +
    ggplot2::scale_y_continuous(breaks = seq(0,8,by=1), limits=c(0, 10)) +
    ggplot2::expand_limits(x=0, y=0) +
    ggplot2::theme_bw()+ ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  print(Volcano_Plot)

}

#' STEPP Analysis
#'
#' This function will take expression data and STEPP (TMT10plex) data directly from ProteomeDiscoverer (v. 2.3) and generate a volcano plot and hit list
#' based on user-defined SD and p-values
#'
#' @param Expression_data Raw ProteomeDiscoverer (v. 2.3) data
#' @param TPP_data Raw STEPP (TMT10plex) ProteomeDiscoverer (v. 2.3) data
#' @param SD_cutoff A number
#' @param p_cutoff A number
#' @return Volcano plot, Hit list, Total proteins assayed and Unique protein hits
#' @return Hit list, volcano plot
#' @export

STEPP.Fun <- function(Expression_data, STEPP_Raw, SD_cutoff, p_cutoff){

  High_Confidence_Exp_Level <-
    Expression_data[Expression_data$`Protein FDR Confidence: Combined` == "High",]

  # This dataframe takes all useful columns from the raw data and will be outputed at the end for reference

  Exp_Level_Concise <- High_Confidence_Exp_Level[, 2:27]

  # Remove all blank rows in Dataframe

  NA_Removed <- tidyr::drop_na(Exp_Level_Concise)

  # Step 1

  Exp_Level_1 <- NA_Removed[, 3:4]
  Exp_Level_2 <- NA_Removed[, 17:26]
  step1 <- cbind(Exp_Level_1, Exp_Level_2)

  # Step 2

  Normalized_Exp_1 <- (step1$`Abundances (Grouped): 126` * mean(step1$`Abundances (Grouped): 129N`)) / (step1$`Abundances (Grouped): 129N` * mean(step1$`Abundances (Grouped): 126`))
  Normalized_Exp_2 <- (step1$`Abundances (Grouped): 127N` * mean(step1$`Abundances (Grouped): 129C`)) / (step1$`Abundances (Grouped): 129C` * mean(step1$`Abundances (Grouped): 127N`))
  Normalized_Exp_3 <- (step1$`Abundances (Grouped): 127C` * mean(step1$`Abundances (Grouped): 130N`)) / (step1$`Abundances (Grouped): 130N` * mean(step1$`Abundances (Grouped): 127C`))
  Normalized_Exp_4 <- (step1$`Abundances (Grouped): 128N` * mean(step1$`Abundances (Grouped): 130C`)) / (step1$`Abundances (Grouped): 130C` * mean(step1$`Abundances (Grouped): 128N`))
  Normalized_Exp_5 <- (step1$`Abundances (Grouped): 128C` * mean(step1$`Abundances (Grouped): 131`)) / (step1$`Abundances (Grouped): 131` * mean(step1$`Abundances (Grouped): 128C`))

  step2 <-
    cbind(step1,
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

  Exp_Log_Avg <- rowMeans(step3[, 18:22], na.rm = T)
  step4 <- cbind(step3, Exp_Log_Avg)

  # Step 5

  c1 <- step4[, 18:22]
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

  Hit_ident1_Exp <- step9[, 1:2]
  Log_Avg_Exp <- step9$Exp_Log_Avg
  Log_10_Exp <- step9$log_10_Exp
  Hit_ident_Exp <- cbind(Hit_ident1_Exp, Log_2_Exp_1, Log_2_Exp_2, Log_2_Exp_3, Log_2_Exp_4, Log_2_Exp_5, z_score_Exp, Log_10_Exp, Log_Avg_Exp)

  n1 <-
    (Hit_ident_Exp$z_score_Exp > SD_cutoff &
       Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff)) |
    (Hit_ident_Exp$Log_Avg < -z_score_Exp & Hit_ident_Exp$Log_10_Exp > -log10(p_cutoff))
  n2 <- as.data.frame(n1)

  Sig_Check_Exp <-
    dplyr::if_else(n2$n1 == "TRUE", "Significant", "Not Significant")
  step10 <- cbind(Hit_ident_Exp, Sig_Check_Exp)


  # Only get Semi-Tryptic Peptides


  STEPP_Confidence <- STEPP_Raw[STEPP_Raw$Confidence == "High",]

  STEPP_Semi_Tryp <- dplyr::filter(STEPP_Confidence, Modifications == "1xTMT6plex [N-Term]")

  STEPP_removed_AA1 <- dplyr::mutate_all(STEPP_Semi_Tryp, ~gsub("\\[[A-Z]].", "", .))
  STEPP_removed_AA2 <- dplyr::mutate_all(STEPP_removed_AA1, ~gsub("\\.\\[[A-Z]]", "", .))
  STEPP_removed_AA3 <- dplyr::mutate_all(STEPP_removed_AA2, ~gsub("\\[[+]].", "", .))
  STEPP_removed_AA4 <- dplyr::mutate_all(STEPP_removed_AA3, ~gsub("\\[[-]].", "", .))
  STEPP_removed_AA5 <- dplyr::mutate_all(STEPP_removed_AA4, ~gsub("\\.\\[[+]]", "", .))
  STEPP_removed_AA6 <- dplyr::mutate_all(STEPP_removed_AA5, ~gsub("\\.\\[[-]]", "", .))

  # Go through TPP analysis

  STEPP_1 <- STEPP_removed_AA6[, 10]
  STEPP_2 <- STEPP_removed_AA6[, 3:4]
  STEPP_3 <- STEPP_removed_AA6[, 15:24]
  step1_STEPP <- cbind(STEPP_1, STEPP_2, STEPP_3)

  # Filter blanks

  NA_Removed_STEPP <- tidyr::drop_na(step1_STEPP)
  NA_Removed_STEPP_Renamed <- magrittr::"%>%" (NA_Removed_STEPP,
    dplyr::rename("Accession" = "Master Protein Accessions"))

  STEPP_merger <- merge(Hit_ident_Exp, NA_Removed_STEPP_Renamed)
  STEPP_merge <-
    dplyr::filter(
      STEPP_merger,
      `Abundances (Grouped): 126` > 0,
      `Abundances (Grouped): 127N` > 0,
      `Abundances (Grouped): 127C` > 0,
      `Abundances (Grouped): 128N` > 0,
      `Abundances (Grouped): 128C` > 0,
      `Abundances (Grouped): 129N` > 0,
      `Abundances (Grouped): 129C` > 0,
      `Abundances (Grouped): 130N` > 0,
      `Abundances (Grouped): 130C` > 0,
      `Abundances (Grouped): 131` > 0
    )

  # Step 11

  Normalized_STEPP_T1 <- (as.numeric(STEPP_merge$`Abundances (Grouped): 126`) * mean(as.numeric(STEPP_merge$`Abundances (Grouped): 129N`))) / (as.numeric(STEPP_merge$`Abundances (Grouped): 129N`) * mean(as.numeric(STEPP_merge$`Abundances (Grouped): 126`)))
  Normalized_STEPP_T2 <- (as.numeric(STEPP_merge$`Abundances (Grouped): 127N`) * mean(as.numeric(STEPP_merge$`Abundances (Grouped): 129C`))) / (as.numeric(STEPP_merge$`Abundances (Grouped): 129C`) * mean(as.numeric(STEPP_merge$`Abundances (Grouped): 127N`)))
  Normalized_STEPP_T3 <- (as.numeric(STEPP_merge$`Abundances (Grouped): 127C`) * mean(as.numeric(STEPP_merge$`Abundances (Grouped): 130N`))) / (as.numeric(STEPP_merge$`Abundances (Grouped): 130N`) * mean(as.numeric(STEPP_merge$`Abundances (Grouped): 127C`)))
  Normalized_STEPP_T4 <- (as.numeric(STEPP_merge$`Abundances (Grouped): 128N`) * mean(as.numeric(STEPP_merge$`Abundances (Grouped): 130C`))) / (as.numeric(STEPP_merge$`Abundances (Grouped): 130C`) * mean(as.numeric(STEPP_merge$`Abundances (Grouped): 128N`)))
  Normalized_STEPP_T5 <- (as.numeric(STEPP_merge$`Abundances (Grouped): 128C`) * mean(as.numeric(STEPP_merge$`Abundances (Grouped): 131`))) / (as.numeric(STEPP_merge$`Abundances (Grouped): 131`) * mean(as.numeric(STEPP_merge$`Abundances (Grouped): 128C`)))

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

  Log_Avg_STEPP <- rowMeans(step12[, 28:32])
  step13 <- cbind(step12, Log_Avg_STEPP)

  # Step  14

  rep1_exp1 <- step13$Log_2_STEPP_T1 - step13$Log_2_Exp_1
  rep2_exp2 <- step13$Log_2_STEPP_T2 - step13$Log_2_Exp_2
  rep3_exp3 <- step13$Log_2_STEPP_T3 - step13$Log_2_Exp_3
  rep4_exp4 <- step13$Log_2_STEPP_T4 - step13$Log_2_Exp_4
  rep5_exp5 <- step13$Log_2_STEPP_T5 - step13$Log_2_Exp_5
  step14 <- cbind(step13, rep1_exp1, rep2_exp2, rep3_exp3, rep4_exp4, rep5_exp5)

  # Step 15

  STEPP_Avg <- rowMeans(step14[, 33:37], na.rm = T)
  step15 <- cbind(step14, STEPP_Avg)

  # Step 16

  k1 <- step15[, 33:37]
  SD_STEPP <- apply(k1, 1, sd)
  step16 <- cbind(step15, SD_STEPP)

  # Step 17

  step17 <- tidyr::drop_na(step16)
  STEPPunique <- dplyr::n_distinct(step17$Accession)

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
  STEPP_p_value <- 2 * pt(hk3, 4, lower.tail = F)
  step20 <- cbind(step19, STEPP_p_value)

  # Step 21

  jk1 <- step20$STEPP_p_value
  STEPP_log_10 <- -log10(jk1)
  step21 <- cbind(step20, STEPP_log_10)

  STEPP_Hit_ident1 <- step21[, 1:2]
  STEPP_Hit_ident2 <- step21[, 11:12]
  STEPP_Log_Avg <- step21[, 39:44]
  STEPP_Hit_ident <- cbind(STEPP_Hit_ident1, STEPP_Hit_ident2, STEPP_Log_Avg)

  # Step 22

  nk1 <-
    (STEPP_Hit_ident$z_score_STEPP > SD_cutoff &
       STEPP_Hit_ident$STEPP_log_10 > -log10(p_cutoff)) |
    (STEPP_Hit_ident$z_score_STEPP < -SD_cutoff & STEPP_Hit_ident$STEPP_log_10 > -log10(p_cutoff))
  nk2 <- as.data.frame(nk1)

  STEPP_Sig_Check <-
    dplyr::if_else(nk2$nk1 == "TRUE", "Significant", "Not Significant")
  step22 <- cbind(STEPP_Hit_ident, STEPP_Sig_Check)

  ddk <- step22$STEPP_Sig_Check

  as.data.frame(ddk)
  abk <- (ddk == "Significant")
  STEPP_Number_Of_Hits <- length(which(abk))

  Hit_List <- step22[ step22$STEPP_Sig_Check == "Significant", ]
  STEPPuniquehits <- dplyr::n_distinct(Hit_List$Accession)
  STEPPExport <- unique(Hit_List$Accession)
  print(step22)
  print(paste("There are", STEPP_Number_Of_Hits, "Hits"))
  print(paste("There are", STEPPuniquehits, "Unique Hits"))
  print(paste(STEPPunique, "Proteins were assayed"))
  print(Hit_List)

  utils::write.table(STEPPExport, file = "STEPPUniqueHits.csv", row.names = FALSE, col.names = FALSE)

  # Visualizing Data

  Volcano_Plot <- ggplot2::ggplot(step22,
                         ggplot2::aes (
                           x = z_score_STEPP ,
                           y = STEPP_log_10 ,
                           label = Accession,
                           colour = STEPP_Sig_Check
                         )) + ggplot2::geom_point(size = 1.5, alpha = 0.5) +
    ggplot2::labs(x="Z Score",y="- log10 (p value)") + ggplot2::scale_colour_manual(breaks = c("Not Significant","Significant"),values = c("grey","red")) +
    ggplot2::scale_x_continuous(breaks = seq(-8,8,by = 1), limits=c(-8, 8)) +
    ggplot2::scale_y_continuous(breaks = seq(0,6,by=1), limits=c(0, 6)) +
    ggplot2::expand_limits(x=0, y=0) +
    ggplot2::theme_bw()+ ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  print(Volcano_Plot)

}


#' Two Method Comparison Venn Diagram
#'
#' This function will make a Venn Diagram showing unique hits for two different proteomics methods
#'
#'
#' @param data1 unqiue hit csv output method 1
#' @param data2 unique hit csv output method 2
#' @param name1 method 1 name
#' @param name2 method 2 name
#' @return Venn Diagram
#' @export


Two_Comparison_Venn_Diagram.Fun <- function(data1, data2, name1, name2){

  z <- utils::read.csv2(file = data1, header = FALSE)
  z1 <- utils::read.csv2(file = data2, header = FALSE)


  VennDiagram::venn.diagram(x = c(z, z1),
               category.names = c(name1, name2),
               filename = 'Two_Comparison_Venn_Diagram.PNG',
               output = TRUE ,
               imagetype="png" ,
               height = 480 ,
               width = 480 ,
               resolution = 300,
               compression = "lzw",
               lwd = 2,
               col= "black",
               lty = "dotted",
               fill = c("cornflowerblue", "darkorchid1"),
               alpha = c(0.3, 0.3),
               cex = 0.65,
               fontfamily = "serif",
               cat.cex = 0.5,
               cat.default.pos = "outer",
               fontface = "bold",
               cat.pos = c(-1, 1),
               cat.dist = c(0.055, 0.055),
               cat.fontfamily = "serif",
               cat.col = c("darkblue", "darkorchid4"),
               scaled = FALSE

  )




}



























