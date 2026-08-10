setwd("/data/chamal/projects/stephanie/SIR-modelling-project/analysis/analysis-thesis-202404/")

#load data
demo_data <- read.csv('demographics-hippocampus.csv')
Jac_data <- read.csv('hipp_averages_new.csv')

# Merge dataframes based on the ID column
data <- merge(demo_data, Jac_data, by = "SubjectID")

####set ref level####
# Convert 'data$Treatment' into a factor
data$Treatment <- as.factor(data$Treatment)
# Check the levels of 'data$Treatment' before releveling
levels(data$Treatment)
# Relevel 'data$Treatment' with "PBS" as the reference level
data$Treatment <- relevel(data$Treatment, ref = "PBS")

# Convert 'data$Sex' into a factor
data$Sex <- as.factor(data$Sex)
# Check the levels of 'data$Sex' before releveling
levels(data$Sex)
# Relevel 'data$Sex' with "PBS" as the reference level
data$Sex <- relevel(data$Sex, ref = "M")

# # Convert 'data$genotype' into a factor
# data$genotype <- as.factor(data$genotype)
# # Check the levels of 'data$genotype' before releveling
# levels(data$genotype)
# # Relevel 'data$genotype' with "PBS" as the reference level
# data$genotype <- relevel(data$genotype, ref = "WT")


#setting up
data$Weight[data$Weight=="??"] <- NA
data$Weight[data$Weight==""] <- NA
# data$litter_size[data$litter_size=="4/8"] <- 8
# data$animal_ID = as.factor(data$animal_ID)
data$Cohort = as.factor(data$Cohort)
# data$batch = as.factor(data$batch)
data$Weight=as.numeric(as.character(data$Weight))
data$dpi=as.numeric(as.character(data$dpi))
# data$litter_size=as.numeric(as.character(data$litter_size))


###run a linear model at each column as of column 42 to last column and save the t-values####
# Create an empty dataframe to store the results
results <- data.frame(matrix(NA, nrow = ncol(Jac_data) - 1, ncol = 4))
colnames(results) <- paste("Time", 1:4)

region_data <- read.csv("merged_label_region_names.csv") ##to be new row name for matching column number in data eg: X23 

for (i in 22:ncol(data)) {
  # Create a vector to store t-values for each Timepoint
  t_values <- c()
  
  # Iterate over each Timepoints (-1,1,2,3)
  for (j in c(-1, 1:3)) {
    
    subset_data_PBS <- subset(data, Timepoint == j &
                                Treatment == "PBS")
    subset_data_PFF <- subset(data, Timepoint == j &
                                Treatment == "PFF")
    if(i == 22) {
      print(paste0("Timepoint: ",j))
      print(paste0("Timepoint: ",j))
      print(paste0("Number of M PBS mice: ", length(which(subset_data_PBS$Sex == "M"))))
      print(paste0("Number of F PBS mice: ", length(which(subset_data_PBS$Sex == "F"))))
      
      print(paste0("Number of M Hu-PFF mice: ", length(which(subset_data_PFF$Sex == "M"))))
      print(paste0("Number of F Hu-PFF mice: ", length(which(subset_data_PFF$Sex == "F"))))
      
    }
    
    # Check if any data points match the conditions
    if (nrow(subset_data_PBS) > 0 && nrow(subset_data_PFF) > 0 ) {
      t_value <- t.test(subset_data_PBS[i],subset_data_PFF[i])$statistic
      # Append the t-value to the vector
      t_values <- c(t_values, t_value)
    } else {
      # Append NA if no data points match the conditions
      t_values <- c(t_values, NA)
    }
  }
  
  # Assign the t-values vector to the corresponding column in the results dataframe
  results[i - 21, ] <- t_values
  
  # Get the column name from data without the "X" in front
  col_name <- substring(names(data)[i], 2)
  
  # Find the corresponding "Structure" in region_data based on label_number
  structure_name <- region_data$Structure[region_data$label_number == col_name]
  
  # Determine if it's one of the first 211 or the last 211
  if (i < 211+22) {
    row_name <- paste0("left ", structure_name)
  } else {
    row_name <- paste0("right ", structure_name)
  }
  
  # Assign the row name in the results dataframe
  rownames(results)[i - 21] <- row_name
}

write.csv(results, '/data/chamal/projects/natvik/sir_extended/analysis/hipp_rgn_t_stats_hemiHuPff_hemiPBS.csv',row.names=TRUE)


#PLOT
# row_index=which(rownames(results)=="right Caudoputamen")
# plot(1:4, results[row_index,], type = "l", xlab = "Time Steps", ylab = "pSyn")
# plot(1:4, results[which(rownames(results)=="right Substantia nigra compact part"),], type = "l", xlab = "Time Steps", ylab = "pSyn")
# plot(1:4, results[which(rownames(results)=="left Globus pallidus internal segment"),], type = "l", xlab = "Time Steps", ylab = "pSyn")

###load the simulated atrophy csv####
#################################sim_atrophy_data <- read.csv('/Users/stephanietullo/Documents/VSC/SIR_mouse-main/model/model/abm_spread_v.730.0.spread_rate.0.07.dt.0.1.seed.40.injection_amount.2.5e-06.clearance_gene.None.k1.0.9.k2.0.9.csv', header = FALSE)
# targets_data <- read.csv('/Users/stephanietullo/Documents/VSC/SIR_mouse-main/model/model/targets.csv', header = FALSE)
ABA_region_names <- read.csv('ABAnewatlas_final_withacros.csv', header = TRUE)
ABA_region_names<-rbind(ABA_region_names,ABA_region_names)

sim_atrophy_data <- cbind(ABA_region_names, sim_atrophy_data) 

###add hemisphere for the simulated atrophy####
alternating_values <- rep(c("right", "left"), each = 213)
sim_atrophy_data$hemisphere <- alternating_values

# Create unique row names based on first and last columns 
rownames(sim_atrophy_data) <- paste(sim_atrophy_data[, ncol(sim_atrophy_data)], sim_atrophy_data[, 1], sep = " ")

# delete first 5 columns and last column
sim_atrophy_data <- sim_atrophy_data[, -(1:5)]
sim_atrophy_data <- sim_atrophy_data[, -ncol(sim_atrophy_data)]

#PLOT
library(ggplot2)
time_steps <- 1:ncol(sim_atrophy_data)  # Assuming the time steps are represented by column indices

# plot(time_steps, sim_atrophy_data["right Caudoputamen",], type = "l", xlab = "Time Steps", ylab = "Atrophy")
# plot(1:4, results["right Caudoputamen",], type = "l", xlab = "Time Steps", ylab = "Atrophy")


#####connectivity data
connectivity <- read.csv('CS-atlas-edited.csv', header = TRUE) ##from Oh et al suppl table 3 right injection sites by contralteral (left) followed by ipsilateral (right) regions
rownames(connectivity) <- connectivity$X   # Set the first column as rownames
connectivity <- connectivity[, -1]  # Remove the original 'Name' column (optional)

col_index_to_delete <- which(colnames(connectivity) == "SUBd") # Find the column index corresponding to "SUBd"
connectivity <- connectivity[, -col_index_to_delete, drop = FALSE] # Delete the column

col_index_to_delete <- which(colnames(connectivity) == "SUBd.1") # Find the column index corresponding to "SUBd"
connectivity <- connectivity[, -col_index_to_delete, drop = FALSE] # Delete the column

col_index_to_delete <- which(colnames(connectivity) == "SUBv") # Find the column index corresponding to "SUBd"
connectivity <- connectivity[, -col_index_to_delete, drop = FALSE] # Delete the column

col_index_to_delete <- which(colnames(connectivity) == "SUBv.1") # Find the column index corresponding to "SUBd"
connectivity <- connectivity[, -col_index_to_delete, drop = FALSE] # Delete the column

seed_connectivity <- connectivity["DG",, drop = FALSE]
seed_connectivity<-t(seed_connectivity[1, ])

# plot(density(as.numeric(seed_connectivity[,1])), main = "Density Plot of Seed Connectivity", xlab = "Seed Connectivity")
# boxplot( as.numeric(seed_connectivity[,1]), main = "Boxplot of Seed Connectivity", ylab = "Seed Connectivity")


###########based on asyn gene expression########
aSyn <- read.csv('asyn-gene-expression.csv', header = TRUE)
aSyn <- rbind(aSyn, aSyn) #double for each hemispher
aSyn$hemisphere <- alternating_values #add hemisphere for the simulated atrophy
rownames(aSyn) <- paste(aSyn[, ncol(aSyn)], aSyn[, 1], sep = " ") # Create unique row names based on first and last columns 
aSyn <- aSyn[, -1] # delete first columns and last column
aSyn <- aSyn[, -ncol(aSyn)]

row_index_to_delete <- which(rownames(aSyn) == "right Subiculum dorsal") # Find the row index corresponding to "SUBd"
aSyn <- aSyn[-row_index_to_delete,, drop = FALSE] # Delete the row

row_index_to_delete <- which(rownames(aSyn) == "left Subiculum dorsal") # Find the row index corresponding to "SUBd"
aSyn <- aSyn[-row_index_to_delete,, drop = FALSE] # Delete the row

row_index_to_delete <- which(rownames(aSyn) == "right Subiculum ventral") # Find the row index corresponding to "SUBd"
aSyn <- aSyn[-row_index_to_delete,, drop = FALSE] # Delete the row

row_index_to_delete <- which(rownames(aSyn) == "left Subiculum ventral") # Find the row index corresponding to "SUBd"
aSyn <- aSyn[-row_index_to_delete,, drop = FALSE] # Delete the row


#########merge aSyn and seed_connectivty datasets#########
expression<- cbind(aSyn,seed_connectivity)

#####Get the common row names between empirical and simulated data#######
common_row_names <- intersect(rownames(results), rownames(sim_atrophy_data))
empirical_common <- results[common_row_names, ]
simulated_common <- sim_atrophy_data[common_row_names, ]

#########merge expression, simualted and empirical datasets#########
# Merge datasets by row names
merged_data <- merge(expression, empirical_common, by = "row.names")

# Reassign column 'Row.names' to be the row names of the merged dataset
rownames(merged_data) <- merged_data$Row.names
merged_data$Row.names <- NULL  # Remove the redundant column

merged_data <- merge(merged_data, simulated_common, by = "row.names")

# Reassign column 'Row.names' to be the row names of the merged dataset
rownames(merged_data) <- merged_data$Row.names
merged_data$Row.names <- NULL  # Remove the redundant column


#############################
# Subset merged_data to keep only Time 1 column
merged_data_T1 <- merged_data[, !grepl("Time 2|Time 3|Time 4", names(merged_data))]
# Subset merged_data to keep only Time 2 column
merged_data_T2 <- merged_data[, !grepl("Time 1|Time 3|Time 4", names(merged_data))]
# Subset merged_data to keep only Time 3 column
merged_data_T3 <- merged_data[, !grepl("Time 2|Time 1|Time 4", names(merged_data))]
# Subset merged_data to keep only Time 4 column
merged_data_T4 <- merged_data[, !grepl("Time 2|Time 3|Time 1", names(merged_data))]


##########Get rid of regions with volume increases (ie PFF > PBS)###################
merged_data_T1<-merged_data_T1[merged_data_T1$`Time 1`>=0,,drop=FALSE]
merged_data_T2<-merged_data_T2[merged_data_T2$`Time 2`>=0,,drop=FALSE]
merged_data_T3<-merged_data_T3[merged_data_T3$`Time 3`>=0,,drop=FALSE]
merged_data_T4<-merged_data_T4[merged_data_T4$`Time 4`>=0,,drop=FALSE]


####Plot correlations #####

# Initialize a vector to store the correlations
correlations_tp1 <- vector("numeric", length = length(time_steps))
correlations_tp2 <- vector("numeric", length = length(time_steps))
correlations_tp3 <- vector("numeric", length = length(time_steps))
correlations_tp4 <- vector("numeric", length = length(time_steps))

# Calculate Spearman correlations
for (i in 1:length(time_steps)) {
  correlations_tp1[i] <- cor(merged_data_T1$`Time 1`, merged_data_T1[[paste0("V", i)]], method = "spearman")
}
for (i in 1:length(time_steps)) {
  correlations_tp2[i] <- cor(merged_data_T2$`Time 2`, merged_data_T2[[paste0("V", i)]], method = "spearman")
}
for (i in 1:length(time_steps)) {
  correlations_tp3[i] <- cor(merged_data_T3$`Time 3`, merged_data_T3[[paste0("V", i)]], method = "spearman")
}
for (i in 1:length(time_steps)) {
  correlations_tp4[i] <- cor(merged_data_T4$`Time 4`, merged_data_T4[[paste0("V", i)]], method = "spearman")
}


# Plot the correlations against the time steps
# Create a dataframe for plotting
correlation_data_t1 <- data.frame(Time_Steps = time_steps, Correlations = correlations_tp1)
correlation_data_t2 <- data.frame(Time_Steps = time_steps, Correlations = correlations_tp2)
correlation_data_t3 <- data.frame(Time_Steps = time_steps, Correlations = correlations_tp3)
correlation_data_t4 <- data.frame(Time_Steps = time_steps, Correlations = correlations_tp4)


# # Find the time step with peak correlation 
print(paste("Time step with peak correlation (tp1):", 
            correlation_data_t1$Time_Steps[which.max(correlation_data_t1$Correlations)],
            "max correlation:", max(correlation_data_t1$Correlations)))
print(paste("Time step with peak correlation (tp2):", 
            correlation_data_t2$Time_Steps[which.max(correlation_data_t2$Correlations)],
            "max correlation:", max(correlation_data_t2$Correlations)))
print(paste("Time step with peak correlation (tp3):", 
            correlation_data_t3$Time_Steps[which.max(correlation_data_t3$Correlations)],
            "max correlation:", max(correlation_data_t3$Correlations)))
print(paste("Time step with peak correlation (tp4):", 
            correlation_data_t4$Time_Steps[which.max(correlation_data_t4$Correlations)],
            "max correlation:", max(correlation_data_t4$Correlations)))


#######Peak correlation and time step########
# Combine all correlation dataframes into a single dataframe
combined_correlation_data <- rbind(correlation_data_t1, correlation_data_t2, correlation_data_t3, correlation_data_t4)

# Find the row with the maximum correlation
max_correlation_row <- combined_correlation_data[which.max(combined_correlation_data$Correlations), ]

# Extract the maximum correlation and its corresponding time step
peak_corr <- max_correlation_row$Correlations
peak_timestep <- max_correlation_row$Time_Steps


#############
library(ggplot2)

# Create the ggplot
ggplot() +
  geom_line(data = correlation_data_t1, aes(x = Time_Steps, y = Correlations, color = "-7 dpi"), linetype = "solid") +
  geom_line(data = correlation_data_t2, aes(x = Time_Steps, y = Correlations, color = "30 dpi"), linetype = "solid") +
  geom_line(data = correlation_data_t3, aes(x = Time_Steps, y = Correlations, color = "90 dpi"), linetype = "solid") +
  geom_line(data = correlation_data_t4, aes(x = Time_Steps, y = Correlations, color = "120 dpi"), linetype = "solid") +
  labs(x = "Simulated Time Steps", y = "Spearman Correlation", title = "Model Fit: Empirical vs Simulated Atrophy") +
  scale_color_manual(
    values = c("#E69F00", "#009E73", "#0072B2", "#CC79A7"),
    breaks = c("-7 dpi", "30 dpi", "90 dpi", "120 dpi"),
    guide = guide_legend(title = "")
  ) +
  geom_point(
    aes(x = peak_timestep, y = peak_corr),
    color = "#D55E00",
    size = 2
  ) +
  annotate(
    "text",
    x = peak_timestep,
    y = peak_corr + 0.02,  # Adjust the value to position the text above the point
    label = paste("Peak Time Step:", peak_timestep, "\nPeak Correlation:", round(peak_corr, 4)),
    hjust = 0,
    vjust = -0.2,  # Make text appear above the point
    color = "#D55E00",
    size = 6
  ) +
  theme_minimal() +
  ylim(-0.2, 0.4) +
  theme(legend.position = "right", text = element_text(size = 28), plot.title = element_text(hjust = 0.25))


###Plot peak correlations
ggplot(data = merged_data_T4, aes(x = merged_data_T4[[paste0("V", peak_timestep)]], y = merged_data_T4$`Time 4`)) +
  geom_point(aes(colour=rownames(merged_data_T4)),size=3) +
  geom_smooth(method="lm",formula=y~x)+
  labs(x = "Simulated Atrophy", y = "Empirical Atrophy", title = "Simulated vs. Empirical Atrophy at Peak Correlation Time Step") +
  theme_minimal() +   theme(legend.position = "none",
                            text = element_text(size = 14), plot.title=element_text(size = 13))

####plot sim vs empr - with color coded regions of interest###############
# Create the ggplot with the highlighted rows
library(dplyr) # Load the dplyr package for data manipulation

# Specify the rownames you want to highlight
highlighted_rownames <- c("right Field CA2", "left Field CA2", 
                          "right Field CA1", "left Field CA1",
                          "right Field CA3", "left Field CA3",
                          "left Dentate gyrus", "right Dentate gyrus",
                          "right Parasubiculum","left Parasubiculum",
                          "right Presubiculum","left Presubiculum",
                          "right Postsubiculum","left Postsubiculum",
                          "right Subiculum","left Subiculum",
                          "left Entorhinal area", "right Entorhinal area")

highlighted_data <- merged_data_T4 %>%
  mutate(
    highlight = ifelse(rownames(merged_data_T4) %in% highlighted_rownames, "highlight", "other")
  )

highlighted_data <- highlighted_data %>%
  rename(empirical = `Time 4`, 
         timestep = !!sym(paste0("V", peak_timestep)))

ggplot(data = highlighted_data, aes(x = timestep, y = empirical)) +
  geom_point(aes(colour = highlight), size = 3) +
  scale_color_manual(values = c("highlight" = "#33CC00", "other" = "grey50")) +
  geom_text(data = data.frame( timestep = highlighted_data$timestep[highlighted_data$highlight == "highlight"],
                               empirical = highlighted_data$empirical[highlighted_data$highlight == "highlight"],
                              label = rownames(highlighted_data)[highlighted_data$highlight == "highlight"]),
            aes(label = label), hjust = 0, vjust = 1, color = "black", size = 4) +
  geom_smooth(method="lm",formula=y~x, colour="#0072B2")+
  labs(x = "Simulated Atrophy", y = "Empirical Atrophy"
       # , title = "Simulated vs. Empirical Atrophy at Peak Correlation Time Step"
  ) +
  theme_minimal() +   
  theme(legend.position = "none",
        text = element_text(size = 25), plot.title=element_text(size = 18))


####plot sim vs empr - with color coded based on connecitivty###############
ggplot(data = highlighted_data, aes(x = timestep, y = empirical, color = DG)) +
    geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  labs(x = "Simulated Atrophy", y = "Empirical Atrophy"
       # , title = "Simulated vs. Empirical Atrophy at Peak Correlation Time Step"
  ) +
  scale_color_gradientn(
    colors = c("grey50", "#56B4E9"),
    values = scales::rescale(c(0, min(highlighted_data$DG[highlighted_data$DG != 0]), max(highlighted_data$DG)), c(0, 1)),
    limits = range(highlighted_data$DG),
    breaks = c(min(highlighted_data$DG[highlighted_data$DG != 0]), max(highlighted_data$DG)) , # Define the desired breaks or levels for the color scale
    name = "Seed Connectivity      "  # Set the legend title here
  ) +
  theme_minimal() +
  theme(
    # legend.position = "none",
    text = element_text(size = 25),
    plot.title = element_text(size = 13)
  )


####plot sim vs empr - with color coded based on aSyn expression###############

ggplot(data = highlighted_data, aes(x = timestep, y = empirical, color = asyn)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  labs(x = "Simulated Atrophy", y = "Empirical Atrophy"
       # , title = "Simulated vs. Empirical Atrophy at Peak Correlation Time Step"
  ) +
  # scale_color_gradientn(
  #   colors = colorRampPalette(c("grey50", "#CC79A7"))(n = 422),  # Adjust the number of colors as needed
  #   limits = range(aSyn_common$asyn),
  #   breaks = custom_breaks,
  #   name = "aSyn gene expression"  # Set the legend title here
  # ) +
  scale_color_gradientn(
    colors = c("grey50", "#CC79A7"),
    values = scales::rescale(c(0, min(highlighted_data$asyn[highlighted_data$asyn != 0]), max(highlighted_data$asyn)), c(0, 1)),
    limits = range(highlighted_data$asyn),
    breaks = c(min(highlighted_data$asyn[highlighted_data$asyn != 0]), max(highlighted_data$asyn)) , # Define the desired breaks or levels for the color scale
    
    name = "aSyn gene expression"  # Set the legend title here
  ) +
  theme_minimal() +
  theme(
    # legend.position = "none",
    text = element_text(size = 25),
    plot.title = element_text(size = 13)
  )


###################check correlation between asyn gene expression and simulated atrophy##########

###sim atrophy vs SNCA expression####
correlations_sim_SNCA <- cor(highlighted_data$asyn, highlighted_data$timestep, method = "spearman")
print(correlations_sim_SNCA)

# Create data frame for plotting
plot_data <- data.frame(
  x = highlighted_data$asyn,
  y = highlighted_data$timestep,
  region = rownames(highlighted_data)
)

# Scatter plot with ggplot
ggplot(data = plot_data, aes(x = x, y = y, color = region)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  geom_text(x = min(plot_data$x), y = max(plot_data$y), 
            label = paste0("Spearman correlation = ", round(correlations_sim_SNCA, 3)), hjust = 0, vjust = 1,color="black",size=6) +
  labs(x = "SNCA expression", y = "Simulated Atrophy", title = "") +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 18),
    plot.title = element_text(size = 13)
  )


###sim atrophy vs connectivity####
correlations_sim_cxn <- cor(highlighted_data$DG, highlighted_data$timestep, method = "spearman")
print(correlations_sim_cxn)

# Create data frame for plotting
plot_data <- data.frame(
  x = highlighted_data$DG,
  y = highlighted_data$timestep,
  region = rownames(highlighted_data)
)

# Scatter plot with ggplot
ggplot(data = plot_data, aes(x = x, y = y, color = region)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  geom_text(x = min(plot_data$x), y = max(highlighted_data$timestep), 
            label = paste0("Spearman correlation = ", round(correlations_sim_cxn, 3)), hjust = 0, vjust = 1,color="black",size=6) +
  labs(x = "Connectivity Strength", y = "Simulated Atrophy", title = "") +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 18),
    plot.title = element_text(size = 13)
  )


####### with labelled top hits#######
ggplot(data = plot_data, aes(x = x, y = y, color = region)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  geom_text(x = min(plot_data$x), y = max(highlighted_data$timestep), 
            label = paste0("Spearman correlation = ", round(correlations_sim_cxn, 3)), hjust = 0, vjust = 1,color="black",size=6) +
  geom_text(data = subset(plot_data, x>0.05), 
            aes(x = x, y = y, label = region),
            hjust = 0, vjust = 1, color = "black", size = 3) +
  labs(x = "Connectivity Strength", y = "Simulated Atrophy", title = "") +
  xlim(0,2) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 18),
    plot.title = element_text(size = 13)
  )


##########check correlation between asyn gene expression and empirical atrophy at 90 dpi##########

###empirical atrophy versus SNCA####
correlations_emp_SNCA <- cor(highlighted_data$asyn, highlighted_data$empirical, method = "spearman")
print(correlations_emp_SNCA)

# Create data frame for plotting
plot_data <- data.frame(
  x = highlighted_data$asyn,
  y = highlighted_data$empirical,
  region = rownames(highlighted_data)
)

# Scatter plot with ggplot
ggplot(data = plot_data, aes(x = x, y = y, color = region)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  geom_text(x = min(plot_data$x), y = max(plot_data$y), 
            label = paste0("Spearman correlation = ", round(correlations_emp_SNCA, 3)), hjust = 0, vjust = 1,color="black",size=6) +
  labs(x = "SNCA expression", y = "Empirical Atrophy", title = "") +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 18),
    plot.title = element_text(size = 13)
  )


###empirical atrophy versus connectivity####

correlations_emp_cxn <- cor(highlighted_data$DG, highlighted_data$empirical, method = "spearman")
print(correlations_emp_cxn)

# Create data frame for plotting
plot_data <- data.frame(
  x = highlighted_data$DG,
  y = highlighted_data$empirical,
  region = rownames(highlighted_data)
)

# Scatter plot with ggplot
ggplot(data = plot_data, aes(x = x, y = y, color = region)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  geom_text(x = min(plot_data$x), y = max(plot_data$y), 
            label = paste0("Spearman correlation = ", round(correlations_emp_cxn, 3)), hjust = 0, vjust = 1,color="black",size=6) +
  labs(x = "Connectivity Strength", y = "Empirical Atrophy", title = "") +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 18),
    plot.title = element_text(size = 13)
  )


####### with labelled top hits#######
ggplot(data = plot_data, aes(x = x, y = y, color = region)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  geom_text(x = 0.5, y = max(plot_data$y), 
            label = paste0("Spearman correlation = ", round(correlations_emp_cxn, 3)), hjust = 0, vjust = 1,color="black",size=6) +
  geom_text(data = subset(plot_data, y > 3), 
            aes(x = x, y = y, label = region),
            hjust = 0, vjust = 1, color = "black", size = 3) +
  labs(x = "Connectivity Strength", y = "Empirical Atrophy", title = "") +
  xlim(0,2)+
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 18),
    plot.title = element_text(size = 13)
  )
