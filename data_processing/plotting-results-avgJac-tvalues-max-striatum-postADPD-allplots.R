###load data
setwd('/data/chamal/projects/natvik/sir_extended/analysis/')
demo_data <- read.csv('/data/chamal/projects/natvik/sir_extended/preprocessed/steph_inputs/tp-data-DBM-analysis202202.csv')
Jac_data <- read.csv('/data/chamal/projects/natvik/sir_extended/preprocessed/steph_inputs/averages_new.csv')

#############TOGGLE BETWEEN Hu-PFF and Ms-PFF in this section of the script to produce the output files##################

# Merge dataframes based on the ID column
data <- merge(demo_data, Jac_data, by = "SubjectID")

####set ref level####s
# Convert 'data$injection' into a factor
data$injection <- as.factor(data$injection)
# Check the levels of 'data$injection' before releveling
levels(data$injection)
# Relevel 'data$injection' with "PBS" as the reference level
data$injection <- relevel(data$injection, ref = "PBS")

# Convert 'data$sex' into a factor
data$sex <- as.factor(data$sex)
# Check the levels of 'data$sex' before releveling
levels(data$sex)
# Relevel 'data$sex' with "PBS" as the reference level
data$sex <- relevel(data$sex, ref = "M")

# Convert 'data$genotype' into a factor
data$genotype <- as.factor(data$genotype)
# Check the levels of 'data$genotype' before releveling
levels(data$genotype)
# Relevel 'data$genotype' with "PBS" as the reference level
data$genotype <- relevel(data$genotype, ref = "WT")


#setting up
data$weight[data$weight=="??"] <- NA
data$weight[data$weight==""] <- NA
data$litter_size[data$litter_size=="4/8"] <- 8
data$animal_ID = as.factor(data$animal_ID)
data$cohort = as.factor(data$cohort)
data$batch = as.factor(data$batch)
data$weight=as.numeric(as.character(data$weight))
data$X.dpi=as.numeric(as.character(data$X.dpi))
data$litter_size=as.numeric(as.character(data$litter_size))


###run a linear model at each column as of column 42 to last column and save the t-values####
# Create an empty dataframe to store the results
results <- data.frame(matrix(NA, nrow = ncol(data) - 42, ncol = 4))
colnames(results) <- paste("Time", 1:4)

region_data <- read.csv("/data/chamal/projects/natvik/sir_extended/preprocessed/steph_inputs/merged_label_region_names.csv") ##to be new row name for matching column number in data eg: X23 

for (i in 43:ncol(data)) {
  # Create a vector to store t-values for each timepoint
  t_values <- c()
  
  # Iterate over each timepoint
  for (j in 1:4) {
    
    #change for looking at WT-PBS (Steph's DBM maps) vs. hemi-PBS (direct comparison with Janice's maps)
    #subset_data_PBS <- subset(data, timepoint == j &
    #                            genotype == "WT" &
    #                            injection == "PBS"#&
    #                          #sex=="M"
    #)
    subset_data_PBS <- subset(data, timepoint == j &
                                 genotype == "hemi" &
                                 injection == "PBS" #&
                                #sex=="M"
    )

    subset_data_PFF <- subset(data, timepoint == j &
                                  genotype == "hemi" &
                                  injection == "Hu-PFF"#&
                                # sex=="M"
    )
    
    if(i == 43) {
      print(paste0("Timepoint: ",j))
      print(paste0("Number of M PBS mice: ", length(which(subset_data_PBS$sex == "M"))))
      print(paste0("Number of F PBS mice: ", length(which(subset_data_PBS$sex == "F"))))
      
      print(paste0("Number of M Hu-PFF mice: ", length(which(subset_data_PFF$sex == "M"))))
      print(paste0("Number of F Hu-PFF mice: ", length(which(subset_data_PFF$sex == "F"))))
      
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
  results[i - 42, ] <- t_values
  
  # Get the column name from data without the "X" in front
  col_name <- substring(names(data)[i], 2)
  
  # Find the corresponding "Structure" in region_data based on label_number
  structure_name <- region_data$Structure[region_data$label_number == col_name]
  
  # Determine if it's one of the first 211 or the last 211
  if (i < 211+43) {
    row_name <- paste0("left ", structure_name)
  } else {
    row_name <- paste0("right ", structure_name)
  }
  
  # Assign the row name in the results dataframe
  rownames(results)[i - 42] <- row_name
}

#write.csv(results, 'rgn_t_stats_full_hemiMsPff_wtPBS.csv',row.names=TRUE)
write.csv(results, 'rgn_t_stats_full_hemiHuPff_hemiPBS.csv',row.names=TRUE)

#Compare different t-stat maps at time point 3########################################
wt_pbs_stats <- as.data.frame(read.csv('rgn_t_stats_full_hemiMsPff_wtPBS.csv',row.names=1))
hemi_pbs_stats <- as.data.frame(read.csv('rgn_t_stats_full_hemiMsPff_hemiPBS.csv',row.names=1))
plot_data <- data.frame(cbind(wt_pbs_stats$Time.3, hemi_pbs_stats$Time.3))
colnames(plot_data) <- c("WT_PBS", "hemi_PBS")
ggplot(data = plot_data, aes(x = WT_PBS, y = hemi_PBS)) +
  geom_point(size=3) +
  geom_smooth(method="lm",formula=y~x)+
  labs(x = "hemi MsPff vs. WT PBS t-stats", y = "hemi MsPff vs. hemi PBS t-stats", title = "T Stat Maps: hemi vs. WT") +
  theme_minimal() +   theme(legend.position = "none",
                            text = element_text(size = 14), plot.title=element_text(size = 13))

#######################################################################################
###load the simulated atrophy csv####
sim_atrophy_data <- read.csv('./SIR_clearance_ADPD_20240909/abm_spread_v.810.2701089378711.spread_rate.0.3445526721851775.dt.0.1.seed.35.injection_amount.30.096005758222525.clearance_gene.Zfhx2os.k1.0.4433975132729524.k2.0.7345818373987805.csv', header = FALSE)
results <- read.csv('rgn_t_stats_full_hemiMsPff_hemiPBS.csv')
rownames(results) <- results[,1]
results <- results[,-1]

# targets_data <- read.csv('/Users/stephanietullo/Documents/VSC/SIR_mouse-main/model/model/targets.csv', header = FALSE)
ABA_region_names <- read.csv('./yohan_source_full.csv', header = TRUE)
ABA_region_names<-rbind(ABA_region_names,ABA_region_names)
ABA_region_names <- ABA_region_names[,-1]


sim_atrophy_data <- cbind(ABA_region_names, sim_atrophy_data) 

###add hemisphere for the simulated atrophy####
nreg <- 209
alternating_values <- rep(c("right", "left"), each = nreg)
sim_atrophy_data$hemisphere <- alternating_values

# Create unique row names based on first and last columns 
rownames(sim_atrophy_data) <- paste(sim_atrophy_data[, ncol(sim_atrophy_data)], sim_atrophy_data[, 1], sep = " ")

# delete first columns and last column
sim_atrophy_data <- sim_atrophy_data[, -1]
sim_atrophy_data <- sim_atrophy_data[, -ncol(sim_atrophy_data)]

#PLOT
library(ggplot2)
time_steps <- 1:ncol(sim_atrophy_data)  # Assuming the time steps are represented by column indices

# plot(time_steps, sim_atrophy_data["right Caudoputamen",], type = "l", xlab = "Time Steps", ylab = "Atrophy")
# plot(1:4, results["right Caudoputamen",], type = "l", xlab = "Time Steps", ylab = "Atrophy")


# Get the common row names
common_row_names <- intersect(rownames(results), rownames(sim_atrophy_data))

# Sort matrices based on common row names
empirical_common <- results[common_row_names, ]
simulated_common <- sim_atrophy_data[common_row_names, ]

# # Sort matrices based on row names
# empirical_common <- empirical_common[order(rownames(empirical_common)), ]
# simulated_common <- simulated_common[order(rownames(simulated_common)), ]

# Find rows not common between the two dataframes
# not_common_rows <- setdiff(rownames(results), rownames(sim_atrophy_data))
not_common_rows <- setdiff(rownames(sim_atrophy_data), rownames(results))


####Plot correlations #####

# Initialize a vector to store the correlations
correlations_tp1 <- vector("numeric", length = ncol(simulated_common))
correlations_tp2 <- vector("numeric", length = ncol(simulated_common))
correlations_tp3 <- vector("numeric", length = ncol(simulated_common))
correlations_tp4 <- vector("numeric", length = ncol(simulated_common))

# Calculate Spearman correlation with each column in sim_atrophy_data
for (i in 1:ncol(simulated_common)) {
  correlations_tp1[i] <- cor(empirical_common[, 1], simulated_common[, i], method = "spearman")
  correlations_tp2[i] <- cor(empirical_common[, 2], simulated_common[, i], method = "spearman")
  correlations_tp3[i] <- cor(empirical_common[, 3], simulated_common[, i], method = "spearman")
  correlations_tp4[i] <- cor(empirical_common[, 4], simulated_common[, i], method = "spearman")
  
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

#############
library(ggplot2)

# Calculate the peak timestep and correlation value
peak_timestep <- correlation_data_t3$Time_Steps[which.max(correlation_data_t3$Correlations)]
peak_corr <- max(correlation_data_t3$Correlations)

# Create the ggplot
ggplot() +
  #geom_line(data = correlation_data_t1, aes(x = Time_Steps, y = Correlations, color = "-7 dpi"), linetype = "solid") +
  #geom_line(data = correlation_data_t2, aes(x = Time_Steps, y = Correlations, color = "30 dpi"), linetype = "solid") +
  geom_line(data = correlation_data_t3, aes(x = Time_Steps, y = Correlations, color = "90 dpi"), linetype = "solid") +
  #geom_line(data = correlation_data_t4, aes(x = Time_Steps, y = Correlations, color = "120 dpi"), linetype = "solid") +
  labs(x = "Simulated Time Steps", y = "Spearman Correlation", title = "Model Fit: Empirical vs Simulated Atrophy") +
  scale_color_manual(
    #values = c("#E69F00", "#009E73", "#0072B2", "#CC79A7"),
    values = c("#0072B2"),
    #breaks = c("-7 dpi", "30 dpi", "90 dpi", "120 dpi"),
    breaks=c("90 dpi"),
    guide = guide_legend(title = "")
  ) +
  geom_point(
    data = correlation_data_t3[correlation_data_t3$Time_Steps == peak_timestep, ],
    aes(x = Time_Steps, y = Correlations),
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
  ylim(0,0.8) +
  theme_minimal() +
  theme(legend.position = "right", text = element_text(size = 28), plot.title = element_text(hjust = 0.25))


###size 900 by 500


#########

# # Find the time step with peak correlation 

peak_timestep=correlation_data_t3$Time_Steps[which.max(correlation_data_t3$Correlations)]
peak_corr=max(correlation_data_t3$Correlations)
simulated_peak <- simulated_common[, peak_timestep]


sim_v_emp <- data.frame(simulated = simulated_peak, empirical = empirical_common)

###Plot peak correlations
ggplot(data = sim_v_emp, aes(x = simulated, y = empirical.Time.3)) +
  geom_point(size=3)+
  #geom_point(aes(colour=rownames(sim_v_emp)),size=3) +
  geom_smooth(method="lm",formula=y~x)+
  labs(x = "Simulated Atrophy", y = "Empirical Atrophy", title = "Simulated vs. Empirical Atrophy at Peak Correlation Time Step") +
  theme_minimal() +   theme(legend.position = "none",
                            text = element_text(size = 14), plot.title=element_text(size = 13))

###########IGNORE###
# Create the ggplot with the highlighted rows
library(dplyr) # Load the dplyr package for data manipulation

# Specify the rownames you want to highlight
highlighted_rownames <- c("right Caudoputamen", "left Caudoputamen", 
                          "right Substantia nigra pars compacta", "left Substantia nigra pars compacta",
                          "left Globus pallidus external segment", "right Globus pallidus external segment",
                          "left Globus pallidus internal segment", "right Globus pallidus internal segment",
                          "right Primary motor area", "left Primary motor area",
                          "right Secondary motor area","left Secondary motor area",
                          "left Primary somatosensory area lower limb", "right Primary somatosensory area lower limb",
                          "left Primary somatosensory area nose","right Primary somatosensory area nose",
                          "left Primary somatosensory area mouth","right Primary somatosensory area mouth",
                          "left Primary somatosensory area barrel field", "right Primary somatosensory area barrel field",
                          "left Primary somatosensory area trunk", "right Primary somatosensory area trunk",
                          "right Nucleus accumbens", "left Nucleus accumbens",
                          "left Primary somatosensory area upper limb", "right Primary somatosensory area upper limb")

####plot empirical atrophy values at 90 dpi####
library(tidyverse)

highlighted_data_emp <- empirical_common %>%
  mutate(
    highlight = ifelse(rownames(empirical_common) %in% highlighted_rownames, "highlight", "other")
  )

ggplot(data = highlighted_data_emp, aes(x = rownames(empirical_common), y = empirical_common[, 3])) +
  # geom_point(aes(colour=rownames(empirical_common)),size=3) 
  geom_point(aes(colour = highlight), size = 3) +
  scale_color_manual(
    values = c("highlight" = "#D55E00", "other" = "grey50")
  ) +
  theme_minimal() +
  theme(legend.position = "none")

####plot simulated atrophy values at peak corr time step####

highlighted_data_sim <- simulated_common %>%
  mutate(
    highlight = ifelse(rownames(simulated_common) %in% highlighted_rownames, "highlight", "other")
  )
ggplot(data = highlighted_data_sim, aes(x = rownames(simulated_common), y = simulated_common[, peak_timestep])) +
  geom_point(aes(colour = highlight), size = 3) +
  scale_color_manual(
    values = c("highlight" = "#D55E00", "other" = "grey50")
  ) +
  theme_minimal() +
  theme(legend.position = "none")


####plot sim vs empr
highlighted_data_both <- sim_v_emp %>%
  mutate(
    highlight = ifelse(rownames(sim_v_emp) %in% highlighted_rownames, "highlight", "other")
  )

ggplot(data = highlighted_data_both, aes(x = simulated, y = empirical.Time.3)) +
  geom_point(aes(colour = highlight), size = 3) +
  scale_color_manual(
    values = c("highlight" = "#D55E00", "other" = "grey50")
  ) +
  geom_smooth(method="lm",formula=y~x, colour="#0072B2")+
  labs(x = "Simulated Atrophy", y = "Empirical Atrophy", title = "Simulated vs. Empirical Atrophy at Peak Correlation Time Step") +
  theme_minimal() +   theme(legend.position = "none",
                            text = element_text(size = 18), plot.title=element_text(size = 16))


# ggplot(data = sim_v_emp, aes(x = simulated, y = empirical.Time.3)) +
#   geom_point(aes(colour=rownames(sim_v_emp)),size=3) +
#   geom_smooth(method="lm",formula=y~x, colour="#0072B2")+
#   labs(x = "Simulated Atrophy", y = "Empirical Atrophy", title = "Simulated vs. Empirical Atrophy at Peak Correlation Time Step") +
#   theme_minimal() +   theme(legend.position = "none",
#                             text = element_text(size = 14), plot.title=element_text(size = 13))
# 


#####################################################################


### START HERE - ADAPT TO INCLUDE CL GENE###
####instead of choosing values; highlight connectivity strength
connectivity <- read.csv('/data/chamal/projects/stephanie/SIR-modelling-project/analysis/analysis-thesis-202404/CS-atlas-edited.csv', header = TRUE) ##from Oh et al suppl table 3 right injection sites by contralteral (left) followed by ipsilateral (right) regions
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

#drop some more regions: [1] "left Median preoptic nucleus"  "left Nucleus raphe magnus"     "right Median preoptic nucleus"
#[4] "right Nucleus raphe magnus"

col_index_to_delete <- which(colnames(connectivity) == "MEPO")
connectivity <- connectivity[, -col_index_to_delete, drop = FALSE] # Delete the column

col_index_to_delete <- which(colnames(connectivity) == "MEPO.1")
connectivity <- connectivity[, -col_index_to_delete, drop = FALSE] # Delete the column

col_index_to_delete <- which(colnames(connectivity) == "RM")
connectivity <- connectivity[, -col_index_to_delete, drop = FALSE] # Delete the column

col_index_to_delete <- which(colnames(connectivity) == "RM.1")
connectivity <- connectivity[, -col_index_to_delete, drop = FALSE] # Delete the column

seed_connectivity <- connectivity["CP",, drop = FALSE]

seed_connectivity<-t(seed_connectivity[1, ])

plot(density(as.numeric(seed_connectivity[,1])), main = "Density Plot of Seed Connectivity", xlab = "Seed Connectivity")
boxplot( as.numeric(seed_connectivity[,1]), main = "Boxplot of Seed Connectivity", ylab = "Seed Connectivity")


# Create the ggplot with a reversed color scale and custom breaks
ggplot(data = sim_v_emp, aes(x = simulated, y = empirical.Time.3, color = seed_connectivity[,1])) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  labs(x = "Simulated Atrophy", y = "Empirical Atrophy", title = "Simulated vs. Empirical Atrophy at Peak Correlation Time Step") +
  scale_color_gradientn(
    colors = colorRampPalette(c("grey50", "#56B4E9"))(n = 418),  # Adjust the number of colors as needed
    limits = range(seed_connectivity[,1]),
    # breaks = custom_breaks,
    name = "Seed Connectivity"  # Set the legend title here
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 18),
    plot.title = element_text(size = 13)
  )


###show regions connected or not -BINARY- connected defined by a not zero connection strength value

# Create a new column for the color based on seed_connectivity
sim_v_emp$color <- ifelse(seed_connectivity[,1] != 0, "Connected", "Not Connected")

# Create the ggplot
ggplot(data = sim_v_emp, aes(x = simulated, y = empirical.Time.3, color = color)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  labs(x = "Simulated Atrophy", y = "Empirical Atrophy", title = "Simulated vs. Empirical Atrophy at Peak Correlation Time Step") +
  scale_color_manual(values = c("Connected" = "#56B4E9", "Not Connected" = "grey50"), name = "Seed Connectivity") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 14),
    plot.title = element_text(size = 13)
  )



######better gradient#########
# Create the ggplot with custom color breaks
ggplot(data = sim_v_emp, aes(x = simulated, y = empirical.Time.3, color = seed_connectivity[,1])) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  labs(x = "Simulated Atrophy", y = "Empirical Atrophy"
       # , title = "Simulated vs. Empirical Atrophy at Peak Correlation Time Step"
  ) +
  scale_color_gradientn(
    colors = c("grey50", "#56B4E9"),
    values = scales::rescale(c(0, min(seed_connectivity[seed_connectivity[, 1] != 0, 1]), max(seed_connectivity[,1])), c(0, 1)),
    limits = range(seed_connectivity[,1]),
    breaks = c(6.1e-7, 0.285) , # Define the desired breaks or levels for the color scale
    name = "Seed Connectivity      "  # Set the legend title here
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 25),
    plot.title = element_text(size = 16)
  )


###########based on asyn gene expression########
aSyn <- read.csv('~/Documents/VSC/SIR_mouse-main/model/model/asyn-gene-expression.csv', header = TRUE)

aSyn <- rbind(aSyn, aSyn) 

###add hemisphere for the simulated atrophy####
aSyn$hemisphere <- alternating_values

# Create unique row names based on first and last columns 
rownames(aSyn) <- paste(aSyn[, ncol(aSyn)], aSyn[, 1], sep = " ")

# delete first columns and last column
aSyn <- aSyn[, -1]
aSyn <- aSyn[, -ncol(aSyn)]


# Get the common row names and not common rows -- double check
common_row_names <- intersect(rownames(sim_v_emp), rownames(aSyn))
not_common_rows <- setdiff(rownames(aSyn), rownames(sim_v_emp))

# Sort matrices based on common row names
sim_v_emp_aSyn <- sim_v_emp[common_row_names, ]
aSyn_common <- aSyn[common_row_names, ]

# Create the ggplot with color scale and custom breaks
ggplot(data = sim_v_emp_aSyn, aes(x = simulated, y = empirical.Time.3, color = aSyn_common$asyn)) +
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
    values = scales::rescale(c(0, min(aSyn_common$asyn), max(aSyn_common$asyn)), c(0, 1)), # Define the desired breaks or levels for the color scale
    limits = range(aSyn_common$asyn),
    breaks = c(0, max(aSyn_common$asyn)),
    name = "aSyn gene expression"  # Set the legend title here
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 25),
    plot.title = element_text(size = 13)
  )


###################check correlation between asyn gene expression and simulated atrophy##########

###sim atrophy vs SNCA expression####
correlations_sim_SNCA <- cor(aSyn_common[,2], sim_v_emp_aSyn[,1], method = "spearman")
print(correlations_sim_SNCA)

# Create data frame for plotting
plot_data <- data.frame(
  x = aSyn_common[,2],
  y = sim_v_emp_aSyn[,1],
  region = rownames(sim_v_emp_aSyn)
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
correlations_sim_cxn <- cor(seed_connectivity[,1], sim_v_emp[,1], method = "spearman")
print(correlations_sim_cxn)

# Create data frame for plotting
plot_data <- data.frame(
  x = seed_connectivity[,1],
  y = sim_v_emp[,1],
  region = rownames(sim_v_emp_aSyn)
)

# Scatter plot with ggplot
ggplot(data = plot_data, aes(x = x, y = y, color = region)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  geom_text(x = 0.1, y = max(plot_data$y), 
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
  geom_text(x = 0.2, y = max(plot_data$y), 
            label = paste0("Spearman correlation = ", round(correlations_sim_cxn, 3)), hjust = 0, vjust = 1,color="black",size=6) +
  geom_text(data = subset(plot_data, x > 0.03), 
            aes(x = x, y = y, label = region),
            hjust = 0, vjust = 1, color = "black", size = 3) +
  labs(x = "Connectivity Strength", y = "Simulated Atrophy", title = "") +
  theme_minimal() +
  xlim(0,0.4)+
  theme(
    legend.position = "none",
    text = element_text(size = 18),
    plot.title = element_text(size = 13)
  )



##########check correlation between asyn gene expression and empirical atrophy at 90 dpi##########

###empirical atrophy versus SNCA####
correlations_emp_SNCA <- cor(aSyn_common[,2], sim_v_emp_aSyn[,4], method = "spearman")
print(correlations_emp_SNCA)

# Create data frame for plotting
plot_data <- data.frame(
  x = aSyn_common[,2],
  y = sim_v_emp_aSyn[,4],
  region = rownames(sim_v_emp_aSyn)
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

correlations_emp_cxn <- cor(seed_connectivity[,1], sim_v_emp[,4], method = "spearman")
print(correlations_emp_cxn)

# Create data frame for plotting
plot_data <- data.frame(
  x = seed_connectivity[,1],
  y = sim_v_emp[,4],
  region = rownames(sim_v_emp_aSyn)
)

# Scatter plot with ggplot
ggplot(data = plot_data, aes(x = x, y = y, color = region)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  geom_text(x = 0.5e-5, y = max(plot_data$y), 
            label = paste0("Spearman correlation = ", round(correlations_emp_cxn, 3)), hjust = 0, vjust = 1,color="black",size=6) +
  labs(x = "Connectivity Strength", y = "Empirical Atrophy", title = "") +
  theme_minimal() +
  xlim(0,0.5e-4)+
  theme(
    legend.position = "none",
    text = element_text(size = 18),
    plot.title = element_text(size = 13)
  )


####### with labelled top hits#######
ggplot(data = plot_data, aes(x = x, y = y, color = region)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  geom_text(x = 0.5e-5, y = max(plot_data$y)+1.25, 
            label = paste0("Spearman correlation = ", round(correlations_emp_cxn, 3)), hjust = 0, vjust = 1,color="black",size=6) +
  geom_text(data = subset(plot_data, x > 0.03), 
            aes(x = x, y = y, label = region),
            hjust = 0, vjust = 1, color = "black", size = 3) +
  labs(x = "Connectivity Strength", y = "Empirical Atrophy", title = "") +
  theme_minimal() +
  xlim(0, 0.39) +
  theme(
    legend.position = "none",
    text = element_text(size = 18),
    plot.title = element_text(size = 13)
  )

####END HERE, ASK STEPH#### 
