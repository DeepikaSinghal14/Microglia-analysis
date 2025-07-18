# Microglia-analysis

##################################################################################################
###Mciroglia analsysis for the dtermination of the characteristics of cell types in each. group###
##################################################################################################

###Load the csv files from each group and sexes###

df1 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/CNN_MALES.csv")
df2<- read.csv("/Users/Deepika/Desktop/Microglia_anal/CNN_females.csv")
df3 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/CSN_FEMALES.csv")
df4 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/CSN_MALES.csv")
df5 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/SSN_FEMALES.csv")
df6 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/SSN_MALES.csv")
df7 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/SSS_FEMALES.csv")
df8 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/SSS_MALES.csv")

##Bind the df together###

df <- rbind(df1, df2, df3,df4,df5,df6,df7,df8)

# Check column names for each data frame
colnames(df1)
colnames(df2)
colnames(df3)
colnames(df4)
colnames(df5)
colnames(df6)
colnames(df7)
colnames(df8)

# List of data frames
data_frames <- list(df1, df2, df3, df4, df5, df6, df7, df8)

# Remove the first column from each data frame
data_frames <- lapply(data_frames, function(x) x[, -1])

# Assign names to the list if desired
names(data_frames) <- c("df1", "df2", "df3", "df4", "df5", "df6", "df7", "df8")

# Check the modified data frames
print(data_frames$df1)
print(data_frames$df2)
print(data_frames$df3)

# Assign back to individual data frames
df1 <- data_frames$df1
df2 <- data_frames$df2
df3 <- data_frames$df3
df4 <- data_frames$df4
df5 <- data_frames$df5
df6 <- data_frames$df6
df7 <- data_frames$df7
df8 <- data_frames$df8

###rbind in one file###

microglia_data <- rbind(df1, df2, df3, df4, df5, df6, df7, df8)
# Remove the first two columns by indexing

# View the updated data frame
head(microglia_data)

library(dplyr)
library(tidyverse)
# Drop rows with missing values
microglia_data <- microglia_data %>% drop_na()

# Clean all columns except the last two by removing the "0..." prefix
microglia_data[, 1:(ncol(microglia_data) - 2)] <- data.frame(lapply(microglia_data[, 1:(ncol(microglia_data) - 2)], function(x) {
  gsub("0\\.\\.\\.", "", x)  # Remove "0..." and keep the numbers
}))

# Convert all columns to numeric where applicable, except for the last three columns
microglia_data[, 1:(ncol(microglia_data) - 3)] <- lapply(microglia_data[, 1:(ncol(microglia_data) - 3)], function(x) {
  as.numeric(as.character(x))  # Convert cleaned data to numeric
})

# Handle missing values (if necessary)
# Option 1: Replace NA with 0
microglia_data[is.na(microglia_data)] <- 0

# Option 2: Remove rows with NA (uncomment if this is preferred)
#df <- na.omit(df)

# View the cleaned dataset
print("Cleaned Data:")
print(head(microglia_data))


write.csv(microglia_data, "/Users/Deepika/Desktop/Microglia_anal/microglia_data_withROI.csv")

df <- microglia_data

# Select Relevant Features
desired_columns <- c("CellTerritoryVol.um3", "CellVolumes", "RamificationIndex", 
                     "NumOfEndpoints", "NumOfBranchpoints", "AvgBranchLength", 
                     "MaxBranchLength", "MinBranchLength")

features <- df[, intersect(desired_columns, colnames(df))]

# Convert Factors and Characters to Numeric
features <- as.data.frame(lapply(features, function(x) {
  if (is.factor(x)) {
    as.numeric(as.character(x))
  } else if (is.character(x)) {
    as.numeric(x)  # Convert character columns to numeric if needed
  } else {
    x
  }
}))

# Handle Missing Values (if any)
features[is.na(features)] <- 0  # Or use other imputation methods, depending on your needs

# Scale Features
scaled_features <- scale(features)

# Check the structure of the scaled features to ensure everything is numeric
str(scaled_features)

# Scale Features
scaled_features <- scale(features)

scaled_features <- as.data.frame(scaled_features)

# Correlation Matrix & Heatmap
install.packages("pheatmap")
library(pheatmap)
cor_matrix <- cor(scaled_features, use = "pairwise.complete.obs")
pheatmap(cor_matrix, color = colorRampPalette(c("blue", "white", "red"))(100), main = "Correlation Heatmap")

# Load required libraries
library(caret)
library(umap)
library(pheatmap)
library(uwot)

# Scale Features
scaled_features <- scale(features)
scaled_features <- as.data.frame(scaled_features)

# Compute Correlation Matrix & Plot Heatmap
cor_matrix <- cor(scaled_features, use = "pairwise.complete.obs")
print(cor_matrix)
pheatmap(cor_matrix, color = colorRampPalette(c("blue", "white", "red"))(100), main = "Correlation Heatmap")

# Identify Highly Correlated Variables (correlation > 0.9)
high_cor_vars <- findCorrelation(cor_matrix, cutoff = 0.9)
high_cor_names <- colnames(scaled_features)[high_cor_vars]  # Get column names

#  Print Highly Correlated Variables
print(high_cor_names)
high_cor_names <- microglia_data

# Remove Highly Correlated Variables
scaled_features_filtered <- scaled_features[, !(colnames(scaled_features) %in% high_cor_names)]
scaled_features_filtered <- as.data.frame(scaled_features_filtered)



# Convert the correlation matrix to absolute values
abs_cor_matrix <- abs(cor_matrix)

# Calculate the mean absolute correlation for each feature
# For each feature (each column), calculate the mean of absolute correlations with all other features
mean_abs_cor_per_feature <- apply(abs_cor_matrix, 2, function(x) mean(x[x != 1]))  # exclude self-correlation (diagonal)

# Print the mean absolute correlation for each feature
cat("Mean Absolute Correlation for Each Feature:\n")
print(mean_abs_cor_per_feature)


# Ensure that all columns are numeric (you can adjust the index if needed)
scaled_features_filtered <- as.data.frame(lapply(scaled_features_filtered, as.numeric))

# Check the structure of the data
str(scaled_features_filtered)

############################
############################
####optimal number of clusters
############################
#########################

install.packages("factoextra")
library(factoextra)

# 1. Elbow method using fviz_nbclust
elbow_plot <- fviz_nbclust(scaled_features_filtered, FUN = kmeans, method = "wss") + 
  ggtitle("Elbow Method for Clustering") +
  theme(
    plot.title = element_text(size = 10, face = "bold"),  # Increase title size
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 14)  # Increase axis text size
  )

print(elbow_plot)
#Add a vertical line at the optimal number of clusters (e.g., k=3)
elbow_plot 


### 2. Gap Statistic ###
install.packages("cluster")
library(cluster)

set.seed(42)
gap_stat <- clusGap(scaled_features_filtered, FUN=kmeans, nstart=25, K.max=4, B=50)

# Plot Gap Statistic
fviz_gap_stat(gap_stat) +
  theme(text = element_text(size = 14))

print(gap_stat)


############################
############################
####UMAP+PAM clustering = Partioning around medoids
############################
############################


# Load required libraries
library(uwot)
library(ggplot2)
library(cluster)  # For PAM clustering
library(RcppHNSW)

# Step 1: Apply UMAP with adjusted parameters
set.seed(123)  # For reproducibility
umap_result <- umap(
  scaled_features_filtered,
  n_neighbors = 7,   # Adjust the number of neighbors
  min_dist = 0.001,   # Set minimum distance to spread clusters
  verbose = TRUE,     # Print progress messages
  nn_method = "fnn",  # Use HNSW for nearest neighbor search (faster)
  metric = "euclidean", # Change metric to Euclidean
  
)


# Extract only the embedding (UMAP coordinates)
umap_coords <- umap_result  # Extracts transformed coordinates

# Convert to a DataFrame
umap_df <- as.data.frame(umap_coords)

# Rename columns
colnames(umap_df) <- c("UMAP1", "UMAP2")

# Verify structure
str(umap_df)  # Should show a data frame with two numeric columns


# Step 2: Perform PAM clustering on UMAP results
set.seed(123)
k <- 3  # Number of clusters
pam_result <- pam(umap_df, k)  # Apply PAM clustering

# Step 3: Assign cluster labels to the UMAP DataFrame
umap_df$cluster_pam <- as.factor(pam_result$clustering)  # Assign cluster labels
# Extract medoids correctly
medoid_indices <- pam_result$id.med  # Fix: Use id.med instead of medoids
medoids <- umap_df[medoid_indices, , drop = FALSE]  # Ensure correct extraction

# Check extracted medoids
print(medoids)

# Step 5: Pull all points closer to their respective medoids
adjust_factor <- 0.8  # Factor to control how much to move points towards medoids

# Create a copy of the original UMAP DataFrame to store adjusted points
adjusted_points <- umap_df

# Loop through each cluster
for (i in 1:nrow(medoids)) {
  # Get points in the current cluster
  cluster_points <- umap_df[umap_df$cluster_pam == i, c("UMAP1", "UMAP2")]
  
  # Get the medoid for the current cluster
  medoid <- medoids[i, c("UMAP1", "UMAP2")]
  
  # Move all points in the cluster towards the medoid
  adjusted_points[umap_df$cluster_pam == i, "UMAP1"] <- 
    cluster_points$UMAP1 + (medoid[1] - cluster_points$UMAP1) * adjust_factor
  adjusted_points[umap_df$cluster_pam == i, "UMAP2"] <- 
    cluster_points$UMAP2 + (medoid[2] - cluster_points$UMAP2) * adjust_factor
}

# Assign custom cluster labels
umap_df$cluster_label <- factor(umap_df$cluster_pam, 
                                levels = c(3, 2, 1), 
                                labels = c("Activated", "Primed", "Surveillant"))


# Plot the UMAP with cluster labels and medoids
ggplot() +
  # Plot all data points with assigned PAM clusters
  geom_point(data = umap_df, aes(x = UMAP1, y = UMAP2, color = cluster_label), 
             size = 1, alpha = 0.6) +
  
  # Highlight medoids as black stars
  geom_point(data = medoids, aes(x = UMAP1, y = UMAP2), 
             shape = 8, size = 1.5, color = "black", stroke = 1.5) +
  
  # Custom cluster colors
  scale_color_manual(
    values = c("Activated" = "red", "Primed" = "darkgreen", "Surveillant" = "blue"),
    name = "PAM Clusters"
  ) +
  
  # Labels and theme settings
  labs(title = "UMAP Visualization of PAM Clusters with Medoids",
       x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black", size = 0.5),  
    legend.position = "right",
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(size = 14),  
    axis.text.y = element_text(size = 14)   
  ) +
  
  # Add legend for medoids separately
  guides(shape = guide_legend(override.aes = list(color = "black", size = 3)))


# Check the class of 'scaled_features_filtered'
# If it's not a data frame, convert it into one
if (!is.data.frame(scaled_features_filtered)) {
  scaled_features_filtered <- as.data.frame(scaled_features_filtered)
}

# Add the 'cluster_pam' column from 'umap_df' to 'scaled_features_filtered'
scaled_features_filtered$cluster_pam <- umap_df$cluster_pam

# Check if the column 'cluster_pam' exists
head(scaled_features_filtered)  # This will print the first few rows to confirm

# Compute the mean of each feature per cluster using dplyr
cluster_means <- scaled_features_filtered %>%
  group_by(cluster_pam) %>%
  summarise(across(everything(), mean))

#  Print the cluster means
print(cluster_means)



#  Check the class of 'scaled_features_filtered'
# If it's not a data frame, convert it into one
if (!is.data.frame(features)) {
  features <- as.data.frame(features)
}

#  Add the 'cluster_pam' column from 'umap_df' to 'scaled_features_filtered'
features$cluster_pam <- umap_df$cluster_pam

# Check if the column 'cluster_pam' exists
head(features)  # This will print the first few rows to confirm

###################################################################################
###################################################################################

#filter the outliers

# Filter out outliers (but keep medoids)
filtered_df <- umap_df %>%
  filter(!(cluster_label == "Surveillant" & (UMAP2 > 10 | UMAP2 < -2))) %>%
  filter(!(cluster_label == "Primed" & (UMAP2 < -10 | UMAP1 < -5)))

# Identify the cells that were removed from Surveillant and Primed clusters

# For Surveillant cluster:
removed_surv_cells <- umap_df %>%
  filter(cluster_label == "Surveillant" & (UMAP2 > 10 | UMAP2 < -2)) %>%
  select(cell_id, cluster_label, UMAP1, UMAP2)

# For Primed cluster:
removed_primed_cells <- umap_df %>%
  filter(cluster_label == "Primed" & (UMAP2 < -10 | UMAP1 < -5)) %>%
  select(cell_id, cluster_label, UMAP1, UMAP2)

# Combine removed cells into one data frame for easier viewing
removed_cells <- bind_rows(removed_surv_cells, removed_primed_cells)

# View removed cells (Surveillant and Primed)
print("Removed cells (Surveillant and Primed clusters):")
print(removed_cells)

# Identify removed cells (as already done above)
removed_cells <- bind_rows(removed_surv_cells, removed_primed_cells)

# Merge the removed cells with the original features in df
removed_features <- inner_join(removed_cells, df, by = "cell_id")

# View the features of the removed cells
print("Features of removed cells:")
head(removed_features)  # Display first few rows to inspect


# Filter out the cells that were included in the final analysis (i.e., not removed)
included_df <- umap_df %>%
  filter(!(cluster_label == "Surveillant" & (UMAP2 > 10 | UMAP2 < -2))) %>%
  filter(!(cluster_label == "Primed" & (UMAP2 < -10 | UMAP1 < -5)))

#  Merge the included cells with the original features in df
final_included_df <- inner_join(included_df, df, by = "cell_id")

# View the final data frame with included cells
print("Final data frame with included cells including cell_id:")
head(final_included_df)  # Display first few rows to inspect

# Optionally: Save to CSV if needed
write.csv(final_included_df, "final_included_cells.csv", row.names = FALSE)

final_included_df <- read.csv("/Users/Deepika/Desktop/Microglia_anal/final_included_cells.csv")

# Load necessary package
library(dplyr)

# Perform inner join by cell_id to get the matching rows between final_included_df and microglia_data
new_df <- final_included_df %>%
  inner_join(microglia_data, by = "cell_id")  # Join by cell_id column

# View the new dataframe
print(new_df)

# Define the columns you want to keep
columns_to_keep <- c("cell_id", "cluster_label", "CellTerritoryVol.um3.x", "CellVolumes.x", 
                     "RamificationIndex.x", "NumOfEndpoints.x", "NumOfBranchpoints.x", 
                     "AvgBranchLength.x", "MaxBranchLength.x", "MinBranchLength.x", "Rat_num", "ROI", "Group", "Sex")

# Filter out the columns you want to keep from new_df
new_df <- new_df %>%
  select(all_of(columns_to_keep))

# View the new dataframe
print(new_df)


#Compute the mean of each feature per cluster using dplyr
columns_to_keep <- c( "cluster_label", "CellTerritoryVol.um3", "CellVolumes", 
                     "RamificationIndex", "NumOfEndpoints", "NumOfBranchpoints", 
                     "AvgBranchLength", "MaxBranchLength", "MinBranchLength")
cluster_final_df <- final_included_df %>%
  select(all_of(columns_to_keep))


cluster_means <- cluster_final_df %>%
  group_by(cluster_label) %>%
  summarise(across(everything(), list(mean = mean, sd = sd)))
cluster_means


# Step 5: Print the cluster means
print(cluster_means)
write.csv(cluster_means,  "/Users/Deepika/Desktop/Microglia_anal/includingROI/cluster_means_new.csv")

####Cell counts####

cluster_counts <- scaled_features_filtered %>%
  group_by(cluster_pam) %>%
  tally()  # Tally the number of rows for each cluster

# Step 4: Print the total cell count for each cluster
print(cluster_counts)

df <- microglia_data

# Load necessary library
library(dplyr)

# Assign cluster names based on cluster numbers
df <- df %>%
  mutate(cluster_name = case_when(
    scaled_features_filtered$cluster_pam == 1 ~ "Surveillant",
    scaled_features_filtered$cluster_pam == 2 ~ "Primed",
    scaled_features_filtered$cluster_pam == 3 ~ "Activated",
    TRUE ~ as.character(scaled_features_filtered$cluster_pam)  # Keep other values unchanged
  ))

# Print the first few rows to check the new cluster names
head(df)
df

umap_df <- umap_df %>% mutate(cell_id = row_number())
df <- inner_join(df %>% select(cell_id), features, by = "cell_id")



##########
##########
library(cluster)

# Run PAM clustering (example)
pam_result <- pam(scaled_features_filtered, k = 3)

# Extract medoids (these are the actual rows in your dataset)
medoid_indices <- pam_result$medoids

# Get feature values of medoids
medoid_cells <- scaled_features_filtered[medoid_indices, ]
library(dplyr)

# Compute median feature values per cluster
median_features <- scaled_features_filtered %>%
  group_by(cluster_pam) %>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))

# Find the closest cell to the median per cluster
closest_to_median <- scaled_features_filtered %>%
  group_by(cluster_pam) %>%
  mutate(distance = rowSums((across(where(is.numeric)) - median_features[match(cluster_pam, median_features$cluster_pam), -1])^2)) %>%
  slice_min(distance) %>%
  select(-distance)

closest_to_median

library(dplyr)

cluster_medians <- features %>%
  group_by(cluster_pam) %>%
  summarise(across(everything(), list(median = median, sd = sd), .names = "{.col}_{.fn}"))

print(cluster_medians)

#################################################################################################
#################################################################################################

head(df)

df <- new_df

df <- df %>% rename(cluster_name = cluster_label)
write.csv(df, "final_data.csv")

df$X <- 1:nrow(df)
library(dplyr)

###subsetting the data
df
CNN_females <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Females",
         Group == "CNN")

print(CNN_females)

write.csv(CNN_females, "/Users/Deepika/Desktop/Microglia_anal/includingROI/CNN_female.csv")


CSN_females <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Females",
         Group == "CSN")
print(CSN_females)
write.csv(CSN_females, "/Users/Deepika/Desktop/Microglia_anal/includingROI/CSN_female.csv")


SSN_females <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Females",
         Group == "SSN")
print(SSN_females)
write.csv(SSN_females, "/Users/Deepika/Desktop/Microglia_anal/includingROI/SSN_female.csv")


SSS_females <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Females",
         Group == "SSS")
print(SSS_females)
write.csv(SSS_females, "/Users/Deepika/Desktop/Microglia_anal/includingROI/SSS_female.csv")

#####Subsetting for males

CNN_males <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Males",
         Group == "CNN")
print(CNN_males)

write.csv(CNN_males, "/Users/Deepika/Desktop/Microglia_anal/includingROI/CNN_male.csv")

CSN_males <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Males",
         Group == "CSN")
print(CSN_males)

write.csv(CSN_males, "/Users/Deepika/Desktop/Microglia_anal/includingROI/CSN_male.csv")

SSN_males <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Males",
         Group == "SSN")
print(SSN_males)

write.csv(SSN_males, "/Users/Deepika/Desktop/Microglia_anal/includingROI/SSN_male.csv")

SSS_males <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Males",
         Group == "SSS")
print(SSS_males)

write.csv(SSS_males, "/Users/Deepika/Desktop/Microglia_anal/includingROI/SSS_male.csv")

#################################################################################################
##############Calculate percentages of cell types################################################
#################################################################################################

CNN_females
###For CNN_females###

# Count the number of cells in each state for each rat and ROI
state_counts_CNN_fem <- aggregate(X ~ Rat_num + ROI + cluster_name, data = CNN_females, FUN = length)

# Rename the columns for clarity
colnames(state_counts_CNN_fem) <- c("Rat_num", "ROI", "State", "Cell_Count")

# Create a data frame of all combinations of Rat_num, ROI, and State (including missing ones)
rat_state_combinations <- expand.grid(
  Rat_num = unique(state_counts_CNN_fem$Rat_num),
  ROI = unique(state_counts_CNN_fem$ROI),
  State = unique(state_counts_CNN_fem$State)
)

# Merge with the existing state counts to include all combinations
state_counts_CNN_fem <- merge(rat_state_combinations, state_counts_CNN_fem, 
                              by = c("Rat_num", "ROI", "State"), all.x = TRUE)

# Replace NAs in Cell_Count with 0, indicating no cells for that state in that rat and ROI
state_counts_CNN_fem$Cell_Count[is.na(state_counts_CNN_fem$Cell_Count)] <- 0

# Calculate total cells per Rat_num and ROI (instead of just Rat_num)
total_cells <- aggregate(Cell_Count ~ Rat_num + ROI, data = state_counts_CNN_fem, FUN = sum)
colnames(total_cells) <- c("Rat_num", "ROI", "Total_Cells")

# Merge total cell counts back into state_counts_CNN_fem
state_counts_CNN_fem <- merge(state_counts_CNN_fem, total_cells, by = c("Rat_num", "ROI"))

# Calculate the percentage of cells in each state within each Rat_num and ROI
state_counts_CNN_fem$Percentage <- (state_counts_CNN_fem$Cell_Count / state_counts_CNN_fem$Total_Cells) * 100

# View the updated data with percentages
print(state_counts_CNN_fem)

#####CSN_females###

# Count the number of cells in each state for each rat and ROI
state_counts_CSN_fem <- aggregate(X ~ Rat_num + ROI + cluster_name, data = CSN_females, FUN = length)

# Rename the columns for clarity
colnames(state_counts_CSN_fem) <- c("Rat_num", "ROI", "State", "Cell_Count")

# Create a data frame of all combinations of Rat_num, ROI, and State (including missing ones)
rat_state_combinations <- expand.grid(
  Rat_num = unique(state_counts_CSN_fem$Rat_num),
  ROI = unique(state_counts_CSN_fem$ROI),
  State = unique(state_counts_CSN_fem$State)
)

# Merge with the existing state counts to include all combinations
state_counts_CSN_fem <- merge(rat_state_combinations, state_counts_CSN_fem, 
                              by = c("Rat_num", "ROI", "State"), all.x = TRUE)

# Replace NAs in Cell_Count with 0, indicating no cells for that state in that rat and ROI
state_counts_CSN_fem$Cell_Count[is.na(state_counts_CSN_fem$Cell_Count)] <- 0

# Calculate total cells per Rat_num and ROI (instead of just Rat_num)
total_cells <- aggregate(Cell_Count ~ Rat_num + ROI, data = state_counts_CSN_fem, FUN = sum)
colnames(total_cells) <- c("Rat_num", "ROI", "Total_Cells")

# Merge total cell counts back into state_counts_CSN_fem
state_counts_CSN_fem <- merge(state_counts_CSN_fem, total_cells, by = c("Rat_num", "ROI"))

# Calculate the percentage of cells in each state within each Rat_num and ROI
state_counts_CSN_fem$Percentage <- (state_counts_CSN_fem$Cell_Count / state_counts_CSN_fem$Total_Cells) * 100

# View the updated data with percentages
print(state_counts_CSN_fem)

#####SSN females


# Count the number of cells in each state for each rat
# Count the number of cells in each state for each rat and ROI
state_counts_SSN_fem <- aggregate(X ~ Rat_num + ROI + cluster_name, data = SSN_females, FUN = length)

# Rename the columns for clarity
colnames(state_counts_SSN_fem) <- c("Rat_num", "ROI", "State", "Cell_Count")

# Create a data frame of all combinations of Rat_num, ROI, and State (including missing ones)
rat_state_combinations <- expand.grid(
  Rat_num = unique(state_counts_SSN_fem$Rat_num),
  ROI = unique(state_counts_SSN_fem$ROI),
  State = unique(state_counts_SSN_fem$State)
)

# Merge with the existing state counts to include all combinations
state_counts_SSN_fem <- merge(rat_state_combinations, state_counts_SSN_fem, 
                              by = c("Rat_num", "ROI", "State"), all.x = TRUE)

# Replace NAs in Cell_Count with 0, indicating no cells for that state in that rat and ROI
state_counts_SSN_fem$Cell_Count[is.na(state_counts_SSN_fem$Cell_Count)] <- 0

# Calculate total cells per Rat_num and ROI (instead of just Rat_num)
total_cells <- aggregate(Cell_Count ~ Rat_num + ROI, data = state_counts_SSN_fem, FUN = sum)
colnames(total_cells) <- c("Rat_num", "ROI", "Total_Cells")

# Merge total cell counts back into state_counts_SSN_fem
state_counts_SSN_fem <- merge(state_counts_SSN_fem, total_cells, by = c("Rat_num", "ROI"))

# Calculate the percentage of cells in each state within each Rat_num and ROI
state_counts_SSN_fem$Percentage <- (state_counts_SSN_fem$Cell_Count / state_counts_SSN_fem$Total_Cells) * 100

# View the updated data with percentages
print(state_counts_SSN_fem)

#####SSS Females

# Count the number of cells in each state for each rat and ROI
state_counts_SSS_fem <- aggregate(X ~ Rat_num + ROI + cluster_name, data = SSS_females, FUN = length)

# Rename the columns for clarity
colnames(state_counts_SSS_fem) <- c("Rat_num", "ROI", "State", "Cell_Count")

# Create a data frame of all combinations of Rat_num, ROI, and State (including missing ones)
rat_state_combinations <- expand.grid(
  Rat_num = unique(state_counts_SSS_fem$Rat_num),
  ROI = unique(state_counts_SSS_fem$ROI),
  State = unique(state_counts_SSS_fem$State)
)

# Merge with the existing state counts to include all combinations
state_counts_SSS_fem <- merge(rat_state_combinations, state_counts_SSS_fem, 
                              by = c("Rat_num", "ROI", "State"), all.x = TRUE)

# Replace NAs in Cell_Count with 0, indicating no cells for that state in that rat and ROI
state_counts_SSS_fem$Cell_Count[is.na(state_counts_SSS_fem$Cell_Count)] <- 0

# Calculate total cells per Rat_num and ROI (instead of just Rat_num)
total_cells <- aggregate(Cell_Count ~ Rat_num + ROI, data = state_counts_SSS_fem, FUN = sum)
colnames(total_cells) <- c("Rat_num", "ROI", "Total_Cells")

# Merge total cell counts back into state_counts_SSS_fem
state_counts_SSS_fem <- merge(state_counts_SSS_fem, total_cells, by = c("Rat_num", "ROI"))

# Calculate the percentage of cells in each state within each Rat_num and ROI
state_counts_SSS_fem$Percentage <- (state_counts_SSS_fem$Cell_Count / state_counts_SSS_fem$Total_Cells) * 100

# View the updated data with percentages
print(state_counts_SSS_fem)


######CSN males
# Count the number of cells in each state for each rat
# Count the number of cells in each state for each rat and ROI
state_counts_CSN_mal <- aggregate(X ~ Rat_num + ROI + cluster_name, data = CSN_males, FUN = length)

# Rename the columns for clarity
colnames(state_counts_CSN_mal) <- c("Rat_num", "ROI", "State", "Cell_Count")

# Create a data frame of all combinations of Rat_num, ROI, and State (including missing ones)
rat_state_combinations <- expand.grid(
  Rat_num = unique(state_counts_CSN_mal$Rat_num),
  ROI = unique(state_counts_CSN_mal$ROI),
  State = unique(state_counts_CSN_mal$State)
)

# Merge with the existing state counts to include all combinations
state_counts_CSN_mal <- merge(rat_state_combinations, state_counts_CSN_mal, 
                              by = c("Rat_num", "ROI", "State"), all.x = TRUE)

# Replace NAs in Cell_Count with 0, indicating no cells for that state in that rat and ROI
state_counts_CSN_mal$Cell_Count[is.na(state_counts_CSN_mal$Cell_Count)] <- 0

# Calculate total cells per Rat_num and ROI
total_cells <- aggregate(Cell_Count ~ Rat_num + ROI, data = state_counts_CSN_mal, FUN = sum)
colnames(total_cells) <- c("Rat_num", "ROI", "Total_Cells")

# Merge total cell counts back into state_counts
state_counts_CSN_mal <- merge(state_counts_CSN_mal, total_cells, by = c("Rat_num", "ROI"))

# Calculate the percentage of cells in each state within each Rat_num and ROI
state_counts_CSN_mal$Percentage <- (state_counts_CSN_mal$Cell_Count / state_counts_CSN_mal$Total_Cells) * 100

# View the updated data with percentages
print(state_counts_CSN_mal)


######SSN males
# Count the number of cells in each state for each rat
# Count the number of cells in each state for each rat and ROI

# Count the number of cells in each state for each rat and ROI
state_counts_SSN_mal <- aggregate(X ~ Rat_num + ROI + cluster_name, data = SSN_males, FUN = length)

# Rename the columns for clarity
colnames(state_counts_SSN_mal) <- c("Rat_num", "ROI", "State", "Cell_Count")

# Create a data frame of all combinations of Rat_num, ROI, and State (even those with no cells)
rat_state_combinations <- expand.grid(
  Rat_num = unique(state_counts_SSN_mal$Rat_num),
  ROI = unique(state_counts_SSN_mal$ROI),
  State = unique(state_counts_SSN_mal$State)
)

# Merge with existing state counts to include all combinations
state_counts_SSN_mal <- merge(rat_state_combinations, state_counts_SSN_mal, 
                              by = c("Rat_num", "ROI", "State"), all.x = TRUE)

# Replace NAs in Cell_Count with 0, indicating no cells for that state in that rat and ROI
state_counts_SSN_mal$Cell_Count[is.na(state_counts_SSN_mal$Cell_Count)] <- 0

# Calculate total cells per Rat_num and ROI
total_cells <- aggregate(Cell_Count ~ Rat_num + ROI, data = state_counts_SSN_mal, FUN = sum)
colnames(total_cells) <- c("Rat_num", "ROI", "Total_Cells")

# Merge total cell counts back into state_counts
state_counts_SSN_mal <- merge(state_counts_SSN_mal, total_cells, by = c("Rat_num", "ROI"))

# Calculate the percentage of cells in each state within each Rat_num and ROI
state_counts_SSN_mal$Percentage <- (state_counts_SSN_mal$Cell_Count / state_counts_SSN_mal$Total_Cells) * 100

# View the updated data with percentages
print(state_counts_SSN_mal)

### SSS Males

# Count the number of cells in each state for each rat
# Count the number of cells in each state for each rat and ROI
state_counts_SSS_mal <- aggregate(X ~ Rat_num + ROI + cluster_name, data = SSS_males, FUN = length)

# Rename the columns for clarity
colnames(state_counts_SSS_mal) <- c("Rat_num", "ROI", "State", "Cell_Count")

# Create a data frame of all combinations of Rat_num, ROI, and State (even those with no cells)
rat_state_combinations <- expand.grid(
  Rat_num = unique(state_counts_SSS_mal$Rat_num),
  ROI = unique(state_counts_SSS_mal$ROI),
  State = unique(state_counts_SSS_mal$State)
)

# Merge with existing state counts to include all combinations
state_counts_SSS_mal <- merge(rat_state_combinations, state_counts_SSS_mal, 
                              by = c("Rat_num", "ROI", "State"), all.x = TRUE)

# Replace NAs in Cell_Count with 0, indicating no cells for that state in that rat and ROI
state_counts_SSS_mal$Cell_Count[is.na(state_counts_SSS_mal$Cell_Count)] <- 0

# Calculate total cells per Rat_num and ROI
total_cells <- aggregate(Cell_Count ~ Rat_num + ROI, data = state_counts_SSS_mal, FUN = sum)
colnames(total_cells) <- c("Rat_num", "ROI", "Total_Cells")

# Merge total cell counts back into state_counts
state_counts_SSS_mal <- merge(state_counts_SSS_mal, total_cells, by = c("Rat_num", "ROI"))

# Calculate the percentage of cells in each state within each Rat_num and ROI
state_counts_SSS_mal$Percentage <- (state_counts_SSS_mal$Cell_Count / state_counts_SSS_mal$Total_Cells) * 100

# View the updated data with percentages
print(state_counts_SSS_mal)


#### CNN Males

# Count the number of cells in each state for each rat
# Count the number of cells in each state for each rat and ROI
state_counts_CNN_mal <- aggregate(X ~ Rat_num + ROI + cluster_name, data = CNN_males, FUN = length)

# Rename the columns for clarity
colnames(state_counts_CNN_mal) <- c("Rat_num", "ROI", "State", "Cell_Count")

# Create a data frame of all combinations of Rat_num, ROI, and State (even those with no cells)
rat_state_combinations <- expand.grid(
  Rat_num = unique(state_counts_CNN_mal$Rat_num),
  ROI = unique(state_counts_CNN_mal$ROI),
  State = unique(state_counts_CNN_mal$State)
)

# Merge with existing state counts to include all combinations
state_counts_CNN_mal <- merge(rat_state_combinations, state_counts_CNN_mal, 
                              by = c("Rat_num", "ROI", "State"), all.x = TRUE)

# Replace NAs in Cell_Count with 0, indicating no cells for that state in that rat and ROI
state_counts_CNN_mal$Cell_Count[is.na(state_counts_CNN_mal$Cell_Count)] <- 0

# Calculate total cells per Rat_num and ROI
total_cells <- aggregate(Cell_Count ~ Rat_num + ROI, data = state_counts_CNN_mal, FUN = sum)
colnames(total_cells) <- c("Rat_num", "ROI", "Total_Cells")

# Merge total cell counts back into state_counts
state_counts_CNN_mal <- merge(state_counts_CNN_mal, total_cells, by = c("Rat_num", "ROI"))

# Calculate the percentage of cells in each state within each Rat_num and ROI
state_counts_CNN_mal$Percentage <- (state_counts_CNN_mal$Cell_Count / state_counts_CNN_mal$Total_Cells) * 100

# View the updated data with percentages
print(state_counts_CNN_mal)

#################Sex difference_Groupwise########
##################################################
library(dplyr)

####Sex wise SSS######
# Step 1: Add a Sex column to the male and female datasets
state_counts_SSS_mal$Sex <- "Male"
state_counts_SSS_fem$Sex <- "Female"

# Step 2: Combine the datasets for both sexes
SSS <- bind_rows(state_counts_SSS_mal, state_counts_SSS_fem)

# Step 3: Calculate the total number of cells for each rat within each sex
total_cells_per_sex_rat <- SSS %>%
  group_by(Rat_num, Sex) %>%
  summarise(Total_Sex_Cells = sum(Cell_Count), .groups = "drop")

# Step 4: Merge the total cells by rat and sex back into the combined dataset
SSS <- SSS %>%
  left_join(total_cells_per_sex_rat, by = c("Rat_num", "Sex"))

# Step 5: Calculate the percentage of cells in each state for each rat
# Replace NA percentages with 0% for missing states
SSS <- SSS %>%
  mutate(Percentage = ifelse(
    Total_Sex_Cells > 0, 
    (Cell_Count / Total_Sex_Cells) * 100, 
    0
  ))

# View the final dataset
print(SSS)

#write.csv(SSS, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/sexdiff_SSS.csv")

#################################################################################################
#################################################################################################
#################################################################################################


state_counts_CNN_fem$Sex <- "Females"
state_counts_CNN_fem$Group <- "CNN"
state_counts_CNN_mal$Sex <- "Males"
state_counts_CNN_mal$Group <- "CNN"
state_counts_CSN_fem$Sex <- "Females"
state_counts_CSN_fem$Group <- "CSN"
state_counts_CSN_mal$Sex <-"Males"
state_counts_CSN_mal$Group <- "CSN"
state_counts_SSN_fem$Sex <- "Females"
state_counts_SSN_fem$Group <- "SSN"
state_counts_SSN_mal$Sex <- "Males"
state_counts_SSN_mal$Group <- "SSN"
state_counts_SSS_mal$Sex <- "Males"
state_counts_SSS_mal$Group <- "SSS"
state_counts_SSS_fem$Sex <- "Females"
state_counts_SSS_fem$Group <- "SSS"
data_percent <- bind_rows(state_counts_CNN_fem,state_counts_CNN_mal,state_counts_CSN_fem,state_counts_CSN_mal,state_counts_SSN_fem,
                          state_counts_SSN_mal,state_counts_SSS_fem, state_counts_SSS_mal)

data_percent

write.csv(data_percent, "/Users/Deepika/Desktop/Microglia_anal/data_percent.csv")
data_percent <- read.csv("/Users/Deepika/Desktop/Microglia_anal/data_percent.csv")
#################################################################################################

data_percent <- data_percent %>%
  mutate(Group_Rat_num = paste(Group, Rat_num, sep = "_"))
data_percent

data_percent_group <- data_percent %>% 
  group_by(Group_Rat_num, State, Sex) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE), .groups = "drop")

data_percent_stress <- data_percent %>%
  filter(Condition == "Stress")
 

#################################################################################################
##Fig2
library(dplyr)

# Assuming your data is in a dataframe called `df`
# Group by 'State' and 'Sex' and calculate the mean percentage for each group
mean_percentage <- data_percent_group %>%
  group_by(State, Sex) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE))  # Calculate mean percentage for each group

# View the result
mean_percentage

library(tidyplots)
sex <- data_percent_group %>% 
  tidyplot(x = State, y = Percentage, color = Sex, fill = Sex) %>% # Map both color and fill to Sex
  add_mean_bar(alpha = 0.4) %>% 
  add_sem_errorbar() %>% 
  add_data_points_beeswarm() %>% 
  adjust_y_axis(limits = c(0, 100)) %>% 
  remove_x_axis_title() %>% 
  reorder_x_axis_labels("Surveillant", "Primed", "Activated") %>% 
  add(ggplot2::labs(
    x = "Morphological State",
    y = "Proportion of Cells (%)"
  )) %>% 
  add(
    ggplot2::scale_color_manual(
      values = c("Males" = "blue", "Females" = "red") # Blue for males, red for females
    )
  ) %>% 
  add(
    ggplot2::scale_fill_manual(
      values = c("Males" = "blue", "Females" = "red") # Ensure bars match the points
    )
  ) %>% 
  add(
    ggplot2::theme(
      aspect.ratio = 1,
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),  # Adjust x-axis text size
      axis.text.y = element_text(size = 11),  # Adjust y-axis text size
      axis.title.x = element_text(size = 9),  # Adjust x-axis title size
      axis.title.y = element_text(size = 12),  # Adjust y-axis title size
      plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
      legend.text = element_text(size = 11),  # Adjust legend text size
      legend.title = element_text(size = 11)  # Adjust legend title size
    )
  )

sex

###P value and cohen's d
data_percent_sex <- data_percent %>% 
  group_by(State, Sex, Group) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE), .groups = "drop")

data_percent_sex

summary_table_sex <- data_percent_sex %>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Sex, 
    values_from = c(Mean, SD, n),
    names_glue = "{Sex}_{.value}"
  )
summary_table_sex

# Modify the sample sizes (Females_n and Males_n) for each state
summary_table_sex <- summary_table_sex %>%
  mutate(
    Females_n = case_when(
      State == "Activated" ~ 16,  # Change to your desired value
      State == "Primed" ~ 16,     # Change to your desired value
      State == "Surveillant" ~ 16, # Change to your desired value
      TRUE ~ Females_n  # Keep the current value for other rows
    ),
    Males_n = case_when(
      State == "Activated" ~ 16,  # Change to your desired value
      State == "Primed" ~ 16,     # Change to your desired value
      State == "Surveillant" ~ 16, # Change to your desired value
      TRUE ~ Males_n  # Keep the current value for other rows
    )
  )

# View the updated table_table

summary_table_sex

# Define pooled standard deviation function
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Load necessary library
library(dplyr)

comparison_table_sexdiff <- summary_table_sex %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Females_SD, Males_SD, 16, 16),  # Use 16 for both Females_n and Males_n
    Cohens_d = (Females_Mean - Males_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = data_percent$Percentage[data_percent$State == State & data_percent$Sex == "Females"],
      y = data_percent$Percentage[data_percent$State == State & data_percent$Sex == "Males"],
      paired = FALSE
    )$p.value
  )

# View the result
print(comparison_table_sexdiff)


# Perform ANOVA on the subsetted data
model <- aov(Percentage ~ Sex * State, data = data_percent)  
summary(model) 

# Calculate Eta-squared
eta_sq_value <- eta_squared(model)
eta_sq_value <- as.numeric(eta_sq_value)

library(effectsize)

# Get Eta-squared values
eta_sq_values <- eta_squared(model)
eta_sq_values
# Extract the numeric Eta-squared value
eta_sq_value <- as.numeric(eta_sq_values$Eta2)  # Extract the first value

# Degrees of freedom (adjust based on your analysis)
df <- 2  # Change this to the correct df

# Convert to Cohen’s d
cohen_d <- sqrt(eta_sq_value / (1 - eta_sq_value)) * sqrt(df)
cohen_d

#################################################################################################
##Sex diffference dorsal horns

# Load necessary libraries
library(dplyr)
library(tidyr)
data_percent
# Assuming your dataset is named 'data_percent'
# Group the data into Superficial Dorsal Horn (ROI 1 + 2) and Deep Dorsal Horn (ROI 3 + 4)
# Step 1: Add a new column to classify ROI into superficial and deep dorsal horn
data_percent <- data_percent %>%
  mutate(Dorsal_Horn = case_when(
    ROI %in% c(1, 2) ~ 'Superficial',
    ROI %in% c(3, 4) ~ 'Deep'
  ))

# Step 2: Group by Dorsal_Horn, State, Sex, and Group, and summarize the data
grouped_data <- data_percent %>%
  group_by(Dorsal_Horn, State, Sex, Group, Rat_num) %>%
  summarize(
    Total_Cell_Count = sum(Cell_Count),
    Total_Cells = sum(Total_Cells),
    Percentage = mean(Percentage),
    .groups = 'drop' # to drop the grouping after summarizing
  )

# View the grouped data
print(grouped_data)

#############
####################

ROI_summary <- data_percent %>%
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  group_by(State, Dorsal_Horn) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE), .groups = 'drop') %>%
  # Adjust Percentage to be a proportion of the total for each 'Sex'
  group_by(Dorsal_Horn) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

# Plot the stacked bar chart with the specified appearance
library(ggplot2)
library(dplyr)

ROI_summary %>% 
  # Create a new x-axis variable that combines Sex and Dorsal_Horn
  mutate(Dorsal = paste(Dorsal_Horn)) %>% 
  ggplot(aes(x = Dorsal, y = Percentage, fill = State)) +  # Use the new combined variable on x-axis
  geom_bar(stat = "identity", width = 0.5, position = "stack") +  # Stacked bars
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 5, fontface = "bold") +  # White text for labels
  scale_fill_manual(values = c("Surveillant" = "blue",  
                               "Primed" = "darkgreen",  
                               "Activated" = "red")) +  # Custom colors
  labs(title = "Proportional Morphological State Distribution by Sex and Region",  
       x = "Sex and Region",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(  
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    axis.text.x = element_text(angle = 0, hjust = 1, size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 12),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 12, hjust = 0.5),  # Adjust title size and center it
    legend.text = element_text(size = 12),  # Adjust legend text size
    legend.title = element_text(size = 12)  # Adjust legend title size
  ) + 
  scale_y_continuous(labels = function(x) paste0(x, "%"),  
                     expand = c(0,0),  
                     breaks = seq(0, 100, by = 25))  # Y axis as percentage



grouped_data_deep <- grouped_data%>%
  filter(Dorsal_Horn == "Deep")
grouped_data_deep

ROI_summary
data_percent_group

##Sex differences in deep dorsal horn region

sex_dh_deep <- grouped_data_deep %>% 
  tidyplot(x = State, y = Percentage, color = Sex, fill = Sex) %>% # Map both color and fill to Sex
  add_mean_bar(alpha = 0.4) %>% 
  add_sem_errorbar() %>% 
  add_data_points_beeswarm() %>% 
  adjust_y_axis(limits = c(0, 100)) %>% 
  remove_x_axis_title() %>% 
  reorder_x_axis_labels("Surveillant", "Primed", "Activated") %>% 
  add(ggplot2::labs(
    x = "Morphological State",
    y = "Proportion of Cells (%)"
  )) %>% 
  add(
    ggplot2::scale_color_manual(
      values = c("Males" = "blue", "Females" = "red") # Blue for males, red for females
    )
  ) %>% 
  add(
    ggplot2::scale_fill_manual(
      values = c("Males" = "blue", "Females" = "red") # Ensure bars match the points
    )
  ) %>% 
  add(
    ggplot2::theme(
      aspect.ratio = 1,
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),  # Adjust x-axis text size
      axis.text.y = element_text(size = 11),  # Adjust y-axis text size
      axis.title.x = element_text(size = 9),  # Adjust x-axis title size
      axis.title.y = element_text(size = 12),  # Adjust y-axis title size
      plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
      legend.text = element_text(size = 11),  # Adjust legend text size
      legend.title = element_text(size = 11)  # Adjust legend title size
    )
  )


sex_dh_deep

data_percent_sex_deep <- grouped_data_deep %>% 
  group_by(State, Sex, Group) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE), .groups = "drop")

data_percent_sex_deep

summary_table_deep <- data_percent_sex_deep %>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 16,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Sex, 
    values_from = c(Mean, SD, n),
    names_glue = "{Sex}_{.value}"
  )
summary_table_deep

# Modify the sample sizes (Females_n and Males_n) for each state
summary_table_deep <- summary_table_deep %>%
  mutate(
    Females_n = case_when(
      State == "Activated" ~ 16,  # Change to your desired value
      State == "Primed" ~ 16,     # Change to your desired value
      State == "Surveillant" ~ 16, # Change to your desired value
      TRUE ~ Females_n  # Keep the current value for other rows
    ),
    Males_n = case_when(
      State == "Activated" ~ 16,  # Change to your desired value
      State == "Primed" ~ 16,     # Change to your desired value
      State == "Surveillant" ~ 16, # Change to your desired value
      TRUE ~ Males_n  # Keep the current value for other rows
    )
  )

# View the updated table_table

summary_table_deep

# Define pooled standard deviation function
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Load necessary library
library(dplyr)

comparison_table_sexdiff_deep <- summary_table_deep %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Females_SD, Males_SD, 16, 16),  # Use 16 for both Females_n and Males_n
    Cohens_d = (Females_Mean - Males_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = grouped_data_deep$Percentage[grouped_data_deep$State == State & grouped_data_deep$Sex == "Females"],
      y = grouped_data_deep$Percentage[grouped_data_deep$State == State & grouped_data_deep$Sex == "Males"],
      paired = FALSE
    )$p.value
  )

# View the result
print(comparison_table_sexdiff_deep)

# Perform ANOVA on the subsetted data
model <- aov(Percentage ~ Sex * State, data = grouped_data_deep)  
summary(model) 

# Calculate Eta-squared
eta_sq_value <- eta_squared(model)
eta_sq_value <- as.numeric(eta_sq_value)

library(effectsize)

# Get Eta-squared values
eta_sq_values <- eta_squared(model)
eta_sq_values
# Extract the numeric Eta-squared value
eta_sq_value <- as.numeric(eta_sq_values$Eta2)  # Extract the first value

# Degrees of freedom (adjust based on your analysis)
df <- 2  # Change this to the correct df

# Convert to Cohen’s d
cohen_d <- sqrt(eta_sq_value / (1 - eta_sq_value)) * sqrt(df)
cohen_d


### Sex differences in superficial dorsal horn region

grouped_data_sup <- grouped_data%>%
  filter(Dorsal_Horn == "Superficial")
grouped_data_sup

sex_dh_sup <- grouped_data_sup %>% 
  tidyplot(x = State, y = Percentage, color = Sex, fill = Sex) %>% # Map both color and fill to Sex
  add_mean_bar(alpha = 0.4) %>% 
  add_sem_errorbar() %>% 
  add_data_points_beeswarm() %>% 
  adjust_y_axis(limits = c(0, 100)) %>% 
  remove_x_axis_title() %>% 
  reorder_x_axis_labels("Surveillant", "Primed", "Activated") %>% 
  add(ggplot2::labs(
    x = "Morphological State",
    y = "Proportion of Cells (%)"
  )) %>% 
  add(
    ggplot2::scale_color_manual(
      values = c("Males" = "blue", "Females" = "red") # Blue for males, red for females
    )
  ) %>% 
  add(
    ggplot2::scale_fill_manual(
      values = c("Males" = "blue", "Females" = "red") # Ensure bars match the points
    )
  ) %>% 
  add(
    ggplot2::theme(
      aspect.ratio = 1,
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),  # Adjust x-axis text size
      axis.text.y = element_text(size = 11),  # Adjust y-axis text size
      axis.title.x = element_text(size = 9),  # Adjust x-axis title size
      axis.title.y = element_text(size = 12),  # Adjust y-axis title size
      plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
      legend.text = element_text(size = 11),  # Adjust legend text size
      legend.title = element_text(size = 11)  # Adjust legend title size
    )
  )

sex_dh_sup

data_percent_sex_sup <- grouped_data_sup %>% 
  group_by(State, Sex, Group) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE), .groups = "drop")

data_percent_sex_sup

summary_table_sup <- data_percent_sex_sup %>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 16,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Sex, 
    values_from = c(Mean, SD, n),
    names_glue = "{Sex}_{.value}"
  )
summary_table_sup

# Modify the sample sizes (Females_n and Males_n) for each state
summary_table_sup <- summary_table_sup %>%
  mutate(
    Females_n = case_when(
      State == "Activated" ~ 16,  # Change to your desired value
      State == "Primed" ~ 16,     # Change to your desired value
      State == "Surveillant" ~ 16, # Change to your desired value
      TRUE ~ Females_n  # Keep the current value for other rows
    ),
    Males_n = case_when(
      State == "Activated" ~ 16,  # Change to your desired value
      State == "Primed" ~ 16,     # Change to your desired value
      State == "Surveillant" ~ 16, # Change to your desired value
      TRUE ~ Males_n  # Keep the current value for other rows
    )
  )

# View the updated table_table

summary_table_sup

# Define pooled standard deviation function
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Load necessary library
library(dplyr)

comparison_table_sexdiff_sup <- summary_table_sup %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Females_SD, Males_SD, 16, 16),  # Use 16 for both Females_n and Males_n
    Cohens_d = (Females_Mean - Males_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = grouped_data_sup$Percentage[grouped_data_sup$State == State & grouped_data_sup$Sex == "Females"],
      y = grouped_data_sup$Percentage[grouped_data_sup$State == State & grouped_data_sup$Sex == "Males"],
      paired = FALSE
    )$p.value
  )

# View the result
print(comparison_table_sexdiff_sup)

# Perform ANOVA on the subsetted data
model <- aov(Percentage ~ Sex * State, data = grouped_data_sup)  
summary(model) 

# Calculate Eta-squared
eta_sq_value <- eta_squared(model)
eta_sq_value <- as.numeric(eta_sq_value)

library(effectsize)

# Get Eta-squared values
eta_sq_values <- eta_squared(model)
eta_sq_values
# Extract the numeric Eta-squared value
eta_sq_value <- as.numeric(eta_sq_values$Eta2)  # Extract the first value

# Degrees of freedom (adjust based on your analysis)
df <- 2  # Change this to the correct df

# Convert to Cohen’s d
cohen_d <- sqrt(eta_sq_value / (1 - eta_sq_value)) * sqrt(df)
cohen_d


### Sex differences in injection groups

#SSN_CSN


library(dplyr)
library(tidyr)
library(broom)

# Step 1: Summary Table for Mean, SD, and Sample Size by State and Group (Separated by Sex)
summary_table_SSN_CSN <- SSN_CSN %>% 
  group_by(State, Group) %>% 
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 8,  # Set n to 8 for both groups (SSN and CSN)
    .groups = "drop"
  ) %>% 
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )

# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
comparison_table_SSN_CSN <- summary_table_SSN_CSN %>% 
  filter(!is.na(SSN_SD) & !is.na(CSN_SD)) %>% 
  rowwise() %>% 
  mutate(
    # Use fixed n1 and n2 values (8) for both groups
    Pooled_SD = pooled_sd(SSN_SD, CSN_SD, 8, 8),  # Set both n1 and n2 to 8
    Cohens_d = (SSN_Mean - CSN_Mean) / Pooled_SD,
    t_test_p_value = wilcox.test(
      x = SSN_CSN$Percentage[SSN_CSN$State == State & SSN_CSN$Group == "SSN"],
      y = SSN_CSN$Percentage[SSN_CSN$State == State & SSN_CSN$Group == "CSN"],
      paired = FALSE
    )$p.value
  ) %>% 
  ungroup()

# View the result
comparison_table_SSN_CSN


# Step 4: Display Results
print(comparison_table_SSN_CSN)

write.csv(comparison_table_SSN_CSN, "/Users/Deepika/Desktop/Microglia_anal/includingROI/d_P_SSN_CSN.csv")

#################################################################################################
#################################################################################################
#Fig3

# SSN vs CSN

SSN <- data_percent %>%
  filter(Group == "SSN")
CSN <- data_percent %>%
  filter(Group == "CSN")

SSN_CSN <- bind_rows(SSN,CSN)
SSN_CSN
summary_table_SSN_CSN
summary_table_SSN_CSN <- SSN_CSN %>%
  group_by(Group, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Group) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

summary_table_SSN_CSN %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  mutate(Group = factor(Group, levels = c("SSN", "CSN"))) %>%
  ggplot(aes(x = Group, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 4, fontface = "bold") +  # White text for labels
  scale_fill_manual(values = c("Surveillant" = "blue", 
                               "Primed" = "darkgreen", 
                               "Activated" = "red")) +  # Custom colors
  labs(title = "Proportional Morphological State Distribution by SSN_CSN",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 12),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 11)  # Adjust legend title size
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),  
                     expand = c(0,0),  
                     breaks = seq(0, 100, by = 25))  # Proper % labels

######

library(dplyr)
library(tidyr)
library(broom)

# Step 1: Summary Table for Mean, SD, and Sample Size by State and Group (Separated by Sex)
summary_table_SSN_CSN <- SSN_CSN %>% 
  group_by(State, Group) %>% 
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 8,  # Set n to 8 for both groups (SSN and CSN)
    .groups = "drop"
  ) %>% 
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )

# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
comparison_table_SSN_CSN <- summary_table_SSN_CSN %>% 
  filter(!is.na(SSN_SD) & !is.na(CSN_SD)) %>% 
  rowwise() %>% 
  mutate(
    # Use fixed n1 and n2 values (8) for both groups
    Pooled_SD = pooled_sd(SSN_SD, CSN_SD, 8, 8),  # Set both n1 and n2 to 8
    Cohens_d = (SSN_Mean - CSN_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSN_CSN$Percentage[SSN_CSN$State == State & SSN_CSN$Group == "SSN"],
      y = SSN_CSN$Percentage[SSN_CSN$State == State & SSN_CSN$Group == "CSN"],
      paired = FALSE
    )$p.value
  ) %>% 
  ungroup()

# View the result
comparison_table_SSN_CSN

#################################################################################################
#SSN vs SSS
SSN
SSS
SSN <- data_percent %>%
  filter(Group == "SSN")
SSS <- data_percent %>%
  filter(Group == "SSS")
SSS_SSN <- bind_rows(SSN,SSS)
SSS_SSN

summary_table_SSS_SSN <- SSS_SSN %>% 
  group_by(State, Group) %>% 
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 8,  # Set n to 8 for both groups (SSN and CSN)
    .groups = "drop"
  ) %>% 
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )
summary_table_SSS_SSN
# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
comparison_table_SSS_SSN <- summary_table_SSS_SSN %>% 
  filter(!is.na(SSS_SD) & !is.na(SSN_SD)) %>% 
  rowwise() %>% 
  mutate(
    # Use fixed n1 and n2 values (8) for both groups
    Pooled_SD = pooled_sd(SSS_SD, SSN_SD, 8, 8),  # Set both n1 and n2 to 8
    Cohens_d = (SSS_Mean - SSN_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSS_SSN$Percentage[SSS_SSN$State == State & SSS_SSN$Group == "SSS"],
      y = SSS_SSN$Percentage[SSS_SSN$State == State & SSS_SSN$Group == "SSN"],
      paired = FALSE
    )$p.value
  ) %>% 
  ungroup()

# View the result
comparison_table_SSS_SSN

summary_table_SSS_SSN
summary_table_SSS_SSN <- SSS_SSN %>%
  group_by(Group, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Group) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

summary_table_SSS_SSN %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  mutate(Group = factor(Group, levels = c("SSN", "SSS"))) %>%
  ggplot(aes(x = Group, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 4, fontface = "bold") +  # White text for labels
  scale_fill_manual(values = c("Surveillant" = "blue", 
                               "Primed" = "darkgreen", 
                               "Activated" = "red")) +  # Custom colors
  labs(title = "Proportional Morphological State Distribution by SSN_CSN",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 12),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 11)  # Adjust legend title size
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),  
                     expand = c(0,0),  
                     breaks = seq(0, 100, by = 25))  # Proper % labels


summary_table_SSS_SSN
SSS_SSN
library(dplyr)
library(tidyr)
library(broom)

# Step 1: Summary Table for Mean, SD, and Sample Size by State and Group (Separated by Sex)
summary_table_SSS_SSN <- SSS_SSN %>%
  group_by(State, Group) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Group) 

mutate(Percentage = Percentage / sum(Percentage) * 100) %>% # Make the total 100%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 8,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )


summary_table_SSS_SSN

#
# Define pooled SD function
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 1: Compute Mean and SD for SSN and SSS
summary_table_SSS_SSN <- SSS_SSN %>%
  group_by(State, Group) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),  # Compute mean
    SD = sd(Percentage, na.rm = TRUE),      # Compute standard deviation
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Group, values_from = c(Mean, SD))  # Reshape data


# Step 3: Compute Cohen’s d and t-test p-values
comparison_table_SSS_SSN <- summary_table_SSS_SSN %>%
  filter(!is.na(SD_SSN) & !is.na(SD_SSS)) %>%  # Ensure no missing SDs
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(SD_SSN, SD_SSS, SSN_n, SSS_n),  # Compute pooled SD
    Cohens_d = (Mean_SSN - Mean_SSS) / Pooled_SD,         # Compute Cohen's d
    t_test_p_value = wilcox.test(
      SSS_SSN$Percentage[SSS_SSN$State == State & SSS_SSN$Group == "SSS"],
      SSS_SSN$Percentage[SSS_SSN$State == State & SSS_SSN$Group == "SSN"],
      paired = FALSE
    )$p.value
  ) %>%
  ungroup()

# View the final table
print(comparison_table_SSS_SSN)

# Step 4: Display Results
print(comparison_table_SSS_SSN)


#################################################################################################
#CSN vs CNN

CSN <- data_percent %>%
  filter(Group == "CSN")
CNN <- data_percent %>%
  filter(Group == "CNN")
CSN_CNN <- bind_rows(CSN,CNN)
CSN_CNN

summary_table_CSN_CNN <- CSN_CNN %>%
  group_by(Group, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Group) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

summary_table_CSN_CNN %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  
  ggplot(aes(x = Group, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 4, fontface = "bold") +  # White text for labels
  scale_fill_manual(values = c("Surveillant" = "blue", 
                               "Primed" = "darkgreen", 
                               "Activated" = "red")) +  # Custom colors
  labs(title = "Proportional Morphological State Distribution by CSN_CNN",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 9),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 11)  # Adjust legend title size
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),  
                     expand = c(0,0),  
                     breaks = seq(0, 100, by = 25))  # Proper % labels



# Step 1: Summary Table for Mean, SD, and Sample Size by State and Group (Separated by Sex)
summary_table_CSN_CNN <- CSN_CNN %>%
  group_by(State, Group) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 8,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )

# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
comparison_table_CSN_CNN <- summary_table_CSN_CNN %>%
  filter(!is.na(CSN_SD) & !is.na(CNN_SD)) %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(CSN_SD,CNN_SD, CSN_n, CNN_n),
    Cohens_d = (CSN_Mean - CNN_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = CSN_CNN$Percentage[CSN_CNN$State == State & CSN_CNN$Group == "CSN"],
      y =CSN_CNN$Percentage[CSN_CNN$State == State & CSN_CNN$Group == "CNN"],
      paired = FALSE
    )$p.value
  ) %>%
  ungroup()

# Step 4: Display Results
print(comparison_table_CSN_CNN)


#################################################################################################
#################################################################################################
#Fig 4

##SSS
library(ggplot2)
library(dplyr)
library(scales)  # For percentage formatting

# Create the 100% stacked bar chart
library(dplyr)

# Aggregate data to get the mean percentage for each State-Sex combination
SSS_summary <- SSS %>%
  group_by(Sex, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE))  # Take the average percentage

# Now use this aggregated data in ggplot
ggplot(SSS_summary, aes(x = Sex, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Reduce bar width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  
            position = position_stack(vjust = 0.5),  
            color = "white", size = 3, fontface = "bold") +  # Add percentage labels
  scale_fill_manual(values = c("Surveillant" = "blue",  
                               "Primed" = "darkgreen",  
                               "Activated" = "red", 
                               "s" = "pink")) +  
  labs(title = "Proportional Morphological State Distribution by SSS",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 11),  # Adjust y-axis text size
    axis.title.x = element_text(size = 9),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 11)  # Adjust legend title size
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),  
                     expand = c(0,0),  
                     breaks = seq(0, 100, by = 25))  # Proper % labels


SSS

# Step 1: Summary Table for Mean, SD, and Sample Size by State and Group (Separated by Sex)
summary_table_SSS <- SSS%>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 4,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Sex, 
    values_from = c(Mean, SD, n),
    names_glue = "{Sex}_{.value}"
  )
summary_table_SSS

# Modify the sample sizes (Females_n and Males_n) for each state
summary_table_SSS <- summary_table_SSS %>%
  mutate(
    Females_n = case_when(
      State == "Activated" ~ 4,  # Change to your desired value
      State == "Primed" ~ 4,     # Change to your desired value
      State == "Surveillant" ~ 4, # Change to your desired value
      TRUE ~ Females_n  # Keep the current value for other rows
    ),
    Males_n = case_when(
      State == "Activated" ~ 4,  # Change to your desired value
      State == "Primed" ~ 4,     # Change to your desired value
      State == "Surveillant" ~ 4, # Change to your desired value
      TRUE ~ Males_n  # Keep the current value for other rows
    )
  )

summary_table_SSS

# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
comparison_table_SSS <- summary_table_SSS %>%
  filter(!is.na(Males_SD) & !is.na(Females_SD)) %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Males_SD, Females_SD, 4, 4),
    Cohens_d = (Males_Mean - Females_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSS$Percentage[SSS$State == State & SSS$Sex == "Males"],
      y =SSS$Percentage[SSS$State == State & SSS$Sex == "Females"],
      paired = FALSE
    )$p.value
  ) %>%
  ungroup()

# Step 4: Display Results
print(comparison_table_SSS)


##SSS


library(ggplot2)
library(dplyr)
library(scales)  # For percentage formatting

SSN_summary <- SSN %>%
  group_by(Sex, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Sex) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

SSN_summary %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  
  ggplot(aes(x = Sex, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3, fontface = "bold") +  # White text for labels
  scale_fill_manual(values = c("Surveillant" = "blue", 
                               "Primed" = "darkgreen", 
                               "Activated" = "red")) +  # Custom colors
  labs(title = "Proportional Morphological State Distribution by SSN",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 11),  # Adjust y-axis text size
    axis.title.x = element_text(size = 9),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 11)  # Adjust legend title size
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),  
                     expand = c(0,0),  
                     breaks = seq(0, 100, by = 25))  # Proper % labels



# Step 1: Summary Table for Mean, SD, and Sample Size by State and Group (Separated by Sex)
summary_table_SSN <- SSN %>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Sex, 
    values_from = c(Mean, SD, n),
    names_glue = "{Sex}_{.value}"
  )
summary_table_SSN

# Modify the sample sizes (Females_n and Males_n) for each state
summary_table_SSN <- summary_table_SSN %>%
  mutate(
    Females_n = case_when(
      State == "Activated" ~ 4,  # Change to your desired value
      State == "Primed" ~ 4,     # Change to your desired value
      State == "Surveillant" ~ 4, # Change to your desired value
      TRUE ~ Females_n  # Keep the current value for other rows
    ),
    Males_n = case_when(
      State == "Activated" ~ 4,  # Change to your desired value
      State == "Primed" ~ 4,     # Change to your desired value
      State == "Surveillant" ~ 4, # Change to your desired value
      TRUE ~ Males_n  # Keep the current value for other rows
    )
  )

summary_table_SSN

# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
library(dplyr)

# Perform t-test for each state
# Adding Cohen's d and p-value (paired t-test)

comparison_table_SSN <- summary_table_SSN %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Males_SD, Females_SD, Males_n, Females_n),
    Cohens_d = (Males_Mean - Females_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSN$Percentage[SSN$State == State & SSN$Sex == "Males"],
      y = SSN$Percentage[SSN$State == State & SSN$Sex == "Females"],
      paired = FALSE
    )$p.value
  )
comparison_table_SSN

write.csv(comparison_table_SSN, "/Users/Deepika/Desktop/Microglia_anal/includingROI/p_D_SSN.csv")

##CNN

# Summary table for mean, SD, and sample size by state and sex
summary_table_CNN <- CNN %>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 4
  ) %>%
  pivot_wider(
    names_from = Sex, 
    values_from = c(Mean, SD, n),
    names_glue = "{Sex}_{.value}"
  )


# Function to calculate pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}
summary_table_CNN


# Modify the sample sizes (Females_n and Males_n) for each state
summary_table_CNN<- summary_table_CNN %>%
  mutate(
    Females_n = case_when(
      State == "Activated" ~ 4,  # Change to your desired value
      State == "Primed" ~ 4,     # Change to your desired value
      State == "Surveillant" ~ 4, # Change to your desired value
      TRUE ~ Females_n  # Keep the current value for other rows
    ),
    Males_n = case_when(
      State == "Activated" ~ 4,  # Change to your desired value
      State == "Primed" ~ 4,     # Change to your desired value
      State == "Surveillant" ~ 4, # Change to your desired value
      TRUE ~ Males_n  # Keep the current value for other rows
    )
  )

summary_table_CNN

# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
library(dplyr)

# Perform t-test for each state
# Adding Cohen's d and p-value (paired t-test)

comparison_table_CNN <- summary_table_CNN %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Males_SD, Females_SD, Males_n, Females_n),
    Cohens_d = (Males_Mean - Females_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = CNN$Percentage[CNN$State == State & CNN$Sex == "Males"],
      y = CNN$Percentage[CNN$State == State & CNN$Sex == "Females"],
      paired = FALSE
    )$p.value
  )
comparison_table_CNN

library(ggplot2)
library(dplyr)
library(scales)  # For percentage formatting

CNN_summary <- CNN %>%
  group_by(Sex, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Sex) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

CNN_summary %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  
  ggplot(aes(x = Sex, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3, fontface = "bold") +  # White text for labels
  scale_fill_manual(values = c("Surveillant" = "blue", 
                               "Primed" = "darkgreen", 
                               "Activated" = "red")) +  # Custom colors
  labs(title = "Proportional Morphological State Distribution by CNN",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 11),  # Adjust y-axis text size
    axis.title.x = element_text(size = 9),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 11)  # Adjust legend title size
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),  
                     expand = c(0,0),  
                     breaks = seq(0, 100, by = 25))  # Proper % labels







#################################################################################################
#################################################################################################



# SSN vs CSN for females

SSN <- data_percent %>%
  filter(Group == "SSN")
CSN <- data_percent %>%
  filter(Group == "CSN")
SSN_CSN

SSN_CSN <- bind_rows(SSN,CSN)

SSN_CSN_females <- SSN_CSN %>%
  filter(Sex == "Females")


summary_table_SSN_CSN_fem <- SSN_CSN_females %>%
  group_by(Group, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Group) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

summary_table_SSN_CSN_fem %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  mutate(Group = factor(Group, levels = c("SSN", "CSN"))) %>%
  ggplot(aes(x = Group, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 4, fontface = "bold") +  # White text for labels
  scale_fill_manual(values = c("Surveillant" = "blue", 
                               "Primed" = "darkgreen", 
                               "Activated" = "red")) +  # Custom colors
  labs(title = "Proportional Morphological State Distribution by SSN_CSN_fem",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 12),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 11)  # Adjust legend title size
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),  
                     expand = c(0,0),  
                     breaks = seq(0, 100, by = 25))  # Proper % labels

######

library(dplyr)
library(tidyr)
library(broom)

# Step 1: Summary Table for Mean, SD, and Sample Size by State and Group (Separated by Sex)
summary_table_SSN_CSN_fem <- SSN_CSN_females %>% 
  group_by(State, Group) %>% 
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 4,  # Set n to 8 for both groups (SSN and CSN)
    .groups = "drop"
  ) %>% 
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )

# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
comparison_table_SSN_CSN_fem <- summary_table_SSN_CSN_fem %>%
  filter(!is.na(SSN_SD) & !is.na(CSN_SD)) %>%
  rowwise() %>%
  mutate(
    # Use fixed n1 and n2 values (8) for both groups
    Pooled_SD = pooled_sd(SSN_SD, CSN_SD, 8, 8),
    Cohens_d_raw = (SSN_Mean - CSN_Mean) / Pooled_SD,
    Cohens_d = round(Cohens_d_raw, 4),
    t_test_p_value_raw = t.test(
      x = SSN_CSN_females$Percentage[SSN_CSN_females$State == State & SSN_CSN_females$Group == "SSN"],
      y = SSN_CSN_females$Percentage[SSN_CSN_females$State == State & SSN_CSN_females$Group == "CSN"],
      paired = FALSE
    )$p.value,
    t_test_p_value = round(t_test_p_value_raw, 4)
  ) %>%
  ungroup()

comparison_table_SSN_CSN_fem



########SSN_CSN_males

# SSN vs CSN for males

SSN <- data_percent %>%
  filter(Group == "SSN")
CSN <- data_percent %>%
  filter(Group == "CSN")
SSN_CSN

SSN_CSN <- bind_rows(SSN,CSN)

SSN_CSN_males <- SSN_CSN %>%
  filter(Sex == "Males")


summary_table_SSN_CSN_mal <- SSN_CSN_males %>%
  group_by(Group, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Group) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

summary_table_SSN_CSN_mal %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  mutate(Group = factor(Group, levels = c("SSN", "CSN"))) %>%
  ggplot(aes(x = Group, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 4, fontface = "bold") +  # White text for labels
  scale_fill_manual(values = c("Surveillant" = "blue", 
                               "Primed" = "darkgreen", 
                               "Activated" = "red")) +  # Custom colors
  labs(title = "Proportional Morphological State Distribution by SSN_CSN_mal",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 12),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 11)  # Adjust legend title size
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),  
                     expand = c(0,0),  
                     breaks = seq(0, 100, by = 25))  # Proper % labels

######

library(dplyr)
library(tidyr)
library(broom)

# Step 1: Summary Table for Mean, SD, and Sample Size by State and Group (Separated by Sex)
summary_table_SSN_CSN_mal <- SSN_CSN_males %>% 
  group_by(State, Group) %>% 
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 8,  # Set n to 8 for both groups (SSN and CSN)
    .groups = "drop"
  ) %>% 
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )

# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
comparison_table_SSN_CSN_mal <- summary_table_SSN_CSN_mal %>% 
  filter(!is.na(SSN_SD) & !is.na(CSN_SD)) %>% 
  rowwise() %>% 
  mutate(
    # Use fixed n1 and n2 values (8) for both groups
    Pooled_SD = pooled_sd(SSN_SD, CSN_SD, 8, 8),  # Set both n1 and n2 to 8
    Cohens_d = (SSN_Mean - CSN_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSN_CSN$Percentage[SSN_CSN$State == State & SSN_CSN$Group == "SSN"],
      y = SSN_CSN$Percentage[SSN_CSN$State == State & SSN_CSN$Group == "CSN"],
      paired = FALSE
    )$p.value
  ) %>% 
  ungroup()

# View the result
comparison_table_SSN_CSN_mal



#################################################################################################
#################################################################################################


###CSN_CNN for females

CSN <- data_percent %>%
  filter(Group == "CSN")
CNN <- data_percent %>%
  filter(Group == "CNN")
CSN_CNN <- bind_rows(CSN,CNN)
CSN_CNN_females <-CSN_CNN %>%
  filter(Sex == "Females")

summary_table_CSN_CNN_fem <- CSN_CNN_females %>%
  group_by(Group, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Group) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

summary_table_CSN_CNN_fem %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  
  ggplot(aes(x = Group, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 4, fontface = "bold") +  # White text for labels
  scale_fill_manual(values = c("Surveillant" = "blue", 
                               "Primed" = "darkgreen", 
                               "Activated" = "red")) +  # Custom colors
  labs(title = "Proportional Morphological State Distribution by CSN_CNN_female",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 9),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 11)  # Adjust legend title size
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),  
                     expand = c(0,0),  
                     breaks = seq(0, 100, by = 25))  # Proper % labels



# Step 1: Summary Table for Mean, SD, and Sample Size by State and Group (Separated by Sex)
summary_table_CSN_CNN_fem <- CSN_CNN_females %>%
  group_by(State, Group) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 8,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )

# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
comparison_table_CSN_CNN_fem <- summary_table_CSN_CNN_fem %>%
  filter(!is.na(CSN_SD) & !is.na(CNN_SD)) %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(CSN_SD,CNN_SD, CSN_n, CNN_n),
    Cohens_d = (CSN_Mean - CNN_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = CSN_CNN_females$Percentage[CSN_CNN_females$State == State & CSN_CNN_females$Group == "CSN"],
      y =CSN_CNN_females$Percentage[CSN_CNN_females$State == State & CSN_CNN_females$Group == "CNN"],
      paired = FALSE
    )$p.value
  ) %>%
  ungroup()

# Step 4: Display Results
print(comparison_table_CSN_CNN_fem)



###CSN_CNN for males

CSN <- data_percent %>%
  filter(Group == "CSN")
CNN <- data_percent %>%
  filter(Group == "CNN")
CSN_CNN <- bind_rows(CSN,CNN)
CSN_CNN_males <-CSN_CNN %>%
  filter(Sex == "Males")

summary_table_CSN_CNN_mal <- CSN_CNN_males %>%
  group_by(Group, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Group) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

summary_table_CSN_CNN_mal %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  
  ggplot(aes(x = Group, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 4, fontface = "bold") +  # White text for labels
  scale_fill_manual(values = c("Surveillant" = "blue", 
                               "Primed" = "darkgreen", 
                               "Activated" = "red")) +  # Custom colors
  labs(title = "Proportional Morphological State Distribution by CSN_CNN_male",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 9),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 11)  # Adjust legend title size
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),  
                     expand = c(0,0),  
                     breaks = seq(0, 100, by = 25))  # Proper % labels



# Step 1: Summary Table for Mean, SD, and Sample Size by State and Group (Separated by Sex)
summary_table_CSN_CNN_mal <- CSN_CNN_males %>%
  group_by(State, Group) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 8,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )

# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
comparison_table_CSN_CNN_mal <- summary_table_CSN_CNN_mal %>%
  filter(!is.na(CSN_SD) & !is.na(CNN_SD)) %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(CSN_SD,CNN_SD, CSN_n, CNN_n),
    Cohens_d = (CSN_Mean - CNN_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = CSN_CNN_males$Percentage[CSN_CNN_males$State == State & CSN_CNN_males$Group == "CSN"],
      y =CSN_CNN_males$Percentage[CSN_CNN_males$State == State & CSN_CNN_males$Group == "CNN"],
      paired = FALSE
    )$p.value
  ) %>%
  ungroup()

# Step 4: Display Results
print(comparison_table_CSN_CNN_mal)













##########################################################
##########################################################
#SSN vs SSS for females

SSN
SSS
SSN <- data_percent %>%
  filter(Group == "SSN")
SSS <- data_percent %>%
  filter(Group == "SSS")
SSS_SSN <- bind_rows(SSN,SSS)

SSS_SSN_females <- SSS_SSN %>%
  filter (Sex == "Females")

summary_table_SSS_SSN_fem <- SSS_SSN_females %>%
  group_by(State, Group) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 8,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )
summary_table_SSS_SSN_fem
# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
comparison_table_SSS_SSN_fem <- summary_table_SSS_SSN_fem %>% 
  filter(!is.na(SSS_SD) & !is.na(SSN_SD)) %>% 
  rowwise() %>% 
  mutate(
    # Use fixed n1 and n2 values (8) for both groups
    Pooled_SD = pooled_sd(SSS_SD, SSN_SD, 8, 8),  # Set both n1 and n2 to 8
    Cohens_d = (SSS_Mean - SSN_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSS_SSN_females$Percentage[SSS_SSN_females$State == State & SSS_SSN_females$Group == "SSS"],
      y = SSS_SSN_females$Percentage[SSS_SSN_females$State == State & SSS_SSN_females$Group == "SSN"],
      paired = FALSE
    )$p.value
  ) %>% 
  ungroup()

# View the result
comparison_table_SSS_SSN_fem


summary_table_SSS_SSN_fem <- SSS_SSN_females %>%
  group_by(Group, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Group) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

summary_table_SSS_SSN_fem %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  mutate(Group = factor(Group, levels = c("SSN", "SSS"))) %>%
  ggplot(aes(x = Group, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 4, fontface = "bold") +  # White text for labels
  scale_fill_manual(values = c("Surveillant" = "blue", 
                               "Primed" = "darkgreen", 
                               "Activated" = "red")) +  # Custom colors
  labs(title = "Proportional Morphological State Distribution by SSN_CSN",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 12),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 11)  # Adjust legend title size
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),  
                     expand = c(0,0),  
                     breaks = seq(0, 100, by = 25))  # Proper % labels


summary_table_SSS_SSN_fem

#################################################################
###################################################################

#SSS_SSN males

SSS_SSN_males <- SSS_SSN %>%
  filter (Sex == "Males")

summary_table_SSS_SSN_mal <- SSS_SSN_males %>%
  group_by(State, Group) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = 8,
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )
summary_table_SSS_SSN_mal
# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
comparison_table_SSS_SSN_mal <- summary_table_SSS_SSN_mal %>% 
  filter(!is.na(SSS_SD) & !is.na(SSN_SD)) %>% 
  rowwise() %>% 
  mutate(
    # Use fixed n1 and n2 values (8) for both groups
    Pooled_SD = pooled_sd(SSS_SD, SSN_SD, 8, 8),  # Set both n1 and n2 to 8
    Cohens_d = (SSS_Mean - SSN_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSS_SSN_males$Percentage[SSS_SSN_males$State == State & SSS_SSN_males$Group == "SSS"],
      y = SSS_SSN_males$Percentage[SSS_SSN_males$State == State & SSS_SSN_males$Group == "SSN"],
      paired = FALSE
    )$p.value
  ) %>% 
  ungroup()

# View the result
comparison_table_SSS_SSN_mal
comparison_table_SSS_SSN_fem

summary_table_SSS_SSN_mal <- SSS_SSN_males %>%
  group_by(Group, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Group) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

summary_table_SSS_SSN_mal %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  mutate(Group = factor(Group, levels = c("SSN", "SSS"))) %>%
  ggplot(aes(x = Group, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 4, fontface = "bold") +  # White text for labels
  scale_fill_manual(values = c("Surveillant" = "blue", 
                               "Primed" = "darkgreen", 
                               "Activated" = "red")) +  # Custom colors
  labs(title = "Proportional Morphological State Distribution by SSN_CSN males",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),  # Adjust x-axis text size
    axis.text.y = element_text(size = 12),  # Adjust y-axis text size
    axis.title.x = element_text(size = 12),  # Adjust x-axis title size
    axis.title.y = element_text(size = 12),  # Adjust y-axis title size
    plot.title = element_text(size = 10, hjust = 0.5),  # Adjust title size and center it
    legend.text = element_text(size = 11),  # Adjust legend text size
    legend.title = element_text(size = 11)  # Adjust legend title size
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),  
                     expand = c(0,0),  
                     breaks = seq(0, 100, by = 25))  # Proper % labels


summary_table_SSS_SSN_mal

