df1 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/CNN_MALES.csv")
df2<- read.csv("/Users/Deepika/Desktop/Microglia_anal/CNN_females.csv")
df3 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/CSN_FEMALES.csv")
df4 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/CSN_MALES.csv")
df5 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/SSN_FEMALES.csv")
df6 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/SSN_MALES.csv")
df7 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/SSS_FEMALES.csv")
df8 <- read.csv("/Users/Deepika/Desktop/Microglia_anal/SSS_MALES.csv")


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

# 1. Scale Features
scaled_features <- scale(features)
scaled_features <- as.data.frame(scaled_features)
# 2. Compute Correlation Matrix & Plot Heatmap
cor_matrix <- cor(scaled_features, use = "pairwise.complete.obs")
pheatmap(cor_matrix, color = colorRampPalette(c("blue", "white", "red"))(100), main = "Correlation Heatmap")

# 3. Identify Highly Correlated Variables (correlation > 0.9)
high_cor_vars <- findCorrelation(cor_matrix, cutoff = 0.9)
high_cor_names <- colnames(scaled_features)[high_cor_vars]  # Get column names

# 4. Print Highly Correlated Variables
print(high_cor_names)
high_cor_names <- microglia_data
# 5. Remove Highly Correlated Variables
scaled_features_filtered <- scaled_features[, !(colnames(scaled_features) %in% high_cor_names)]
scaled_features_filtered <- as.data.frame(scaled_features_filtered)



# Ensure that all columns are numeric (you can adjust the index if needed)
scaled_features_filtered <- as.data.frame(lapply(scaled_features_filtered, as.numeric))

# Check the structure of the data
str(scaled_features_filtered)

# 1. Compute the correlation matrix
cor_matrix <- cor(scaled_features, use = "pairwise.complete.obs")

# 2. Take the absolute value
abs_cor_matrix <- abs(cor_matrix)

# 3. Set the diagonal to NA to exclude self-correlation
diag(abs_cor_matrix) <- NA

# 4. Compute the mean absolute correlation for each feature
mean_abs_cor <- rowMeans(abs_cor_matrix, na.rm = TRUE)

# 5. Convert to a data frame for easier viewing
mean_abs_cor_df <- data.frame(
  Feature = names(mean_abs_cor),
  Mean_Absolute_Correlation = round(mean_abs_cor, 4)
)

# 6. View the result
print(mean_abs_cor_df)

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
gap_stat <- clusGap(scaled_features_filtered, FUN=kmeans, nstart=25, K.max=10, B=50)

# Plot Gap Statistic
fviz_gap_stat(gap_stat) +
  theme(text = element_text(size = 14))

print(gap_stat)


### 3. Silhouette score##
# Silhouette Method using fviz_nbclust
silhouette_plot <- fviz_nbclust(scaled_features, FUN = kmeans, method = "silhouette") + 
  ggtitle("Silhouette Method for Optimal Clusters")
silhouette_plot


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

# Step 1: Assign custom cluster labels
umap_df$cluster_label <- factor(umap_df$cluster_pam, 
                                levels = c(3, 2, 1), 
                                labels = c("Activated", "Primed", "Surveillant"))

# Step 3: Plot the UMAP with cluster labels and medoids
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

# Load the necessary libraries
library(dplyr)

# Step 1: Check the class of 'scaled_features_filtered'
# If it's not a data frame, convert it into one
if (!is.data.frame(scaled_features_filtered)) {
  scaled_features_filtered <- as.data.frame(scaled_features_filtered)
}

# Step 2: Add the 'cluster_pam' column from 'umap_df' to 'scaled_features_filtered'
scaled_features_filtered$cluster_pam <- umap_df$cluster_pam

# Step 3: Check if the column 'cluster_pam' exists
head(scaled_features_filtered)  # This will print the first few rows to confirm

# Step 4: Compute the mean of each feature per cluster using dplyr
cluster_means <- scaled_features_filtered %>%
  group_by(cluster_pam) %>%
  summarise(across(everything(), mean))

# Step 5: Print the cluster means
print(cluster_means)



# Step 1: Check the class of 'scaled_features_filtered'
# If it's not a data frame, convert it into one
if (!is.data.frame(features)) {
  features <- as.data.frame(features)
}

# Step 2: Add the 'cluster_pam' column from 'umap_df' to 'scaled_features_filtered'
features$cluster_pam <- umap_df$cluster_pam

# Step 3: Check if the column 'cluster_pam' exists
head(features)  # This will print the first few rows to confirm

# Step 4: Compute the mean of each feature per cluster using dplyr
cluster_means <- features %>%
  group_by(cluster_pam) %>%
  summarise(across(everything(), list(mean = mean, sd = sd)))
cluster_means
library(dplyr)

# Step 4: Compute the mean and standard error of each feature per cluster
cluster_stats <- scaled_features_filtered %>%
  group_by(cluster_pam) %>%
  summarise(across(everything(), list(mean = mean, sd = sd)))

# Step 5: Print the cluster statistics
print(cluster_stats)

# Step 5: Print the cluster means
print(cluster_means)

#write.csv(cluster_means,  "/Users/Deepika/Desktop/Microglia_anal/newfigs/cluster_means.csv")

####Cell counts####

cluster_counts <- scaled_features_filtered %>%
  group_by(cluster_pam) %>%
  tally()  # Tally the number of rows for each cluster

# Step 4: Print the total cell count for each cluster
print(cluster_counts)


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

########
##########g
###############

df$X <- 1:nrow(df)
library(dplyr)

###subsetting the data

CNN_females <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Females",
         Group == "CNN")

print(CNN_females)

#write.csv(CNN_females, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/CNN_female.csv")


CSN_females <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Females",
         Group == "CSN")
print(CSN_females)
#write.csv(CSN_females, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/CSN_female.csv")


SSN_females <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Females",
         Group == "SSN")
print(SSN_females)
#write.csv(SSN_females, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/SSN_female.csv")


SSS_females <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Females",
         Group == "SSS")
print(SSS_females)
#write.csv(SSS_females, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/SSS_female.csv")

#####Subsetting for males

CNN_males <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Males",
         Group == "CNN")
print(CNN_males)

#write.csv(CNN_males, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/CNN_male.csv")

CSN_males <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Males",
         Group == "CSN")
print(CSN_males)

#write.csv(CSN_males, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/CSN_male.csv")

SSN_males <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Males",
         Group == "SSN")
print(SSN_males)

#write.csv(SSN_males, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/SSN_male.csv")

SSS_males <- df %>%
  filter(Rat_num %in%  c("R1", "R2", "R3", "R4"),
         Sex == "Males",
         Group == "SSS")
print(SSS_males)

#write.csv(SSS_males, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/SSS_male.csv")


###########Calculatepercentages of cell types#########
########################################################

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
####################
library(tidyplots)

SSS %>%
  tidyplot(x = State, y = Percentage, color = Sex, fill = Sex) %>% # Map both color and fill to Sex
  add_mean_bar(alpha = 0.4) %>%
  add_sem_errorbar() %>%
  add_data_points_beeswarm() %>%
  adjust_y_axis(limits = c(0, 100)) %>%
  remove_x_axis_title() %>%
  reorder_x_axis_labels("Surveillant", "Primed", "Transitional", "Amoeboid") %>%
  add(ggplot2::labs(
    title = "Morphological State Percentages by Sex",
    x = "Morphological State",
    y = "Percentage of Cells (%)"
  )) %>%
  add(
    ggplot2::scale_color_manual(
      values = c("Male" = "blue", "Female" = "red") # Blue for males, red for females
    )
  ) %>%
  add(
    ggplot2::scale_fill_manual(
      values = c("Male" = "blue", "Female" = "red") # Ensure bars match the points
    )
  ) %>%
  add(
    ggplot2::theme(
      aspect.ratio = 1, # Adjust the height-to-width ratio
      
    )
  )

###########
###########
############

# Load necessary libraries
library(dplyr)


# Summary table for mean, SD, and sample size by state and sex
summary_stats <- SSS %>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n()
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

# Adding Cohen's d and p-value (paired t-test)
comparison_table_SSS <- summary_stats %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Male_SD, Female_SD, Male_n, Female_n),
    Cohens_d = (Male_Mean - Female_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSS$Percentage[SSS$State == State & SSS$Sex == "Male"],
      y = SSS$Percentage[SSS$State == State & SSS$Sex == "Female"],
      paired = FALSE
    )$p.value
  )

# Print the final table
print(comparison_table_SSS)
#write.csv(comparison_table_SSS, "/Users/Deepika/Desktop/Microglia_anal/newfigs/D_P_SSS.csv")
###############
###Sex wise SSN###
# Step 1: Add a Sex column to the male and female datasets
state_counts_SSN_mal$Sex <- "Male"
state_counts_SSN_fem$Sex <- "Female"

# Step 2: Combine the datasets for both sexes
SSN <- bind_rows(state_counts_SSN_mal, state_counts_SSN_fem)

# Step 3: Calculate the total number of cells for each rat within each sex
total_cells_per_sex_rat <- SSN %>%
  group_by(Rat_num, Sex) %>%
  summarise(Total_Sex_Cells = sum(Cell_Count), .groups = "drop")

# Step 4: Merge the total cells by rat and sex back into the combined dataset
SSN <- SSN %>%
  left_join(total_cells_per_sex_rat, by = c("Rat_num", "Sex"))


# View the final dataset
print(SSN)

#write.csv(SSN, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/SSNmalesvsfemales.csv")

library(tidyplots)

SSN %>%
  tidyplot(x = State, y = Percentage, color = Sex, fill = Sex) %>% # Map both color and fill to Sex
  add_mean_bar(alpha = 0.4) %>%
  add_sem_errorbar() %>%
  add_data_points_beeswarm() %>%
  adjust_y_axis(limits = c(0, 100)) %>%
  remove_x_axis_title() %>%
  reorder_x_axis_labels("Surveillant", "Primed", "Transitional", "Amoeboid") %>%
  add(ggplot2::labs(
    title = "Morphological State Percentages by Sex",
    x = "Morphological State",
    y = "Percentage of Cells (%)"
  )) %>%
  add(
    ggplot2::scale_color_manual(
      values = c("Male" = "blue", "Female" = "red") # Blue for males, red for females
    )
  ) %>%
  add(
    ggplot2::scale_fill_manual(
      values = c("Male" = "blue", "Female" = "red") # Ensure bars match the points
    )
  ) %>%
  add(
    ggplot2::theme(
      aspect.ratio = 1, # Adjust the height-to-width ratio
      
    )
  )

# Load necessary libraries
library(dplyr)


# Summary table for mean, SD, and sample size by state and sex
summary_stats <- SSN %>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n()
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

# Adding Cohen's d and p-value (paired t-test)
comparison_table_SSN <- summary_stats %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Male_SD, Female_SD, Male_n, Female_n),
    Cohens_d = (Male_Mean - Female_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSN$Percentage[SSN$State == State & SSN$Sex == "Male"],
      y = SSN$Percentage[SSN$State == State & SSN$Sex == "Female"],
      paired = FALSE
    )$p.value
  )

# Print the final table
print(comparison_table_SSN)
#write.csv(comparison_table_SSN, "/Users/Deepika/Desktop/Microglia_anal/newfigs/d_P_SSN.csv")

###Sex wise CSN######
###################

# Step 1: Add a Sex column to the male and female datasets
state_counts_CSN_mal$Sex <- "Male"
state_counts_CSN_fem$Sex <- "Female"

# Step 2: Combine the datasets for both sexes
CSN <- bind_rows(state_counts_CSN_mal, state_counts_CSN_fem)

# Step 3: Calculate the total number of cells for each rat within each sex
total_cells_per_sex_rat <- CSN%>%
  group_by(Rat_num, Sex) %>%
  summarise(Total_Sex_Cells = sum(Cell_Count), .groups = "drop")

# Step 4: Merge the total cells by rat and sex back into the combined dataset
CSN <- CSN %>%
  left_join(total_cells_per_sex_rat, by = c("Rat_num", "Sex"))
# View the final dataset
print(CSN)

#write.csv(CSN, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/CSN_malesvsfemales.csv")

library(tidyplots)

CSN %>%
  tidyplot(x = State, y = Percentage, color = Sex, fill = Sex) %>% # Map both color and fill to Sex
  add_mean_bar(alpha = 0.4) %>%
  add_sem_errorbar() %>%
  add_data_points_beeswarm() %>%
  adjust_y_axis(limits = c(0, 100)) %>%
  remove_x_axis_title() %>%
  reorder_x_axis_labels("Surveillant", "Primed", "Transitional", "Amoeboid") %>%
  add(ggplot2::labs(
    title = "Morphological State Percentages by Sex",
    x = "Morphological State",
    y = "Percentage of Cells (%)"
  )) %>%
  add(
    ggplot2::scale_color_manual(
      values = c("Male" = "blue", "Female" = "red") # Blue for males, red for females
    )
  ) %>%
  add(
    ggplot2::scale_fill_manual(
      values = c("Male" = "blue", "Female" = "red") # Ensure bars match the points
    )
  ) %>%
  add(
    ggplot2::theme(
      aspect.ratio = 1, # Adjust the height-to-width ratio
      
    )
  )

# Load necessary libraries
library(dplyr)


# Summary table for mean, SD, and sample size by state and sex
summary_stats <- CSN %>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n()
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

# Adding Cohen's d and p-value (paired t-test)
comparison_table_CSN <- summary_stats %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Male_SD, Female_SD, Male_n, Female_n),
    Cohens_d = (Male_Mean - Female_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = CSN$Percentage[CSN$State == State & CSN$Sex == "Male"],
      y = CSN$Percentage[CSN$State == State & CSN$Sex == "Female"],
      paired = FALSE
    )$p.value
  )

# Print the final table
print(comparison_table_CSN)

#write.csv(comparison_table_CSN, "/Users/Deepika/Desktop/Microglia_anal/newfigs/d_P_CSN.csv")

####Sex wise CNN#######
########################
state_counts_CNN_mal$Group <- "CNN"
# Step 1: Add a Sex column to the male and female datasets
state_counts_CNN_mal$Sex <- "Male"
state_counts_CNN_fem$Sex <- "Female"

# Step 2: Combine the datasets for both sexes
CNN <- bind_rows(state_counts_CNN_mal, state_counts_CNN_fem)

# Step 3: Calculate the total number of cells for each rat within each sex
total_cells_per_sex_rat <- CNN %>%
  group_by(Rat_num, Sex) %>%
  summarise(Total_Sex_Cells = sum(Cell_Count), .groups = "drop")

# Step 4: Merge the total cells by rat and sex back into the combined dataset
CNN <- CNN %>%
  left_join(total_cells_per_sex_rat, by = c("Rat_num", "Sex"))


# View the final dataset
print(CNN)

#write.csv(CNN, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/sexdiff_CNN.csv")

library(tidyplots)

CNN %>%
  tidyplot(x = State, y = Percentage, color = Sex, fill = Sex) %>% # Map both color and fill to Sex
  add_mean_bar(alpha = 0.4) %>%
  add_sem_errorbar() %>%
  add_data_points_beeswarm() %>%
  adjust_y_axis(limits = c(0, 100)) %>%
  remove_x_axis_title() %>%
  reorder_x_axis_labels("Surveillant", "Primed", "Transitional", "Amoeboid") %>%
  add(ggplot2::labs(
    title = "Morphological State Percentages by Sex",
    x = "Morphological State",
    y = "Percentage of Cells (%)"
  )) %>%
  add(
    ggplot2::scale_color_manual(
      values = c("Male" = "blue", "Female" = "red") # Blue for males, red for females
    )
  ) %>%
  add(
    ggplot2::scale_fill_manual(
      values = c("Male" = "blue", "Female" = "red") # Ensure bars match the points
    )
  ) %>%
  add(
    ggplot2::theme(
      aspect.ratio = 1, # Adjust the height-to-width ratio
      
    )
  )

# Load necessary libraries
library(dplyr)


# Summary table for mean, SD, and sample size by state and sex
summary_stats <- CNN %>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n()
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

# Adding Cohen's d and p-value (paired t-test)
comparison_table_CNN <- summary_stats %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Male_SD, Female_SD, Male_n, Female_n),
    Cohens_d = (Male_Mean - Female_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = CNN$Percentage[CNN$State == State & CNN$Sex == "Male"],
      y = CNN$Percentage[CNN$State == State & CNN$Sex == "Female"],
      paired = FALSE
    )$p.value
  )

# Print the final table
print(comparison_table_CNN)
#write.csv(comparison_table_CNN, "/Users/Deepika/Desktop/Microglia_anal/newfigs/d_P_CNN.csv")

#####################################
###injection differences#############
#####################################

##############
################
################
##CSN+SSN

# Load necessary libraries
library(dplyr)

# Add a "Group" column to distinguish datasets
CSN$Group <- "CSN"
SSN$Group <- "SSN"

# Combine the datasets
SSN_CSN <- bind_rows(CSN, SSN)
SSN_CSN
# Summary table for mean, SD, and sample size by state and sex
summary_stats <- SSN_CSN %>%
  group_by(State, Group) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n()
  ) %>%
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )

# Function to calculate pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Adding Cohen's d and p-value (paired t-test)
comparison_table_SSN_CSN <- summary_stats %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(SSN_SD, CSN_SD, SSN_n, CSN_n),
    Cohens_d = (SSN_Mean - CSN_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSN_CSN$Percentage[SSN_CSN$State == State & SSN_CSN$Group == "SSN"],
      y = SSN_CSN$Percentage[SSN_CSN$State == State & SSN_CSN$Group == "CSN"],
      paired = FALSE
    )$p.value
  )
comparison_table_SSN_CSN
#write.csv(comparison_table_SSN_CSN, "/Users/Deepika/Desktop/Microglia_anal/newfigs/SSNvsCSN.csv")

################
##SSS+SSN

# Load necessary libraries
library(dplyr)

# Add a "Group" column to distinguish datasets
SSS$Group <- "SSS"
SSN$Group <- "SSN"

# Combine the datasets
SSS_SSN<- bind_rows(SSS, SSN)
SSS_SSN
# Summary table for mean, SD, and sample size by state and sex
summary_stats <- SSS_SSN %>%
  group_by(State, Group) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n()
  ) %>%
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )

# Function to calculate pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Adding Cohen's d and p-value (paired t-test)
comparison_table_SSS_SSN <- summary_stats %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(SSS_SD, SSN_SD, SSS_n, SSN_n),
    Cohens_d = (SSS_Mean - SSN_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSS_SSN$Percentage[SSS_SSN$State == State & SSS_SSN$Group == "SSS"],
      y = SSS_SSN$Percentage[SSS_SSN$State == State & SSS_SSN$Group == "SSN"],
      paired = FALSE
    )$p.value
  )
comparison_table_SSS_SSN
#write.csv(comparison_table_SSS_SSN, "/Users/Deepika/Desktop/Microglia_anal/newfigs/SSSvsSSN.csv")

###Tiplot
library(tidyplots)
SSS_SSN %>%
  tidyplot(x = State, y = Percentage, color = Group, fill = Group) %>% # Map both color and fill to Sex
  add_mean_bar(alpha = 0.4) %>%
  add_sem_errorbar() %>%
  add_data_points_beeswarm() %>%
  adjust_y_axis(limits = c(0, 100)) %>%
  remove_x_axis_title() %>%
  reorder_x_axis_labels("Surveillant", "Primed", "Transitional", "Amoeboid") %>%
  add(ggplot2::labs(
    title = "Morphological State Percentages by Sex",
    x = "Morphological State",
    y = "Percentage of Cells (%)"
  )) %>%
  add(
    ggplot2::scale_color_manual(
      values = c("SSN" = "red", "SSS" = "black") 
    )
  ) %>%
  add(
    ggplot2::scale_fill_manual(
      values = c("SSN" = "red", "SSS" = "black") # Ensure bars match the points
    )
  ) %>%
  add(
    ggplot2::theme(
      aspect.ratio = 1, # Adjust the height-to-width ratio
      
    )
  )

library(tidyplots)
library(ggplot2)

SSS_SSN %>%
  tidyplot(x = State, y = Percentage, fill = Group) %>%  # Fill bars by Group
  add_column_bar(position = "stack") %>%  # Stack bars to represent parts of whole
  adjust_y_axis(limits = c(0, 100), labels = scales::percent_format(scale = 1)) %>%  # Convert y-axis to percent format
  remove_x_axis_title() %>% 
  reorder_x_axis_labels("Surveillant", "Primed", "Transitional", "Amoeboid") %>%
  add(ggplot2::labs(
    title = "Morphological State Composition by Sex",
    x = "Morphological State",
    y = "Percentage of Cells (%)"
  )) %>%
  add(
    ggplot2::scale_fill_manual(
      values = c("SSN" = "red", "SSS" = "black") # Ensure colors remain distinct
    )
  ) %>%
  add(
    ggplot2::theme(
      aspect.ratio = 1,  # Keep square aspect ratio
      legend.position = "right" # Move legend to the side for clarity
    )
  )



##CSNvsCNN



# Load necessary libraries
library(dplyr)

# Add a "Group" column to distinguish datasets
CSN$Group <- "CSN"
CNN$Group <- "CNN"

# Combine the datasets
CSN_CNN<- bind_rows(CSN, CNN)
CSN_CNN
# Summary table for mean, SD, and sample size by state and sex
summary_stats <- CSN_CNN %>%
  group_by(State, Group) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n()
  ) %>%
  pivot_wider(
    names_from = Group, 
    values_from = c(Mean, SD, n),
    names_glue = "{Group}_{.value}"
  )

# Function to calculate pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Adding Cohen's d and p-value (paired t-test)
comparison_table_CSN_CNN <- summary_stats %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(CSN_SD, CNN_SD, CSN_n, CNN_n),
    Cohens_d = (CSN_Mean - CNN_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = CSN_CNN$Percentage[CSN_CNN$State == State & CSN_CNN$Group == "CSN"],
      y = CSN_CNN$Percentage[CSN_CNN$State == State & CSN_CNN$Group == "CNN"],
      paired = FALSE
    )$p.value
  )
comparison_table_CSN_CNN
write.csv(comparison_table_CSN_CNN, "/Users/Deepika/Desktop/Microglia_anal/newfigs/CSNvsCNN.csv")

###Tiplot

library(tidyplots)
CSN_CNN %>%
  tidyplot(x = State, y = Percentage, color = Group, fill = Group) %>% # Map both color and fill to Sex
  add_mean_bar(alpha = 0.4) %>%
  add_sem_errorbar() %>%
  add_data_points_beeswarm() %>%
  adjust_y_axis(limits = c(0, 100)) %>%
  remove_x_axis_title() %>%
  reorder_x_axis_labels("Surveillant", "Primed", "Transitional", "Amoeboid") %>%
  add(ggplot2::labs(
    title = "Morphological State Percentages by Sex",
    x = "Morphological State",
    y = "Percentage of Cells (%)"
  )) %>%
  add(
    ggplot2::scale_color_manual(
      values = c("CNN" = "red", "CSN" = "black") 
    )
  ) %>%
  add(
    ggplot2::scale_fill_manual(
      values = c("CNN" = "red", "CSN" = "black") # Ensure bars match the points
    )
  ) %>%
  add(
    ggplot2::theme(
      aspect.ratio = 1, # Adjust the height-to-width ratio
      
    )
  )
##################
#################
################

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

write.csv(data_percent, "/Users/Deepika/Desktop/Microglia_anal/includingROI/data_ROI_all.csv")
# 3-way Anova Sex difference

# Perform a three-way ANOVA
anova_result <- aov( Percentage ~ Sex * Group * State, data = data_percent)

# Summary of the ANOVA
summary(anova_result)

# Post-hoc analysis for multiple comparisons
TukeyHSD(anova_result)

library(ggplot2)


data_percent <- data_percent %>%
  mutate(Group_Rat_num = paste(Group, Rat_num, sep = "_"))
data_percent

data_percent_group <- data_percent %>% 
  group_by(Group_Rat_num, State, Sex) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE), .groups = "drop")

library(dplyr)

# Assuming your data is in a dataframe called `df`
# Group by 'State' and 'Sex' and calculate the mean percentage for each group
mean_percentage <- data_percent_group %>%
  group_by(State, Sex) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE))  # Calculate mean percentage for each group

# View the result
mean_percentage


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
data_percent_group

library(dplyr)
library(ggplot2)
library(tidyverse)
library(effsize)




# Perform t-test and compute Cohen's d for each state
sex_comparison <- data_percent_group %>%
  group_by(State) %>%
  summarise(
    p_value = t.test(Percentage ~ Sex)$p.value,  # Extract p-value directly from the t-test result
    cohen_d_value = cohen.d(Percentage ~ Sex, data = .)$estimate  # Extract Cohen's d from cohen.d result
  )

# View results
sex_comparison


write.csv(comparison_table_sexdiff
          , "/Users/Deepika/Desktop/Microglia_anal/newfigs/sexdiff.csv")




data_percent%>%
  tidyplot(x = State, y = Percentage, color = Sex, fill = Sex) %>%
  add(ggplot2::facet_wrap(~ Group)) %>%
  add_data_points_beeswarm() %>%
  add_mean_bar(alpha = 0.4) %>%
  add_sem_errorbar() %>%
  add_data_points_beeswarm() %>%
  add(ggplot2::labs(
    title = "Interaction plot Sex*Group*State",
    x = "Morphological State",
    y = "Percentage of Cells (%)"
    
  ))

# Cohen's f for a factor in ANOVA
f_square <- sqrt((F_value * df_effect) / (df_effect * F_value + df_residual))

# For significant interactions, you can run post-hoc tests:
library(emmeans)
emm <- emmeans(anova_result, pairwise ~ Sex * Group * State)
summary(emm)


###SSS
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
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),  # Adjust x-axis text size
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
summary_table_SSS <- SSS %>%
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
summary_table_SSS

# Step 2: Function to Calculate Pooled SD
pooled_sd <- function(sd1, sd2, n1, n2) {
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

# Step 3: Compute Cohen's d and p-values for **Sex Differences** (Males vs. Females)
library(dplyr)

# Perform t-test for each state
# Adding Cohen's d and p-value (paired t-test)

comparison_table_SSS <- summary_table_SSS %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Male_SD, Female_SD, Male_n, Female_n),
    Cohens_d = (Male_Mean - Female_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSS$Percentage[SSS$State == State & SSS$Sex == "Male"],
      y = SSS$Percentage[SSS$State == State & SSS$Sex == "Female"],
      paired = FALSE
    )$p.value
  )


comparison_table_SSS
#####Parts of whole bar graphs for the males vs females
##SSN##

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
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),  # Adjust x-axis text size
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


SSN_summary
# Load necessary libraries
library(ggplot2)
library(dplyr)


# Example dataset
SSN_summary <- SSN %>%
  group_by(Sex, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE))  # Take the average percentage


library(ggplot2)

ggplot(SSN, aes(x = Sex, y = Percentage, fill = State)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Transparent boxplots, hide outliers
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.6) +  # Individual points
  stat_summary(fun = mean, geom = "point", shape = 4, size = 4, color = "black") +  # Mean as "X"
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") +  # SE error bars
  facet_wrap(~ State) +  # Separate plots for each state
  labs(title = "Cell Population Distribution by Sex",
       x = "Sex",
       y = "Percentage",
       fill = "State") +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove grid lines
        axis.line = element_line(color = "black", linewidth = 1))  # Solid axes


library(ggplot2)
library(dplyr)

# Reorder the 'State' factor levels
SSN$State <- factor(SSN$State, levels = c("Surveillant", "Primed", "Activated"))

# Create the ggplot with no gridlines and individual points
ggplot(SSN, aes(x = Sex, y = Percentage, fill = State)) +
  # Bar plot to show the mean percentage of each group
  geom_bar(stat = "identity", width = 0.6, alpha = 0.7) +
  
  # Add individual data points on top of the bars with jitter
  geom_jitter(position = position_jitter(width = 0.2), size = 2, aes(color = State)) +
  
  # Separate facets for each State
  facet_wrap(~ State) +
  
  # Customize plot labels and titles
  labs(title = "Mean Cell Population Percentage by Sex",
       x = "Sex",
       y = "Mean Percentage",
       fill = "State") +
  
  # Apply a minimal theme, removing grid lines and customizing axes
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove gridlines
        axis.line = element_line(color = "black", linewidth = 0.5)) +  # Solid axes
  
  # Set y-axis range from 0 to 100
  ylim(0, 100)





library(ggplot2)
library(dplyr)

ggplot(SSN, aes(x = Sex, y = Percentage, color = State)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +  # Individual data points
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.2, color = "black", linewidth = 1.2) +  # SE bars
  stat_summary(fun = mean, geom = "point", 
               size = 4, shape = 21, fill = "black") +  # Mean points
  facet_wrap(~ State) +  # Separate by cell population state
  labs(title = "Cell Population Distribution by Sex",
       x = "Sex",
       y = "Percentage",
       color = "State") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 1))






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
    Pooled_SD = pooled_sd(Male_SD, Female_SD, Male_n, Female_n),
    Cohens_d = (Male_Mean - Female_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSN$Percentage[SSN$State == State & SSN$Sex == "Male"],
      y = SSN$Percentage[SSN$State == State & SSN$Sex == "Female"],
      paired = FALSE
    )$p.value
  )


comparison_table_SSN


##CSN###
CSN
CSN_summary <- CSN %>%
  group_by(Sex, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Sex) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

CSN_summary %>%
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
  labs(title = "Proportional Morphological State Distribution by CSN",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),  # Adjust x-axis text size
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

####CNN##

library(ggplot2)
library(dplyr)
library(scales)  # For percentage formatting
CNN_summary
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
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),  # Adjust x-axis text size
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



data_percent

CNN <- CNN %>%
  mutate(Group = "CNN")
#######
#####CSN+SSN###



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


# Step 4: Display Results
print(comparison_table_SSN_CSN)

write.csv(comparison_table_SSN_CSN, "/Users/Deepika/Desktop/Microglia_anal/includingROI/d_P_SSN_CSN.csv")



summary_table_SSN_CSN <- SSN_CSN %>%
  group_by(Group, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Group) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

summary_table_SSN_CSN %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  
  ggplot(aes(x = Group, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3, fontface = "bold") +  # White text for labels
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
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),  # Adjust x-axis text size
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






####SSS+SSN###

summary_table_SSS_SSN <- SSS_SSN %>%
  group_by(Group, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Group) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%

summary_table_SSS_SSN %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  
  ggplot(aes(x = Group, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3, fontface = "bold") +  # White text for labels
  scale_fill_manual(values = c("Surveillant" = "blue", 
                               "Primed" = "darkgreen", 
                               "Activated" = "red")) +  # Custom colors
  labs(title = "Proportional Morphological State Distribution by SSS_SSN",  
       x = "Sex",  
       y = "Proportion of Cells (%)",  
       fill = "State") +  
  theme_minimal() +  
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black"),  # Solid axes
    aspect.ratio = 1,  # Set aspect ratio
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),  # Adjust x-axis text size
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


library(dplyr)
library(tidyr)
library(broom)

# Step 1: Summary Table for Mean, SD, and Sample Size by State and Group (Separated by Sex)
summary_table_SSS_SSN <- SSS_SSN %>%
  group_by(State, Group) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n(),
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
comparison_table_SSS_SSN <- summary_table_SSS_SSN %>%
  filter(!is.na(SSN_SD) & !is.na(SSN_SD)) %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(SSN_SD, SSS_SD, SSN_n, SSS_n),
    Cohens_d = (SSN_Mean - SSS_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSS_SSN$Percentage[SSS_SSN$State == State & SSS_SSN$Group == "SSS"],
      y = SSS_SSN$Percentage[SSS_SSN$State == State & SSS_SSN$Group == "SSN"],
      paired = FALSE
    )$p.value
  ) %>%
  ungroup()

# Step 4: Display Results
print(comparison_table_SSS_SSN)

write.csv(comparison_table_SSN_CSN, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/d_P_SSS_SSN.csv")
#



###CSN_CNN###
CSN_CNN
library(ggplot2)

library(dplyr)

summary_table_CSN_CNN <- CSN_CNN %>%
  group_by(Group, State) %>%
  summarise(Percentage = mean(Percentage, na.rm = TRUE)) %>%  # Average percentage
  # Convert Percentage to proportion of the total for each 'Sex'
  group_by(Sex) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100)  # Make the total 100%
summary_table_CSN_CNN

summary_table_CSN_CNN %>%
  # Reorder the levels of the 'State' factor
  mutate(State = factor(State, levels = c("Activated", "Primed", "Surveillant"))) %>%
  
  ggplot(aes(x = Group, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Bar chart with reduced width
  geom_text(aes(label = paste0(round(Percentage, 1), "%")),  # Add percentage labels
            position = position_stack(vjust = 0.5), 
            color = "white", size = 3, fontface = "bold") +  # White text for labels
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
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),  # Adjust x-axis text size
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
summary_table_CSN_CNN <- CSN_CNN %>%
  group_by(State, Group) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n(),
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

###############
############
#################
#Sex difference#
# Step 2: Aggregate total percentages for each state within each sex
sex_diff <- data_percent %>%
  group_by(Sex, State) %>%
  summarise(Total_Cells = sum(Cell_Count), .groups = "drop") %>%
  group_by(Sex) %>%
  mutate(Percentage = (Total_Cells / sum(Total_Cells)) * 100)  # Normalize within sex to 100%

# Step 3: Create a stacked bar chart
ggplot(sex_diff, aes(x = Sex, y = Percentage, fill = State)) +
  geom_bar(stat = "identity", width = 0.5) +  # Reduce bar width
  scale_fill_manual(values = c("Activated" = "red", "Primed" = "green", "Surveillant" = "blue")) +  
  labs(title = "Proportional Morphological State Distribution by SSS_SSN",
       x = "Sex",
       y = "Proportion of Cells (%)",
       fill = "State") +  
  theme_minimal() +  
  theme(panel.grid = element_blank(),  # Remove grid lines
        axis.line = element_line(color = "black", linewidth = 1)) +  # Solid axes
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0,0), breaks = seq(0, 100, by = 25))  # Proper % labels
###
#######
##########
#########
###########
#################


# Summary table for mean, SD, and sample size by state and sex
summary_stats <- SSS %>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n()
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

# Adding Cohen's d and p-value (paired t-test)
comparison_table_SSS <- summary_stats %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Male_SD, Female_SD, Male_n, Female_n),
    Cohens_d = (Male_Mean - Female_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSS$Percentage[SSS$State == State & SSS$Sex == "Male"],
      y = SSS$Percentage[SSS$State == State & SSS$Sex == "Female"],
      paired = FALSE
    )$p.value
  )

# Print the final table
print(comparison_table_SSS)
write.csv(comparison_table_SSS, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/d_P_SSS.csv")


#########
#########

# Summary table for mean, SD, and sample size by state and sex
summary_stats <- SSN %>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n()
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

# Adding Cohen's d and p-value (paired t-test)
comparison_table_SSN <- summary_stats %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Male_SD, Female_SD, Male_n, Female_n),
    Cohens_d = (Male_Mean - Female_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = SSN$Percentage[SSN$State == State & SSN$Sex == "Male"],
      y = SSN$Percentage[SSN$State == State & SSN$Sex == "Female"],
      paired = FALSE
    )$p.value
  )

# Print the final table
print(comparison_table_SSN)
write.csv(comparison_table_SSN, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/d_P_SSN.csv")


######
######
########
###########

# Summary table for mean, SD, and sample size by state and sex
summary_stats <- CSN %>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n()
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

# Adding Cohen's d and p-value (paired t-test)
comparison_table_CSN <- summary_stats %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Male_SD, Female_SD, Male_n, Female_n),
    Cohens_d = (Male_Mean - Female_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = CSN$Percentage[CSN$State == State & CSN$Sex == "Male"],
      y = CSN$Percentage[CSN$State == State & CSN$Sex == "Female"],
      paired = FALSE
    )$p.value
  )

# Print the final table
print(comparison_table_CSN)
write.csv(comparison_table_CSN, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/d_P_CSN.csv")

#############
#############
#############
#CNN
# Summary table for mean, SD, and sample size by state and sex
summary_stats <- CNN %>%
  group_by(State, Sex) %>%
  summarise(
    Mean = mean(Percentage, na.rm = TRUE),
    SD = sd(Percentage, na.rm = TRUE),
    n = n()
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

# Adding Cohen's d and p-value (paired t-test)
comparison_table_CNN <- summary_stats %>%
  rowwise() %>%
  mutate(
    Pooled_SD = pooled_sd(Male_SD, Female_SD, Male_n, Female_n),
    Cohens_d = (Male_Mean - Female_Mean) / Pooled_SD,
    t_test_p_value = t.test(
      x = CNN$Percentage[CNN$State == State & CNN$Sex == "Male"],
      y = CNN$Percentage[CNN$State == State & CNN$Sex == "Female"],
      paired = FALSE
    )$p.value
  )

# Print the final table
print(comparison_table_CNN)
write.csv(comparison_table_CNN, "/Users/Deepika/Desktop/Microglia_anal/Graphsandtables/d_P_CNN.csv")

##########
##############
##################
write.csv(cluster_stats, "/Users/Deepika/Desktop/Microglia_anal/cluster_stats.csv")
write.csv(cluster_means, "/Users/Deepika/Desktop/Microglia_anal/cluster_means_features.csv")

# Data: Means, SDs, and sample sizes for each feature
# Means and SDs for the features
# Create a data frame with the given data
# Create a data frame from your new dataset
data <- data.frame(
  cluster_name = c("Activated", "Primed", "Surveillant"),
  feature = c("Cell Territorial Vol (μm³)", "Cell Volume (μm³)", "Ramification index",
              "Number of endpoints", "Number of branchpoints", "Average Branch length (μm)",
              "Max Branch length (μm)", "Min Branch length (μm)"),
  CellTerritoryVol_um3_mean = c(55651.91396, 61994.92032, 82803.43109),
  CellTerritoryVol_um3_sd = c(60361.25523, 64084.57902, 79887.5512),
  CellVolumes_mean = c(9157.99278, 10346.95189, 13374.0155),
  CellVolumes_sd = c(8447.214827, 9036.293386, 9840.951477),
  RamificationIndex_mean = c(5.656583585, 5.800060495, 5.911876459),
  RamificationIndex_sd = c(3.30079233, 3.152904238, 3.003499307),
  NumOfEndpoints_mean = c(9.78700361, 11.55237515, 14.62015504),
  NumOfEndpoints_sd = c(10.12859942, 13.50229619, 15.51625399),
  NumOfBranchpoints_mean = c(7.241877256, 8.425700365, 10.79069767),
  NumOfBranchpoints_sd = c(7.593092357, 9.838409088, 11.34668996),
  AvgBranchLength_mean = c(63.76876485, 70.39982501, 79.04487162),
  AvgBranchLength_sd = c(40.75224258, 43.18174435, 44.62292084),
  MaxBranchLength_mean = c(122.616716, 137.020853, 159.3826151),
  MaxBranchLength_sd = c(90.52319411, 104.8782026, 130.4129168),
  MinBranchLength_mean = c(19.11850527, 20.38407964, 21.49710823),
  MinBranchLength_sd = c(15.94972, 15.646236, 15.3257215),
  Surveillant_n = 1592,
  Primed_n = 1647,
  Activated_n = 277
)
# Print the data frame
print(data)
# Cohen's d Calculation Function
cohens_d <- function(m1, m2, sd1, sd2, n1, n2) {
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  d <- (m1 - m2) / pooled_sd
  return(d)
}

# Loop through each feature and calculate Cohen's d and p-value
results <- data.frame(feature = data$feature, 
                      cohen_Surveillant_Primed = NA, 
                      cohen_Primed_Activated = NA, 
                      p_Surveillant_Primed = NA, 
                      p_Primed_Activated = NA)

for (i in 1:nrow(data)) {
  # Cohen's d for Surveillant to Primed
  results$cohen_surveillant_primed[i] <- cohens_d(data$surveillant_mean[i], data$primed_mean[i], 
                                                  data$surveillant_sd[i], data$primed_sd[i], 
                                                  data$surveillant_n, data$primed_n)
  
  # Cohen's d for Primed to Activated
  results$cohen_primed_activated[i] <- cohens_d(data$primed_mean[i], data$activated_mean[i], 
                                                data$primed_sd[i], data$activated_sd[i], 
                                                data$primed_n, data$activated_n)
  
  # t-test to calculate p-values for Surveillant vs Primed
  t_test_surveillant_primed <- t.test(rnorm(data$surveillant_n, data$surveillant_mean[i], data$surveillant_sd[i]), 
                                      rnorm(data$primed_n, data$primed_mean[i], data$primed_sd[i]))
  
  # t-test to calculate p-values for Primed vs Activated
  t_test_primed_activated <- t.test(rnorm(data$primed_n, data$primed_mean[i], data$primed_sd[i]), 
                                    rnorm(data$activated_n, data$activated_mean[i], data$activated_sd[i]))
  
  # Storing the p-values in the results dataframe
  results$p_surveillant_primed[i] <- t_test_surveillant_primed$p.value
  results$p_primed_activated[i] <- t_test_primed_activated$p.value
}

# Display results
print(results)


#############
###########


# Data: Means, SDs, and sample sizes for each feature
data <- data.frame(
  feature = c("Cell Territorial Vol (μm³)", "Cell Volume (μm³)", "Ramification index",
              "Number of endpoints", "Number of branchpoints", "Average Branch length (μm)",
              "Max Branch length (μm)", "Min Branch length (μm)"),
  surveillant_mean = c(0.684553706, 0.618984941, 0.530450387, 
                       0.572861312, 0.568244302, 0.623604289, 
                       0.594886405, 0.013777315),
  primed_mean = c(-0.52610636, -0.476842092, -0.340322378, 
                  -0.436606218, -0.399771783, -0.321339783, 
                  -0.365577085, 0.210972366),
  activated_mean = c(-0.806181683, -0.722256679, -1.025148226, 
                     -0.696407102, -0.888883759, -1.673398579, 
                     -1.245320209, -1.33359196),
  surveillant_sd = c(1.111323612, 1.146750765, 1.049184999, 
                     1.238752512, 1.228219586, 0.950420004, 
                     1.141126557, 0.790652079),
  primed_sd = c(0.314500854, 0.40233421, 0.696677354, 
                0.248352883, 0.268652826, 0.5582754, 
                0.373939514, 1.089598584),
  activated_sd = c(0.138715885, 0.314963721, 0.355740163, 
                   0.069361581, 0, 0, 
                   0, 0),
  surveillant_n = 1445,  # Adjusted from your original numbers
  primed_n = 1673,  
  activated_n = 376  
)

# Cohen's d Calculation Function
cohens_d <- function(m1, m2, sd1, sd2, n1, n2) {
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  d <- (m1 - m2) / pooled_sd
  return(d)
}

# Loop through each feature and calculate Cohen's d and p-value
results <- data.frame(feature = data$feature, 
                      cohen_surveillant_primed = NA, 
                      cohen_primed_activated = NA, 
                      p_surveillant_primed = NA, 
                      p_primed_activated = NA)

for (i in 1:nrow(data)) {
  # Cohen's d for Surveillant to Primed
  results$cohen_surveillant_primed[i] <- cohens_d(data$surveillant_mean[i], data$primed_mean[i], 
                                                  data$surveillant_sd[i], data$primed_sd[i], 
                                                  data$surveillant_n, data$primed_n)
  
  # Cohen's d for Primed to Activated
  results$cohen_primed_activated[i] <- cohens_d(data$primed_mean[i], data$activated_mean[i], 
                                                data$primed_sd[i], data$activated_sd[i], 
                                                data$primed_n, data$activated_n)
  
  # t-test to calculate p-values for Surveillant vs Primed
  t_test_surveillant_primed <- wilcox.test(rnorm(data$surveillant_n, data$surveillant_mean[i], data$surveillant_sd[i]), 
                                      rnorm(data$primed_n, data$primed_mean[i], data$primed_sd[i]))
  
  # t-test to calculate p-values for Primed vs Activated
  t_test_primed_activated <- wilcox.test(rnorm(data$primed_n, data$primed_mean[i], data$primed_sd[i]), 
                                    rnorm(data$activated_n, data$activated_mean[i], data$activated_sd[i]))
  
  # t-test to calculate p-values for Activated vs Surveillant
  t_test_activated_surveillant <- wilcox.test(rnorm(data$activated_n, data$activated_mean[i], data$activated_sd[i]), 
                                         rnorm(data$surveillant_n, data$surveillant_mean[i], data$surveillant_sd[i]))
  
  
  # Storing the p-values in the results dataframe
  results$p_surveillant_primed[i] <- t_test_surveillant_primed$p.value
  results$p_primed_activated[i] <- t_test_primed_activated$p.value
  results$p_activated_surveillant[i] <- t_test_activated_surveillant$p.value
}

# Display results
print(results)






# Convert all columns (except the first, which is likely categorical) to numeric
features_numeric <- features
features_numeric[, -1] <- lapply(features[, -1], as.numeric)

# Apply Shapiro-Wilk test to each numeric column
normality_results <- apply(features_numeric[, -1], 2, function(x) shapiro.test(x)$p.value)

# Print results
print(normality_results)



# Check if the data is normalized by checking the mean and standard deviation
check_normalization <- function(data) {
  result <- data.frame(
    feature = data$feature,
    mean_surveillant = data$surveillant_mean,
    mean_primed = data$primed_mean,
    mean_activated = data$activated_mean,
    sd_surveillant = data$surveillant_sd,
    sd_primed = data$primed_sd,
    sd_activated = data$activated_sd,
    is_surveillant_normalized = abs(data$surveillant_mean) < 0.1 & abs(data$surveillant_sd - 1) < 0.1,
    is_primed_normalized = abs(data$primed_mean) < 0.1 & abs(data$primed_sd - 1) < 0.1,
    is_activated_normalized = abs(data$activated_mean) < 0.1 & abs(data$activated_sd - 1) < 0.1
  )
  return(result)
}

# Call the function to check if the data is normalized
check_normalization(data)








library(dplyr)

# Assuming your data is in a data frame called df
features_surveillant <- features %>% filter(cluster_pam == 1)
features_surveillant

# Assuming your data is in a data frame called df
features_primed <- features %>% filter(cluster_pam == 2)
features_primed

# Assuming your data is in a data frame called df
features_activated <- features %>% filter(cluster_pam == 3)
features_activated





# Convert all columns (except the first, which is likely categorical) to numeric
features_numeric <- features_surveillant
features_numeric[, -1] <- lapply(features[, -1], as.numeric)

# Apply Shapiro-Wilk test to each numeric column
normality_results <- apply(features_numeric[, -1], 2, function(x) shapiro.test(x)$p.value)

# Print results
print(normality_results)




# Define the dataset (as provided)
data <- data.frame(
  cluster_pam = c("Surveillant", "Primed", "Activated"),
  CellTerritoryVol_um3_mean = c(119701.5098, 32669.58774, 12535.5349),
  CellTerritoryVol_um3_sd = c(79890.8238, 22608.83514, 9972.006528),
  CellVolumes_mean = c(17392.82035, 7021.816636, 4699.191336),
  CellVolumes_sd = c(10852.9504, 3807.726456, 2980.844442),
  RamificationIndex_mean = c(7.487929119, 4.788488832, 2.665493537),
  RamificationIndex_sd = c(3.252527375, 2.159735574, 1.102812772),
  NumOfEndpoints_mean = c(21.47675879, 6.58955677, 2.758122744),
  NumOfEndpoints_sd = c(18.26860038, 3.66260373, 1.022915384),
  NumOfBranchpoints_mean = c(15.67525126, 5.261687917, 0),
  NumOfBranchpoints_sd = c(13.21273753, 2.890068941, 0),
  AvgBranchLength_mean = c(101.265866, 59.6069804, 0),
  AvgBranchLength_sd = c(41.90029805, 24.61217732, 0),
  MaxBranchLength_mean = c(216.3742645, 103.4415211, 0),
  MaxBranchLength_sd = c(134.1753787, 43.96837108, 0),
  MinBranchLength_mean = c(20.88398784, 23.94047661, 0),
  MinBranchLength_sd = c(12.25496878, 16.88858727, 0),
  Surveillant_n = 1592,
  Primed_n = 1647,
  Activated_n = 277
)

# Define the Cohen's d function
cohen_d <- function(group1, group2) {
  mean_diff <- mean(group1) - mean(group2)
  pooled_sd <- sqrt(((length(group1) - 1) * sd(group1)^2 + (length(group2) - 1) * sd(group2)^2) / (length(group1) + length(group2) - 2))
  return(mean_diff / pooled_sd)
}

# Extract data for Surveillant and Primed cell types
cell_type_1 <- data[data$cluster_pam == "Surveillant", ]
cell_type_2 <- data[data$cluster_pam == "Primed", ]

# Calculate Cohen's d for each feature (mean columns only)
cohen_d_results <- data.frame(feature = character(), effect_size = numeric(), stringsAsFactors = FALSE)

# Loop through each feature column (mean columns only)
for (feature in colnames(data)[grepl("_mean", colnames(data))]) {
  group1_values <- cell_type_1[[feature]]
  group2_values <- cell_type_2[[feature]]
  
  # Calculate Cohen's d
  effect_size <- cohen_d(group1_values, group2_values)
  
  # Store results
  cohen_d_results <- rbind(cohen_d_results, data.frame(feature = feature, effect_size = effect_size))
}

# View Cohen's d results
print(cohen_d_results)

# Define the Cohen's d calculation function from means and SDs
cohen_d_from_stats <- function(mean1, mean2, sd1, sd2) {
  return((mean1 - mean2) / sqrt((sd1^2 + sd2^2) / 2))
}

# Initialize a results dataframe for Cohen's d and p-values
cohen_d_results <- data.frame(feature = character(), cohen_d = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each feature to calculate Cohen's d and p-value
for (feature in colnames(features_surveillant)) {
  # Extract data for the current feature
  group_surveillant <- features_surveillant[[feature]]
  group_primed <- features_primed[[feature]]
  
  # Calculate the means and standard deviations for both groups
  mean_surveillant <- mean(group_surveillant)
  mean_primed <- mean(group_primed)
  sd_surveillant <- sd(group_surveillant)
  sd_primed <- sd(group_primed)
  
  # Calculate Cohen's d for the current feature
  effect_size <- cohen_d_from_stats(mean_surveillant, mean_primed, sd_surveillant, sd_primed)
  
  # Perform a t-test to calculate the p-value
  t_test_result <- t.test(group_surveillant, group_primed)
  p_value <- t_test_result$p.value
  
  # Store the results
  cohen_d_results <- rbind(cohen_d_results, data.frame(feature = feature, cohen_d = effect_size, p_value = p_value))
}

# Print the results
print(cohen_d_results)






# Define the Cohen's d calculation function from means and SDs
cohen_d_from_stats <- function(mean1, mean2, sd1, sd2) {
  return((mean1 - mean2) / sqrt((sd1^2 + sd2^2) / 2))
}

# Initialize a results dataframe for Cohen's d and p-values
cohen_d_results_primed_activated <- data.frame(feature = character(), cohen_d = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each feature to calculate Cohen's d and p-value for Primed vs Activated
for (feature in colnames(features_primed)) {
  # Extract data for the current feature
  group_primed <- features_primed[[feature]]
  group_activated <- features_activated[[feature]]
  
  # Calculate the means and standard deviations for both groups
  mean_primed <- mean(group_primed)
  mean_activated <- mean(group_activated)
  sd_primed <- sd(group_primed)
  sd_activated <- sd(group_activated)
  
  # Calculate Cohen's d for the current feature
  effect_size <- cohen_d_from_stats(mean_primed, mean_activated, sd_primed, sd_activated)
  
  # Perform a t-test to calculate the p-value
  t_test_result <- wilcox.test(group_primed, group_activated)
  p_value <- t_test_result$p.value
  
  # Store the results
  cohen_d_results_primed_activated <- rbind(cohen_d_results_primed_activated, data.frame(feature = feature, cohen_d = effect_size, p_value = p_value))
}

# Print the results
print(cohen_d_results_primed_activated)





























# Define the Cohen's d calculation function from means and SDs
cohen_d_from_stats <- function(mean1, mean2, sd1, sd2) {
  return((mean1 - mean2) / sqrt((sd1^2 + sd2^2) / 2))
}

# Initialize a results dataframe for Cohen's d, p-values, means and SDs
cohen_d_results_primed_activated <- data.frame(
  feature = character(),
  cohen_d = numeric(),
  p_value = numeric(),
  mean_primed = numeric(),
  sd_primed = numeric(),
  mean_activated = numeric(),
  sd_activated = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each feature to calculate Cohen's d, p-value, means, and SDs for Primed vs Activated
for (feature in colnames(features_primed)) {
  # Extract data for the current feature
  group_primed <- features_primed[[feature]]
  group_activated <- features_activated[[feature]]
  
  # Calculate the means and standard deviations for both groups
  mean_primed <- mean(group_primed)
  mean_activated <- mean(group_activated)
  sd_primed <- sd(group_primed)
  sd_activated <- sd(group_activated)
  
  # Calculate Cohen's d for the current feature
  effect_size <- cohen_d_from_stats(mean_primed, mean_activated, sd_primed, sd_activated)
  
  # Perform Mann-Whitney U test to calculate the p-value
  mann_whitney_result <- wilcox.test(group_primed, group_activated)
  p_value <- mann_whitney_result$p.value
  
  # Store the results with mean ± SD information
  cohen_d_results_primed_activated <- rbind(cohen_d_results_primed_activated, data.frame(
    feature = feature,
    cohen_d = effect_size,
    p_value = p_value,
    mean_primed = mean_primed,
    sd_primed = sd_primed,
    mean_activated = mean_activated,
    sd_activated = sd_activated
  ))
}

# Print the results
print(cohen_d_results_primed_activated)






# Define the Cohen's d calculation function from means and SDs
cohen_d_from_stats <- function(mean1, mean2, sd1, sd2) {
  return((mean1 - mean2) / sqrt((sd1^2 + sd2^2) / 2))
}

# Initialize a results dataframe for Cohen's d, p-values, means, and SDs
cohen_d_results_surveillant_primed <- data.frame(
  feature = character(),
  cohen_d = numeric(),
  p_value = numeric(),
  mean_surveillant = numeric(),
  sd_surveillant = numeric(),
  mean_primed = numeric(),
  sd_primed = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each feature to calculate Cohen's d, p-value, means, and SDs for Surveillant vs Primed
for (feature in colnames(features_surveillant)) {
  # Extract data for the current feature
  group_surveillant <- features_surveillant[[feature]]
  group_primed <- features_primed[[feature]]
  
  # Calculate the means and standard deviations for both groups
  mean_surveillant <- mean(group_surveillant)
  mean_primed <- mean(group_primed)
  sd_surveillant <- sd(group_surveillant)
  sd_primed <- sd(group_primed)
  
  # Calculate Cohen's d for the current feature
  effect_size <- cohen_d_from_stats(mean_surveillant, mean_primed, sd_surveillant, sd_primed)
  
  # Perform Mann-Whitney U test to calculate the p-value
  mann_whitney_result <- wilcox.test(group_surveillant, group_primed)
  p_value <- mann_whitney_result$p.value
  
  # Store the results with mean ± SD information
  cohen_d_results_surveillant_primed <- rbind(cohen_d_results_surveillant_primed, data.frame(
    feature = feature,
    cohen_d = effect_size,
    p_value = p_value,
    mean_surveillant = mean_surveillant,
    sd_surveillant = sd_surveillant,
    mean_primed = mean_primed,
    sd_primed = sd_primed
  ))
}

# Print the results
print(cohen_d_results_surveillant_primed)

write.csv(cohen_d_results_surveillant_primed, "/Users/Deepika/Desktop/Microglia_anal/cohen'sd_sur_pri.csv")

write.csv(cohen_d_results_primed_activated, "/Users/Deepika/Desktop/Microglia_anal/cohen_d_results_primed_activated.csv")






# Define the function for Cohen's d (mean-based calculation)
cohen_d_from_stats <- function(mean1, mean2, sd1, sd2) {
  return((mean1 - mean2) / sqrt((sd1^2 + sd2^2) / 2))
}

# Define the Mann-Whitney U test function for p-value
mann_whitney_p_value <- function(group1, group2) {
  result <- wilcox.test(group1, group2)
  return(result$p.value)
}

# Initialize a results dataframe for Cohen's d, p-values, medians, and SDs
results <- data.frame(
  feature = character(),
  cohen_d = numeric(),
  p_value = numeric(),
  median_group1 = numeric(),
  median_group2 = numeric(),
  mean_group1 = numeric(),
  mean_group2 = numeric(),
  sd_group1 = numeric(),
  sd_group2 = numeric(),
  stringsAsFactors = FALSE
)

# Loop through features to calculate Cohen's d and p-value for median comparisons
for (feature in colnames(features_surveillant)) {
  # Extract data for the current feature
  group_surveillant <- features_surveillant[[feature]]
  group_primed <- features_primed[[feature]]
  
  # Calculate the means and standard deviations for both groups
  mean_surveillant <- mean(group_surveillant)
  mean_primed <- mean(group_primed)
  sd_surveillant <- sd(group_surveillant)
  sd_primed <- sd(group_primed)
  
  # Calculate Cohen's d for the current feature
  cohen_d <- cohen_d_from_stats(mean_surveillant, mean_primed, sd_surveillant, sd_primed)
  
  # Perform Mann-Whitney U test to calculate the p-value
  p_value <- mann_whitney_p_value(group_surveillant, group_primed)
  
  # Calculate the medians for both groups
  median_surveillant <- median(group_surveillant)
  median_primed <- median(group_primed)
  
  # Store the results with mean, median, SD, and p-value information
  results <- rbind(results, data.frame(
    feature = feature,
    cohen_d = cohen_d,
    p_value = p_value,
    median_group1 = median_surveillant,
    median_group2 = median_primed,
    mean_group1 = mean_surveillant,
    mean_group2 = mean_primed,
    sd_group1 = sd_surveillant,
    sd_group2 = sd_primed
  ))
}

# Print the results
print(results)



# Calculate the median for each feature in each dataset
medians_surveillant <- apply(features_surveillant, 2, median)
medians_activated <- apply(features_activated, 2, median)
medians_primed <- apply(features_primed, 2, median)

# Print the medians for each feature
print("Medians for Surveillant group:")
print(medians_surveillant)

print("Medians for Activated group:")
print(medians_activated)

print("Medians for Primed group:")
medians_primed



















# Your data manually input as a data frame

data <- data.frame(
  feature = c("Cell Territorial Vol (μm³)", "Cell Volume (μm³)", "Ramification index",
              "Number of endpoints", "Number of branchpoints", "Average Branch length (μm)",
              "Max Branch length (μm)", "Min Branch length (μm)"),
  surveillant_mean = c(82803.43109, 13374.0155, 5.911876459, 
                       14.62015504, 10.79069767, 79.04487162, 
                       159.3826151, 21.49710823),
  primed_mean = c(61994.92032, 10346.95189, 5.800060495, 
                  11.55237515, 8.425700365, 70.39982501, 
                  137.020853, 20.38407964),
  activated_mean = c(55651.91396, 9157.99278, 5.656583585, 
                     9.78700361, 7.241877256, 63.76876485, 
                     122.616716, 19.11850527),
  surveillant_sd = c(79887.5512, 9840.951477, 3.003499307, 
                     15.51625399, 11.34668996, 44.62292084, 
                     130.4129168, 15.3257215),
  primed_sd = c(64084.57902, 9036.293386, 3.152904238, 
                13.50229619, 9.838409088, 43.18174435, 
                104.8782026, 15.646236),
  activated_sd = c(60361.25523, 8447.214827, 3.30079233, 
                   10.12859942, 7.593092357, 40.75224258, 
                   90.52319411, 15.94972),
  surveillant_n = 1445,
  primed_n = 1673,
  activated_n = 376
)

# Function to calculate Cohen's d
cohens_d <- function(m1, m2, sd1, sd2, n1, n2) {
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  (m1 - m2) / pooled_sd
}

# Initialize results dataframe
results <- data.frame(
  feature = data$feature,
  cohen_surveillant_primed = NA,
  cohen_primed_activated = NA,
  cohen_activated_surveillant = NA,
  p_surveillant_primed = NA,
  p_primed_activated = NA,
  p_activated_surveillant = NA
)

# Loop through each row/feature
for (i in 1:nrow(data)) {
  # Simulate data
  surveillant_vals <- rnorm(data$surveillant_n[i], data$surveillant_mean[i], data$surveillant_sd[i])
  primed_vals <- rnorm(data$primed_n[i], data$primed_mean[i], data$primed_sd[i])
  activated_vals <- rnorm(data$activated_n[i], data$activated_mean[i], data$activated_sd[i])
  
  # Calculate Cohen's d
  results$cohen_surveillant_primed[i] <- cohens_d(data$surveillant_mean[i], data$primed_mean[i],
                                                  data$surveillant_sd[i], data$primed_sd[i],
                                                  data$surveillant_n[i], data$primed_n[i])
  
  results$cohen_primed_activated[i] <- cohens_d(data$primed_mean[i], data$activated_mean[i],
                                                data$primed_sd[i], data$activated_sd[i],
                                                data$primed_n[i], data$activated_n[i])
  
  results$cohen_activated_surveillant[i] <- cohens_d(data$activated_mean[i], data$surveillant_mean[i],
                                                     data$activated_sd[i], data$surveillant_sd[i],
                                                     data$activated_n[i], data$surveillant_n[i])
  
  # Calculate Wilcoxon p-values (non-parametric)
  results$p_surveillant_primed[i] <- wilcox.test(surveillant_vals, primed_vals)$p.value
  results$p_primed_activated[i] <- wilcox.test(primed_vals, activated_vals)$p.value
  results$p_activated_surveillant[i] <- wilcox.test(activated_vals, surveillant_vals)$p.value
}

# Display the results
print(results)

