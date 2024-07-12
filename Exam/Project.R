library(dplyr)
library(reshape2)
library(viridis)
library(caret)
library(cluster)
library(fpc)
library(data.table)
library(factoextra)
library(dbscan)
library(rrcov)
library(NbClust)
library(plotly)
library(gridExtra)
library(ClusterStability)
library(kernlab)
library(ggplot2)
library(motifcluster)
library(Spectrum)
library(grid)

df_metrics <- data.frame(Name = "", 
                         DB = as.numeric(-1), 
                         CH = as.numeric(-1), 
                         SI = as.numeric(-1))


data <- read.csv("C:\\Users\\mirko\\Desktop\\Progetto Esame\\CC GENERAL.csv")
color_palette_k <- c("#FF6347", "#4682B4", "#33FF57", "#FFD700", "#9D4EDD")

# Remove CUST_ID column
data <- data[, !names(data) %in% "CUST_ID"]

numeric_cols <- sapply(data, is.numeric)
numeric_data <- data[, numeric_cols & names(data) != "TENURE"]
colors <- rainbow(ncol(numeric_data))

par(mfrow = c(2, 4), mar = c(1, 1, 1, 1) + 0.1) 
for (col in names(numeric_data)) {
  boxplot(numeric_data[[col]], 
          main = col, 
          col = "#4682B4",
          cex.main = 0.8)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# Remove rows with missing values in CREDIT_LIMIT column
data <- data[!is.na(data$CREDIT_LIMIT), ]

#Fix problem MINIMUM_PAYMENTS column, NA value, imputation
payments_mean <- mean(data$PAYMENTS, na.rm = TRUE)
na_count <- sum(is.na(data$MINIMUM_PAYMENTS))
na_count
data$MINIMUM_PAYMENTS[is.na(data$MINIMUM_PAYMENTS) & data$PAYMENTS == 0] <- 0
na_count <- sum(is.na(data$MINIMUM_PAYMENTS))
na_count
data$MINIMUM_PAYMENTS[is.na(data$MINIMUM_PAYMENTS) & data$PAYMENTS > 0 & data$PAYMENTS <= payments_mean] <- data$PAYMENTS
na_count <- sum(is.na(data$MINIMUM_PAYMENTS))
na_count
data$MINIMUM_PAYMENTS[is.na(data$MINIMUM_PAYMENTS)] <- payments_mean
na_count <- sum(is.na(data$MINIMUM_PAYMENTS))
na_count

#Noise detection
Q1 <- apply(data, 2, quantile, probs = 0.25)
Q3 <- apply(data, 2, quantile, probs = 0.75)
IQR <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

count_outliers <- function(x, lower, upper) {
  sum(x < lower | x > upper, na.rm = TRUE)
}

count_out <- sapply(seq_along(data), function(i) {
  count_outliers(data[[i]], lower_bound[i], upper_bound[i])
})

# Create a dataframe with counts of outliers and print it 
df_count_out <- data.frame(count_out)
names(df_count_out) <- "count_out"
print(df_count_out)


#Remove the outliers from the dataset
no_outliers <- function(x, lower, upper) {
  (x >= lower) & (x <= upper)
}

non_outlier_matrix <- mapply(no_outliers, data, lower_bound, upper_bound)
rows_to_keep <- rowSums(non_outlier_matrix) == ncol(data)
cleaned_data <- data[rows_to_keep, ]

# New box-plot without outliers, problem to much outliers, we can't kill all it
numeric_cols <- sapply(cleaned_data, is.numeric)
for (col in names(cleaned_data)[numeric_cols]) {
  boxplot(cleaned_data[[col]], 
          main = "", 
          ylab = col)
}

#Then Soft Noise elimination
data <- data[data$BALANCE < 15000, ]
data <- data[data$PURCHASES < 40000, ]
data <- data[data$ONEOFF_PURCHASES < 30000, ]
data <- data[data$INSTALLMENTS_PURCHASES < 20000, ]
data <- data[data$CASH_ADVANCE < 40000, ]
data <- data[data$MINIMUM_PAYMENTS < 60000, ]

#Density Plot
n_columns=c("BALANCE","PURCHASES","ONEOFF_PURCHASES","INSTALLMENTS_PURCHASES",
           "PURCHASES_TRX","CREDIT_LIMIT","PAYMENTS","TENURE")
plots <- list()

# Create plots
for (name in n_columns) {
  column_min <- min(data[[name]], na.rm = TRUE)
  column_max <- max(data[[name]], na.rm = TRUE)
  column_range <- column_max - column_min
  line_positions <- c(column_min, (column_range / 4) + column_min, 
                      (column_range / 2) + column_min, 
                      (3 * column_range / 4) + column_min, column_max)
  p <- ggplot(data, aes_string(x = name)) +
    geom_density(fill = "#4682B4", color = "#e9ecef", alpha = 1) +
    geom_vline(xintercept = line_positions, linetype = "dashed", color = "DimGray", alpha = 0.5)+
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(size = 5))
  plots[[name]] <- p
}

# Arrange plots in a 3x3 grid
grid.arrange(
  plots[[1]], plots[[2]], plots[[3]],plots[[4]],
  plots[[5]], plots[[6]], plots[[7]],plots[[8]],
  ncol = 4
)


#Heatmap correlation
cor_matrix <- cor(data, use = "complete.obs")
cor_melt <- melt(cor_matrix)

ggplot(cor_melt, aes(Var1, Var2, fill= value)) + 
  geom_tile() +
  geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 2.5)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 6),axis.text.y = element_text(size = 6),     axis.title.x = element_blank(),
        axis.title.y = element_blank()) 

# Normalization/Scaling

df_model <- data
# Standardize the data using scale() function from caret
df_Standardize <- as.data.frame(scale(df_model))
colnames(df_Standardize) <- colnames(df_model)
print(df_Standardize)

# PCA Standard, outliers problem
pca_prcomp <- prcomp(df_Standardize)
summary(pca_prcomp)

plot(pca_prcomp)
biplot(pca_prcomp, cex=0.6, col=c("#35354f", "#c9658f"), pc.biplot=TRUE)
screeplot(pca_prcomp, type="lines", col="#008080") 
pca_prcomp <- prcomp(df_Standardize, rank. = 3)
summary(pca_prcomp)
PCA_PR <- pca_prcomp[["x"]]

#Robust PCA WITH MCD - Good decision
Rpca <- PcaCov(df_Standardize)
summary(Rpca)
PcaClassic(df_Standardize)

biplot(Rpca, cex=0.6, col=c("#35354f", "#c9658f"))
screeplot(Rpca, type="lines", col="#008080") 

#2PCA Correct
Rpca2 <- PcaCov(df_Standardize, k=2) 
summary(Rpca2)
X_pca2 <- Rpca2@scores

plot(Rpca2)
biplot(Rpca2, cex=0.6, col=c("#35354f", "#c9658f"))
screeplot(Rpca2, type="lines", col="#008080") 

#3PCA Correct
Rpca3 <- PcaCov(df_Standardize, k=3) 
summary(Rpca3)
X_pca3 <- Rpca3@scores

plot(Rpca3)
biplot(Rpca3, cex=0.6, col=c("#35354f", "#c9658f"))
screeplot(Rpca3, type="lines", col="#008080") 

contains_na <- any(is.na(Rpca3@scores))
print(contains_na)

#OPTIMAL NUMBER OF CLUSTER FOR KMEANS

p1 <- fviz_nbclust(df_Standardize, kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")
set.seed(32)
p2 <- fviz_nbclust(df_Standardize, kmeans, method = "wss") + 
  labs(subtitle = "Elbow method")
set.seed(32)
library(NbClust)
gap_stat <- clusGap(df_Standardize, FUN = kmeans, K.max = 8, B = 300)
p3 <- fviz_gap_stat(gap_stat, maxSE = list(method = "Tibs2001SEmax")) + 
  labs(subtitle = "Gap statistic method")
grid.arrange(p1, p2, p3, nrow = 1)

#K-MEANS

km <- kmeanspp(df_Standardize, k=4)
cluster_countskm <- table(km$cluster)
cluster_percentageskm <- cluster_countskm / length(km$cluster) * 100 
cluster_percentageskm

sil <- silhouette(km$cluster, dist(df_Standardize))
summary_info <- summary(sil)
mean_sil_width_km <- summary_info$avg.width
mean_sil_width_km

fviz_silhouette(sil, 
                main = "Silhouette Plot of Clusters", 
                palette = color_palette_k, 
                ggtheme = theme_minimal(), 
                geom = "bar") +
  labs(x = "Silhouette Width", y = "Cluster Labels") +
  theme(
    plot.title = element_text(size = 12, face = "bold", family = "serif"),
    legend.text = element_text(size = 9)
  )

#PLOT CLUSTER
color_palette_k <- c("#FF6347", "#4682B4", "#33FF57", "#FFD700")
hullplot(X_pca3, km$cluster, main = "KMeans Clustering")

X_pca3_df <- as.data.frame(X_pca3)
X_pca3_df$cluster =km$cluster
p <- plot_ly(X_pca3_df, x=~PC1, y=~PC2, 
             z=~PC3, color=~cluster, colors=color_palette_k) %>%
  add_markers(size=1.5)%>%
  layout(legend = list(title = list(text = 'Cluster')))
print(p)

X_pca3_df$cluster =as.factor(km$cluster)

cluster_counts <- count(X_pca3_df, cluster)
cluster_percentages <- cluster_counts %>%
  mutate(percentage = n / sum(n) * 100)
X_pca3_df <- left_join(X_pca3_df, cluster_percentages, by = "cluster")


mainp <- ggplot(X_pca3_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75) +
  labs(title = "Scatter Plot Clusters Distributions") +
  scale_color_manual(values = color_palette_k, 
                     labels = paste0('Cluster ', 1:5, ' (', round(cluster_percentages$percentage, 1), '%)'), 
                     breaks = c(1, 2, 3, 4, 5),
                     drop = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

subset_df <- subset(X_pca3_df, cluster == 1)
sub1 <- ggplot(subset_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75, color=color_palette_k[1]) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

subset_df <- subset(X_pca3_df, cluster == 2)
sub2 <- ggplot(subset_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75, color=color_palette_k[2]) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

subset_df <- subset(X_pca3_df, cluster == 3)
sub3 <- ggplot(subset_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75, color=color_palette_k[3]) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

subset_df <- subset(X_pca3_df, cluster == 4)
sub4 <- ggplot(subset_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75, color=color_palette_k[4]) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

print(mainp)

box <- rectGrob(gp=gpar(fill="white", col=NA))
grid.arrange(
  arrangeGrob(
    mainp, 
    ncol = 1 
  ),
  arrangeGrob(
    sub1,
    sub2,
    box,
    ncol = 1
  ),
  arrangeGrob(
    sub3,
    sub4,
    box,
    ncol = 1
  ),
  widths = c(3,1,1)
)

#Metriche KMEANS

DB_SCORE_KMEANS <- davies_bouldin_score(df_Standardize, km$cluster)
CH_SCORE_KMEANS <- calinski_harabasz_score(df_Standardize, km$cluster)
cat(".: Davies-Bouldin Index: \033[1;32m", DB_SCORE_KMEANS, "\033[0m\n")
cat(".: Silhouette Score: \033[1;32m", mean_sil_width_km, "\033[0m\n")
cat(".: Calinski Harabasz Index: \033[1;32m", CH_SCORE_KMEANS, "\033[0m\n")
df_metrics[nrow(df_metrics)+1, ] <- c(Name = "KMeans", DB = DB_SCORE_KMEANS, CH = CH_SCORE_KMEANS, SI = mean_sil_width_km)

#Hierarchical clustering best K

set.seed(32)
p1 <- fviz_nbclust(df_Standardize, hcut, method = "silhouette", hc_method="ward.D2") +
  labs(subtitle = "Silhouette method")
set.seed(32)
p2 <- fviz_nbclust(df_Standardize, hcut, method = "wss", hc_method="ward.D2") + 
  labs(subtitle = "Elbow method")
grid.arrange(p1, p2, nrow = 1)

#Hierarchical clustering

hc = hclust(dist(df_Standardize,method="euclidean"), method="ward.D2")
dend <- as.dendrogram(hc)
plot(dend, ylab = "Height") # ,horiz = TRUE
abline(h = 115, col = "red", lty = 2)

hc_clusters <- cutree(hc, k = 5) # K = 5
cluster_countshc <- table(hc_clusters)
cluster_percentageshc <- cluster_countshc / length(hc_clusters) * 100 
cluster_percentageshc

color_palette_k <- c("#FF6347", "#4682B4", "#33FF57", "#FFD700", "#9D4EDD")
sil <- silhouette(hc_clusters, dist(df_Standardize, method="euclidean"))
summary_info <- summary(sil)
mean_sil_width_HC <- summary_info$avg.width

fviz_silhouette(sil, 
                main = "Silhouette Plot of Clusters", 
                palette = color_palette_k, 
                ggtheme = theme_minimal(), 
                geom = "bar") +
  labs(x = "Silhouette Width", y = "Cluster Labels") +
  theme(
    plot.title = element_text(size = 12, face = "bold", family = "serif"),
    legend.text = element_text(size = 9)
  )

#PLOT CLUSTER

hullplot(X_pca3, hc_clusters, main = "Hierarchical Clustering")

X_pca3_df <- as.data.frame(X_pca3)
X_pca3_df$cluster = hc_clusters
p <- plot_ly(X_pca3_df, x=~PC1, y=~PC2, 
             z=~PC3, color=~cluster, colors=color_palette_k) %>%
  add_markers(size=1.5)%>%
  layout(legend = list(title = list(text = 'Cluster')))
print(p)

X_pca3_df$cluster =as.factor(hc_clusters)

cluster_counts <- count(X_pca3_df, cluster)
cluster_percentages <- cluster_counts %>%
  mutate(percentage = n / sum(n) * 100)
X_pca3_df <- left_join(X_pca3_df, cluster_percentages, by = "cluster")


mainp <- ggplot(X_pca3_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75) +
  labs(title = "Scatter Plot Clusters Distributions") +
  scale_color_manual(values = color_palette_k, 
                     labels = paste0('Cluster ', 1:5, ' (', round(cluster_percentages$percentage, 1), '%)'), 
                     breaks = c(1, 2, 3, 4, 5),
                     drop = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

subset_df <- subset(X_pca3_df, cluster == 1)
sub1 <- ggplot(subset_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75, color=color_palette_k[1]) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

subset_df <- subset(X_pca3_df, cluster == 2)
sub2 <- ggplot(subset_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75, color=color_palette_k[2]) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

subset_df <- subset(X_pca3_df, cluster == 3)
sub3 <- ggplot(subset_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75, color=color_palette_k[3]) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

subset_df <- subset(X_pca3_df, cluster == 4)
sub4 <- ggplot(subset_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75, color=color_palette_k[4]) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

subset_df <- subset(X_pca3_df, cluster == 5)
sub5 <- ggplot(subset_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75, color=color_palette_k[5]) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
box <- rectGrob(gp=gpar(fill="white", col=NA))
grid.arrange(
  arrangeGrob(
    mainp,
    ncol = 1
  ),
  arrangeGrob(
    sub1,
    sub2,
    sub3,
    ncol = 1
  ),
  arrangeGrob(
    sub4,
    sub5,
    box,
    ncol = 1
  ),
  widths = c(3, 1, 1)
)

#METRICHE Hierarchical clustering
DB_SCORE_HC <- davies_bouldin_score(df_Standardize, hc_clusters)
CH_SCORE_HC <- calinski_harabasz_score(df_Standardize, hc_clusters)
cat(".: Davies-Bouldin Index: \033[1;32m", DB_SCORE_HC, "\033[0m\n")
cat(".: Silhouette Score: \033[1;32m", mean_sil_width_HC, "\033[0m\n")
cat(".: Calinski Harabasz Index: \033[1;32m", CH_SCORE_HC, "\033[0m\n")
df_metrics[nrow(df_metrics)+1, ] <- c(Name = "Hierarchical", DB = DB_SCORE_HC, CH = CH_SCORE_HC, SI = mean_sil_width_HC)

#SPECTRAL CLUSTERING
X_pca3 <- Rpca3@scores
r <- Spectrum(df_Standardize)
sc <- specc(X_pca3, centers = 3, nystrom.red=TRUE)

cluster_labels <- sc@.Data
sc

sil <- silhouette(cluster_labels, dist(df_Standardize))
summary_info <- summary(sil)
mean_sil_width_sc <- summary_info$avg.width
mean_sil_width_sc

fviz_silhouette(sil, 
                main = "Silhouette Plot of Clusters", 
                palette = color_palette_k, 
                ggtheme = theme_minimal(), 
                geom = "bar") +
  labs(x = "Silhouette Width", y = "Cluster Labels") +
  theme(
    plot.title = element_text(size = 12, face = "bold", family = "serif"),
    legend.text = element_text(size = 9)
  )

#PLOT Spectral Clustering

hullplot(X_pca3, cluster_labels, main = "Spectral Clustering")
X_pca3_df <- as.data.frame(X_pca3)

color_palette_k3 <- c("#FF6347", "#4682B4", "#33FF57")
X_pca3_df$cluster = cluster_labels
p <- plot_ly(X_pca3_df, x=~PC1, y=~PC2, 
             z=~PC3, color=~cluster, colors=color_palette_k3) %>%
  add_markers(size=1.5)%>%
  layout(legend = list(title = list(text = 'Cluster')))
print(p)

X_pca3_df$cluster =as.factor(cluster_labels)

cluster_counts <- count(X_pca3_df, cluster)
cluster_percentages <- cluster_counts %>%
  mutate(percentage = n / sum(n) * 100)
X_pca3_df <- left_join(X_pca3_df, cluster_percentages, by = "cluster")

mainp <- ggplot(X_pca3_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75) +
  scale_color_manual(values = color_palette_k, 
                     labels = paste0('Cluster ', 1:3, ' (', round(cluster_percentages$percentage, 1), '%)'), 
                     breaks = c(1, 2, 3),
                     drop = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

subset_df <- subset(X_pca3_df, cluster == 1)
sub1 <- ggplot(subset_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75, color=color_palette_k[1]) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

subset_df <- subset(X_pca3_df, cluster == 2)
sub2 <- ggplot(subset_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75, color=color_palette_k[2]) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

subset_df <- subset(X_pca3_df, cluster == 3)
sub3 <- ggplot(subset_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2.5, alpha = 0.75, color=color_palette_k[3]) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

grid.arrange(
  arrangeGrob(
    mainp,
    ncol = 1
  ),
  arrangeGrob(
    sub1,
    sub2,
    sub3,
    ncol = 1 
  ),
  widths = c(3, 1)
)

#METRICHE Spectral clustering
DB_SCORE_SC <- davies_bouldin_score(df_Standardize, cluster_labels)
CH_SCORE_SC <- calinski_harabasz_score(df_Standardize, cluster_labels)
cat(".: Davies-Bouldin Index: \033[1;32m", DB_SCORE_SC, "\033[0m\n")
cat(".: Silhouette Score: \033[1;32m", mean_sil_width_sc, "\033[0m\n")
cat(".: Calinski Harabasz Index: \033[1;32m", CH_SCORE_SC, "\033[0m\n")
df_metrics[nrow(df_metrics)+1, ] <- c(Name = "Spectral", DB = DB_SCORE_SC, CH = CH_SCORE_SC, SI = mean_sil_width_sc)


#PLOT FINAL RESULT CLUSTERING PROCESS


#--Scatter Plot Kmeans

eda <- data[]
eda$CLUSTER <- as.character(km$cluster)

mainp <- ggplot(eda, aes(x = CREDIT_LIMIT, y = BALANCE, color = CLUSTER)) +
  geom_point(size = 2.5, alpha = 0.75) +
  labs(
    title = "Scatter Plot Clusters Distributions",
    subtitle = "Clusters 1 and 4 have the highest balance and credit limit"
  ) +
  scale_color_manual(values = color_palette_k, 
                     labels = paste0('Cluster ', 1:5), 
                     breaks = c(1, 2, 3, 4, 5),
                     drop = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

eda$point_color <- ifelse(eda$CLUSTER == 1, color_palette_k[1], "grey")
eda$point_alpha <- ifelse(eda$CLUSTER == 1, 1, 0.1)
sub1 <- ggplot(eda, aes(x = CREDIT_LIMIT, y = BALANCE)) +
  geom_point(aes(color = point_color, alpha = point_alpha), size = 2.5) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_color_identity() +
  scale_alpha_identity()

eda$point_color <- ifelse(eda$CLUSTER == 2, color_palette_k[2], "grey")
eda$point_alpha <- ifelse(eda$CLUSTER == 2, 1, 0.1)
sub2 <- ggplot(eda, aes(x = CREDIT_LIMIT, y = BALANCE)) +
  geom_point(aes(color = point_color, alpha = point_alpha), size = 2.5) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_color_identity() +
  scale_alpha_identity()

eda$point_color <- ifelse(eda$CLUSTER == 3, color_palette_k[3], "grey")
eda$point_alpha <- ifelse(eda$CLUSTER == 3, 1, 0.1)
sub3 <- ggplot(eda, aes(x = CREDIT_LIMIT, y = BALANCE)) +
  geom_point(aes(color = point_color, alpha = point_alpha), size = 2.5) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_color_identity() +
  scale_alpha_identity()

eda$point_color <- ifelse(eda$CLUSTER == 4, color_palette_k[4], "grey")
eda$point_alpha <- ifelse(eda$CLUSTER == 4, 1, 0.1)
sub4 <- ggplot(eda, aes(x = CREDIT_LIMIT, y = BALANCE)) +
  geom_point(aes(color = point_color, alpha = point_alpha), size = 2.5) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_color_identity() + 
  scale_alpha_identity()

grid.arrange(
  arrangeGrob(
    mainp, 
    ncol = 1 
  ),
  arrangeGrob(
    sub1,
    sub2,
    sub3,
    sub4,
    ncol = 1 
  ),
  widths = c(3, 1)
)

eda <- data[]
eda$CLUSTER <- as.character(km$cluster)

mainp <- ggplot(eda, aes(x = TENURE, y = PAYMENTS, color = CLUSTER)) +
  geom_point(size = 2.5, alpha = 0.75, position = position_jitter(width = 0.25, height = 0)) +
  labs(
    title = "Scatter Plot Clusters Distributions",
    subtitle = "Most customers in clusters 2 and 3 have zero payments"
  ) +
  scale_color_manual(values = color_palette_k, 
                     labels = paste0('Cluster ', 1:5), 
                     breaks = c(1, 2, 3, 4, 5),
                     drop = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

eda$point_color <- ifelse(eda$CLUSTER == 1, color_palette_k[1], "grey")
eda$point_alpha <- ifelse(eda$CLUSTER == 1, 1, 0.1)
sub1 <- ggplot(eda, aes(x = TENURE, y = PAYMENTS)) +
  geom_point(aes(color = point_color, alpha = point_alpha), size = 2.5, position = position_jitter(width = 0.25, height = 0)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_color_identity() +
  scale_alpha_identity()

eda$point_color <- ifelse(eda$CLUSTER == 2, color_palette_k[2], "grey")
eda$point_alpha <- ifelse(eda$CLUSTER == 2, 1, 0.1)
sub2 <- ggplot(eda, aes(x = TENURE, y = PAYMENTS)) +
  geom_point(aes(color = point_color, alpha = point_alpha), size = 2.5, position = position_jitter(width = 0.25, height = 0)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_color_identity() +
  scale_alpha_identity()

eda$point_color <- ifelse(eda$CLUSTER == 3, color_palette_k[3], "grey")
eda$point_alpha <- ifelse(eda$CLUSTER == 3, 1, 0.1)
sub3 <- ggplot(eda, aes(x = TENURE, y = PAYMENTS)) +
  geom_point(aes(color = point_color, alpha = point_alpha), size = 2.5, position = position_jitter(width = 0.25, height = 0)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_color_identity() +
  scale_alpha_identity()

eda$point_color <- ifelse(eda$CLUSTER == 4, color_palette_k[4], "grey")
eda$point_alpha <- ifelse(eda$CLUSTER == 4, 1, 0.1)
sub4 <- ggplot(eda, aes(x = TENURE, y = PAYMENTS)) +
  geom_point(aes(color = point_color, alpha = point_alpha), size = 2.5, position = position_jitter(width = 0.25, height = 0)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_color_identity() +
  scale_alpha_identity()

grid.arrange(
  arrangeGrob(
    mainp,
    ncol = 1 
  ),
  arrangeGrob(
    sub1,
    sub2,
    sub3,
    sub4,
    ncol = 1
  ),
  widths = c(3, 1)
)

eda <- data[]
eda$CLUSTER <- as.character(km$cluster)

mainp <- ggplot(eda, aes(x = CREDIT_LIMIT, y = INSTALLMENTS_PURCHASES, color = CLUSTER)) +
  geom_point(size = 2.5, alpha = 0.75) +
  labs(
    title = "Scatter Plot Clusters Distributions",
    subtitle = "Clusters 1 and 3 are more active in making installment purchases"
  ) +
  scale_color_manual(values = color_palette_k, 
                     labels = paste0('Cluster ', 1:5), 
                     breaks = c(1, 2, 3, 4, 5),
                     drop = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

eda$point_color <- ifelse(eda$CLUSTER == 1, color_palette_k[1], "grey")
eda$point_alpha <- ifelse(eda$CLUSTER == 1, 1, 0.1)
sub1 <- ggplot(eda, aes(x = CREDIT_LIMIT, y = INSTALLMENTS_PURCHASES)) +
  geom_point(aes(color = point_color, alpha = point_alpha), size = 2.5) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_color_identity() +
  scale_alpha_identity()

eda$point_color <- ifelse(eda$CLUSTER == 2, color_palette_k[2], "grey")
eda$point_alpha <- ifelse(eda$CLUSTER == 2, 1, 0.1)
sub2 <- ggplot(eda, aes(x = CREDIT_LIMIT, y = INSTALLMENTS_PURCHASES)) +
  geom_point(aes(color = point_color, alpha = point_alpha), size = 2.5) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_color_identity() +
  scale_alpha_identity()

eda$point_color <- ifelse(eda$CLUSTER == 3, color_palette_k[3], "grey")
eda$point_alpha <- ifelse(eda$CLUSTER == 3, 1, 0.1)
sub3 <- ggplot(eda, aes(x = CREDIT_LIMIT, y = INSTALLMENTS_PURCHASES)) +
  geom_point(aes(color = point_color, alpha = point_alpha), size = 2.5) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_color_identity() +
  scale_alpha_identity()

eda$point_color <- ifelse(eda$CLUSTER == 4, color_palette_k[4], "grey")
eda$point_alpha <- ifelse(eda$CLUSTER == 4, 1, 0.1)
sub4 <- ggplot(eda, aes(x = CREDIT_LIMIT, y = INSTALLMENTS_PURCHASES)) +
  geom_point(aes(color = point_color, alpha = point_alpha), size = 2.5) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_color_identity() +
  scale_alpha_identity()

grid.arrange(
  arrangeGrob(
    mainp,
    ncol = 1
  ),
  arrangeGrob(
    sub1,
    sub2,
    sub3,
    sub4,
    ncol = 1
  ),
  widths = c(3, 1)
)


#--Densità dei cluster
eda <- data[]
color_palette_k <- c("#FF6347", "#4682B4", "#33FF57", "#FFD700")
eda$CLUSTER <- as.character(km$cluster)
plots <- list()
for (var in colnames(data)) {
  p <- ggplot(eda, aes_string(x = var, fill = "factor(km$cluster)", color = "factor(km$cluster)")) + 
    geom_density(alpha = 0.1, size= 0.7, adjust = 3) +
    scale_fill_manual(values = color_palette_k) +
    scale_color_manual(values = color_palette_k) +
    labs(y = "",
         x = var,
         fill = NULL,
         color = NULL)+
    theme(axis.text = element_text(size = 5),
          axis.title.x = element_text(size = 6.5),
          legend.position = "top",
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.5, "lines"))
  print(p)
  plots[[var]] <- p
}

grid.arrange(
  plots[[1]], plots[[2]], plots[[3]],plots[[4]],
  plots[[5]], plots[[6]], plots[[7]],plots[[8]],
  ncol = 4
)

grid.arrange(
  plots[[9]], plots[[10]], plots[[11]],plots[[12]],plots[[13]],
  plots[[14]], plots[[15]],plots[[16]],plots[[17]],
  ncol = 5
)
