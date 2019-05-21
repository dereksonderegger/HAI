library(tidyverse)
library(devtools)
library(readxl)

# ## Testing out how the clustering functions work
# ## Create a distance matrix
# n <- 5
# data <- data.frame( x=runif(n, 0, 10), y=runif(n, 0, 10) ) %>%
#   mutate_all(round)
# D <- dist(data, diag=TRUE, method = 'manhattan')
# ## Now create clusters
# clust <- hclust(D, method='complete')
# plot(clust)
# data$cluster <- cutree(clust, h=4)   # group with distances <= 4


# Now do this with the Oxfordshire data
Oxford.Diss <- read_excel( 'data-raw/eyre_pairwise_snps.xlsx') %>% as.data.frame()
rownames(Oxford.Diss) <- Oxford.Diss$X__1
Oxford.Diss <- Oxford.Diss %>% select( -X__1 )

# Oxford.D <- dist(Oxford.Diss)
# clust <- hclust(Oxford.D, method='complete')
SNP_Thresholds = c(2, 4, 10)
Oxford <- NULL
for( SNP_Threshold in SNP_Thresholds ){
  clust <- dist(Oxford.Diss) %>% hclust(method='complete')
  Oxford <- data.frame( Patient_ID = rownames(Oxford.Diss),
                        Cluster_ID = cutree(clust, h=SNP_Threshold),
                        SNP_Threshold = SNP_Threshold) %>%
    group_by(Cluster_ID) %>% mutate( Cluster_Size = n() ) %>%
    rbind(Oxford, .)
}
Oxford <- Oxford %>% group_by(SNP_Threshold) %>%
    mutate( True_HAI_rate = calc_HAI(Cluster_ID) )


# Save the Oxfordshire data to the package data repository
use_data(Oxford, overwrite = TRUE)
use_data(Oxford.Diss, overwrite = TRUE)


# Just check how the distribution of Oxford Data cluster sizes
ggplot(Oxford, aes(x=Cluster_Size) ) + geom_bar() +
  labs(x='Cluster Size', y='Number of Individuals') +
  ggtitle('Oxfordshire Cluster Sizes')



