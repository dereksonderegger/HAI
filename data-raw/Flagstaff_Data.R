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
Flagstaff.Diss <- read_excel( 'data-raw/flinn_pairwise_snps.xlsx') %>% as.data.frame()
Flagstaff.Diss <- Flagstaff.Diss %>% select( -...1 )

# Oxford.D <- dist(Oxford.Diss)
# clust <- hclust(Oxford.D, method='complete')
SNP_Thresholds = c(2, 4, 10)
Flagstaff <- NULL
for( SNP_Threshold in SNP_Thresholds ){
  clust <- dist(Flagstaff.Diss) %>% hclust(method='complete')
  Flagstaff <- data.frame( Patient_ID = rownames(Flagstaff.Diss),
                           Cluster_ID = cutree(clust, h=SNP_Threshold),
                           SNP_Threshold = SNP_Threshold) %>%
    group_by(Cluster_ID) %>% mutate( Cluster_Size = n() ) %>%
    rbind(Flagstaff, .)
}

Flagstaff <-
  Flagstaff %>% group_by(SNP_Threshold) %>%
    do({
      calc_HAI(.$Cluster_ID, alpha = 1, bias = FALSE, CI = FALSE) %>%
        as.data.frame() %>% t() %>% as.data.frame() %>%
        select(-lwr1, -lwr2, -upr1, -upr2) %>%
        rename(True_HAI_rate = est)
      }) %>%
  left_join( Flagstaff, .)


# Save the Oxfordshire data to the package data repository
use_data(Flagstaff, overwrite = TRUE)
use_data(Flagstaff.Diss, overwrite = TRUE)


# Just check how the distribution of Oxford Data cluster sizes
ggplot(Flagstaff, aes(x=Cluster_Size) ) + geom_bar() +
  facet_grid( SNP_Threshold ~ . ) +
  labs(x='Cluster Size', y='Number of Individuals') +
  ggtitle('Flagstaff Cluster Sizes')



