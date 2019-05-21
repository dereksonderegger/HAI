# If the population acutally is actually Poisson distributed
# then we really need 30% sampling or more to actually get any accuracy

n <- 1000 # Number of patients
N <- 2   # Number of times to repeat the simulation
lambda_range = seq(.5,8, by=.5)
mixing_range = seq(0, .20, by = .05)
sampling_range = c(.1, .2, .3, .4, .6, .8, 1)
#sampling_range = c(.1, .2)

out <- NULL
for( mixing_frac in mixing_range ){
  for( sample_frac in  sampling_range ){
    for( lambda in lambda_range ){
      for( i in 1:N ){
        Full.Data <- generate_population(n, lambda, mixing_frac)
        Sample.Data <- Full.Data %>% sample_frac(sample_frac, replace = FALSE)

        out <- rbind(out, data.frame(
          lambda=lambda, n=n, sample_frac=sample_frac, mixing_frac=mixing_frac,
          Naive.HAI.rate = calc_HAI( Sample.Data$clusterID ),
          HyperTPoisson.HAI.rate = calc_HAI(clusters=Sample.Data$clusterID, method = 'HypergeometricPoisson',
                                            alpha=sample_frac, N=nrow(Full.Data) ),
          True.HAI.rate = calc_HAI( Full.Data$clusterID ) ))
      }
    }
  }
}
Poisson_Simulation <- out

Poisson_Simulation_CIs <-
  Poisson_Simulation %>%
  group_by(mixing_frac, sample_frac) %>%
  mutate( Obs.HAI.rate.group = as.character( cut(Obs.HAI.rate, breaks = 10) ) ) %>%
  group_by(mixing_frac, sample_frac, Obs.HAI.rate.group) %>%
  summarize( Obs.HAI.rate = mean(Obs.HAI.rate),
             lwr = quantile(True.HAI.rate, 0.025),
             upr = quantile(True.HAI.rate, 0.975) )

usethis::use_data(Poisson_Simulation, overwrite=TRUE)
usethis::use_data(Poisson_Simulation_CIs, overwrite=TRUE)




## Oxford Simulation
N <- 100   # Number of times to repeat the simulation
sampling_range = c(.1, .2, .3, .4, .6, .8, 1)
SNP_Thresholds = c(2, 4, 10)
Oxford_Simulation <- NULL
for( sample_frac in  sampling_range ){
  for( SNP_Threshold in SNP_Thresholds ){
    for( i in 1:N ){
      index <- 1:811 %>% sample( sample_frac * 811 )
      Sampled.Diss <- Oxford.Diss[index, index]

      D <- dist(Sampled.Diss)
      clust <- hclust(D, method='complete')
      Sampled.Data <-
        data.frame( Patient_ID = rownames(Sampled.Diss)) %>%
        mutate( Cluster_ID = cutree(clust, h=SNP_Threshold) )  %>%
        group_by(Cluster_ID) %>% mutate( Cluster_Size = n() )

      Oxford_Simulation <- rbind(Oxford_Simulation, data.frame(
        sample_frac=sample_frac, SNP_Threshold = SNP_Threshold,
        Obs_HAI_rate = calc_HAI( Sampled.Data, 'Cluster_ID' ) ))
    }
  }
}
True_HAIs <- NULL
for( SNP_Threshold in SNP_Thresholds ){
  Oxford.D <- dist(Oxford.Diss)
  clust <- hclust(Oxford.D, method='complete')
  Oxford <- data.frame( Patient_ID = rownames(Oxford.Diss),
                        Cluster_ID = cutree(clust, h=SNP_Threshold) ) %>%
    group_by(Cluster_ID) %>% mutate( Cluster_Size = n() )
  True_HAIs <- rbind( True_HAIs, data.frame(
    SNP_Threshold = SNP_Threshold,
    True_HAI_rate = calc_HAI(Oxford, 'Cluster_ID') ))
}

Oxford_Simulation <- Oxford_Simulation %>%
  left_join(True_HAIs, by='SNP_Threshold' )


Oxford_Simulation_CIs <- Oxford_Simulation %>%
  group_by(sample_frac, SNP_Threshold) %>%
  summarize( lwr = quantile(Obs_HAI_rate, 0.025),
             upr = quantile(Obs_HAI_rate, 0.975) )

usethis::use_data(Oxford_Simulation, overwrite=TRUE)
usethis::use_data(Oxford_Simulation_CIs, overwrite=TRUE)

