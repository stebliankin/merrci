# Compute Correlations for Causal Graph

# Input:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data_file_name <- "../results/analysis-out/4-Select_major_ABR/PTR_species_filtered_metadata_major_AMR.csv"

bootfile <- "boot_strength.csv"
filename_dir <- "directed.csv"
filename_undir <- "undirected.csv"
network_csv_file <- "disease_bn.csv"

# Output:
bn_correlation <- "bn_correlation.csv"
con_file_name <- "output_file.xgmml" #output file


#### Drop out of redundant edges
data <- read.csv(filename_undir,header = TRUE, colClasses=c("from"="character","to"="character"))
bn_df <- data.frame(data)


redundant_undi_edge <- c()

for (i in 1:((nrow(bn_df))))
{
  for (j in 1:nrow(bn_df))
  {
    if(i!=j & !(i %in% redundant_undi_edge))
    {
      if(bn_df[i,1]==bn_df[j,2] & bn_df[i,2]==bn_df[j,1])
      {
        redundant_undi_edge <- append(redundant_undi_edge,j)
      }
      
    }
  }
  
}

#### Boot strength

boot_data <- read.csv(bootfile,header = TRUE, colClasses=c("from"="character","to"="character","strength"="double","direction"="double"))


#### Undirected Edge
undirected_edge <- bn_df [c(redundant_undi_edge),]
undirected_edge$directed <- rep(FALSE,nrow(undirected_edge))

undirected_weights <- c()

for (i in 1:nrow(undirected_edge))
{
  for(j in 1:nrow(boot_data))
  {
    if (undirected_edge[i,1]==boot_data[j,1] & undirected_edge [i,2] == boot_data [j,2])
    {
      undirected_weights <- append(undirected_weights,boot_data[j,3])
    }
    
  }
  
}

undirected_edge$weight <- undirected_weights

#### Directed Edge

directed_edge <- read.csv(filename_dir,header = TRUE, colClasses=c("from"="character","to"="character"))
directed_edge <- data.frame(directed_edge)

directed_edge_weight <- c()



for (i in 1:nrow(directed_edge))
{
  for( j in 1:nrow(boot_data))
  {
    if(directed_edge[i,1]== boot_data[j,1] & directed_edge[i,2]== boot_data[j,2])
    {
      directed_edge_weight <- append(directed_edge_weight,boot_data[j,3])
    }
  }
  
}

directed_edge$directed <- rep(TRUE,nrow(directed_edge))
directed_edge$weight <- directed_edge_weight


#### Writing Network file
kera_gingiva_net_file <- rbind(directed_edge,undirected_edge)

#### Writing xgmml file for cytosacpe-backend

data_file <- read.csv(data_file_name,header=TRUE)
#data_file <- data.frame(data_file)
#data_file <- lapply(data_file, as.numeric)
#### Correlation Data
cor_value <- cor(data_file, method = "pearson")
cv_col <- colnames(cor_value)
cv_row <- rownames(cor_value)

correlation <- c()

for(i in 1:nrow(kera_gingiva_net_file))
{
  for (j in 1: length(cv_col))
    
  {
    if (kera_gingiva_net_file[i,1]==cv_col[j])
    {
      x = j
    }
    
  }  
  
  for (k in 1: length(cv_row))
    
  {
    if (kera_gingiva_net_file[i,2]==cv_row[k])
    {
      y = k
    }
    
  } 
  
  correlation <- append(correlation,cor_value[x,y])
  
}

comp_bn_cor <- kera_gingiva_net_file[,1:4]
comp_bn_cor$pearson <- correlation

###
spearman_cor <- cor(data_file, method = "spearman")
spearman <- c()

for(i in 1:nrow(kera_gingiva_net_file))
{
  for (j in 1: length(cv_col))
    
  {
    if (kera_gingiva_net_file[i,1]==cv_col[j])
    {
      x = j
    }
    
  }  
  
  for (k in 1: length(cv_row))
    
  {
    if (kera_gingiva_net_file[i,2]==cv_row[k])
    {
      y = k
    }
    
  } 
  
  spearman <- append(spearman,spearman_cor[x,y])
  
}
comp_bn_cor$spearman <- spearman

###
kendall_cor <- cor(data_file, method = "kendall")
kendall <- c()

for(i in 1:nrow(kera_gingiva_net_file))
{
  for (j in 1: length(cv_col))
    
  {
    if (kera_gingiva_net_file[i,1]==cv_col[j])
    {
      x = j
    }
    
  }  
  
  for (k in 1: length(cv_row))
    
  {
    if (kera_gingiva_net_file[i,2]==cv_row[k])
    {
      y = k
    }
    
  } 
  
  kendall <- append(kendall,kendall_cor[x,y])
  
}
comp_bn_cor$kendall <- kendall
write.csv(comp_bn_cor,"bn_correlation.csv",row.names = F)
