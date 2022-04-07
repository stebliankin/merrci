#rm(list = ls())

# Cohort description:
#   Control=0
#   Antibiotic=1


print(Sys.time())
print("Starting Causal Graph Learning...")
# Castalia:
#setwd("/lclhome/vsteb002/RAPToR/B-Causality")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(bnlearn)
library(parallel)

N_TASKS<-48

cl=makeCluster(N_TASKS)
clusterSetRNGStream(cl, N_TASKS)
# Input:
data_file_name <- "../results/analysis-out/4-Select_major_ABR/PTR_species_filtered_metadata_major_AMR_RPKM.csv"
blacklist_file <- "../results/B-out/blacklist.csv"

# Intermediate:
bootfile <- "../results/B-out-RPKM/boot_strength.csv"
filename_dir <- "../results/B-out-RPKM/directed.csv"
filename_undir <- "../results/B-out-RPKM/undirected.csv"
network_csv_file <- "../results/B-out-RPKM/disease_bn.csv"

# Output:
bn_correlation <- "../results/B-out-RPKM/bn_correlation.csv"
con_file_name <- "../results/B-out-RPKM/output_file.xgmml" #output file


graph_id <-  1001

#### Pc stable Algorithm
df <- read.csv(data_file_name)
dades <- lapply(df, as.numeric)
bn_df <- data.frame(dades)
bn_df <- na.omit(bn_df)

blacklist_df <- read.csv(blacklist_file)
bkacklist_df <- data.frame(blacklist_df)

resGS<- pc.stable(bn_df, cluster = cl,whitelist = NULL, blacklist = blacklist_df, test = "zf",
                  alpha = 0.05, B = NULL, debug = F, undirected = FALSE)

####
plot(resGS)

sum_of_col <- data.frame(colnames(bn_df),colSums(bn_df))
col_name <- colnames(bn_df)

#### V-structures
vstruct <- vstructs(resGS, arcs = FALSE, moral = TRUE, debug = FALSE)
write.csv(vstruct,"kera_gingiva_vstructure.csv")

#### End of V-structures
di_arcs_kera_gingiva <- directed.arcs(resGS)
write.csv(di_arcs_kera_gingiva,filename_dir,row.names=F)

undi_arcs_kera_gingiva <- undirected.arcs(resGS)
write.csv(undi_arcs_kera_gingiva,filename_undir,row.names=F)
#stop()
temp_dat <- na.omit(bn_df)

print(Sys.time())
print("Running bootsrap...")

boot_strength <- boot.strength(temp_dat, R = 100, algorithm = "pc.stable", algorithm.args=list(blacklist=blacklist_df),cluster = cl) #Change bootstrap
write.csv(boot_strength,bootfile,row.names=F)
print("Done with bootsrap...")
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
      directed_edge_weight <- append(directed_edge_weight,boot_data[j,4])
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
data_file <- subset(data_file, select = -c(Cohort))

# remove cohort.

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
write.csv(comp_bn_cor,bn_correlation,row.names = F)

data_file_colname <- colnames(data_file)
data_colname <- data.frame(data_file_colname)

node_id <- c()

for (i in 1:nrow(data_colname))
{
  node_id <- append(node_id,i+1000)
}

node_data <- data.frame(data_colname,node_id)

edge_id <- c()
source_id <- c()
target_id <- c()

for (i in 1:nrow(kera_gingiva_net_file))
{
  edge_id <- append(edge_id,i+100)
  for(j in 1:nrow(node_data))
  {
    if(kera_gingiva_net_file[i,1]==node_data[j,1])
    {
      source_id <- append(source_id,node_data[j,2])
      
    }
  }
}

for (i in 1:nrow(kera_gingiva_net_file))
{
  for(j in 1:nrow(node_data))
  {
    if(kera_gingiva_net_file[i,2]==node_data[j,1])
    {
      target_id <- append(target_id,node_data[j,2])
      
    }
  }
}

kera_gingiva_net_file$edgeid <- edge_id
kera_gingiva_net_file$sourceid<- source_id
kera_gingiva_net_file$targetid<- target_id

kera_gingiva_net_file <- na.omit(kera_gingiva_net_file)

##
## xml file writing
late <- c("Selen1","Selen2","Porph1","Porph2","Trepo2")
early <- c("Veill1","Actin2","Strep1","Strep2","Actin1","Strep3","Actin3")

name <- network_csv_file

#### -------------------CoN--------------------------------------------

sink(con_file_name)
cat("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n")
cat(sprintf("<graph id=\"%d\" label=\"%s\" directed=\"1\" cy:documentVersion=\"3.0\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:cy=\"http://www.cytoscape.org\" xmlns=\"http://www.cs.rpi.edu/XGMML\">",graph_id,name))
cat("\n<att name=\"networkMetadata\">")
cat("\n<rdf:RDF>")
cat("\n<rdf:Description rdf:about=\"http://www.cytoscape.org/\">")
cat("\n<dc:type>Protein-Protein Interaction</dc:type>")
cat("\n<dc:description>N/A</dc:description>")
cat("\n<dc:identifier>N/A</dc:identifier>")
cat("\n<dc:date>2018-02-28 19:57:24</dc:date>\n")
cat(sprintf("<dc:title>%s</dc:title>",name))
cat("\n<dc:source>http://www.cytoscape.org/</dc:source>")
cat("\n<dc:format>Cytoscape-XGMML</dc:format>")
cat("\n</rdf:Description>")
cat("\n</rdf:RDF>")
cat("\n </att>")
cat(sprintf("\n<att name=\"shared name\" value=\"%s\" type=\"string\"/>",name))
cat(sprintf("\n<att name=\"name\" value=\"%s\" type=\"string\"/>",name))
cat("\n<att name=\"selected\" value=\"1\" type=\"boolean\"/>")
cat("\n<att name=\"__Annotations\" type=\"list\">")
cat("\n</att>")
cat("\n<att name=\"layoutAlgorithm\" value=\"Hierarchical Layout\" type=\"string\" cy:hidden=\"1\"/>")
cat("\n
    <graphics>
    \n\t<att name=\"NETWORK_CENTER_Y_LOCATION\" value=\"200.0\" type=\"string\"/>
    \n\t<att name=\"NETWORK_SCALE_FACTOR\" value=\"0.3942098543284541\" type=\"string\"/>
    \n\t<att name=\"NETWORK_DEPTH\" value=\"0.0\" type=\"string\"/>
    \n\t<att name=\"NETWORK_WIDTH\" value=\"833.0\" type=\"string\"/>
    \n\t<att name=\"NETWORK_EDGE_SELECTION\" value=\"true\" type=\"string\"/>
    \n\t<att name=\"NETWORK_NODE_SELECTION\" value=\"true\" type=\"string\"/>
    \n\t<att name=\"NETWORK_HEIGHT\" value=\"600.0\" type=\"string\"/>
    \n\t<att name=\"NETWORK_CENTER_X_LOCATION\" value=\"400\" type=\"string\"/>
    \n\t<att name=\"NETWORK_CENTER_Z_LOCATION\" value=\"0.0\" type=\"string\"/>
    \n\t<att name=\"NETWORK_BACKGROUND_PAINT\" value=\"#FFFFFF\" type=\"string\"/>
    \n\t<att name=\"NETWORK_TITLE\" value=\"\" type=\"string\"/>
    \n</graphics>
    ")

#### Writing node information
for(i in 1:nrow(node_data))
{
  cat(sprintf("\n<node id=\"%d\" label=\"%s\">",node_data[i,2],node_data[i,1]))
  cat(sprintf("\n<att name=\"shared name\" value=\"%s\" type=\"string\"/>",node_data[i,2]))
  cat(sprintf("\n<att name=\"name\" value=\"%s\" type=\"string\"/>",node_data[i,1]))
  cat("\n<att name=\"selected\" value=\"0\" type=\"boolean\"/>")
  
  ## early
  for (j in 1:length(early))
  {
    if (node_data[i,1] == early[j])
    {
      for(k in 1:nrow(sum_of_col))
      {
        if (sum_of_col[k:k,1:1]==node_data[i,1])
        {
          abundance <- sum_of_col[k:k,2:2] 
          
        }
      }
      
      if (abundance >= 0)
        abun <- log(abundance,2)*12
      if(abun<=10)
        abun <- 10
      
      cat(sprintf("\n<graphics  outline=\"#9400D3\"  h=\"%f\" w=\"%f\" fill=\"#FFFFFF\" type=\"ELLIPSE\" width=\"5.0\">",abun,abun))
      
      cat(sprintf("
                  
                  \n\t<att name=\"NODE_NESTED_NETWORK_IMAGE_VISIBLE\" value=\"true\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_TRANSPARENCY\" value=\"255\" type=\"string\"/>
                  \n\t<att name=\"NODE_VISIBLE\" value=\"true\" type=\"string\"/>
                  \n\t<att name=\"NODE_DEPTH\" value=\"0.0\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_WIDTH\" value=\"200.0\" type=\"string\"/>
                  \n\t<att name=\"NODE_SELECTED_PAINT\" value=\"#FFFF00\" type=\"string\"/>
                  \n\t<att name=\"NODE_BORDER_TRANSPARENCY\" value=\"255\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_COLOR\" value=\"#9400D3\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_FONT_SIZE\" value=\"20\" type=\"string\"/>
                  \n\t<att name=\"NODE_BORDER_STROKE\" value=\"SOLID\" type=\"string\"/>
                  \n\t<att name=\"NODE_TRANSPARENCY\" value=\"0\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL\" value=\"%s\" type=\"string\"/>
                  
                  ", node_data[i,1]))
      
      cat("\n</graphics>")
      
    }
    
  }
  
  ## late
  for (j in 1:length(late))
  {
    if (node_data[i,1] == late[j])
    {
      for(k in 1:nrow(sum_of_col))
      {
        if (sum_of_col[k:k,1:1]==node_data[i,1])
        {
          abundance <- sum_of_col[k:k,2:2] 
          
        }
      }
      
      if (abundance >= 0)
        abun <- log(abundance,2)*12
      if(abun<=10)
        abun <- 10
      
      cat(sprintf("\n<graphics  outline=\"#FF0033\"  h=\"%f\" w=\"%f\" fill=\"#FFFFFF\" type=\"ELLIPSE\" width=\"5.0\">",abun,abun))
      
      cat(sprintf("
                  
                  \n\t<att name=\"NODE_NESTED_NETWORK_IMAGE_VISIBLE\" value=\"true\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_TRANSPARENCY\" value=\"255\" type=\"string\"/>
                  \n\t<att name=\"NODE_VISIBLE\" value=\"true\" type=\"string\"/>
                  \n\t<att name=\"NODE_DEPTH\" value=\"0.0\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_WIDTH\" value=\"200.0\" type=\"string\"/>
                  \n\t<att name=\"NODE_SELECTED_PAINT\" value=\"#FFFF00\" type=\"string\"/>
                  \n\t<att name=\"NODE_BORDER_TRANSPARENCY\" value=\"255\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_COLOR\" value=\"#FF0033\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_FONT_SIZE\" value=\"20\" type=\"string\"/>
                  \n\t<att name=\"NODE_BORDER_STROKE\" value=\"SOLID\" type=\"string\"/>
                  \n\t<att name=\"NODE_TRANSPARENCY\" value=\"0\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL\" value=\"%s\" type=\"string\"/>
                  
                  ", node_data[i,1]))
      
      cat("\n</graphics>")
      
    }
    
  }
  
  
  # setdiff(x, y)
  
  rest_set <- setdiff(col_name,early)
  rest_set <- setdiff(rest_set,late)
  
  
  ## Rest
  for (j in 1:length(rest_set))
  {
    if (node_data[i,1] == rest_set[j])
    {
      for(k in 1:nrow(sum_of_col))
      {
        if (sum_of_col[k:k,1:1]==node_data[i,1])
        {
          abundance <- sum_of_col[k:k,2:2] 
          
        }
      }
      
      # if (abundance >= 0)
      #   abun <- 25
      # if(abun<=10)
      #   abun <- 25
      
      abun <- 25
      
      
      cat(sprintf("\n<graphics  outline=\"#000000\"  h=\"%f\" w=\"%f\" fill=\"#FFFFFF\" type=\"ELLIPSE\" width=\"5.0\">",abun,abun))
      
      cat(sprintf("
                  
                  \n\t<att name=\"NODE_NESTED_NETWORK_IMAGE_VISIBLE\" value=\"true\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_TRANSPARENCY\" value=\"255\" type=\"string\"/>
                  \n\t<att name=\"NODE_VISIBLE\" value=\"true\" type=\"string\"/>
                  \n\t<att name=\"NODE_DEPTH\" value=\"0.0\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_WIDTH\" value=\"200.0\" type=\"string\"/>
                  \n\t<att name=\"NODE_SELECTED_PAINT\" value=\"#FFFF00\" type=\"string\"/>
                  \n\t<att name=\"NODE_BORDER_TRANSPARENCY\" value=\"255\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_COLOR\" value=\"#000000\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL_FONT_SIZE\" value=\"20\" type=\"string\"/>
                  \n\t<att name=\"NODE_BORDER_STROKE\" value=\"SOLID\" type=\"string\"/>
                  \n\t<att name=\"NODE_TRANSPARENCY\" value=\"0\" type=\"string\"/>
                  \n\t<att name=\"NODE_LABEL\" value=\"%s\" type=\"string\"/>
                  
                  ", node_data[i,1]))
      
      cat("\n</graphics>")
      
    }
    
  }
  
  
  ### End of node graphics
  
  cat("\n</node>")
}

#### Writing edge information

cor_wt <- comp_bn_cor$spearman
kera_gingiva_net_file$spearman <- cor_wt
kera_gingiva_net_file <- na.omit(kera_gingiva_net_file)


for (i in 1:nrow(kera_gingiva_net_file))
{
  cat(sprintf("\n<edge id=\"%d\" label=\"%s (pp) %s\" source=\"%d\" target=\"%d\" cy:directed=\"1\">",kera_gingiva_net_file[i,5],kera_gingiva_net_file[i,1],kera_gingiva_net_file[i,2],kera_gingiva_net_file[i,6],kera_gingiva_net_file[i,7]))
  cat(sprintf("\n\t<att name=\"shared name\" value=\"%s (pp) %s\" type=\"string\"/>",kera_gingiva_net_file[i,1],kera_gingiva_net_file[i,2]))   
  cat("\n\t<att name=\"shared interaction\" value=\"pp\" type=\"string\"/>")
  cat(sprintf("\n\t<att name=\"name\" value=\"%s (pp) %s\" type=\"string\"/>",kera_gingiva_net_file[i,1],kera_gingiva_net_file[i,2]))
  cat("\n\t<att name=\"selected\" value=\"0\" type=\"boolean\"/>")
  cat("\n\t<att name=\"interaction\" value=\"pp\" type=\"string\"/>")
  
  pos_flag <- 0
  neg_flag <- 0
  
  if(kera_gingiva_net_file[i,3]==TRUE)
  {
    cat("\n\t<att name=\"Directed\" value=\"1\" type=\"boolean\"/>")
    
  }
  if(kera_gingiva_net_file[i,3]==FALSE)
  {
    cat("\n\t<att name=\"Directed\" value=\"0\" type=\"boolean\"/>")
  }
  
  
  if (kera_gingiva_net_file[i,8]>0)
  {
    edge_wt <- kera_gingiva_net_file[i,8]
    pos_flag <-1
    
  }
  
  else if (kera_gingiva_net_file[i,8]<0)
  {
    edge_wt <- kera_gingiva_net_file[i,8]*(-1)
    neg_flag <- 1
  }
  
  
  cat(sprintf("\n\t<att name=\"weight\" value=\"%f\" type=\"real\"/>",cor_wt[i]))
  
  
  if (pos_flag == 1)
    #cat(sprintf("\n<graphics width=\"%f\" fill=\"#00CC00\">", edge_wt*12))
    cat(sprintf("\n<graphics width=\"%f\" fill=\"#00CC00\">", 2.75))
  
  if (neg_flag==1 )
    #cat(sprintf("\n<graphics width=\"%f\" fill=\"#FF0000\">", edge_wt*12))
    cat(sprintf("\n<graphics width=\"%f\" fill=\"#FF0000\">", 2.75))
  
  if(kera_gingiva_net_file[i,3]==TRUE)
  {
    
    cat("\n\t<att name=\"EDGE_TARGET_ARROW_SHAPE\" value=\"DELTA\" type=\"string\"/>")
    cat("\n\t<att name=\"EDGE_TARGET_ARROW_UNSELECTED_PAINT\" value=\"#000000\" type=\"string\"/>")
    
    
  }
  
  #<att name="EDGE_TRANSPARENCY" value="100" type="string"/>
  #cat(sprintf("\n\t<att name=\"EDGE_TRANSPARENCY\" value=\"%f\" type=\"string\"/>",kera_gingiva_net_file[i,4]*255))
  
  
  # cat(sprintf("\n\t<att Transparency=\"%f\"/>",kera_gingiva_net_file[i,4]*255))
  
  cat("\n</graphics>")
  
  cat("\n</edge>")
  
}

cat("\n</graph>")

sink()
print(Sys.time())
print("Done...")


