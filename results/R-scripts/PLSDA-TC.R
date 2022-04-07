setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(mixOmics)
library(snow)
library(base)
library(dplyr)
library(ggpubr)
library(ggplot2)

out_folder <- "."
in_file <- "../analysis-out/4-Select_major_ABR/abr_tim.csv"
out_features <- paste(out_folder, "/features-tc.csv", sep="")
out_features_significant <- paste(out_folder, "/features-tc_significant.csv", sep="")


dataMatrix <- read.csv(in_file, sep=",", header=TRUE)

original_df <- as.data.frame(dataMatrix)
Antibiotics <- original_df[,1]
df <- original_df[,-1]


#----------------------------------------------------
# Separate X and Y
Y <- Antibiotics
X <- df

# Check if data contains missing values
list_na <- colnames(X)[ apply(X, 2, anyNA) ]
list_na
#----------------------------------------------------
# 1. Build general PLSDA model
# Start with 10 components

plsda_model <- plsda(X, Y, ncomp=20, scale=TRUE)
plotIndiv(plsda_model, ind.names = TRUE, ellipse = TRUE, legend=TRUE, scale=TRUE, title="PLS-DA visualization")



#----------------------------------------------------
# 2. Finding optimal number of components
# Cross validate the model using different number of components
# 5-fold cross-validation repeated 10 times
#set.seed(2543) # for reproducibility here, only when the `cpus' argument is not used
#perf.plsda <- perf(plsda_model, validation = "Mfold", folds = 10, 
#                   progressBar = FALSE, auc = TRUE, nrepeat = 10, scale=TRUE) 

# Plot error rates of different number of components using different distance metrics
#plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")
#ggsave(paste(out_folder,"/",common_name,"error_rate.png",sep=""))
# From this plot we can see that optimal number of components is 5

#----------------------------------------------------
# 3. Variable selection

# List of how many features to keep
list.keepX <- c(10,seq(20, 200, 10))

list.keepX

# Use max.dist because it shows the lowest classificational error
tune.splsda.srbct <- tune.splsda(X, Y, ncomp = 5, validation = 'Mfold', folds = 10, 
                                 progressBar = TRUE, dist = 'max.dist', measure = "BER",
                                 test.keepX = list.keepX, nrepeat = 10, cpus = 6, scale=TRUE)


error <- tune.splsda.srbct$error.rate  # error rate per component for the keepX grid

select.keepX <- tune.splsda.srbct$choice.keepX # minimum value of error for selected components

# optimal number of components based on t-tests
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp 
ncomp
# Our optimal number of components is 3

plot(tune.splsda.srbct, col = color.jet(5))

#----------------------------------------------------
# 3. Run sepparation with selected important features:

splsda.srbct <- splsda(X, Y, ncomp = ncomp, keepX = select.keepX, scale=TRUE) 



# Performance of final PLSDA
perf.splsda <- perf(splsda.srbct, validation = "Mfold", folds = 10, 
                    progressBar = FALSE, auc = TRUE, nrepeat = 10, scale=TRUE) 

error_rates <- perf.splsda$error.rate
error_all_comp <- error_rates$overall[ncomp,1]

# Final Model
plotIndiv(splsda.srbct, comp = c(1,2), point.lwd = 1,style = 'ggplot2',
          group = Y, ind.names = TRUE,centroid=FALSE,
          ellipse = TRUE, legend = TRUE,
          ellipse.level = 0.8,
          title = paste("PLS-DA; CV error rate:", error_all_comp, sep=""), scale=TRUE)
ggsave(paste(out_folder,"/1-finalModel.png",sep=""))


# here we match the selected variables to the stable features
#ind.match = match(selectVar(splsda.srbct, comp = 1)$name, 
#                 names(perf.srbct$features$stable[[1]]))
#extract the frequency of selection of those selected variables
#Freq = as.numeric(perf.srbct$features$stable[[1]][ind.match])

important_features1<-data.frame(selectVar(splsda.srbct, comp = 1)$value)
important_features1$component <- 1
important_features1$name <- selectVar(splsda.srbct, comp = 1)$name
#rownames(important_features1) <- NULL

important_features2<-data.frame(selectVar(splsda.srbct, comp = 2)$value)
important_features2$component <- 2
important_features2$name <- selectVar(splsda.srbct, comp = 2)$name

#important_features3<-data.frame(selectVar(splsda.srbct, comp = 3)$value)
#important_features3$component <- 3
#important_features3$name <- selectVar(splsda.srbct, comp = 3)$name

#start_index <- dim(important_features1)[1] +1
#end_index <- dim(important_features2)[1]+dim(important_features1)[1]
#rownames(important_features2) <- c(start_index:end_index)



#important_features3<-data.frame(selectVar(splsda.srbct, comp = 3)$value)
#jpeg(file=paste("PLSDA/", common_name,"-features1.jpeg", sep=""), quality = 3000,width=1300, height=800, res=130)
#plotLoadings(splsda.srbct, comp = 1, title = paste(common_name,': Features component 1', sep=""), 
#             contrib = 'max', method = 'mean')
#dev.off()

#ggsave(paste("PLSDA/", common_name,"-features1.jpeg", sep=""))
#plotLoadings(splsda.srbct, comp = 2, title = paste(common_name, ': Features component 2', sep=""), 
#             contrib = 'max', method = 'mean')

#plotLoadings(splsda.srbct, comp = 3, title = paste(label, ': Features PLSDA component 3', sep=""), 
#             contrib = 'max', method = 'mean')

#all_features <- rbind(important_features1, important_features2, important_features3)
all_features <- rbind(important_features1, important_features2)

rownames(all_features) <- NULL
#all_features <- all_features[all_features[,"value.var"]>0,]

# Perform Wilcoxon t-test on each important gene

# Wilcoxon test
factor_vector <- Antibiotics
p_value <- c()
mean_TC <- c()
mean_Amp.Mero <- c()
diff <- c()
for (species in all_features$name){
  #print(species)
  curr_mean_TC <- mean(original_df[original_df$Treatment=="TIM",species])
  curr_mean_Amp <-  mean(original_df[original_df$Treatment=="AMP/MEM",species])
  diff <- c(diff,curr_mean_TC - curr_mean_Amp)
  test <- wilcox.test(original_df[,species]~factor_vector)
  p_value <- c(p_value, test$p.value)
  mean_TC <- c(mean_TC, curr_mean_TC)
  mean_Amp.Mero <- c(mean_Amp.Mero,curr_mean_Amp)
  }
  
  
p_value_adjusted <- p.adjust(p_value, method="fdr")
all_features$pvalue <- p_value_adjusted
all_features$mean_TC <- mean_TC
all_features$mean_Amp.Mero <- mean_Amp.Mero
all_features$diff <- diff

write.csv(all_features, out_features, row.names = FALSE)

all_features_sign <- all_features[all_features$diff>0,]
all_features_sign <- all_features_sign[all_features_sign$pvalue<0.05,]

write.csv(all_features_sign, out_features_significant, row.names = FALSE)

glm(formula = Treatment ~ wt + disp, family = binomial, data = mtcars)

