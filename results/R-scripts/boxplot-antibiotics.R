library(ggplot2)
library(ggpubr)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
out_fig <- "../out_figures/Fig2C.png"
out_fig_violin <- "../out_figures/Fig2B.png"

data_file <- read.csv("../analysis-out/3-Divide_High-Low/high-low.csv")
data_all <- as.data.frame(data_file)



#### Compare different antibiotics
ant_to_compare<-c("Gentamicin","Ampicillin","Meropenem","Vancomycin","Ticarcillin.Clavulanate")


# Select values of recently applied antibiotics
ant_data <- c()

for (ant in ant_to_compare){
  for (row in 1:nrow(data_all)){
    abundance <- data_all[row, "abundance"]
    PTR <- as.numeric(data_all[row, "PTR"])
    dosage <- data_all[row, paste("r_",ant,sep="")]
    ab_level <- as.character(data_all[row, "Abundance.Level"])
    if (dosage>0){
      if (ant=="Ticarcillin.Clavulanate"){
        antibiotics<-"Ticarcillin.Clavulanate"
      }else{
        antibiotics<-ant
      }
      ant_data<-rbind(ant_data,c(abundance, PTR, antibiotics, ab_level))
    }
  }
}

# Select control:
for (row in 1:nrow(data_all)){
  abundance <- data_all[row, "abundance"]
  PTR <- data_all[row, "PTR"]
  if (data_all[row, "Cohort"]=="Control"){
    antibiotics<-ant
    ab_level <- as.character(data_all[row, "Abundance.Level"])
    ant_data<-rbind(ant_data,c(abundance, PTR, "Control", ab_level))
  }
}


colnames(ant_data)<-c("abundance","PTR", "Antibiotics", "Abundance.Level")
ant_data <- data.frame(ant_data)
ant_data$PTR <- as.numeric(as.character(ant_data$PTR))

# Remove medium
ant_data_dropped <- ant_data[ant_data[,"Abundance.Level"]!="medium",]
ant_data_dropped$Antibiotics <- factor(ant_data_dropped$Antibiotics, 
                                       levels=c("Control","Gentamicin", "Vancomycin","Ampicillin","Meropenem","Ticarcillin.Clavulanate" ))

ant_data_dropped$Abundance.Level <- factor(ant_data_dropped$Abundance.Level, levels=c("low","high"))


high_low_plot <- ggboxplot(ant_data_dropped, x="Abundance.Level" , y = "PTR",
                           color = "Abundance.Level", palette = "jco", short.panel.labs = FALSE, facet.by ="Antibiotics", ncol=6) +
 theme(axis.title.x = element_text(size=0, face="bold"), axis.text.x =element_text(size=12, face="bold")) +
  stat_compare_means() + ylim(1,3.5)
high_low_plot
ggsave(out_fig,width=12,height = 4)


high_low_plot <- ggboxplot(ant_data_dropped, x="Abundance.Level" , y = "PTR",
                           color = "Abundance.Level", palette = "jco", short.panel.labs = FALSE, ncol=6) +
  theme(axis.title.x = element_text(size=0, face="bold"), axis.text.x =element_text(size=12, face="bold")) +
  stat_compare_means() + ylim(1,3.5) + facet_grid(. ~ Antibiotics)
high_low_plot
ggsave(out_fig,width=12,height = 4)

# Violin Plot
ant_data_dropped$Antibiotics <- factor(ant_data_dropped$Antibiotics, 
                                       levels=c("Control","Ampicillin","Meropenem","Ticarcillin.Clavulanate","Gentamicin", "Vancomycin" ))


antibiotics_plot <- ggviolin(ant_data_dropped, x="Antibiotics" , y = "PTR",
                              fill = "Antibiotics", palette = "jco", short.panel.labs = FALSE) +
  xlab("") + theme(axis.title.x = element_text(size=0, face="bold"), axis.text.x =element_text(size=10, face="bold")) + 
  stat_compare_means(comparisons=my_comparisons) + geom_boxplot(width=0.1) + theme_minimal(base_size = 12)

antibiotics_plot
ggsave(out_fig_violin,width=10,height = 5)

