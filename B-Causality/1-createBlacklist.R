# Objective:
#   Create a set of blacklist. Restrict edges coming to clinical variables

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data_file_name <- "../results/analysis-out/4-Select_major_ABR/PTR_species_filtered_metadata_major_AMR.csv"

out_blacklist <- "blacklist.csv"

df <- read.csv(data_file_name)

clinical_var_list <- c("Day_of_Life", "PostMenst_Age", "Gestational_Age",
                       "Birthweight", "Gentamicin", "Cefazolin","Ampicillin", "Trimethoprim.Sulfamathoxazole", "Meropenem",
                       "Vancomycin", "Ticarcillin.Clavulanate", "Clindamycin", "Cefotaxime", "Total_abx", "r_Gentamicin",
                       "r_Meropenem", "r_Ticarcillin.Clavulanate", "r_Vancomycin", "r_Ampicillin",
                       "r_Cefotaxime","r_TOTAL","Human_Milk","Maternal_Milk", "Donor_Milk", "Formula","Fortification","Vitamin_A",
                       "Caffeine","Iron","Furosemide_Lasix","m_ampicillin","m_ceftriaxone","m_azithromycin",
                       "m_amoxicillin", "m_cefazolin","m_erythromycin","m_gentamicin","m_penicillin","m_vancomycin",
                       "m_clindamycin","m_cefotaxime", "dur_membrane_rupture","Total.Antibiotic.Days", "Cohort", "CRIB.II.Score")
#clinical_var_list <- c("Day_of_Life","PostMenst_Age	Gestational_Age	Birthweight	Gentamicin	Cefazolin	Ampicillin	Trimethoprim-Sulfamathoxazole	Meropenem	Vancomycin	Ticarcillin-Clavulanate	Clindamycin	Cefotaxime	Total_abx	r_Gentamicin	r_Meropenem	r_Ticarcillin-Clavulanate	r_Vancomycin	r_Ampicillin	r_Cefotaxime	r_TOTAL	Human_Milk	Maternal_Milk	Donor_Milk	Formula	Fortification	Vitamin_A	Caffeine	Iron	Furosemide_Lasix	m_ampicillin	m_ceftriaxone	m_azithromycin	m_amoxicillin	m_cefazolin	m_erythromycin	m_gentamicin	m_penicillin	m_vancomycin	m_clindamycin	m_cefotaxime	dur_membrane_rupture	CRIB II Score	Total Antibiotic Days")
cols <- colnames(df)

# Restrict edges from all columns to Clinical Variables
restricted_vec <- c()
for(from in cols){
  for(to in cols){
    if(to %in% clinical_var_list){
      restricted_vec <- rbind(restricted_vec, c(from, to))
    }
  }
}

# Restrict edges between averagePTR and PTR of any taxa
#for(col in cols){
#  print(col)
#  if(grepl(".PTR", col) | grepl(".abundance", col)){
#    restricted_vec <- rbind(restricted_vec, c(col, "AveragePTR"))
#    restricted_vec <- rbind(restricted_vec, c("AveragePTR", col))
#  }
#}

# Restrict edges between PTR-PTR
for (col1 in cols){
  if (grepl("PTR",col1)){
    taxa1 <- gsub(".PTR","", col1)
    for (col2 in cols){
      # PTR-PTR
      if (grepl("PTR",col2)){
        restricted_vec <- rbind(restricted_vec, c(col1, col2))
      # PTR - abundace
      }else if(grepl("abundance",col2)){
        taxa2 <- gsub(".abundance","",col2)
        if (taxa1!=taxa2){
          restricted_vec <- rbind(restricted_vec, c(col1, col2))
        }
      }else if(grepl("ARO",col2)){
        restricted_vec <- rbind(restricted_vec, c(col1, col2))
      }
      
    }
  }
}


for (col1 in cols){
  if (grepl("abundance",col1)){
    for (col2 in cols){
      if (grepl("PTR",col2)){
        restricted_vec <- rbind(restricted_vec, c(col1, col2))
      }
      
    }
  }
}

# Restrinct gene-gene
for (col1 in cols){
  if (grepl("ARO.",col1)){
    for (col2 in cols){
      if (grepl("ARO.",col2)){
        restricted_vec <- rbind(restricted_vec, c(col1, col2))
      }
      
    }
  }
}
df <- data.frame(restricted_vec)
colnames(df)<- c("from", "to")

write.csv(df, out_blacklist, row.names = FALSE)

