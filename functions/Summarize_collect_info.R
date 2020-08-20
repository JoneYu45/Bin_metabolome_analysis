#Environment setup
library(ggplot2)
couunt_element <- function(data){
  #Make output table
  elements <- array(NA, dim = c(nrow(data), 7))
  colnames(elements) <- c('C', 'H', 'O', 'N', 'P', 'S', 'Cl')
  
  #Collect element number
  for (i in 1:nrow(data)) {
    target <- unlist(strsplit(as.character(data[i,3]), split = ' '))
    
    #Searching C
    if (T %in% grepl('C[^a-zA-Z]', target)) {
      C <- target[grep('C[^a-zA-Z]', target)]
      if (grepl('[0-9]+', C)){
        C = as.numeric(sub('C', '', C))
      }else{
        C = 1
      }
    }else{
      C = 0
    }
    
    #Searching H
    if (T %in% grepl('H', target)) {
      H <- target[grep('H', target)]
      if (grepl('[0-9]+', H)){
        H = as.numeric(sub('H', '', H))
      }else{
        H = 1
      }
    }else{
      H = 0
    }
    
    #Searching O
    if (T %in% grepl('O', target)) {
      O <- target[grep('O', target)]
      if (grepl('[0-9]+', O)){
        O = as.numeric(sub('O', '', O))
      }else{
        O = 1
      }
    }else{
      O = 0
    }
    
    #Searching N
    if (T %in% grepl('N[^a-zA-Z]', target)) {
      N <- target[grep('N[^a-zA-Z]', target)]
      if (grepl('[0-9]+', N)){
        N = as.numeric(sub('N', '', N))
      }else{
        N = 1
      }
    }else{
      N = 0
    }
    
    #Searching P
    if (T %in% grepl('P', target)) {
      P <- target[grep('P', target)]
      if (grepl('[0-9]+', P)){
        P = as.numeric(sub('P', '', P))
      }else{
        P = 1
      }
    }else{
      P = 0
    }
    
    #Searching S
    if (T %in% grepl('S', target)) {
      S <- target[grep('S', target)]
      if (grepl('[0-9]+', S)){
        S = as.numeric(sub('S', '', S))
      }else{
        S = 1
      }
    }else{
      S = 0
    }
    
    #Searching Cl
    if (T %in% grepl('Cl', target)) {
      Cl <- target[grep('Cl', target)]
      if (grepl('[0-9]+', Cl)){
        Cl = as.numeric(sub('Cl', '', Cl))
      }else{
        Cl = 1
      }
    }else{
      Cl = 0
    }
    
    
    elements[i, 1] = C
    elements[i, 2] = H
    elements[i, 3] = O
    elements[i, 4] = N
    elements[i, 5] = P
    elements[i, 6] = S
    elements[i, 7] = Cl
  }
  
  #Output
  return(elements)
}

#Data input
workbook <- 'D:/PythonProject/Bin_metabolome_analysis/outputs/YT9.csv'
rawdata <- read.csv(workbook, row.names = 1)
remove(workbook)

#Summarize all useful compounds
data <- rawdata[which(rawdata[,7] != '__'), c(1,4:7)]
collect_table <- array(NA, dim = c(1,5))
colnames(collect_table) <- c('orf', 'compound', 'formula', 'direct_parent', 'super_class')

for (i in 1:nrow(data)) {
  #Process formula info
  C_num <- sub('___', '', as.character(data[i,2]))
  all_C_num <- strsplit(C_num, split = '_')[[1]]
  #Process formula info
  formula <- sub('___', '', as.character(data[i,3]))
  all_formula <- strsplit(formula, split = '_')[[1]]
  #Process direct parent info
  direct_parent <- sub('___', '', as.character(data[i,4]))
  all_direct_parent <- strsplit(direct_parent, split = '_')[[1]]
  #Process super class info
  super_class <- sub('___', '', as.character(data[i,5]))
  all_super_class <- strsplit(super_class, split = '_')[[1]]
  
  for (j in 1:length(all_formula)) {
    new_info <- cbind(as.character(data[i,1]), all_C_num[j], all_formula[j], 
                      all_direct_parent[j], all_super_class[j])
    colnames(new_info) <- c('orf', 'compound', 'formula', 'direct_parent', 'super_class')
    
    collect_table <- rbind(collect_table, new_info)
  }
}
collect_table <- collect_table[-1,]

#Analyze compound formula
#Preprocess collect_table
for (i in 1:nrow(collect_table)) {
  collect_table[i,3] <- sub('H', ' H', collect_table[i,3])
  collect_table[i,3] <- sub('O', ' O', collect_table[i,3])
  collect_table[i,3] <- sub('N', ' N', collect_table[i,3])
  collect_table[i,3] <- sub('P', ' P', collect_table[i,3])
  collect_table[i,3] <- sub('S', ' S', collect_table[i,3])
  collect_table[i,3] <- sub('Cl', ' Cl', collect_table[i,3])
  collect_table[i,3] <- sub('\\+', '', collect_table[i,3])
  collect_table[i,3] <- sub('Ca', ' Ca', collect_table[i,3])
  collect_table[i,3] <- sub('D', ' D', collect_table[i,3])
  collect_table[i,3] <- sub('F', ' F', collect_table[i,3])
  collect_table[i,3] <- sub('Na', ' Na', collect_table[i,3])
  collect_table[i,3] <- sub('I', ' I', collect_table[i,3])
  collect_table[i,3] <- sub('Mg', ' Mg', collect_table[i,3])
  collect_table[i,3] <- sub('\\-', ' \\-', collect_table[i,3])
  collect_table[i,3] <- sub('Br', ' Br', collect_table[i,3])
}

#Collect elements info
elements <- as.data.frame(couunt_element(collect_table))
##For ACE
ACE <- array(c('ACE', '_', 'C4 H5 N O4 S'), dim = c(1,3))
ACE_element <- couunt_element(ACE)
ACE_element <- as.data.frame(ACE_element)
colnames(ACE_element) <- c('C', 'H', 'O', 'N', 'P', 'S', 'Cl')

#Calculate molecular property
molecular_property <- array(NA, dim = c(nrow(collect_table), 3))
colnames(molecular_property) <- c('HC', 'OC', 'AI')

for (i in 1:nrow(molecular_property)) {
  molecular_property[i,1] <- elements$H[i]/elements$C[i]
  molecular_property[i,2] <- elements$O[i]/elements$C[i]
  DBE = (1+elements$C[i]-elements$O[i]-elements$S[i]-0.5*elements$H[i])
  CAI = (elements$C[i]-elements$O[i]-elements$S[i]-elements$N[i]-elements$P[i])
  if ({DBE <= 0} | {CAI <= 0}) {
    molecular_property[i,3] = 0
  }else{
    molecular_property[i,3] = DBE/CAI
  }                                     
}

##For ACE
ACE_molecular_property <- array(NA, dim = c(1, 3))
colnames(ACE_molecular_property) <- c('HC', 'OC', 'AI')

for (i in 1:nrow(ACE_molecular_property)) {
  ACE_molecular_property[i,1] <- ACE_element$H[i]/ACE_element$C[i]
  ACE_molecular_property[i,2] <- ACE_element$O[i]/ACE_element$C[i]
  DBE = (1+ACE_element$C[i]-ACE_element$O[i]-ACE_element$S[i]-0.5*ACE_element$H[i])
  CAI = (ACE_element$C[i]-ACE_element$O[i]-ACE_element$S[i]-ACE_element$N[i]-ACE_element$P[i])
  if ({DBE <= 0} | {CAI <= 0}) {
    ACE_molecular_property[i,3] = 0
  }else{
    ACE_molecular_property[i,3] = DBE/CAI
  }                                     
}

#Plot the metabolome map
input <- as.data.frame(molecular_property)
ACE_molecular_property <- as.data.frame(ACE_molecular_property)
#HC/OC
ggplot(input, aes(x=OC, y=HC))+
  geom_point(shape=21, size=3, fill='red')+
  geom_vline(xintercept=.9, linetype="dashed",color = "green4", size=1)+
  geom_hline(yintercept=1.5, linetype="dashed",color = "green4", size=1)+
  geom_hline(yintercept=2, linetype="dashed",color = "green4", size=1)+
  geom_point(data=ACE_molecular_property, mapping=aes(x=OC, y=HC), 
             shape=21, size=5, fill='blue', alpha=.5)

#HC/AI
ggplot(input, aes(x=AI, y=HC))+
  geom_point(shape=21, size=3, fill='red')+
  geom_vline(xintercept=.9, linetype="dashed",color = "green4", size=1)+
  geom_hline(yintercept=1.5, linetype="dashed",color = "green4", size=1)+
  geom_hline(yintercept=2, linetype="dashed",color = "green4", size=1)+
  geom_point(data=ACE_molecular_property, mapping=aes(x=AI, y=HC), 
             shape=21, size=5, fill='blue', alpha=.5)

#Summarize chemical taxonomy
direct_parent_unique <- unique(collect_table[,4])
chemical_taxonomy <- array(NA, dim = c(length(direct_parent_unique), 3))
colnames(chemical_taxonomy) <- c('Direct_parent', 'Super_class', 'Substrate_number')

chemical_taxonomy[,1] <- direct_parent_unique
for(i in 1:nrow(chemical_taxonomy)) {
  loc <- which(collect_table[,4] == chemical_taxonomy[i,1])
  chemical_taxonomy[i,2] <- collect_table[loc[1], 5]
  chemical_taxonomy[i,3] <- length(loc)
}

#Plot
input <- as.data.frame(chemical_taxonomy)
input[,3] <- as.numeric(as.character(input[,3]))
input <- input[order(input$Super_class, input$Direct_parent),]

ggplot(input, aes(x=Direct_parent, y=Substrate_number, fill=Super_class))+
  geom_bar(stat="identity")+
  scale_x_discrete(limits = unique(as.character(input$Direct_parent)))+
  coord_flip()
