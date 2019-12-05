library(magrittr)

PBS1_05mM_r1 = read.delim("Data/PBS1_05mM_r1.txt",
               sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame("index" = as.numeric(colnames(PBS1_05mM_r1)),
           "means" = colMeans(PBS1_05mM_r1),
           "strain" = "PBS1",
           "conc" = "0.5",
           "rep" = "1") %>% 
  write.table("Data_means/PBS1_05mM_r1_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_05mM_r2 = read.delim("Data/PBS1_05mM_r2.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame("index" = as.numeric(colnames(PBS1_05mM_r2)),
           "means" = colMeans(PBS1_05mM_r2),
           "strain" = "PBS1",
           "conc" = "0.5",
           "rep" = "2") %>% 
  write.table("Data_means/PBS1_05mM_r2_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_05mM_r3 = read.delim("Data/PBS1_05mM_r3.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame("index" = as.numeric(colnames(PBS1_05mM_r3)),
           "means" = colMeans(PBS1_05mM_r3),
           "strain" = "PBS1",
           "conc" = "0.5",
           "rep" = "3") %>% 
  write.table("Data_means/PBS1_05mM_r3_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_1mM_r1 = read.delim("Data/PBS1_1mM_r1.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame("index" = as.numeric(colnames(PBS1_1mM_r1)),
           "means" = colMeans(PBS1_1mM_r1),
           "strain" = "PBS1",
           "conc" = "1",
           "rep" = "1") %>% 
  write.table("Data_means/PBS1_1mM_r1_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_1mM_r2 = read.delim("Data/PBS1_1mM_r2.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame("index" = as.numeric(colnames(PBS1_1mM_r2)),
           "means" = colMeans(PBS1_1mM_r2),
           "strain" = "PBS1",
           "conc" = "1",
           "rep" = "2") %>% 
  write.table("Data_means/PBS1_1mM_r2_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_1mM_r3 = read.delim("Data/PBS1_1mM_r3.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame("index" = as.numeric(colnames(PBS1_1mM_r3)),
           "means" = colMeans(PBS1_1mM_r3),
           "strain" = "PBS1",
           "conc" = "1",
           "rep" = "3") %>%
  write.table("Data_means/PBS1_1mM_r3_means.txt",
              sep = "\t",
              quote = F,row.names = F)

PBS1_5mM_r1 = read.delim("Data/PBS1_5mM_r1.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame("index" = as.numeric(colnames(PBS1_5mM_r1)),
           "means" = colMeans(PBS1_5mM_r1),
           "strain" = "PBS1",
           "conc" = "5",
           "rep" = "1") %>%
  write.table("Data_means/PBS1_5mM_r1_means.txt",
              sep = "\t",
              quote = F,row.names = F)

PBS1_5mM_r2 = read.delim("Data/PBS1_5mM_r2.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame("index" = as.numeric(colnames(PBS1_5mM_r2)),
           "means" = colMeans(PBS1_5mM_r2),
           "strain" = "PBS1",
           "conc" = "5",
           "rep" = "2") %>%
  write.table("Data_means/PBS1_5mM_r2_means.txt",
              sep = "\t",
              quote = F,row.names = F)

PBS1_5mM_r3 = read.delim("Data/PBS1_5mM_r3.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame("index" = as.numeric(colnames(PBS1_5mM_r3)),
           "means" = colMeans(PBS1_5mM_r3),
           "strain" = "PBS1",
           "conc" = "5",
           "rep" = "3") %>%
  write.table("Data_means/PBS1_5mM_r3_means.txt",
              sep = "\t",
              quote = F,row.names = F)

PBS1_10mM_r1 = read.delim("Data/PBS1_10mM_r1.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame("index" = as.numeric(colnames(PBS1_10mM_r1)),
           "means" = colMeans(PBS1_10mM_r1),
           "strain" = "PBS1",
           "conc" = "10",
           "rep" = "1") %>%
  write.table("Data_means/PBS1_10mM_r1_means.txt",
              sep = "\t",
              quote = F,row.names = F)

PBS1_10mM_r2 = read.delim("Data/PBS1_10mM_r2.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame("index" = as.numeric(colnames(PBS1_10mM_r2)),
           "means" = colMeans(PBS1_10mM_r2),
           "strain" = "PBS1",
           "conc" = "10",
           "rep" = "2") %>%
  write.table("Data_means/PBS1_10mM_r2_means.txt",
              sep = "\t",
              quote = F,row.names = F)

PBS1_10mM_r3 = read.delim("Data/PBS1_10mM_r3.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame("index" = as.numeric(colnames(PBS1_10mM_r3)),
           "means" = colMeans(PBS1_10mM_r3),
           "strain" = "PBS1",
           "conc" = "10",
           "rep" = "3") %>%
  write.table("Data_means/PBS1_10mM_r3_means.txt",
              sep = "\t",
              quote = F,row.names = F)