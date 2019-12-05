library(magrittr)

PBS1_05mM_r1 = read.delim("Data/PBS1_05mM_r1.txt",
               sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBS1_05mM_r1)),
           colMeans(PBS1_05mM_r1)) %>% 
  write.table("Data_means/PBS1_05mM_r1_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_05mM_r2 = read.delim("Data/PBS1_05mM_r2.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBS1_05mM_r2)),
           colMeans(PBS1_05mM_r2)) %>% 
  write.table("Data_means/PBS1_05mM_r2_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_05mM_r3 = read.delim("Data/PBS1_05mM_r3.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBS1_05mM_r3)),
           colMeans(PBS1_05mM_r3)) %>% 
  write.table("Data_means/PBS1_05mM_r3_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_1mM_r1 = read.delim("Data/PBS1_1mM_r1.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBS1_1mM_r1)),
           colMeans(PBS1_1mM_r1)) %>% 
  write.table("Data_means/PBS1_1mM_r1_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_1mM_r2 = read.delim("Data/PBS1_1mM_r2.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBS1_1mM_r2)),
           colMeans(PBS1_1mM_r2)) %>% 
  write.table("Data_means/PBS1_1mM_r2_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_1mM_r3 = read.delim("Data/PBS1_1mM_r3.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBS1_1mM_r3)),
           colMeans(PBS1_1mM_r3)) %>%
  write.table("Data_means/PBS1_1mM_r3_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_5mM_r1 = read.delim("Data/PBS1_5mM_r1.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBS1_5mM_r1)),
           colMeans(PBS1_5mM_r1)) %>%
  write.table("Data_means/PBS1_5mM_r1_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_5mM_r2 = read.delim("Data/PBS1_5mM_r2.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBS1_5mM_r2)),
           colMeans(PBS1_5mM_r2)) %>%
  write.table("Data_means/PBS1_5mM_r2_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_5mM_r3 = read.delim("Data/PBS1_5mM_r3.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBS1_5mM_r3)),
           colMeans(PBS1_5mM_r3)) %>%
  write.table("Data_means/PBS1_5mM_r3_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_10mM_r1 = read.delim("Data/PBS1_10mM_r1.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBS1_10mM_r1)),
           colMeans(PBS1_10mM_r1)) %>%
  write.table("Data_means/PBS1_10mM_r1_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_10mM_r2 = read.delim("Data/PBS1_10mM_r2.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBS1_10mM_r2)),
           colMeans(PBS1_10mM_r2)) %>%
  write.table("Data_means/PBS1_10mM_r2_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBS1_10mM_r3 = read.delim("Data/PBS1_10mM_r3.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBS1_10mM_r3)),
           colMeans(PBS1_10mM_r3)) %>%
  write.table("Data_means/PBS1_10mM_r3_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBH481_05mM_r1 = read.delim("Data/PBH481_05mM_r1.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBH481_05mM_r1)),
           colMeans(PBH481_05mM_r1)) %>% 
  write.table("Data_means/PBH481_05mM_r1_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBH481_05mM_r2 = read.delim("Data/PBH481_05mM_r2.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBH481_05mM_r2)),
           colMeans(PBH481_05mM_r2)) %>% 
  write.table("Data_means/PBH481_05mM_r2_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBH481_05mM_r3 = read.delim("Data/PBH481_05mM_r3.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBH481_05mM_r3)),
           colMeans(PBH481_05mM_r3)) %>% 
  write.table("Data_means/PBH481_05mM_r3_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBH481_1mM_r1 = read.delim("Data/PBH481_1mM_r1.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBH481_1mM_r1)),
           colMeans(PBH481_1mM_r1)) %>% 
  write.table("Data_means/PBH481_1mM_r1_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBH481_1mM_r2 = read.delim("Data/PBH481_1mM_r2.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBH481_1mM_r2)),
           colMeans(PBH481_1mM_r2)) %>% 
  write.table("Data_means/PBH481_1mM_r2_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBH481_1mM_r3 = read.delim("Data/PBH481_1mM_r3.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBH481_1mM_r3)),
           colMeans(PBH481_1mM_r3)) %>%
  write.table("Data_means/PBH481_1mM_r3_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBH481_5mM_r1 = read.delim("Data/PBH481_5mM_r1.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBH481_5mM_r1)),
           colMeans(PBH481_5mM_r1)) %>%
  write.table("Data_means/PBH481_5mM_r1_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBH481_5mM_r2 = read.delim("Data/PBH481_5mM_r2.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBH481_5mM_r2)),
           colMeans(PBH481_5mM_r2)) %>%
  write.table("Data_means/PBH481_5mM_r2_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBH481_5mM_r3 = read.delim("Data/PBH481_5mM_r3.txt",
                         sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBH481_5mM_r3)),
           colMeans(PBH481_5mM_r3)) %>%
  write.table("Data_means/PBH481_5mM_r3_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBH481_10mM_r1 = read.delim("Data/PBH481_10mM_r1.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBH481_10mM_r1)),
           colMeans(PBH481_10mM_r1)) %>%
  write.table("Data_means/PBH481_10mM_r1_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBH481_10mM_r2 = read.delim("Data/PBH481_10mM_r2.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBH481_10mM_r2)),
           colMeans(PBH481_10mM_r2)) %>%
  write.table("Data_means/PBH481_10mM_r2_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)

PBH481_10mM_r3 = read.delim("Data/PBH481_10mM_r3.txt",
                          sep = "\t", header = T, check.names = F)[ ,3:307]
data.frame(as.numeric(colnames(PBH481_10mM_r3)),
           colMeans(PBH481_10mM_r3)) %>%
  write.table("Data_means/PBH481_10mM_r3_means.txt",
              sep = "\t",
              quote = F,row.names = F,col.names = F)