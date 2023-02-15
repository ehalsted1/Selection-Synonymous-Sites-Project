
####Frequencies

##all
mut_count_all<-read.delim("/home/user/Desktop/syn_proyect_trial3/Output/output_syn_3_100_m1,m2,m3,m4,m5_count.sfs")
mut_count_all[,1]<-NULL #deleting the column names of the samples

#Sum of each column
for (i in 1:99){
        mut_count_all$sites[i]<-sum(mut_count_all[,i])
}

#Site frequency spectrum
barplot(mut_count_all$sites[1:98], main = "All mutations (m1-m5)",
        xlab = "Sites", ylab = "Frequencies",
        col = colors()[10:60])


##Non-coding mutations (m1)
mut_count_m1<-read.delim("/home/user/Desktop/lab/syn_proyect_trial2/Output/output_syn_2_100_m1_count.sfs")
mut_count_m1[,1]<-NULL #deleting the column names of the samples

#Sum of each column
for (i in 1:99){
        mut_count_m1$sites[i]<-sum(mut_count_m1[,i])
        }

#Site frequency spectrum
barplot(mut_count_m1$sites[1:98], main = "Non-coding mutations (m1)",
        xlab = "Sites", ylab = "Frequencies",
        col = colors()[10:60])


sites<-mut_count_m1$sites
write.table(sites, "sitesm1.txt")

##Synonimous neutral mutations (m2)
mut_count_m2<-read.delim("/home/user/Desktop/syn_proyect_trial2/Output/output_syn_2_100_m2_count.sfs")
mut_count_m2[,1]<-NULL #deleting the column names of the samples

#Sum of each column
for (i in 1:99){
        mut_count_m2$sites[i]<-sum(mut_count_m2[,i])
}

#Site frequency spectrum
barplot(mut_count_m2$sites[1:98], main = "Synonimous neutral mutations (m2)",
        xlab = "Sites", ylab = "Frequencies",
        col = colors()[10:60])


##Deleteripos mutations (m3)
mut_count_m3<-read.delim("/home/user/Desktop/lab/syn_proyect_trial2/Output/output_syn_2_100_m3_count.sfs")
mut_count_m3[,1]<-NULL #deleting the column names of the samples

#Sum of each column
for (i in 1:99){
        mut_count_m3$sites[i]<-sum(mut_count_m3[,i])
}

sitesm3<-mut_count_m3$sites
write.table(sitesm3, "sitesm3.txt")

#Site frequency spectrum
barplot(mut_count_m3$sites[1:98], main = "Deleteripos mutations (m3)",
        xlab = "Sites", ylab = "Frequencies",
        col = colors()[10:60])

##Beneficial mutations (m4)
mut_count_m4<-read.delim("/home/user/Desktop/syn_proyect_trial2/Output/output_syn_2_100_m4_count.sfs")
mut_count_m4[,1]<-NULL #deleting the column names of the samples

#Sum of each column
for (i in 1:99){
        mut_count_m4$sites[i]<-sum(mut_count_m4[,i])
}

#Site frequency spectrum
barplot(mut_count_m4$sites[1:98], main = "Beneficial mutations (m4)",
        xlab = "Sites", ylab = "Frequencies",
        col = colors()[10:60])

##Synonymous deleterious mutations (m5)
mut_count_m5<-read.delim("/home/user/Desktop/lab/syn_proyect_trial2/Output/output_syn_2_100_m5_count.sfs")
mut_count_m5[,1]<-NULL #deleting the column names of the samples

#Sum of each column
for (i in 1:99){
        mut_count_m5$sites[i]<-sum(mut_count_m5[,i])
}


sitesm5<-mut_count_m5$sites
write.table(sitesm5, "sites_m5.txt")


#Site frequency spectrum
barplot(mut_count_m5$sites[1:98], main = "Synonymous deleterious mutations (m5)",
        xlab = "Sites", ylab = "Frequencies",
        col = colors()[10:60])












## Binding datasets
mut_freq<-cbind(mut_freq_all, mut_freq_m1, mut_freq_m2, mut_freq_m3, mut_freq_m4, mut_freq_m5)

mut_freq_table <- table(mut_freq$all, mut_freq$m1, mut_freq$m2, mut_freq$m3, mut_freq$m4, mut_freq$m5)

## Barplot
library(dplyr)
library(ggplot2)



mut_freq<-t(mut_freq)

spineplot(t(mut_freq))


frequencies <- barplot(as.matrix(mut_freq), xlab='Individuals',  col = colors()[10:60],
              main = "Frequencies of all types of mutations", ylab='Frequencies', las=1, beside=T)

text(x=frequencies, y=genotype, pos=3, cex=0.8, col="red",
     label=round(genotype, 4))

mut_freq_table

#### SFS colours
library(viridis)


#Cutting the last 600 individuals
first_ind <- mut_freq[, 1:335]

op <- par(cex = 0.4)
barplot(first_ind, main="Site frequency Spectrum", ylab="Frequencies", xlab="Individuals", col=ViridisColors, legend = c("all", "m1", "m2", "m3", "m4", "m5"), beside=TRUE, cex.lab=2, cex.axis=2,cex.names=2, args.legend = list(x = "top",cex=1.5, bty = "n"), cex.main = 2, ylim = c(0,0.5), yaxt="n",cex = 0.5, border=NA)
barplot(log10(counts) - log10(0.0000000001))










#Counts
mut_counts<- read.delim("/home/user/Desktop/syn_proyect_trial1/Output_slim_1/output_syn_2_1000_m1,m2,m3,m4,m5_count.sfs")
mut_counts<-as.data.frame(t(mut_counts))
sum(mut_counts$V1)



as.numeric(mut_counts$V1)

barplot(mut_counts$V1, main = "Counts of all types of mutations",
        xlab = "Individuals", ylab = "Counts",
        col = colors()[10:60])


mut1_counts<- read.delim("/home/user/Desktop/syn_proyect_trial1/Output_slim_1/output_syn_2_1000_m1_count.sfs")
mut1_counts<-as.data.frame(t(mut1_counts))
sum(mut1_counts$V1)

mut2_counts<- read.delim("/home/user/Desktop/syn_proyect_trial1/Output_slim_1/output_syn_2_1000_m2_count.sfs")
mut2_counts<-as.data.frame(t(mut2_counts))
sum(mut2_counts$V1)

mut3_counts<- read.delim("/home/user/Desktop/syn_proyect_trial1/Output_slim_1/output_syn_2_1000_m3_count.sfs")
mut3_counts<-as.data.frame(t(mut3_counts))
sum(mut3_counts$V1)

mut4_counts<- read.delim("/home/user/Desktop/syn_proyect_trial1/Output_slim_1/output_syn_2_1000_m4_count.sfs")
mut4_counts<-as.data.frame(t(mut4_counts))
sum(mut4_counts$V1)


#### SFS colours
library(viridis)

ViridisColors <- viridis(5)
counts <- table(mtcars$vs, mtcars$gear)
counts <- rbind(counts,c(1,2,3))
counts <- rbind(counts,c(1,2,3))
counts <- rbind(counts,c(1,2,3))


for (i in 4:20){
counts <- cbind(counts,c(1,1,1,1,1))
}

counts <- cbind(counts,c(1,1,1,1,1))
counts[1,1] <- 1
counts[1,2] <- 1
counts[1,3] <- 1
counts[4,1] <- 1
counts[4,2] <- 1
counts[4,3] <- 1
counts[5,1] <- 1
counts[5,2] <- 1
counts[5,3] <- 1
counts[2,1] <- 1
counts[2,2] <- 1
counts[2,3] <- 1
counts[3,1] <- 1
counts[3,2] <- 1
counts[3,3] <- 1


BoykoParams <- 1
colnames(counts) <- c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40", "40-45", "45-50","50-55", "55-60", "60-65", "65-70", "70-75", "75-80", "80-85", "85-90", "90-95", "95-100", ">100")
rownames(counts) <- c("Real P(4Ns | DFE, D)", "Inferred P(4Ns| DFE, D)", "Inferred P(4Ns| DFE, D) with dadi. Full SFS with between 300-500 1% frequency variants", "Inferred P(4Ns| DFE, D) with dadi. Full SFS with 300-500 SNPs after taking the segregating sites with variants at all frequencies.", "Inferred P(4Ns| 1%, DFE, D)")
op <- par(cex = 0.4)
barplot(log10(counts) - log10(0.0000000001), main="A) Constant Size - Human DFE", ylab="Probability", xlab=expression(bolditalic(s[j]) * " (4Ns intervals)"), col=ViridisColors, legend = c(expression("Real" ~ 'P'[psi] * "(" * bolditalic(s[j]) * ")"), expression("Inferred" ~ 'P'[psi] * "(" * bolditalic(s[j]) * ") with fitdadi. Full SFS with between 300-500 1% frequency variants"), expression("Inferred" ~ 'P'[psi] * "(" * bolditalic(s[j]) * ") with fitdadi. Full SFS with 300-500 SNPs after taking the segregating sites with variants at all frequencies"), expression("Inferred" ~ 'P'[psi] * "(" * bolditalic(s[j]) * ") with our method using data from 300 1% frequency variants"), expression('P'[psi] * "(" * bolditalic(s[j]) * " | "* italic(f) * ", "* italic(D) *") with our method using data from 300 1% frequency variants")), beside=TRUE, cex.lab=2, cex.axis=2,cex.names=2, args.legend = list(x = "top",cex=1.5, bty = "n"), cex.main = 2, ylim = c(0,12), yaxt="n",cex = 0.5)
17:28
barplot(log10(counts) - log10(0.0000000001))





