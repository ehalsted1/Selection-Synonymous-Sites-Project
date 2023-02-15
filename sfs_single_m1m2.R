
##Non-coding mutations (m1)
mut_count_m1<-read.delim("/Users/arletlaar/Desktop/Data-UNAMlab/syn_project_trial4/Output/output_single_pop_m1_100_m1_count.sfs")
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
mut_count_m2<-read.delim("/Users/arletlaar/Desktop/Data-UNAMlab/syn_project_trial4/Output/output_single_pop_m2_100_m2_count.sfs")
mut_count_m2[,1]<-NULL #deleting the column names of the samples

#Sum of each column
for (i in 1:99){
  mut_count_m2$sites[i]<-sum(mut_count_m2[,i])
}

#Site frequency spectrum
barplot(mut_count_m2$sites[1:98], main = "Synonimous neutral mutations (m2)",
        xlab = "Sites", ylab = "Frequencies",
        col = colors()[10:60])


sites<-mut_count_m2$sites
write.table(sites, "sitesm2.txt")