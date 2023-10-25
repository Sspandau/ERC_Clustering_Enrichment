library("dplyr")
library("ape")
library("ctc")

args = commandArgs(trailingOnly=TRUE)

#load in distance matrix
distance_matrix <-read.csv(args[1])

#load clusters csv
clusters<-read.csv(args[2])

#clusters > 200
clu <-subset(clusters, count > args[5]) # which clusters to split up

for(i in clu$cluster){
  print(i)
  clusterdf <- subset(clusters, cluster == i)
  cluster_prot <- clusterdf$proteins
  matrix_clust <- distance_matrix[(rownames(distance_matrix) %in% cluster_prot),(rownames(distance_matrix) %in% cluster_prot)]
  rho.dist <- as.dist(nrow(distance_matrix)-matrix_clust)
  rho.tree <- hclust(rho.dist, method="complete")
  rho.dend <-as.dendrogram(rho.tree)
  subsize <- subset(clu, cluster == i)
  size <- subsize$count
  size2 <- size%/%100
  print(size2)
  if(i == 1){
    voorhees_subcluster1<-cutree(rho.tree, k = size2) 
    voorhees_subcluster1 <- cbind(matrix_clust, voorhees_subcluster1)
    rownames(voorhees_subcluster1)<-cluster_prot
    voorhees_subcluster1<-data.frame(proteins = c(rownames(voorhees_subcluster1)), subcluster = c(voorhees_subcluster1$voorhees_subcluster1))
    voorhees_subcluster1$cluster<-c(i)
    b1 <- voorhees_subcluster1 %>% group_by(subcluster) %>% summarize(count=n())
    b1$cluster<-c(i)
  }
  else{
    voorhees_subcluster<-cutree(rho.tree, k = size2) 
    voorhees_subcluster <- cbind(matrix_clust, voorhees_subcluster)
    rownames(voorhees_subcluster)<-cluster_prot
    voorhees_subcluster<-data.frame(proteins = c(rownames(voorhees_subcluster)), subcluster = c(voorhees_subcluster$voorhees_subcluster))
    voorhees_subcluster$cluster<-c(i)
    voorhees_subcluster1<-rbind(voorhees_subcluster1, voorhees_subcluster)
    b <- voorhees_subcluster %>% group_by(subcluster) %>% summarize(count=n())
    b$cluster<-c(i)
    b1 <- rbind(b1, b)
  }
}

write.csv(voorhees_subcluster1, args[3])
write.csv(b1, args[4])