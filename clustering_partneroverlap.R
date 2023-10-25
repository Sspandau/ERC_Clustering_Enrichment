library("dplyr")
library("ape")
library("ctc")

args = commandArgs(trailingOnly=TRUE)

#load in distance matrix
distance_matrix <-read.csv(args[1])
distance_matrix<-distance_matrix[, -c(1)]

#save vector of protein names and make them the row and col names
protein_names <- colnames(distance_matrix)
rownames(distance_matrix) <- protein_names

#clustering for partner overlap
rho.dist <- as.distance(nrow(distance_matrix)-distance_matrix)
rho.tree <- hclust(rho.dist, method="complete")
rho.dend <-as.dendrogram(rho.tree)

#save cluster tree to newick file
newick<-hc2Newick(rho.tree)
write(hc2Newick(rho.tree),file=args[2])

#plot tree
plot(rho.tree, cex = 0.01, main = "Coevoling Protein Partner Overlap Dendrogram")

# cut tree at predetermined threshold (either number of clusters or max distance)
if(args[5]=="k"){
  clusters<-cutree(rho.tree, k = args[6]) #k is number of clusters, h is distance 
}else{
  clusters<-cutree(rho.tree, h = args[6]) #k is number of clusters, h is distance 
}

#add clusters to df
clusters_df <- cbind(distance_matrix, clusters)
rownames(clusters_df)<-protein_names
clusters_df<-data.frame(proteins = c(rownames(clusters_df)), cluster = c(clusters_df$clusters))

#count cluster sizes
clusters_count <- clusters_df %>% group_by(cluster) %>% summarize(count=n())

#save clusters
write.csv(clusters_df, args[3])
write.csv(clusters_count, args[4])
