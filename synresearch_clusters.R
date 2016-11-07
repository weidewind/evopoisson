library(ggplot2)
install.packages("ggdendro")
library(ggdendro)
protein = "h1"
path = file.path(dirname(getwd()), "input", "synresearch", "syn", paste(protein, "_reversals_list", sep = ""), fsep = .Platform$file.sep)
sdat = read.csv(path,comment.char = "#", stringsAsFactors=FALSE)
path = file.path(dirname(getwd()), "input", "synresearch", "nsyn", paste(protein, "_reversals_list", sep = ""), fsep = .Platform$file.sep)
ndat = read.csv(path,comment.char = "#", stringsAsFactors=FALSE)
dat <- rbind(sdat, ndat)

distances <- sapply(dat$list, function(e){
  vect <- sapply(dat$list, function(d){
          1-jaccard_index(e,d)
          })
} )
m <- as.matrix(distances)
rownames(m) <- dat$anc
colnames(m) <- dat$anc
dm<- as.dist(m)
clusters <- hclust(dm)
ggdendrogram(clusters, rotate = FALSE, size = 2)
#clusterCut <- cutree(clusters, 11)
clusterCut <- cutree(clusters, h=0.75) # for h1 all clusters at height 0.25 contain identical elements
length(clusterCut)
tab <- table(clusterCut, dat$anc)
identicals <- tab[rowSums(tab)>1,]
for (i in c(1:nrow(identicals))){
  ancs <- identicals[i,identicals[i,]>0]
  sapply(names(ancs), function(anc){
    print(dat[dat$ancestor == anc, ])
  })
  print ("-------")
  #print(identicals[i,identicals[i,]>0])
}


jaccard_index<- function(list1, list2){
  list1 = strsplit(list1, ";")[[1]]
  list2 = strsplit(list2, ";")[[1]]
  overlap <- length(intersect(list1, list2))
  union <-length(list1)+length(list2)-overlap
  jaccard <- overlap/union
  }

