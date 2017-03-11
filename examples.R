library(data.table)
library(IceT)
ercc = fread("/Users/rongma/Box Sync/Rong-Mingyao/ercc.txt")
ercccount = as.matrix(as.data.frame(ercc)[,c(-1,-2)])

mrna = fread("~/Box Sync/Rong-Mingyao/bio_mrna_read_counts.txt")
mrna = as.matrix(as.data.frame(mrna)[,-1])

#estimate TASC model parameters
alphabeta = ab.est(ercc$concentration, ercccount)

#recover expression data using TASC model
r.data = recover(mrna, alphabeta)

#cell type clustering via IceT
icet0 = icet(r.data)

#plot cell clusters
icet.plot(icet0)


#reset some parameters (e.g. eps) in icet.plot
icet.plot(icet0, reset = T, recover.data = r.data, eps=1.3)

cluster5 = which(icet0$cluster==5)
r2.data = list(recover = r.data$recover[,icet0$cluster == 5], recover.log=r.data$recover.log[,icet0$cluster == 5],
               reduced=r.data$reduced[,icet0$cluster == 5], reduced.log=r.data$reduced.log[,icet0$cluster == 5])
icet1 = icet(r2.data, fdr = 0.05)
icet.plot(icet1)
icet.plot(icet1, reset = T, recover.data = r2.data, perplexity = 30, eps = 2.5)


cluster.z = c()
cluster.z[which(cluster.z == 0)] = NA
cluster.z[icet0$cluster != 5] = icet0$cluster[icet0$cluster != 5]
icet1$cluster[which(icet1$cluster == 0)] = NA
cluster.z[icet0$cluster == 5] = icet1$cluster+4
table(cluster.z)

names(cluster.z) = colnames(mrna)[-1]

write.table(cluster.z, "/Users/rongma/Box Sync/Rong-Mingyao/cluster_index.txt")

data = matrix(nrow = 19972, ncol = 3005)
data = r.data$recover
rownames(data) = mrna$cell_id
colnames(data) = colnames(mrna)[-1]

write.table(data, "/Users/rongma/Box Sync/Rong-Mingyao/recovered_count.txt")




