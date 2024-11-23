# Set working directory (replace with your actual path)
# setwd("/path/to/your/working/directory")


# Simulation and data generation
set.seed(20)
u1 = rnorm(100)
set.seed(10)
u2 = rnorm(100)

# Pattern vectors
v1 = c(1, 0.86, 0.66, 0.9, rep(0, 8))
v2 = c(rep(0,4), 0.2, -0.55, -0.35, 0.17, rep(0, 4))

d1 = 10
d2 = 5

X1 = d1*u1%*%t(v1)
X2 = d2*u2%*%t(v2)

set.seed(30)
X = X1 + X2
write.csv(X1, file="X1.csv", row.names=FALSE)

# Noise generation
set.seed(30)
n_rows <- 100
n_cols <- 12
noise_genes <- matrix(0, nrow=n_rows, ncol=n_cols)
noise_values <- matrix(runif(n_rows*4, min=200, max=300), nrow=n_rows, ncol=4)
noise_genes[,9:12] <- noise_values
write.csv(noise_genes, file="noise_genes.csv", row.names=FALSE)

set.seed(30)
random_matrix <- matrix(rnorm(n_rows*n_cols), nrow=n_rows, ncol=n_cols)
X = X1 + X2 + noise_genes + 5*matrix(rnorm(100*12), ncol=12)

# Network edges
y_1 <- list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4), c(5,9), c(6,10))
y_2 <- list(c(5,6), c(5,7), c(5,8), c(6,7), c(6,8), c(7,8), c(1,11), c(3,12))
edges = c(y_1, y_2)

# Weight calculation
num_genes = 12
weights = matrix(1, nrow=num_genes, ncol=1)
connections = rep(0, num_genes)

for (edge in edges) {
    connections[edge[1]] = connections[edge[1]] + 1
    connections[edge[2]] = connections[edge[2]] + 1
}

max_conn = max(connections)
if (max_conn > 0) {
    for (i in 1:num_genes) {
        if (connections[i] > 0) {
            weights[i,1] = 1 + (connections[i]/max_conn)
        }
    }
}

write.table(weights, file="power.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# T-test analysis
X1 = t(X)
num_cols = ncol(X1)
half_cols = num_cols/2
t_test_results = numeric(nrow(X1))

for (i in 1:nrow(X1)) {
    t_test = t.test(X1[i,1:half_cols], X1[i,(half_cols+1):num_cols])
    t_test_results[i] = abs(t_test$statistic)
}

min_val = min(t_test_results)
max_val = max(t_test_results)
normalized_t_test_results = 1 + (t_test_results-min_val)/(max_val-min_val)
write.table(normalized_t_test_results, file="P_val.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# PCA methods comparison
out1 = svd(X, nu=2, nv=2)

source('functions/fun_SPCA.R')
out2 = SPCA(X, k=2, kv=c(4,4))

source('functions/fun_ESPCA.R')
out3 = ESPCA(X, k=2, edges, k.group=3)

source('functions/fun_DM_ESPCA.R')
out4 = DM(X, k=2, edges, k.group=3)

source('functions/fun_AWGE_ESPCA.R')
out5 = AWGE(X, k=2, edges, k.group=3)

# Results compilation
PC.dat = cbind(out1$v, out2$V, out3$V, out4$V, out5$V)
colnames(PC.dat) = c("PCA.PC1", "PCA.PC2", "SPCA.PC1", "SPCA.PC2", "ESPCA.PC1", "ESPCA.PC2",
                     "DM.PC1", "DM.PC2", "AWGE.PC1", "AWGE.PC2")
row.names(PC.dat) = paste0("Var", 1:num_genes)

write.csv(round(PC.dat, 3), "PC_data.csv", row.names=TRUE)
file.remove(c("power.txt", "P_val.txt", "noise_genes.csv", "X1.csv"))
