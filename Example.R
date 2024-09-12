# Set working directory (replace with your actual path)
# setwd("/path/to/your/working/directory")

# Generate gene and pathway data
set.seed(20)
u1 = rnorm(100)
set.seed(10)
u2 = rnorm(100)

v1 = c(1, 0.86, 0.66, 0.9, rep(0, 8))
v2 = c(rep(0,4), 0.2, -0.55, -0.35, 0.17, rep(0, 4))

X1 = 10 * u1 %*% t(v1)
X2 = 5 * u2 %*% t(v2)

# Generate noise and combine data
set.seed(30)
n_rows <- 100
n_cols <- 12
noise_genes <- cbind(matrix(0, nrow = n_rows, ncol = 8), 
                     matrix(runif(n_rows * 4, min = 200, max = 300), nrow = n_rows, ncol = 4))
X = X1 + X2 + noise_genes + 5 * matrix(rnorm(n_rows * n_cols), ncol = n_cols)

# Define and process edges
edges = c(
  list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4), c(5,9), c(6,10)),
  list(c(5,6), c(5,7), c(5,8), c(6,7), c(6,8), c(7,8), c(1,11), c(3,12))
)

num_genes = n_cols
weights = matrix(1, nrow = num_genes, ncol = 1)
connections = sapply(1:num_genes, function(i) sum(sapply(edges, function(e) sum(e == i))))

if (max(connections) > 0) {
  weights[connections > 0, 1] = 1 + connections[connections > 0] / max(connections)
}

write.table(weights, file = "power.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Perform t-test and normalize results
X1 = t(X)
t_test_results = sapply(1:nrow(X1), function(i) {
  abs(t.test(X1[i, 1:(ncol(X1)/2)], X1[i, (ncol(X1)/2 + 1):ncol(X1)])$statistic)
})

normalized_t_test_results = 1 + (t_test_results - min(t_test_results)) / 
  (max(t_test_results) - min(t_test_results))

write.table(normalized_t_test_results, file = "P_val.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Perform various PCA methods
out1 = svd(X, nu = 2, nv = 2)

source('functions/fun_SPCA.R')
out2 = SPCA(X, k = 2, kv = c(4,4))

source('functions/fun_ESPCA.R')
out3 = ESPCA(X, k = 2, edges, k.group = 3)

source('functions/fun_DM_ESPCA.R')
out4 = DM(X, k = 2, edges, k.group = 3)

source('functions/fun_AWGE_ESPCA.R')
out5 = AWGE(X, k = 2, edges, k.group = 3)

# Combine and save results
PC.dat = cbind(out1$v, out2$V, out3$V, out4$V, out5$V)
colnames(PC.dat) = c("PCA.PC1", "PCA.PC2", "SPCA.PC1", "SPCA.PC2", "ESPCA.PC1", "ESPCA.PC2",
                     "DM.PC1", "DM.PC2", "AWGE.PC1", "AWGE.PC2")
row.names(PC.dat) = paste0("Var", 1:num_genes)

write.csv(round(PC.dat, 3), "PC_data.csv", row.names = TRUE)

file.remove(c("power.txt", "P_val.txt"))

