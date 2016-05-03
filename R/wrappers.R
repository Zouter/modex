agglom <- function(E, k=c(5, 10, 15, 20), method=c("ward.D2", "single", "complete", "average")) {
  # although this can look quite cumbersome, this kind of setup is easy to expand to multiple parameters and multiple commands which have to be run progressively
  treeparams = expand_grid(method=method)

  gene_cors = as.dist(1-cor(t(E)))

  trees = mcdapply(treeparams, function(params) {
    stats::hclust(gene_cors, method = params$method)
  })

  cutparams = expand_grid(k=k, method=method)
  modules = mcdapply(cutparams, function(params) {
    tree = trees[[matchparams(treeparams, params)]]
    stats::cutree(tree, params$k)
  })
  list(params=cutparams, modules=modules)
}
