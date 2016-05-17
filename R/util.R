mcdapply <- function(df, f, ..., mc.cores=getOption("mc.cores", default=1)) {
  f2 = function(...) {
    tryCatch(f(...),
             error=function(e) {
               warning(e)
               NULL
             }
    )
  }

  params = lapply(1:nrow(df), function(i) setNames(as.list(df[i,]), colnames(df)))
  parallel::mclapply(params, f2, ..., mc.cores=mc.cores, mc.silent=F)
}

matchparams <- function(df, params) {
  matchparams = intersect(colnames(df), names(params))
  match(interaction(params[matchparams]), interaction(df[,matchparams,drop=F]))[[1]]
}

# expand grid directly without dumb factors and attribute names, why isn't this the default??
expand_grid <- function(...) {
  df = expand.grid(..., stringsAsFactors = F)
  attr(df, "out.attrs") <- NULL
  df
}

labels2modules = function(labels, genes=NULL, ignore_label=NULL) {
  if(!is.null(genes)) names(labels) = genes
  if(!is.null(ignore_label)) labels = labels[labels != ignore_label]
  labels = factor(labels)
  lapply(levels(labels), function(label) {
    names(labels)[labels == label]
  })
}

aacast = function(df, formula, value.var) {
  formula = gsetid~moduleid
  df = reshape2::dcast(scores, formula, value.var = value.var)
  rownames(df) = df[,1]
  as.matrix(df[,2:ncol(df)])
}
