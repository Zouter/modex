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
  parallel::mclapply(params, f2,..., mc.cores=mc.cores)
}

cutparams = expand.grid(k=k, method=method)
lapply(1:nrow(cutparams), function(i) as.list(cutparams[i,]))

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
