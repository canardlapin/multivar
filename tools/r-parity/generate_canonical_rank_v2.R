# Independent base-R score geometry and Gaussian reference calculations.
root <- commandArgs(trailingOnly=TRUE)[[1]]
out <- file.path(root,'modules/inference/fixtures/canonical-rank-v2')
options(digits=17)
orth <- function(x) qr.Q(qr(x))[,seq_len(ncol(x)),drop=FALSE]
roots <- function(x,y) svd(crossprod(orth(x),orth(y)),nu=0,nv=0)$d
wilks <- function(r) -sum(log1p(-r*r))
complete <- function(v) if(nrow(v)==ncol(v)) v else
  cbind(v,qr.Q(qr(v),complete=TRUE)[,(ncol(v)+1L):nrow(v),drop=FALSE])
permutations <- function(x) if(length(x)==1L) matrix(x,1L) else
  do.call(rbind,lapply(seq_along(x),function(i) cbind(x[[i]],permutations(x[-i]))))
input <- as.matrix(read.delim(file.path(out,'affine-input.tsv')))
x <- input[,1:2];y <- input[,3:5]
qx <- orth(x);qy <- orth(y);fit <- svd(crossprod(qx,qy))
u <- qx %*% complete(fit$u);v <- qy %*% complete(fit$v)
actions <- permutations(1:6)
statistics <- t(apply(actions,1L,function(action) vapply(1:2,function(k)
  wilks(roots(u[action,k:ncol(u),drop=FALSE],v[,k:ncol(v),drop=FALSE])),0.0)))
observed <- vapply(1:2,function(k) wilks(fit$d[k:2]),0.0)
counts <- colSums(sweep(statistics,2,observed,`-`) >= -1e-12)
stopifnot(all(counts==c(460,642)))
write.table(cbind(actions,statistics),file.path(out,'score-actions.tsv'),sep='\t',row.names=FALSE,quote=FALSE)

draws <- read.delim(file.path(out,'gaussian-input.tsv'),colClasses=c(rep('integer',5),'character','numeric'))
nulls <- matrix(0,2,39)
for(k in 0:1) for(b in 0:38) {
  a <- subset(draws,step==k & replicate==b & block==0)
  z <- subset(draws,step==k & replicate==b & block==1)
  nulls[k+1,b+1] <- wilks(roots(matrix(a$value,6,2-k,byrow=TRUE),matrix(z$value,6,3,byrow=TRUE)))
}
write.table(nulls,file.path(out,'gaussian-reference.tsv'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
fmt <- function(values) paste(sprintf('%.17g',values),collapse=',')
vector <- function(values) paste0('Vector(',fmt(values),')')
lines <- c('package multivar.inference','',
  '// Generated independently by tools/r-parity/generate_canonical_rank_v2.py and base R.',
  'private[inference] object CanonicalRankV2Fixtures:',
  paste0('  val input = gale.linalg.DMat.dense(6,5,',vector(as.vector(t(input))),')'),
  paste0('  val observed = ',vector(observed)),
  '  val actions = Vector(',
  vapply(seq_len(720),function(i) paste0('    (Vector(',paste(actions[i,]-1L,collapse=','),'),',vector(statistics[i,]),')',if(i<720) ',' else ')'),''),
  paste0('  val gaussian = Vector(',vector(nulls[1,]),',',vector(nulls[2,]),')'),
  paste0('  val hypothesisSeeds = Vector(',paste(paste0(unique(draws$hypothesis_seed),'L'),collapse=','),')'),
  paste0('  val firstNormalBlock = ',vector(subset(draws,step==0 & replicate==0 & block==0)$value)))
writeLines(lines,file.path(root,'modules/inference/shared/src/test/scala/multivar/inference/CanonicalRankV2Fixtures.scala'))
cat('PASS: 720 score actions, 78 independently generated Gaussian references; counts ',counts,'\n')
print(sessionInfo())
