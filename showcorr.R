library( rhdf5 )

cc <- h5read( "~/w/metdense_autocorr/corr_counts.h5", "corr_counts")
cc <- aperm( cc, 4:1 )

ccs <- apply( cc, c(1,3,4), sum )

ccs[4,,]

corr_coef_from_2x2_tbl <- function(tbl) {
  storage.mode(tbl) <- "numeric"
  ( tbl[2,2]*tbl[1,1] - tbl[2,1]*tbl[1,2] ) /
    sqrt( prod(rowSums(tbl)) * prod(colSums(tbl)) )
}

plot(
  2^(1:17),
  apply( ccs, 1, corr_coef_from_2x2_tbl ), log="x" )

corr_coef_from_2x2_tbl( ccs[1,,] )

x <- as.numeric( runif(1000)<.3 )
y <- as.numeric( runif(1000) < .3+x/3 )
cor(x,y)

tbl <- table(x,y)

( tbl[2,2]*tbl[1,1] - tbl[2,1]*tbl[1,2] ) /
sqrt( prod(rowSums(tbl)) * prod(colSums(tbl)) )

