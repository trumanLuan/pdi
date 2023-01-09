#####TRIM#####

trimEdge <- function(x, threshold=1)
{
  res <- list()
  for (i in seq(length(x)))
  {

if(all(x2ic(x[[i]]/sum(x[[i]][,1]))<threshold)){res[[i]]=NA}else{

    res[[i]]=x[[i]]
    for (k in 1:2)
    {
      j=1
      ic  <- x2ic(res[[i]]/sum(res[[i]][,1]))

      while (ic[j] < threshold)
      {
        res[[i]] <- res[[i]][,-1]
        j=j+1
      }
      res[[i]] <- res[[i]][, ncol(res[[i]]):1]
    } # end of k
   } # end of if
  } # end of i
  names(res) <- names(x)
  return (res)
}

x2ic <- function(x)
{
  npos <- ncol(x)
  ic <- numeric(length=npos)
  for (i in 1:npos)
  {
    ic[i] <- 2 + sum(sapply(x[, i], function(x) { 
      if (x > 0) { x*log2(x) } else { 0 }
      }))
    }    
    ic
  }

trimNA<-function(x){
lapply(x,function(x)x[,!apply(is.na(x),2,all)])
}

