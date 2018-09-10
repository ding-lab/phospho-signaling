# Yige Wu @ WashU March 2018
# analyze cohort level mRNA/protein/phosphoprotein expression data and convert to different matrices in a sample-gene format



# functions ---------------------------------------------------------------
unfactorize = function(df){
  for(i in which(sapply(df, class) == "factor")) df[[i]] = as.numeric(as.character(df[[i]]))
  return(df)
}

ecdf_fun = function(x,perc) ecdf(x)(perc)

expression_effect = function(m){ 
  cat("##### EXPRESSION ANALYSIS #####\n")
  minNum = 5
  m = as.matrix(m)
  num = nrow(m)
  m2 = as.matrix(m[rowSums(!is.na(m)) >= minNum, ])
  num_NA= nrow(m2)
  cat(paste("Original number of markers:", num, "; NA filtered:", num_NA, "\n", sep=" "))
  
  # initiate tables
  outlier = matrix(data = NA,nrow=dim(m2)[1],ncol=dim(m2)[2])
  row.names(outlier) = row.names(m2)
  colnames(outlier) = colnames(m2)
  exp_score = outlier
  exp_quantile = outlier
  
  # gene-wise expression score and quantile score
  for (i in 1:nrow(m2)){
    #IQR = quantile(m2[i,], probs=0.75, na.rm=T) - quantile(m2[i,], probs=0.25, na.rm=T) 
    exp_score[i,] = m2[i,]#(m2[i,] - quantile(m2[i,], probs=0.50, na.rm=T))/IQR
    exp_quantile[i,] = ecdf_fun(m2[i,],m2[i,])
  }
  
  return(list("exp_score"=exp_score, "exp_quantile"=exp_quantile))
}


