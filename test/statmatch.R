# --> c.neti (only in data1) and c.neti.bis (only in data2) coded voluntarily in 2 distinct encodings but represent the same information
# --> c.neti is a factor variable corresponding to the person's net income initially categorized in 7 classes of thousand of Euros

 
library(StatMatch)

data(samp.A)
samp.A = samp.A[,c(1:11,13,12)]
c.neti = as.numeric(samp.A$c.neti)

samp.A$c.neti.bis = as.factor(ifelse(c.neti %in% c(1,2),1,
                              ifelse(c.neti %in% c(3,4),2,
                              ifelse(c.neti %in% c(5,6),3,4))))

write.csv2(samp.A, "statmatch.csv")

## Add random NA in covariates (10% by covariates):
#add_NA = function(DB,tx){
#  DB_NA = DB
#  for (j in 1:ncol(DB)){
#     NA_indic = sample(1:nrow(DB),round(nrow(DB)*tx/100),replace=FALSE)
#     DB_NA[NA_indic,j] = rep(NA,length(NA_indic))
#  }
#return(DB_NA)
#}
#
#set.seed(4036)
#df = add_NA(data,10)
