Install and load mbDecoda
```
devtools::install_github('YuxuanZ0ng/mbDecoda_package')
library(mbDecoda)
help(mbDecoda)
```

An example
```
library(mbDecoda)
data = ZINB_sim(seed=1,K=50,n1=20,n2=20,p=0.1,bias ="small",zi=0.3,confounder =F)
group = data[["grp"]]
count = data[["count"]]
out = mbDecoda(count, x=group)
```
