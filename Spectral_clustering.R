#  spectral clustering on data
rm(list=ls())
# for spectral clustering
library(kernlab) 
compounds=read.csv('Liion_comp_528.csv')
set.seed(7)
# the mxrd spectram for 528 compounds
data_to_pro=as.matrix(compounds[,5:4505])
rownames(data_to_pro)=as.character(compounds$formula_id)
# note the input of pca matrix rows are for different samples
# scale should be changed to false to keep original format in the pca calculation
pcaout=prcomp(data_to_pro,center=TRUE,scale=FALSE) 
pcaprop=cumsum(pcaout$sdev^2 / sum(pcaout$sdev^2))
plot(pcaprop)

# clean up variables
rm(first_sc,second_sc,third_sc,fourth_sc,fifth_sc,sixth_sc,seventh_sc)
rm(no1,no2,knext)
rm(ind,ind_2nd,ind_3rd,ind_4th,ind_5th,ind_6th,ind_7th,subdata,zoomindex)
rm(subdata,ggframe)

# build the data frame for plotting using PC
ggframe=data.frame(pcaout$x[,c(1:3)])
ggframe$name=rownames(ggframe)
ggframe$class=0
ggframe$col=0



first_sc=specc(data_to_pro,centers=2)
names(first_sc) # this gives the lists of spectral clustering compound names
for (ind in names(first_sc)){
  ggframe[ind,]$class=first_sc[ind]
}
#check number of elements in different first order clusters
no1=sum(first_sc==1) 
no2=sum(first_sc==2) 
knext=1+as.integer(no1<no2)
write.csv(ggframe[,c(4,5)],"color_structure/level_1_structure_color.csv")


#Now divide the cluster deeper
ind_2nd=names(first_sc[first_sc==knext])
second_sc=specc(data_to_pro[ind_2nd,],centers=2)
for (ind in names(second_sc)){
  ggframe[ind,]$class=second_sc[ind]+2
}
no1=sum(second_sc==1) 
no2=sum(second_sc==2) 
knext=1+as.integer(no1<no2)
write.csv(ggframe[,c(4,5)],"color_structure/level_2_structure_color.csv")

# divide again
ind_3rd=names(second_sc[second_sc==knext])
third_sc=specc(data_to_pro[ind_3rd,],centers=2)
for (ind in names(third_sc)){
  ggframe[ind,]$class=third_sc[ind]+4
}
no1=sum(third_sc==1) 
no2=sum(third_sc==2) 
knext=1+as.integer(no1<no2)
write.csv(ggframe[,c(4,5)],"color_structure/level_3_structure_color.csv")


#### save the data at level 3.
ggframe$col[ggframe$class==5]='cornflowerblue'
save(list = ls(all = TRUE), file= "all.RData")



# divide further
ind_4th=names(third_sc[third_sc==knext])
fourth_sc=specc(data_to_pro[ind_4th,],centers=2)
for (ind in names(fourth_sc)){
  ggframe[ind,]$class=fourth_sc[ind]+6
}
no1=sum(fourth_sc==1) 
no2=sum(fourth_sc==2) 
knext=1+as.integer(no1<no2)
subdata=subset(ggframe,col=="cornflowerblue")
write.csv(subdata[,c(4,5)],"color_structure/level_4_structure_color.csv")


ind_5th=names(fourth_sc[fourth_sc==knext])
fifth_sc=specc(data_to_pro[ind_5th,],centers=2)
for (ind in names(fifth_sc)){
  ggframe[ind,]$class=fifth_sc[ind]+8
}
no1=sum(fifth_sc==1) 
no2=sum(fifth_sc==2) 
knext=1+as.integer(no1<no2)
subdata=subset(ggframe,col=="cornflowerblue")
write.csv(subdata[,c(4,5)],"color_structure/level_5_structure_color.csv")


# divide further
ind_6th=names(fifth_sc[fifth_sc==knext])
sixth_sc=specc(data_to_pro[ind_6th,],centers=2)
for (ind in names(sixth_sc)){
  ggframe[ind,]$class=sixth_sc[ind]+10
}
no1=sum(sixth_sc==1) 
no2=sum(sixth_sc==2) 
knext=1+as.integer(no1<no2)
subdata=subset(ggframe,col=="cornflowerblue")
write.csv(subdata[,c(4,5)],"color_structure/level_6_structure_color.csv")


# divide further
ind_7th=names(sixth_sc[sixth_sc==knext])
seventh_sc=specc(data_to_pro[ind_7th,],centers=2)
for (ind in names(seventh_sc)){
  ggframe[ind,]$class=seventh_sc[ind]+12
}
no1=sum(seventh_sc==1) 
no2=sum(seventh_sc==2) 
knext=1+as.integer(no1<no2)
subdata=subset(ggframe,col=="cornflowerblue")
write.csv(subdata[,c(4,5)],"color_structure/level_7_structure_color.csv")


#### save the clustering at 7th level
save(subdata,file="subdat.Rda")