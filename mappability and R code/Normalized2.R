Normalized2<-function(dat,nci,nri,classid,MappSS, Biotype,Genename,file,classname,output){
	if(nci==0&nri==0){
data<-dat		
	}else if(nci==0&nri>=1){
data<-dat[-c(1:nri),]		
	}else if(nci>=1&nri==0){
data<-dat[,-c(1:nci)]		
	}else if(nci>=1&nri>=1){
data<-dat[-c(1:nri),-c(1:nci)]
}
#data<-as.data.frame(data)
data<-apply(data,2,as.numeric)
allbiotypes = unique(Biotype)
#print(allbiotypes)
allbiotypeslength <- rep(NA,nrow(allbiotypes))
allbiotypescount_data <- matrix(0,nrow(allbiotypes),ncol(data))
allbiotypespercent_data<- matrix(0,nrow(allbiotypes),ncol(data))
#print(allbiotypescount_data)
for (i in 1:nrow(allbiotypes)){

temp <- is.element(Biotype[,1],allbiotypes[,1])

allbiotypeslength[i] <- length(which(temp==T))

if (allbiotypeslength[i]>1){
ind = is.element(Biotype[,1],allbiotypes[i,1])	
subdata<-data[ind,]
if(!is.null(dim(subdata))){
allbiotypescount_data[i,] = apply(subdata,2,sum,na.rm=TRUE) 
allbiotypespercent_data[i,] = allbiotypescount_data[i,]/apply(data,2,sum,na.rm=TRUE)*100
}
}
}
#print(dim(allbiotypespercent_data))
#print(allbiotypespercent_data[c(2,3),])
#print(dim(allbiotypespercent_data))
row.names(allbiotypespercent_data)<-allbiotypes[,1]
colnames(allbiotypespercent_data)<-dat[nri,-c(1:nci)]

#print(allbiotypespercent_data[c(2,3),])
min<-apply(allbiotypespercent_data,1,min)
max<-apply(allbiotypespercent_data,1,max)
average<-apply(allbiotypespercent_data,1,mean)
allbiotypespercent_data<-cbind(min,max,average,allbiotypespercent_data)
#write.csv(allbiotypespercent_data, file =paste("Biotype_",file,".csv"))

proteincodingindx = is.element(Biotype[,1],allbiotypes[3,1])
lincrnaindx = is.element(Biotype,allbiotypes[2,1])
biotypeindx = which(proteincodingindx | lincrnaindx)
length(biotypeindx)
#proteincoding_data<-IGF2BP1ripTPM2[proteincodingindx,]
#install.packages("stringr")              # Install stringr package

library("stringr")

indMT<-str_detect(Genename[,1],'^MT-')
indH1<-str_detect(Genename[,1], '^H1')
indH2<-str_detect(Genename[,1], '^H2')
indH3<-str_detect(Genename[,1], '^H3')
indH4<-str_detect(Genename[,1], '^H4')
indRPL<-str_detect(Genename[,1], '^RPL')
indRPS<-str_detect(Genename[,1], '^RPS')
Geneid<-dat[-c(1:nri),1]
#additgeneid<-Gencode_33_Selected_Geneid[ind,1]
MTgeneid<-Geneid[indMT]
H1geneid<-Geneid[indH1]
H2geneid<-Geneid[indH2]
H3geneid<-Geneid[indH3]
H4geneid<-Geneid[indH4]
RPLgeneid<-Geneid[indRPL]
RPSgeneid<-Geneid[indRPS]

#additionalgenes<-c(MTgeneid,H1geneid,H2geneid,H3geneid,H4geneid,RPLgeneid,RPSgeneid)
#additional_Geneid_GMask<-unique(additionalgenes)
#length(additional_Geneid_GMask)

additionalgenes<-indMT|indH1|indH2|indH3|indH4|indRPL|indRPS
length(additionalgenes)

nonadditionalgenes<-seq(length(Genename[,1]))
nonadditionalgenes<-nonadditionalgenes[!additionalgenes]

mappableindx<-which(MappSS[,1]>50)
length(mappableindx)


finalIndexGeneric = intersect(biotypeindx,intersect(nonadditionalgenes,mappableindx));

sm<-apply(data,1,sum,na.rm=TRUE)
datindex <- which(sm>5)

smrip<-apply(data,1,sum,na.rm=TRUE)
datindexrip <- which(smrip>6)

countindx <- unique(c(datindex,datindexrip))
length(countindx)

finalIndex<-intersect(finalIndexGeneric,countindx)
data1 <- data[finalIndex,]
geneid <- Geneid[finalIndex]
genename <- Genename[finalIndex,1]
mappss<- MappSS[finalIndex,1]
#mappus<- MappUS[finalIndex,1]
biotype<-Biotype[finalIndex,1]

data2 <- data1/mappss*1000
#data3<-apply(data2,2,as.numeric)
data2[is.nan(data2)]<-0
data2[is.infinite(data2)]<-0
data3<-data2
small<-apply(data3,2,sum,na.rm=TRUE)
for(i in 1:length(small)){
data3[,i] <- data2[,i]/small[i]*1000000
}
data4<-cbind(as.data.frame(geneid),as.data.frame(genename),data3)
colnames(data4)[1:2]<-c("x1","x2")
TMP<-rbind(dat[c(1:nri),],data4)
#write.csv(TMP, file=paste("TPM_",file,".csv"))

treatment<-factor(dat[classid,-c(1:nci)])

#colData<-matrix(NA,ncol(data3),1)
#print(treatment[1:20])
#row.names(colData)<- sample[1,]


#colData<-t(treatment)

#print(dim(colData))
#print(colData[1:5,])

#colnames(data) <- sample[1,]

library(DESeq2)
library(IHW)

GMaskFactor <- DESeqDataSetFromMatrix(data1, colData=as.data.frame(treatment),design= ~ treatment)
GMaskFactor <- DESeq(GMaskFactor)

GMaskFactor_vsd<- varianceStabilizingTransformation(GMaskFactor,blind=FALSE)
vsd<-assay(GMaskFactor_vsd)
#print(names(vsd))
data5<-cbind(as.data.frame(geneid),as.data.frame(genename),vsd)
colnames(data5)[1:2]<-c("x1","x2")
VSD<-rbind(dat[c(1:nri),],data5)
#write.csv(VSD,file=paste("VSD_",file,classname,".csv"))

if(output=="TMP"|output=="tmp"){
data6<-TMP	
}else if(output=="VSD"|output=="vsd"){
data6<-VSD	
}

return(data6)
}