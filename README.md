# Perform differential expression analysis on fibrotic and non-fibrotic patients under 4 different treatments on HPC 

<br> 

## Author 

Sam (Cheng-Hsiang) Lu  

Email: Cheng-Hsiang.Lu@cshs.org 

<br> 

## Mentor 

David Casero 

Email: David.Casero@cshs.org 

<br> 

## Background 

### Inflammatory bowel diseases(IBD) 

Inflammatory bowel diseases (IBD) are a group of chronic conditions that cause inflammation and damage to the digestive tract. The two main types of IBD are Crohn's Disease and Ulcerative Colitis. 


Crohn's disease can affect any part of the digestive tract, like small or large intestine, and can cause symptoms such as abdominal pain, diarrhea, weight loss, and fatigue. It can also cause complications such as fistulas (abnormal connections between different parts of the intestine) and strictures (narrowing of the intestine). 


Ulcerative colitis, on the other hand, affects only the colon and rectum and causes symptoms such as bloody diarrhea, abdominal pain, and a frequent need to pass stools. It can also lead to complications such as inflammation of the skin, eyes, and joints. 


Both Crohn's disease and ulcerative colitis are chronic conditions, meaning they can last for a lifetime and require ongoing treatment to manage symptoms and prevent complications. 

### Induced Pluripotent Stem Cells(iPSCs) 

Induced pluripotent stem cells (iPSCs) are a type of stem cell that are generated in the laboratory by reprogramming adult cells, such as skin or blood cells, to a pluripotent state. A pluripotent state means that the cells have the potential to develop into any type of cell in the body, just like embryonic stem cells. 


iPSCs offer several advantages as they can be generated from the patient's own cells, avoiding issues with immune rejection, ethical concerns and the need for embryos. 


iPSCs can be used to study the underlying causes of diseases, test new drugs and therapies, and potentially generate replacement tissues or organs for transplantation.

<br> 

## Aim 

In this project, we will be analyzing RNA-seq data from 19 samples, comprising of 10 samples with fibrotic complications and 9 non-fibrotic samples. Each sample has undergone two runs and 4 different treatments(untreated, TGF-B, TNF-A, and TGF-B+TNF-A), resulting in a total of 151 samples(1 library failed). We used induced pluripotent stem cells (iPSC) to differentiate into myofibroblasts and stimulated the system with different signals to observe its development. The objective is to investigate the effect of four different treatments: untreated, TGF-B, TNF-A, and TGF-B+TNF-A on the development of the system. In the end, we will perform differential expression analysis to identify the genes that are differentially expressed in fibrotic and non-fibrotic samples under 4 treatments. 

<br> 

## Pipelines 

### In HPC

#### Convert 151 Fastq to Fasta files 

 

First, I put all fasq.gz files in one folder and list all fastq files’ name in fastqfiles.txt 

```
ls *q.gz > fastqfiles.txt
``` 

Cut redundant suffix “\_R1_trimmed” and list all fastq files’ name in libraryname.txt and preffix.txt 

```
ls *q.gz | cut -f 1 -d '.' | sed 's/_R1_trimmed//g'  >libraryname.txt 

ls *q.gz | cut -f 1 -d '.' | sed 's/_R1_trimmed//g' > preffix.txt
``` 

Form a table with 3 columns: fastqfiles.txt libraryname.txt preffix.txt 

```
paste fastqfiles.txt libraryname.txt preffix.txt > tofastatable.txt
``` 

Create small-sized fasta-formatted files. To submit this job to the cluster on HPC, you need to read the file, library, and prefix. Once you have done that, run the script "generatefastaFromFastaqz" which will combine the script "DCfastaqTofastaLibraryId.pl". This results in small-sized fasta-formatted files contain only one header and one sequence per read. You can find all scripts in the "scripts" folder. 

``` 
cat tofastatable.txt | awk {print}' | while read file library preffix ; do qsub -cwd -o $PWD -e $PWD -l h_data=2048M,h_rt=8:00:00 $HOME/scripts/generatefastaFromFastaqz $file $library $preffix  

done 
``` 
This is what each fasta-formatted file would look like

![](/Pics/fasta_formatted.JPG)

<br>

#### Generate auxiliary files and directories for each sample 

Put a list of names of all fasta files in the directory and save them in a text file named "fastafiles.txt" 

```
ls *fasta.gz > fastafiles.txt
``` 

Cut the redundant suffix ".fasta.gz" from the names of all fasta files and generate a new list of file names with the suffix removed in a text file named “targetdirectories.GTF.txt” 

```
cat fastafiles.txt |sed 's/.fasta.gz//g' > targetdirectories.GTF.txt
``` 

Create a separate directory for each sample listed in "fastafiles.txt" 

```
cat fastafiles.txt |while read line ; do mkdir ${line/.fasta\.gz/GTFpass1/} ; done
``` 

<br> 

#### Form the submission script called “sendmyof” 

Add a shebang line at the beginning of your script file named "sendmyof" to indicate the interpreter that should be used to execute the script 

```
echo '#!/bin/bash/' > sendmyof
``` 

The command below runs the "generatesendscriptSingleGTFParam" script with several input parameters to map the RNA-seq data with STAR. The input parameters include the list of target directories containing the input data ("targetdirectories.GTF.txt"), the subdirectory name ("GTFpass1"), a parameter file containing settings for STAR alignment ("Parameters.txt"), a prefix for output files ("myof"), the path to the STAR index directory ("/home/luc/RNASEQ_MASTER/Hsapiens/GRC38/INDEXES/GRCh38.primary.33.basicselected.STAR2.7.3a/"), the path to the input data directory ("/home/luc/iPSC/MYOFIBROBLAST/"), the amount of free memory to use ("mem_free=32G"), and the number of threads to use ("8"). 
 
``` 
./generatesendscriptSingleGTFParam targetdirectories.GTF.txt GTFpass1 Parameters.txt myof /home/luc/RNASEQ_MASTER/Hsapiens/GRC38/INDEXES/GRCh38.primary.33.basicselected.STAR2.7.3a/ /home/luc/iPSC/MYOFIBROBLAST/ mem_free=32G 8 >> sendmyof 
``` 

Change sendmyof into executable mode and run sendmyof 

``` 
chmod a+x sendmyof 
. sendmyof 
``` 

It will take less than one day to run through 151 samples and generate each sample a folder which contain every output from STAR

<br> 

#### Create a table summarizing the mapping statistics for each sample 

Change directory into one sample file which ends with "GTFpass1". Extract the first column from the mapping statistics file and store it in “temp2.txt” 

``` 
grep "|" 008iP22TGFbM_S71GTFpass1/008iP22TGFbM_S71GTFpass1Log.final.out | cut -f 1 -d "|" | sed 's/^  *//g' | awk 'NR>3 {print}' > temp2.txt 
```

The first column from the mapping statistics file

![](/Pics/mapping_statistics.jpg)

<br>

Create an empty temporary file for storing intermediate results 

```
rm tempprev.txt 
touch tempprev.txt
```

Extract the total mapped reads from each subsequent mapping statistics file and combine with previous results 

```
ls *pass1/*final.out | while read line ; do  
grep "|" $line | cut -f 2  > temp.txt 
paste tempprev.txt temp.txt > tempnew.txt 
mv tempnew.txt tempprev.txt 
done 
```

Remove the first column and write the final results to a file called “mappingstatsFirstpass.txt” 

```
cut -f 2- tempprev.txt | awk 'NR>3 {print}'  > tempnew.txt 
mv tempnew.txt tempprev.txt 
paste temp2.txt tempprev.txt > mappingstatsFirstpass.txt 
``` 

The mappingstatsFirstpass.txt would look like this

![](/Pics/mappingstatsFirstpass.jpg)

<br> 

#### Compile Counts 

Generate a directory called "COUNTS" and copy all gene count files to this folder. 

``` 
mkdir COUNTS 
cp *pass1/*PerGene* COUNTS/ 
``` 

Clean filenames 

``` 
ls *tab|while read line ; do mv $line ${line/GTFpass1ReadsPerGene.out/} ; done 
```

For each count file, extract and create the five count tables

```
ls *.tab | while read line ; do 
echo $line 
cat $line | awk 'NR==3{print}' | cut -f 2- > ${line/tab/nofeature\.tab} 
cat $line | awk 'NR==4{print}' | cut -f 2- > ${line/tab/ambiguous\.tab} 
cat $line | awk 'NR>4{print}' | cut -f 2 > ${line/tab/nostrand\.tab} 
cat $line | awk 'NR>4{print}' | cut -f 3 > ${line/tab/sense\.tab} 
cat $line | awk 'NR>4{print}' | cut -f 4 > ${line/tab/antisense\.tab} 
done 
```

Make a Geneid list from one of the count tables as "countsannot_GRCh38.primary.Selected.Geneid.txt"

```
ls 008iP22TGFbM_S71.tab | head -1 | while read line; do  
cut -f 1 $line | awk 'NR>4{print}' > countsannot_GRCh38.primary.Selected.Geneid.txt 
done 
```

Create a file listing the names of all samples as "RBarretTNFATGFBsamples.txt"

```
ls *.sense.tab | sed 's/.sense.tab//g' | tr -s " " "\n" | sed 's/_1//g' > RBarretTNFATGFBsamples.txt 
```

Make count tables for sense, anti-sense, nostrand, ambiguous, and nofeature reads

```  
# Combine all sense counts into RBarretTNFATGFB_sense.ALL.cnt 
paste *.sense.tab > RBarretTNFATGFB_sense.ALL.cnt   

# Combine all antisense counts into RBarretTNFATGFB_antisense.ALL.cnt 
paste *.antisense.tab > RBarretTNFATGFB_antisense.ALL.cnt   

# Combine all nostrand counts into RBarretTNFATGFB_nostrand.ALL.cnt 
paste *.nostrand.tab > RBarretTNFATGFB_nostrand.ALL.cnt 

# Combine all ambiguous counts into RBarretTNFATGFB_ambiguous.cnt 
cat *ambiguous.tab > RBarretTNFATGFB_ambiguous.cnt   

# Combine all nofeature counts into RBarretTNFATGFB_nofeature.cnt 
cat *nofeature.tab > RBarretTNFATGFB_nofeature.cnt 
``` 

I am going to use "RBarretTNFATGFB_antisense.ALL.cnt" file for the further analysis

<br> 

### In MATLAB

#### Transfer counts, annotation, and mappability to your local laptop

```
RBarretTNFATGFBCnt = textread('RBarretTNFATGFB_antisense.ALL.cnt','');
RBarretsamplesTNFATGFB = textread('RBarretTNFATGFBsamples.txt','%s');
RBarretsampleskeysTNFATGFB = textread('samplekeys.txt','%s');
RBarretTNFATGFBmeta_seqdepth=round(sum(RBarretTNFATGFBCnt)/1000000);

Gencode_33_Selected_MappSS=textread('mappability and R code/gencode.v33.Selected.ReadsPerGene.out.MappSS.txt','');
Gencode_33_Selected_MappUS=textread('mappability and R code/gencode.v33.Selected.ReadsPerGene.out.MappUS.txt','');
Gencode_33_Selected_Geneid=textread('mappability and R code/gencode.v33.annotation.Selected.geneid.txt','%s\n');
Gencode_33_Selected_Biotype=textread('mappability and R code/gencode.v33.annotation.Selected.biotype.txt','%s\n');
Gencode_33_Selected_Genename=textread('mappability and R code/gencode.v33.annotation.Selected.genename.txt','%s\n');
```

#### Compile counts

```
RBarretTNFATGFBTPM = RBarretTNFATGFBCnt;
for i=1:size(RBarretTNFATGFBCnt,1)
RBarretTNFATGFBTPM(i,:) = RBarretTNFATGFBCnt(i,:)/Gencode_33_Selected_MappSS(i)*1000;
end
RBarretTNFATGFBTPM(isnan(RBarretTNFATGFBTPM)) = 0;
RBarretTNFATGFBTPM(isinf(RBarretTNFATGFBTPM)) = 0;
for i=1:size(RBarretTNFATGFBTPM,2)
RBarretTNFATGFBTPM(:,i) = RBarretTNFATGFBTPM(:,i)/sum(RBarretTNFATGFBTPM(:,i))*1000000;
end
```
#### Make first dendrogram

```
thisrand = unique(randi([1 size(RBarretTNFATGFBTPM,1)],1,1000));
thisdist = pdist(RBarretTNFATGFBTPM(thisrand,:)');
for i=1:9999
thisrand = unique(randi([1 size(RBarretTNFATGFBTPM,1)],1,1000));
thisdist = thisdist+pdist(RBarretTNFATGFBTPM(thisrand,:)');
end
thisdistmat = squareform(thisdist/10000);
thistree = seqlinkage(thisdistmat,'average', RBarretsamplesTNFATGFB)
plot(thistree)
```

#### Biotypes

```
allbiotypes = unique(Gencode_33_Selected_Biotype);

allbiotypeslength = cell(length(allbiotypes),1);
allbiotypescounts = zeros(length(allbiotypes),size(RBarretTNFATGFBCnt,2));
allbiotypescountspercents = zeros(length(allbiotypes),size(RBarretTNFATGFBCnt,2));


for i=1:length(allbiotypes)
temp = strmatch(allbiotypes{i}, Gencode_33_Selected_Biotype);
allbiotypeslength{i} = length(temp);
if length(temp)>1
allbiotypescounts(i,:) = sum(RBarretTNFATGFBCnt(strmatch(allbiotypes{i}, Gencode_33_Selected_Biotype),:)); 
allbiotypescountspercents(i,:) = allbiotypescounts(i,:)./sum(RBarretTNFATGFBCnt)*100; 
end
end

dlmwrite('allbiotypescountspercents.txt', allbiotypescountspercents,'delimiter','\t')
writetable(cell2table(allbiotypes),'allbiotypes.txt','WriteVariableNames',0)
```
#### Keep only protein coding genes

```
allbiotypes=unique(Gencode_33_Selected_Biotype);
proteincodingindx = strmatch(allbiotypes{15}, Gencode_33_Selected_Biotype);
biotypeindx = proteincodingindx;

additionalgenes = [strmatch('MT-',Gencode_33_Selected_Genename) ; strmatch('H1',Gencode_33_Selected_Genename); strmatch('H2',Gencode_33_Selected_Genename); strmatch('H3',Gencode_33_Selected_Genename); strmatch('H4',Gencode_33_Selected_Genename) ; strmatch('RPL',Gencode_33_Selected_Genename) ; strmatch('RPS',Gencode_33_Selected_Genename)];																
nonadditionalgenes = 1:length(Gencode_33_Selected_Genename);
nonadditionalgenes(additionalgenes) = [];

mappableindx = find(Gencode_33_Selected_MappSS>50);
```

```
finalIndexGeneric = intersect(biotypeindx,intersect(nonadditionalgenes,mappableindx));			

countindx = find(sum(RBarretTNFATGFBCnt')'>150);	
finalIndexGeneric=intersect(finalIndexGeneric,countindx);						

	
RBarretTNFATGFBCnt_GMask = RBarretTNFATGFBCnt(finalIndexGeneric,:);

Gencode_33_Selected_Geneid_GMask = Gencode_33_Selected_Geneid(finalIndexGeneric);
Gencode_33_Selected_Genename_GMask = Gencode_33_Selected_Genename(finalIndexGeneric);
Gencode_33_Selected_MappSS_GMask = Gencode_33_Selected_MappSS(finalIndexGeneric);
Gencode_33_Selected_MappUS_GMask = Gencode_33_Selected_MappUS(finalIndexGeneric);


RBarretTNFATGFBExpression_GMask = RBarretTNFATGFBCnt_GMask;
for i=1:size(RBarretTNFATGFBExpression_GMask,2)
RBarretTNFATGFBExpression_GMask(:,i) = RBarretTNFATGFBCnt_GMask(:,i)/sum(RBarretTNFATGFBCnt_GMask(:,i))*1000000;
end
for i=1:size(RBarretTNFATGFBExpression_GMask)
RBarretTNFATGFBExpression_GMask(i,:) = RBarretTNFATGFBExpression_GMask(i,:)/Gencode_33_Selected_MappSS_GMask(i)*1000;
end
RBarretTNFATGFBExpression_GMask(isnan(RBarretTNFATGFBExpression_GMask)) = 0;
RBarretTNFATGFBExpression_GMask(isinf(RBarretTNFATGFBExpression_GMask)) = 0;

RBarretTNFATGFBCPM_GMask = zeros(size(RBarretTNFATGFBCnt_GMask));
for i=1:size(RBarretTNFATGFBCnt_GMask,2)
RBarretTNFATGFBCPM_GMask(:,i) = RBarretTNFATGFBCnt_GMask(:,i)/sum(RBarretTNFATGFBCnt_GMask(:,i))*1000000;
end

RBarretTNFATGFBTPM_GMask = RBarretTNFATGFBCnt_GMask;
for i=1:size(RBarretTNFATGFBCnt_GMask,1)
RBarretTNFATGFBTPM_GMask(i,:) = RBarretTNFATGFBCnt_GMask(i,:)/Gencode_33_Selected_MappSS_GMask(i)*1000;
end
RBarretTNFATGFBTPM_GMask(isnan(RBarretTNFATGFBTPM_GMask)) = 0;
RBarretTNFATGFBTPM_GMask(isinf(RBarretTNFATGFBTPM_GMask)) = 0;
for i=1:size(RBarretTNFATGFBTPM_GMask,2)
RBarretTNFATGFBTPM_GMask(:,i) = RBarretTNFATGFBTPM_GMask(:,i)/sum(RBarretTNFATGFBTPM_GMask(:,i))*1000000;
end
```

#### Make second dendrogram(the clustering is by treatment)

```
thisrand = unique(randi([1 size(RBarretTNFATGFBTPM_GMask,1)],1,1000));
thisdist = pdist(RBarretTNFATGFBTPM_GMask(thisrand,:)');
for i=1:9999
thisrand = unique(randi([1 size(RBarretTNFATGFBTPM_GMask,1)],1,1000));
thisdist = thisdist+pdist(RBarretTNFATGFBTPM_GMask(thisrand,:)');
end
thisdistmat = squareform(thisdist/10000);
thistree = seqlinkage(thisdistmat,'average', RBarretsampleskeysTNFATGFB)
plot(thistree,'ORIENTATION','top')
```
#### What is the percent of the top 100 genes

```
yall=[];
for i=1:151
[x y]=sort(RBarretTNFATGFBTPM_GMask(:,i),'descend');
yall=unique([y(1:100); yall]);
top100percent(i)=sum(RBarretTNFATGFBTPM_GMask(y(1:100),i))/1000000;
end
```

```
writetable(cell2table(Gencode_33_Selected_Geneid_GMask),'Gencode_33_Selected_Geneid_GMask.txt','WriteVariableNames',0)
writetable(cell2table(Gencode_33_Selected_Genename_GMask),'Gencode_33_Selected_Genename_GMask.txt','WriteVariableNames',0)
dlmwrite('Gencode_33_Selected_MappSS_GMask.txt', Gencode_33_Selected_MappSS_GMask,'delimiter','\t')
dlmwrite('RBarretTNFATGFBTPM_GMask.txt', RBarretTNFATGFBTPM_GMask,'delimiter','\t')
dlmwrite('RBarretTNFATGFBCnt_GMask.txt', RBarretTNFATGFBCnt_GMask,'delimiter','\t')
%writetable(cell2table(RBarretsamplesTNFATGFB),'RBarretsamplesTNFATGFB.txt','WriteVariableNames',0)
```

### In R

```
setwd("/Users/LuC/Desktop/Cedars-Sinai/PROJECTS/IBD_RNASeq/RBARRETTNFATGFB/")#setwd("/Users/samuellu/Desktop/Cedars-Sinai/PROJECTS/IBD_RNASeq/RBARRETTNFATGFB/")#if (!requireNamespace("BiocManager", quietly = TRUE))  #install.packages("BiocManager")#BiocManager::install("BiocLite")#BiocManager::install("IHW")#BiocManager::install("DESeq2")#install.packages("ggplot2")library(DESeq2)library(IHW)library(ggplot2)library(ggrepel)RBarretTNFATGFBCntGMask = as.matrix(read.table("RBarretTNFATGFBCnt_GMask.txt"))sampleNameTNFATGFB = as.matrix(read.table("RBarretTNFATGFBsamples.txt"))sampleKeyTNFATGFB = as.matrix(read.table("samplekeys_Sam.txt"))genenames = as.matrix(read.table("Gencode_33_Selected_Genename_GMask.txt"))#have to separate samplekeys_Sam.txt by "_" to get samplekeys_Sam.tabsampleTableTNFATGFB = read.table("samplekeys_Sam.tab")rownames(sampleTableTNFATGFB)<-sampleKeyTNFATGFBcolnames(sampleTableTNFATGFB)<- c("Treatment","Line","Pheno","Sex","Pass")sampleTableTNFATGFB$Factor <- paste(sampleTableTNFATGFB$Treatment,sampleTableTNFATGFB$Pheno,sep="_")#concatenating the "Line" and "Pass" columns with an underscore separatorsampleTableTNFATGFB$Batch <- paste(sampleTableTNFATGFB$Line,sampleTableTNFATGFB$Pass,sep="_")colnames(RBarretTNFATGFBCntGMask) <- sampleKeyTNFATGFBwrite.table(sampleTableTNFATGFB,file="sampleTableTNFATGFB.txt", sep = "\t", col.names = FALSE)RBarretTNFATGFBCntGMaskBatch <- DESeqDataSetFromMatrix(RBarretTNFATGFBCntGMask, colData= sampleTableTNFATGFB,design= ~Batch)RBarretTNFATGFBCntGMaskBatch <- DESeq(RBarretTNFATGFBCntGMaskBatch)RBarretTNFATGFBCntGMaskFactor <- DESeqDataSetFromMatrix(RBarretTNFATGFBCntGMask, colData= sampleTableTNFATGFB,design= ~Factor)RBarretTNFATGFBCntGMaskFactor <- DESeq(RBarretTNFATGFBCntGMaskFactor)RBarretTNFATGFBCntGMaskBatch_vsd <- varianceStabilizingTransformation(RBarretTNFATGFBCntGMaskBatch,blind=FALSE)RBarretTNFATGFBCntGMaskFactor_vsd <- varianceStabilizingTransformation(RBarretTNFATGFBCntGMaskFactor,blind=FALSE)#PCApcabatch <- prcomp(t(assay(RBarretTNFATGFBCntGMaskBatch_vsd)))#gives the percentage of variance explained by each principal componentpercentVarbatch <- round(100*pcabatch$sdev^2/sum(pcabatch$sdev^2))#pcabatch$rotation is a matrix containing the loadings of the principal components.aloadbatch <- abs(pcabatch$rotation)aloadrelativebatch <- sweep(aloadbatch, 2, colSums(aloadbatch), "/")#pcabatch$x is a matrix containing the scores of each sample (i.e., observation) on each principal componentpcabatchALL <- pcabatch$xpcabatchR<- cbind(pcabatchALL,sampleTableTNFATGFB)#subtract the mean from each value, why only pc1?pcabatchR$PC1 <- scale(pcabatchR$PC1, center = TRUE)RBarretTNFATGFBCntGMaskPC1 <- DESeqDataSetFromMatrix(RBarretTNFATGFBCntGMask, colData= pcabatchR,design= ~PC1)RBarretTNFATGFBCntGMaskPC1 <- DESeq(RBarretTNFATGFBCntGMaskPC1)RBarretTNFATGFBCntGMaskPC1_vsd <- varianceStabilizingTransformation(RBarretTNFATGFBCntGMaskPC1,blind=FALSE)#Plotggplot(pcabatchR, aes(PC1, PC2, color= Pheno)) +  geom_point(aes(size= Treatment),alpha=0.6,stroke = 3)+geom_point(aes(size= Pheno),color="black",alpha=0.2) +  xlab(paste0("PC1: ",percentVarbatch[1],"% variance")) +  ylab(paste0("PC2: ",percentVarbatch[2],"% variance")) +  geom_text_repel(aes(label = sampleKeyTNFATGFB),size=4,box.padding   = 0.35, point.padding = 0.5,segment.color = 'grey50')+ theme_bw()ggplot(pcabatchR, aes(PC1, PC2, color= Sex)) +  geom_point(aes(size= Treatment),alpha=0.6,stroke = 3)+geom_point(aes(size= Treatment),color="black",alpha=0.2) +  xlab(paste0("PC1: ",percentVarbatch[1],"% variance")) +  ylab(paste0("PC2: ",percentVarbatch[2],"% variance")) + theme_bw()ggplot(pcabatchR, aes(PC2, PC9, color= Pheno)) +  geom_point(aes(size= Treatment),alpha=0.6,stroke = 3)+geom_point(aes(size= Treatment),color="black",alpha=0.2) +  xlab(paste0("PC2: ",percentVarbatch[2],"% variance")) +  ylab(paste0("PC9: ",percentVarbatch[9],"% variance")) +  geom_text_repel(aes(label = sampleKeyTNFATGFB),size=4,box.padding   = 0.35, point.padding = 0.5,segment.color = 'grey50')+ theme_bw()# The resulting vector coul will contain 12 colors from the "Set3" palette.library(RColorBrewer)coul <- brewer.pal(12, "Set3")# generates colors for a plot based on the batch variablecolors=pcabatchR$Batchallbatches<-unique(pcabatchR$Batch)for (i in 1:38){  colors[pcabatchR$Batch==allbatches[i]]<-coul[i%%12+1]}thinlines=c(seq(4,72,8),75,seq(83,151,8))thicklines=c(seq(8,72,8),79,seq(87,151,8))# first half of the barplot would be the non-fibrotic group and the second part would be the fibrotic group# The order would be CC, TG, TN, TTsamplesorder=c(4,1,3,2,8,5,7,6,28,25,27,26,32,29,31,30,36,33,35,34,40,37,39,38,52,49,51,50,55,53,55,54,68,65,67,66,72,69,71,70,76,73,75,74,80,77,79,78,84,81,83,82,88,85,87,86,116,113,115,114,120,117,119,118,139,136,138,137,143,140,142,141,12,9,11,10,16,13,15,14,20,17,19,18,24,21,23,22,44,41,43,42,48,45,47,46,60,57,59,58,64,61,63,62,92,89,91,90,96,93,95,94,100,97,99,98,104,101,103,102,108,105,107,106,112,109,111,110,124,121,123,122,128,125,127,126,132,129,131,130,135,133,134,147,144,146,145,151,148,150,149)#create 38 barplots and saving each of them as a PNG filefor (i in 1:38) {  filename = paste("PC_",i,".png", sep = "")  png(filename)  barplot(pcabatchALL[samplesorder,i],col=colors[samplesorder],las=2,xaxt='n',space=0)  for (i in 1:length(thinlines)) {    abline(v = thinlines[i], col = "black",lty = 3)  }  for (i in 1:length(thicklines)) {    abline(v = thicklines[i], col = "black",lty = 1)  }  abline(v = 72, col = "red",lty = 1)  dev.off()}write.csv(aloadrelativebatch,file="aloadrelativeMask_batchmodel_filtered.csv")write.csv(pcabatch$x,file="pca_batchmodel_x.csv")# make a spreadsheet## sheeet2sheet2_1 <- list("Biotypes")sheet2_2 <- sampleKeyTNFATGFBcombined_sheet2 <- c(sheet2_1, sheet2_2)combined_spreadsheet2 <- as.matrix(read.table("combine_allbiotypes_percents.txt"))colnames(combined_spreadsheet2) <- combined_sheet2write.table(combined_spreadsheet2,file="combined_spreadsheet2.txt", sep = "\t", row.names = FALSE)## sheeet3sheet3_1 <- list("Genename","Genename","Geneid","Mapp","PC1","PC2","PC3")sheet3_2<- sampleKeyTNFATGFBcombined_headers <- c(headers_1, headers_2)combined_spreadsheet <- as.matrix(read.table("combine_test.txt"))colnames(combined_spreadsheet) <- combined_headerscombined_spreadsheet <- combined_spreadsheet[order(combined_spreadsheet[,1]),] #sort by the first columnwrite.table(combined_spreadsheet,file="combined_spreadsheet.txt", sep = "\t", row.names = FALSE)#clean these two csv files in terminal and in matlab
```



 <br> 

## Future works 



<br> 



## References 

<br> # Myofibroblast
