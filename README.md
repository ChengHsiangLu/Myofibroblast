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

Inflammatory bowel diseases (IBD) are a group of chronic conditions that cause inflammation and damage to the digestive tract. The two main types of IBD are Crohn's Disease (CD) and Ulcerative Colitis (UC). Both Crohn's disease and ulcerative colitis are chronic conditions, meaning they can last for a lifetime and require ongoing treatment to manage symptoms and prevent complications. 

UC affects only the colon and rectum and causes symptoms such as bloody diarrhea, abdominal pain, and a frequent need to pass stools. It can also lead to complications such as inflammation of the skin, eyes, and joints. On the other hand, CD can affect any part of the digestive tract, like small or large intestine, and can cause symptoms such as abdominal pain, diarrhea, weight loss, and fatigue. It can also cause complications such as fistulas (abnormal connections between different parts of the intestine) 

On of the most prevalent complication of CD is the onset of fibrotic complications and strictures (narrowing of the intestine)[1]. The molecular mechanisms involved in these phenotypes remain largely unknown, and as a result, there are currently no effective drugs to prevent or treat stricturing CD.


<br>

### Induced Pluripotent Stem Cells(iPSCs) 

Induced pluripotent stem cells (iPSCs) are a type of stem cell that are generated in the laboratory by reprogramming adult cells, such as skin or blood cells, to a pluripotent state. A pluripotent state means that the cells have the potential to develop into any type of cell in the body, just like embryonic stem cells. 


iPSCs offer several advantages as they can be generated from the patient's own cells, avoiding issues with immune rejection, ethical concerns and the need for embryos. 


iPSCs can be differentiated into multiple cell and tissue types, and can therefore be used to study the underlying causes of diseases, test new drugs and therapies, and potentially generate replacement tissues or organs for transplantation. As iPSCs possess the same genetic background as the patient they are derived from, they are considered an instrumental tool in the field of personalized and precision medicine[2].

<br> 

## Aim 

In this project, we aim to take advantage of iPSC lines to unveil specific signaling pathways specifically affected in patients with fibrotic CD[3]. We will be analyzing RNA-seq data from 19 iPSC lines that were differentiated into gut mesenchymal organoids. This panel comprises 10 iPSC lines derived from Crohn's disease patients that suffered fibrotic complications, and 9 lines from patients with non-fibrotic disease. Each iPSC line was differentiated into mesenchymal organoids in two independent replicates, and each was subjected to 4 different treatments:

* untreated 
* TGFb (a pro-fibrotic cytokine)[4]
* TNFa (a pro-inflammatory cytokine) 
* and the combination of TGF-b+TNF-a

The final RNA-Seq dataset comprised a total of 151 samples(1 library failed). The objective is to 

* investigate if the effect of the four different treatments in iPSC-derived organoids recapitulate the expected responses observed in-vivo[5].
* perform differential expression analysis to identify genes that show differential responses between fibrotic and non-fibrotic patients[6].

<br> 

## Pipelines 

### In HPC

#### Convert 151 Fastq to Fasta files 

First, put all fasq.gz files into one folder and list all fastq files’ name in fastqfiles.txt.

```

ls *q.gz > fastqfiles.txt

``` 

Cut redundant suffix “\_R1_trimmed” and list all fastq files’ name in libraryname.txt and preffix.txt.

```

ls *q.gz | cut -f 1 -d '.' | sed 's/_R1_trimmed//g'  >libraryname.txt 

ls *q.gz | cut -f 1 -d '.' | sed 's/_R1_trimmed//g' > preffix.txt

``` 

Form a table with 3 columns: fastqfiles.txt libraryname.txt preffix.txt.

```

paste fastqfiles.txt libraryname.txt preffix.txt > tofastatable.txt

``` 

Create small-sized fasta-formatted files. To submit this job to the cluster on HPC, you need to read the file, library, and prefix. Once you have done that, run the script "generatefastaFromFastaqz" which combine the script "DCfastaqTofastaLibraryId.pl". This results in small-sized fasta-formatted files contain only one header and one sequence per read. You can find all scripts in the "scripts" folder. 

``` 

cat tofastatable.txt | awk {print}' | while read file library preffix ; do qsub -cwd -o $PWD -e $PWD -l h_data=2048M,h_rt=8:00:00 $HOME/scripts/generatefastaFromFastaqz $file $library $preffix  

done 

``` 
Figure 1 shows what each fasta-formatted file look like.

![](/Pics/Figure_1.png)

<br>

#### Generate auxiliary files and directories for each sample 

Put a list of names of all fasta files in the directory and save them in a text file named "fastafiles.txt".

```

ls *fasta.gz > fastafiles.txt

``` 

Cut the redundant suffix ".fasta.gz" from the names of all fasta files and generate a new list of file names with the suffix removed in a text file named “targetdirectories.GTF.txt”.

```

cat fastafiles.txt |sed 's/.fasta.gz//g' > targetdirectories.GTF.txt

``` 

Create a separate directory for each sample listed in "fastafiles.txt".

```

cat fastafiles.txt |while read line ; do mkdir ${line/.fasta\.gz/GTFpass1/} ; done

``` 

<br> 

#### Form the submission script called “sendmyof” 

Add a shebang line at the beginning of your script file named "sendmyof" to indicate the interpreter that should be used to execute the script.

```

echo '#!/bin/bash/' > sendmyof

``` 

The command below runs the "generatesendscriptSingleGTFParam" script with several input parameters to map the RNA-seq data with STAR. The input parameters include the list of target directories containing the input data ("targetdirectories.GTF.txt"), the directory prefix for pass-1 alignments ("GTFpass1"), a parameter file containing settings for STAR alignment ("Parameters.txt"), a prefix for individual submission scripts to HPC ("myof"), the path to the STAR index directory ("/home/luc/RNASEQ_MASTER/Hsapiens/GRC38/INDEXES/GRCh38.primary.33.basicselected.STAR2.7.3a/"), the path to the input data directory ("/home/luc/iPSC/MYOFIBROBLAST/"), the amount of free memory to use ("mem_free=32G"), and the number of threads to use ("8"). In the end, it will generate a sample-specific sumission script called "processLaneSingleGTFParam" in each sample's folder:
 
``` 

./generatesendscriptSingleGTFParam targetdirectories.GTF.txt GTFpass1 Parameters.txt myof /home/luc/RNASEQ_MASTER/Hsapiens/GRC38/INDEXES/GRCh38.primary.33.basicselected.STAR2.7.3a/ /home/luc/iPSC/MYOFIBROBLAST/ mem_free=32G 8 >> sendmyof 

``` 

Change sendmyof into executable mode and run sendmyof.

``` 

chmod a+x sendmyof 
. sendmyof 

``` 

It will take less than one day to run through 151 samples and generate each sample a folder which contain every output from STAR.

<br> 

#### Create a table summarizing the mapping statistics for each sample 

Change directory into one sample file which ends with "GTFpass1". Extract the first column from the mapping statistics file and store it in “temp2.txt”.

``` 

grep "|" 008iP22TGFbM_S71GTFpass1/008iP22TGFbM_S71GTFpass1Log.final.out | cut -f 1 -d "|" | sed 's/^  *//g' | awk 'NR>3 {print}' > temp2.txt 

```

The first column from the mapping statistics file in Figure 2.

![](/Pics/Figure_2.png)

<br>

Create an empty temporary file for storing intermediate results.

```

rm tempprev.txt 
touch tempprev.txt

```

Extract the total mapped reads from each subsequent mapping statistics file and combine with previous results.

```

ls *pass1/*final.out | while read line ; do  
grep "|" $line | cut -f 2  > temp.txt 
paste tempprev.txt temp.txt > tempnew.txt 
mv tempnew.txt tempprev.txt 
done 

```

Remove the first column and write the final results to a file called “mappingstatsFirstpass.txt”.

```

cut -f 2- tempprev.txt | awk 'NR>3 {print}'  > tempnew.txt 
mv tempnew.txt tempprev.txt 
paste temp2.txt tempprev.txt > mappingstatsFirstpass.txt 

``` 

Figure 3 shows what mappingstatsFirstpass.txt look like.

![](/Pics/Figure_3.png)

The summary statistics show good rates of unique alignments for all samples.

<br> 

#### Counts 

Generate a directory called "COUNTS" and copy all gene count files to this folder and then clean all file names.

``` 

mkdir COUNTS 
cp *pass1/*PerGene* COUNTS/ 

``` 

``` 

ls *tab|while read line ; do mv $line ${line/GTFpass1ReadsPerGene.out/} ; done 

```

For each count file, extract and create the five count tables.

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

Make a Geneid list from one of the count tables as "countsannot_GRCh38.primary.Selected.Geneid.txt".

```

ls 008iP22TGFbM_S71.tab | head -1 | while read line; do  
cut -f 1 $line | awk 'NR>4{print}' > countsannot_GRCh38.primary.Selected.Geneid.txt 
done 

```

Create a file listing the names of all samples as "RBarretTNFATGFBsamples.txt".

```

ls *.sense.tab | sed 's/.sense.tab//g' | tr -s " " "\n" | sed 's/_1//g' > RBarretTNFATGFBsamples.txt 

```

Make count tables for sense, anti-sense, nostrand, ambiguous, and nofeature reads.

```  

# combine all sense counts into RBarretTNFATGFB_sense.ALL.cnt 
paste *.sense.tab > RBarretTNFATGFB_sense.ALL.cnt   

# combine all antisense counts into RBarretTNFATGFB_antisense.ALL.cnt 
paste *.antisense.tab > RBarretTNFATGFB_antisense.ALL.cnt   

# combine all nostrand counts into RBarretTNFATGFB_nostrand.ALL.cnt 
paste *.nostrand.tab > RBarretTNFATGFB_nostrand.ALL.cnt 

# combine all ambiguous counts into RBarretTNFATGFB_ambiguous.cnt 
cat *ambiguous.tab > RBarretTNFATGFB_ambiguous.cnt   

# combine all nofeature counts into RBarretTNFATGFB_nofeature.cnt 
cat *nofeature.tab > RBarretTNFATGFB_nofeature.cnt 

``` 

Next, I am going to use "RBarretTNFATGFB_antisense.ALL.cnt" file for the further analysis, as this matrix contains the counts matching the strand-specificity of the RNA-Seq libraries generated in this study.

<br> 

### In MATLAB

#### Transfer data and import annotation

Transfer the counts, annotation, and mappability data to your local laptop.

```

RBarretTNFATGFBCnt = textread('RBarretTNFATGFB_antisense.ALL.cnt','');
RBarretsamplesTNFATGFB = textread('RBarretTNFATGFBsamples.txt','%s');
RBarretsampleskeysTNFATGFB = textread('samplekeys_Sam.txt','%s');

% calculate the sum of the counts in RBarretTNFATGFBCnt, divides the result by 1000000, and rounds the result to the nearest integer. 
RBarretTNFATGFBmeta_seqdepth=round(sum(RBarretTNFATGFBCnt)/1000000);

```

The following files contain the annotation and gene effective lenghts (mappabilities) for the human gene annoation used for alignment, and can be found in the  "mappability and R code" folder.

```

Gencode_33_Selected_MappSS=textread('mappability and R code/gencode.v33.Selected.ReadsPerGene.out.MappSS.txt','');
Gencode_33_Selected_MappUS=textread('mappability and R code/gencode.v33.Selected.ReadsPerGene.out.MappUS.txt','');
Gencode_33_Selected_Geneid=textread('mappability and R code/gencode.v33.annotation.Selected.geneid.txt','%s\n');
Gencode_33_Selected_Biotype=textread('mappability and R code/gencode.v33.annotation.Selected.biotype.txt','%s\n');
Gencode_33_Selected_Genename=textread('mappability and R code/gencode.v33.annotation.Selected.genename.txt','%s\n');

```

<br>

#### Compile counts

First, initialize a new variable called RBarretTNFATGFBTPM with the same count data as RBarretTNFATGFBCnt.
Then, iterates over each gene in the count data matrix. For each gene, the corresponding row in RBarretTNFATGFBTPM is updated by dividing the count data by the gene effective length from the "Gencode\_33\_Selected\_MappSS", multiplying by 1000, and storing the result in RBarretTNFATGFBTPM.

Finally, iterates over each sample in the TPM data matrix. For each sample, the corresponding column in RBarretTNFATGFBTPM is updated by dividing the values in the column by the sum of the values in the column, multiplying by 1,000,000, and storing the result in RBarretTNFATGFBTPM. This step **normalizes the TPM values** across samples and scales the resulting values to TPM.

```

RBarretTNFATGFBTPM = RBarretTNFATGFBCnt;
for i=1:size(RBarretTNFATGFBCnt,1)
% divid the gene count matrix RBarretTNFATGFBCnt by Gencode_33_Selected_MappSS matrix, which is the sum of the transcript length of each gene
RBarretTNFATGFBTPM(i,:) = RBarretTNFATGFBCnt(i,:)/Gencode_33_Selected_MappSS(i)*1000;
end
% set any NaN or Inf values resulting from the normalization process to 0
RBarretTNFATGFBTPM(isnan(RBarretTNFATGFBTPM)) = 0;
RBarretTNFATGFBTPM(isinf(RBarretTNFATGFBTPM)) = 0;
for i=1:size(RBarretTNFATGFBTPM,2)
% scale the TPM values so that the sum of expression values across each sample of the matrix is equal to 1,000,000. This ensures that the expression values are comparable across different samples and allows meaningful comparisons of gene expression levels between different samples.
RBarretTNFATGFBTPM(:,i) = RBarretTNFATGFBTPM(:,i)/sum(RBarretTNFATGFBTPM(:,i))*1000000;
end

```
<br>

#### Make the first dendrogram

Make a dendrogram to visualize the relationships among samples in the RBarretTNFATGFB dataset based on their gene expression profiles. To downgrade the effect of potential expression outliers in the dendrogram and compute a more robust sample clustering:

1. I generates a random selection of 1,000 genes from the TPM data matrix.
2. Calculates the pairwise distances between the selected genes.
3. Creates a for loop that iterates 9,999 times. For each iteration, a new random selection of 1,000 genes is generated, and the pairwise distances between these genes are added to the previous 'thisdist' calculation.
4. Converts the one-dimensional distance vector 'thisdist' into a distance matrix 'thisdistmat' using the 'squareform' function.
5. Generates a hierarchical clustering tree based on the distance matrix 'thisdistmat'.

Overall, I perform a clustering analysis on a subset of genes in the RBarretTNFATGFB dataset to visualize the relationships among samples based on their gene expression profiles in Figure 4.

```

thisrand = unique(randi([1 size(RBarretTNFATGFBTPM,1)],1,1000));
thisdist = pdist(RBarretTNFATGFBTPM(thisrand,:)');
for i=1:9999
thisrand = unique(randi([1 size(RBarretTNFATGFBTPM,1)],1,1000));
thisdist = thisdist+pdist(RBarretTNFATGFBTPM(thisrand,:)');
end
thisdistmat = squareform(thisdist/10000);
thistree = seqlinkage(thisdistmat,'average', RBarretsamplesTNFATGFB)
plot(thistree, 'ORIENTATION', 'top')

```

![](/Pics/Figure_4.png)

A first observation is that the major clusters are formed by samples from the same treatment, with exceptions. Therefore, the Treatment factor seems to be the dominant source of gene expression variation in this experiment. 
<br>

#### All biotypes counts percents

As part of the preliminary quality control, the following estimates the relative contribution of each gene biotype to the expression matrix:

```

% use the unique function and stored the allbiotypes variable.
allbiotypes = unique(Gencode_33_Selected_Biotype);

% create A cell array allbiotypeslength to store the lengths of each biotype name.
allbiotypeslength = cell(length(allbiotypes),1);

% two new matrices, allbiotypescounts and allbiotypescountspercents, are initialized with zeros. These matrices have dimensions (number of unique biotypes) x (number of samples in the TPM data). They will be used to store the number of reads (counts) and the percentage of total reads (%TPM) for each biotype in each sample.
allbiotypescounts = zeros(length(allbiotypes),size(RBarretTNFATGFBCnt,2));
allbiotypescountspercents = zeros(length(allbiotypes),size(RBarretTNFATGFBCnt,2));

for i=1:length(allbiotypes)
%  finds all the indices of Gencode_33_Selected_Biotype that match the current biotype. Then, returned a vector of **indices** where the biotype occurs in Gencode_33_Selected_Biotype.
temp = strmatch(allbiotypes{i}, Gencode_33_Selected_Biotype);
% Stored the length of the temp vector represents the number of genes with the current biotype in the Gencode_33_Selected_Biotype.
allbiotypeslength{i} = length(temp);
if length(temp)>1
% Sum the expression values for all genes with the current biotype across all samples. The resulting sums are stored in the corresponding row of the allbiotypescounts matrix.
allbiotypescounts(i,:) = sum(RBarretTNFATGFBCnt(strmatch(allbiotypes{i}, Gencode_33_Selected_Biotype),:)); 
allbiotypescountspercents(i,:) = allbiotypescounts(i,:)./sum(RBarretTNFATGFBCnt)*100; 
end
end

dlmwrite('allbiotypescountspercents.txt', allbiotypescountspercents,'delimiter','\t')
writetable(cell2table(allbiotypes),'allbiotypes.txt','WriteVariableNames',0)

```

Using Excel, create a spreadsheet using "allbiotypes.txt" and "allbiotypescountspercents.txt", and calculate the minimum, maximum, and average values for each biotype in Figure 5. You can access my completed spreadsheet [here](/spreadsheet/Barret_Myofibroblast_TGFTNF_MASTER.xlsx). Notably, protein_coding genes exhibit an average of 98.65% among the various biotypes, consistent with my expectations. Moreover, I found no samples with excessive contributions from other biotypes (e.g. mitochondrial and non-coding RNAs), and therefore no library quality issues were found in this step. 

![](/Pics/Figure_5.png)

<br>

#### Protein coding genes

The "allbiotypes.txt" file contains multiple biotypes. For the next step, I will only retain the "protein_coding" biotype. I will also remove some gene classes that tipically show very noisy or variable gene expression across different samples (e.g histone and ribosomal genes, among others). 

```

allbiotypes=unique(Gencode_33_Selected_Biotype);
% finds the 15th unique value, which is protein_coding, of Gencode_33_Selected_Biotype in the array proteincodingindx.
proteincodingindx = strmatch(allbiotypes{15}, Gencode_33_Selected_Biotype);
biotypeindx = proteincodingindx;

% creates an array additionalgenes contains the indices of genes that have certain prefixes such as 'MT-', 'H1', 'H2', 'H3', 'H4', 'RPL', or 'RPS' in their names.
additionalgenes = [strmatch('MT-',Gencode_33_Selected_Genename) ; strmatch('H1',Gencode_33_Selected_Genename); strmatch('H2',Gencode_33_Selected_Genename); strmatch('H3',Gencode_33_Selected_Genename); strmatch('H4',Gencode_33_Selected_Genename) ; strmatch('RPL',Gencode_33_Selected_Genename) ; strmatch('RPS',Gencode_33_Selected_Genename)];

% creates an array nonadditionalgenes with the same length as the Gencode_33_Selected_Genename array.	
nonadditionalgenes = 1:length(Gencode_33_Selected_Genename);
% remove the indices of genes in additionalgenes from the nonadditionalgenes array.
nonadditionalgenes(additionalgenes) = [];

% mappableindx contains the indices of elements in the Gencode_33_Selected_MappSS array that are greater than 50.
mappableindx = find(Gencode_33_Selected_MappSS>50);

```

```

% a new variable finalIndexGeneric which is the intersection of three other variables: biotypeindx, nonadditionalgenes, and mappableindx.
finalIndexGeneric = intersect(biotypeindx,intersect(nonadditionalgenes,mappableindx));			
% find the indices of rows in RBarretTNFATGFBCnt that have a sum greater than 150 (an average of >1 per sample). 
countindx = find(sum(RBarretTNFATGFBCnt')'>150);

% update finalIndexGeneric to be the intersection of finalIndexGeneric and countindx.
finalIndexGeneric=intersect(finalIndexGeneric,countindx);

```						

```

% create a new variable RBarretTNFATGFBCnt_GMask which is a subset of RBarretTNFATGFBCnt corresponding to the rows indexed by finalIndexGeneric.
RBarretTNFATGFBCnt_GMask = RBarretTNFATGFBCnt(finalIndexGeneric,:);

Gencode_33_Selected_Geneid_GMask = Gencode_33_Selected_Geneid(finalIndexGeneric);
Gencode_33_Selected_Genename_GMask = Gencode_33_Selected_Genename(finalIndexGeneric);
Gencode_33_Selected_MappSS_GMask = Gencode_33_Selected_MappSS(finalIndexGeneric);
Gencode_33_Selected_MappUS_GMask = Gencode_33_Selected_MappUS(finalIndexGeneric);

```

```

% normalize the expression data like we did previously, keeping only the filtered set of genes above:
RBarretTNFATGFBExpression_GMask = RBarretTNFATGFBCnt_GMask;
for i=1:size(RBarretTNFATGFBExpression_GMask,2)
RBarretTNFATGFBExpression_GMask(:,i) = RBarretTNFATGFBCnt_GMask(:,i)/sum(RBarretTNFATGFBCnt_GMask(:,i))*1000000;
end
for i=1:size(RBarretTNFATGFBExpression_GMask)
RBarretTNFATGFBExpression_GMask(i,:) = RBarretTNFATGFBExpression_GMask(i,:)/Gencode_33_Selected_MappSS_GMask(i)*1000;
end
RBarretTNFATGFBExpression_GMask(isnan(RBarretTNFATGFBExpression_GMask)) = 0;
RBarretTNFATGFBExpression_GMask(isinf(RBarretTNFATGFBExpression_GMask)) = 0;

% RBarretTNFATGFBCPM_GMask contains the expression data normalized only by CPM, using the same normalization method as the code above.
RBarretTNFATGFBCPM_GMask = zeros(size(RBarretTNFATGFBCnt_GMask));
for i=1:size(RBarretTNFATGFBCnt_GMask,2)
RBarretTNFATGFBCPM_GMask(:,i) = RBarretTNFATGFBCnt_GMask(:,i)/sum(RBarretTNFATGFBCnt_GMask(:,i))*1000000;
end

% RBarretTNFATGFBTPM_GMask contains the expression data normalized only by TPM.
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

<br>

#### Dendrogram with only protein coding genes


Perform hierarchical clustering on the filtered gene expression data stored in the variable RBarretTNFATGFBTPM_GMask:

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

Agains, the samples are clustered largely by their treatment status but in a more consistent fashion as compared to the unfiltered dataset in Figure 6.

![](/Pics/Figure_6.png)

<br>

#### The percent of the top 100 genes

Another item for quality control is achieved by calculating the percent of signal attributed to the top 100 expressed genes in each sample based on their transcript per million (TPM) values in the RBarretTNFATGFBTPM_GMask matrix.

```

% iterates over 151 samples, it first sorts the TPM values of all genes in descending order and stores the indices of the sorted genes in y. The top 100 expressed genes in the sample are obtained by selecting the first 100 indices in y, and these indices are appended to a running list of all top 100 indices yall.
yall=[];
for i=1:151
[x y]=sort(RBarretTNFATGFBTPM_GMask(:,i),'descend');
yall=unique([y(1:100); yall]);
top100percent(i)=sum(RBarretTNFATGFBTPM_GMask(y(1:100),i))/1000000;
end

```

I find that, for some samples, the top 100 most-expressed genes accumulate ~40% of the total TPMs for the sample, while the average is ~25%. I will keep track of these number in case those samples show outlier behaviour in downstream analyses. 

In the end, we store the Gencode\_33\_Selected\_Geneid\_GMask.txt, Gencode\_33\_Selected\_Genename\_GMask.txt, Gencode\_33\_Selected\_MappSS\_GMask.txt, RBarretTNFATGFBTPM\_GMask.txt, and RBarretTNFATGFBCnt\_GMask.txt for ours further analysis in R.

```

writetable(cell2table(Gencode_33_Selected_Geneid_GMask),'Gencode_33_Selected_Geneid_GMask.txt','WriteVariableNames',0)
writetable(cell2table(Gencode_33_Selected_Genename_GMask),'Gencode_33_Selected_Genename_GMask.txt','WriteVariableNames',0)
dlmwrite('Gencode_33_Selected_MappSS_GMask.txt', Gencode_33_Selected_MappSS_GMask,'delimiter','\t')
dlmwrite('RBarretTNFATGFBTPM_GMask.txt', RBarretTNFATGFBTPM_GMask,'delimiter','\t')
dlmwrite('RBarretTNFATGFBCnt_GMask.txt', RBarretTNFATGFBCnt_GMask,'delimiter','\t')

```

<br>

### Differential expression analysis in R

#### Install packages

First, install packages BiocManager, BiocLite, IHW, DESeq2[7], and ggplot2. Then, read in RBarretTNFATGFBCnt\_GMask.txt, RBarretTNFATGFBsamples.txt, samplekeys\_Sam.txt, and Gencode\_33\_Selected\_Genename\_GMask.txt.

```

setwd("/Users/LuC/Desktop/Cedars-Sinai/PROJECTS/IBD_RNASeq/RBARRETTNFATGFB/")
#setwd("/Users/samuellu/Desktop/Cedars-Sinai/PROJECTS/IBD_RNASeq/RBARRETTNFATGFB/")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("BiocLite")
#BiocManager::install("IHW")
#BiocManager::install("DESeq2")
#install.packages("ggplot2")

library(DESeq2)
library(IHW)
library(ggplot2)
library(ggrepel)

RBarretTNFATGFBCntGMask = as.matrix(read.table("RBarretTNFATGFBCnt_GMask.txt"))
sampleNameTNFATGFB = as.matrix(read.table("RBarretTNFATGFBsamples.txt"))
sampleKeyTNFATGFB = as.matrix(read.table("samplekeys_Sam.txt"))
genenames = as.matrix(read.table("Gencode_33_Selected_Genename_GMask.txt"))

```

<br>

#### Form samplekeys_Sam.tab

Separate samplekeys\_Sam.txt by "\_" to get samplekeys\_Sam.tab before next step. Here are my code in terminal.

```

#In terminal
#Create an empty file to store the output
touch samplekeys_Sam.tab

#Loop over the sample names and split them by "_"
for sample in $(cat samplekeys_Sam.txt); do
    IFS=_ read -r col1 col2 col3 col4 col5 <<< "$sample"
    echo -e "$col1\t$col2\t$col3\t$col4\t$col5" >> samplekeys_Sam.tab 
done

```

<br>

#### Generate a sampleTableTNFATGFB 

The sampleTableTNFATGFB contains the experimental factors: Treatment, iPSc line (Line), fibrotic phenotype (Pheno), Sex, number of iPSC passages (Pass), the combination of phenotype and treatment (Factor), and the combination of line and passages (Batch) in Figure 7.

```

sampleTableTNFATGFB = read.table("samplekeys_Sam.tab")
rownames(sampleTableTNFATGFB)<-sampleKeyTNFATGFB
colnames(sampleTableTNFATGFB)<- c("Treatment","Line","Pheno","Sex","Pass")
sampleTableTNFATGFB$Factor <- paste(sampleTableTNFATGFB$Treatment,sampleTableTNFATGFB$Pheno,sep="_")
#concatenating the "Line" and "Pass" columns with an underscore separator
sampleTableTNFATGFB$Batch <- paste(sampleTableTNFATGFB$Line,sampleTableTNFATGFB$Pass,sep="_")
colnames(RBarretTNFATGFBCntGMask) <- sampleKeyTNFATGFB
write.table(sampleTableTNFATGFB,file="sampleTableTNFATGFB.txt", sep = "\t", col.names = FALSE)

```

![](/Pics/Figure_7.png)

<br>

#### DESeq2 package

Different experimental factors are tested while creating the DESEq object for differential expression, to check if there are significant differences. The RBarretTNFATGFBCntGMaskBatch is created with the **Batch** information specified in the design formula, while the RBarretTNFATGFBCntGMaskFactor is created with the **treatment and phenotype** information specified in the design formula. I also tested if fitting the data to the first principal component (PC1, see below) makes a difference in the first steps.

The DESeq function is used to estimate size factors and dispersion values for the DESeqDataSet objects. Using this object, we first use the varianceStabilizingTransformation function  to perform variance stabilizing transformation. This transformation is important for reducing the effect of noise and impose heteroscedasticity in the data, making it more suitable for downstream analyses such as linear modeling and clustering.

```

RBarretTNFATGFBCntGMaskBatch <- DESeqDataSetFromMatrix(RBarretTNFATGFBCntGMask, colData= sampleTableTNFATGFB,design= ~Batch)
RBarretTNFATGFBCntGMaskBatch <- DESeq(RBarretTNFATGFBCntGMaskBatch)

RBarretTNFATGFBCntGMaskFactor <- DESeqDataSetFromMatrix(RBarretTNFATGFBCntGMask, colData= sampleTableTNFATGFB,design= ~Factor)
RBarretTNFATGFBCntGMaskFactor <- DESeq(RBarretTNFATGFBCntGMaskFactor)

RBarretTNFATGFBCntGMaskPC1 <- DESeqDataSetFromMatrix(RBarretTNFATGFBCntGMask, colData= pcabatchR,design= ~PC1)
RBarretTNFATGFBCntGMaskPC1 <- DESeq(RBarretTNFATGFBCntGMaskPC1)

RBarretTNFATGFBCntGMaskBatch_vsd <- varianceStabilizingTransformation(RBarretTNFATGFBCntGMaskBatch,blind=FALSE)
RBarretTNFATGFBCntGMaskFactor_vsd <- varianceStabilizingTransformation(RBarretTNFATGFBCntGMaskFactor,blind=FALSE)
RBarretTNFATGFBCntGMaskPC1_vsd <- varianceStabilizingTransformation(RBarretTNFATGFBCntGMaskPC1,blind=FALSE)

```

<br>

#### Principal Component Analysis (PCA)

```

#perform principal component analysis (PCA) on the variance-stabilized counts data 
pcabatch <- prcomp(t(assay(RBarretTNFATGFBCntGMaskBatch_vsd)))

#give the percentage of variance explained by each principal component
percentVarbatch <- round(100*pcabatch$sdev^2/sum(pcabatch$sdev^2))

#pcabatch$rotation is a matrix containing the loadings of the principal components.
aloadbatch <- abs(pcabatch$rotation)

#normalize the loadings in aloadbatch so that each column (i.e., PC) sums to 1. 
aloadrelativebatch <- sweep(aloadbatch, 2, colSums(aloadbatch), "/")

#pcabatch$x is a matrix containing each sample's coordinate on each principal component
pcabatchALL <- pcabatch$x
pcabatchR<- cbind(pcabatchALL,sampleTableTNFATGFB)






```

<br>

#### PCA plots

```

ggplot(pcabatchR, aes(PC1, PC2, color= Treatment)) +
  geom_point(aes(size= Pheno),alpha=0.6,stroke = 3)+
  xlab(paste0("PC1: ",percentVarbatch[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarbatch[2],"% variance")) +
  geom_text_repel(aes(label = sampleKeyTNFATGFB),size=4,box.padding   = 0.35, point.padding = 0.5,segment.color = 'grey50')+ theme_bw()



```

The plot shows the clustering of all samples using the first two principal components, **PC1 and PC2**, colored by **Pheno** variable, with the point size indicating the **Treatment** variable.

In Figure 8, four distinct groups were formed based on their treatment: the CC group (untreated) is located in the right corner, the TG group (treated with TGF-b) is located at the bottom, the TN group (treated with TNF-a) is located in the right corner, and the TT group (treated with both TGF-b and TNF-a) is located at the top. These 4 groups were differentiated based on the combination of both PC1 and PC2, which accounted for 19% and 17% of the variance, respectively.

![](/Pics/Figure_8.jpeg)

<br>


#### A series of bar plots (one for each principal component)

Each bar plot represents the coordinates of all samples on a given principal component. These plots provide a quick visual evaluation of the potential association of each component with specific experimental factors.

```

#the resulting vector coul will contain 12 colors from the "Set3" palette.
library(RColorBrewer)
coul <- brewer.pal(12, "Set3")
#generates colors for a plot based on the batch variable
colors=pcabatchR$Batch
allbatches<-unique(pcabatchR$Batch)
for (i in 1:38){
  colors[pcabatchR$Batch==allbatches[i]]<-coul[i%%12+1]
}

thinlines=c(seq(4,72,8),75,seq(83,151,8))
thicklines=c(seq(8,72,8),79,seq(87,151,8))

#first half of the barplot would be the non-fibrotic group and the second part would be the fibrotic group
#the order would be CC, TG, TN, TT
samplesorder=c(4,1,3,2,8,5,7,6,28,25,27,26,32,29,31,30,36,33,35,34,40,37,39,38,52,49,51,50,55,53,55,54,68,65,67,66,72,69,71,70,76,73,75,74,80,77,79,78,84,81,83,82,88,85,87,86,116,113,115,114,120,117,119,118,139,136,138,137,143,140,142,141,12,9,11,10,16,13,15,14,20,17,19,18,24,21,23,22,44,41,43,42,48,45,47,46,60,57,59,58,64,61,63,62,92,89,91,90,96,93,95,94,100,97,99,98,104,101,103,102,108,105,107,106,112,109,111,110,124,121,123,122,128,125,127,126,132,129,131,130,135,133,134,147,144,146,145,151,148,150,149)

#create 38 barplots and saving each of them as a PNG file
for (i in 1:38) {
  filename = paste("PC_",i,".png", sep = "")
  png(filename)
  barplot(pcabatchALL[samplesorder,i],col=colors[samplesorder],las=2,xaxt='n',space=0)
  for (i in 1:length(thinlines)) {
    abline(v = thinlines[i], col = "black",lty = 3)
  }
  for (i in 1:length(thicklines)) {
    abline(v = thicklines[i], col = "black",lty = 1)
  }
  abline(v = 72, col = "red",lty = 1)
  dev.off()
}

write.csv(aloadrelativebatch,file="aloadrelativeMask_batchmodel_filtered.csv")
write.csv(pcabatch$x,file="pca_batchmodel_x.csv")

```

To facilitate visualization, a red line is drawn at position 72 in order to separate the non-fibrotic group from the fibrotic group. Within each patient, the treatment order would be CC, TG, TN, TT. The color of each bar represents the batch of the sample, with a unique color assigned to each batch. The vertical lines on the plot separate the data for individual iPSC lines and their two batches.  

The Figure 9 below represents PC1 for all samples. From this plot, one can see that PC1 corresponds to extreme expression after treatment with TNF-a in all cases, even more that after its combination with TGF-b. Therefore, it seems to indicate that PC1 is associated with an interaction between TNF-a and TGF-b in iPSC mesenchymal organoids, an unexpected finding that warrants further analysis. 

![](/Pics/PCs/Figure_9.png)

<br>

#### PCA rank matrix

For easier visualization, I next clean the PCA results for exporting into spreadsheets. Take csv files and converts it to the txt files with the second column onwards. It does this by first removing the first row using awk, replacing multiple commas with tabs using tr, and removing the first column using cut.

```

#In terminal
cat aloadrelativeMask_batchmodel_filtered.csv | awk 'NR>1{print}' | tr -s "," "\t" | cut -f 2- > aloadrelativeMask_batchmodel_filtered.clean.txt
cat pca_batchmodel_x.csv | awk 'NR>1{print}' | tr -s "," "\t" | cut -f 2- > pca_batchmodel_x.clean.txt

```

Read in the preprocessed data files created in the previous steps and store them in variables pcabatch\_samples and pca\_loadings, respectively.

```

#In Matlab
pcabatch_samples = textread('pca_batchmodel_x.clean.txt','');
pca_loadings = textread('aloadrelativeMask_batchmodel_filtered.clean.txt','');

```

Sort the three columns of pca\_loadings in descending order and store the sorted values in variables x1, x2, and x3, and the corresponding indices in y1, y2, and y3.

```

%x = pca_loading number, y = its index

[x1 y1]=sort(pca_loadings(:,1),'descend');
[x2 y2]=sort(pca_loadings(:,2),'descend');
[x3 y3]=sort(pca_loadings(:,3),'descend');

```

Determine the rank of each row in the original order for the first three principal components and store the ranks in a matrix pcarankmatrix.

```

[x y z]=intersect(1:length(y1),y1);
pcarankmatrix(:,1)=z;
[x y z]=intersect(1:length(y2),y2);
pcarankmatrix(:,2)=z;
[x y z]=intersect(1:length(y3),y3);
pcarankmatrix(:,3)=z;

%contains the rank of each feature in the original order for the first three principal components

dlmwrite('pcarankmatrix.txt', pcarankmatrix,'delimiter','\t')

```

To efficiently manage our data with a single glance, I have organized it into an Excel spreadsheet using a combination of command line, Excel, and R.

<br>

#### Spreadsheet

The Figure 10 (patients) built on excel contains patients order, patients id, phenotypes, and sex. You can visit the sheet by clicking [here](/spreadsheet/Barret_Myofibroblast_TGFTNF_MASTER.xlsx).

![](/Pics/Figure_10.png)

<br>

The Figure 11 (allbiotypes_percents) includes the names and percentages of all biotypes, along with their respective minimum, maximum, and average values, providing us with a comprehensive overview. You can visit the sheet by clicking [here](/spreadsheet/Barret_Myofibroblast_TGFTNF_MASTER.xlsx).

```

#In terminal
paste allbiotypes allbiotypescountspercents > combine_allbiotypes_percents.txt

#In R
#allbiotypes_percents
sheet2_1 <- list("Biotypes")
sheet2_2 <- sampleKeyTNFATGFB
combined_sheet2 <- c(sheet2_1, sheet2_2)
combined_spreadsheet2 <- as.matrix(read.table("combine_allbiotypes_percents.txt"))
colnames(combined_spreadsheet2) <- combined_sheet2
write.table(combined_spreadsheet2,file="combined_spreadsheet2.txt", sep = "\t", row.names = FALSE)
#add their respective minimum, maximum, and average values on Excel

```
![](/Pics/Figure_11.png)

<br>

The Figure 12 is the main sheet that includes Genename, Geneid, Mapp, PC1, PC2, PC3, and  patient's TPM values.

```

#In terminal

paste Gencode_33_Selected_Genename_GMask.txt Gencode_33_Selected_Genename_GMask.txt Gencode_33_Selected_Geneid_GMask.txt Gencode_33_Selected_MappSS_GMask.txt pcarankmatrix.txt  > combine_test.txt

#In R

#spreadsheet

sheet3_1 <- list("Genename","Genename","Geneid","Mapp","PC1","PC2","PC3")
sheet3_2<- sampleKeyTNFATGFB
combined_headers <- c(sheet3_1, sheet3_2)
combined_spreadsheet <- as.matrix(read.table("combine_test.txt"))
colnames(combined_spreadsheet) <- combined_headers
combined_spreadsheet <- combined_spreadsheet[order(combined_spreadsheet[,1]),] #sort by the first column
write.table(combined_spreadsheet,file="combined_spreadsheet.txt", sep = "\t", row.names = FALSE)

```

You can sort this sheet with PC1, PC2, and so on to see the correlation between the experimental factors and the expression level of each gene, which allows a quick identification of genes and gene classes more associated with the dominant sources of gene expression variability in this experiment.

![](/Pics/Figure_12.png)

 <br> 
 
#### DESed analysis for pairwise comparisons
 
Reorder data by phenotypes and patients (lines)

```

samplesorder=c(4,1,3,2,8,5,7,6,28,25,27,26,32,29,31,30,36,33,35,34,40,37,39,38,52,49,51,50,56,53,55,54,68,65,67,66,72,69,71,70,76,73,75,74,80,77,79,78,84,81,83,82,88,85,87,86,116,113,115,114,120,117,119,118,139,136,138,137,143,140,142,141,12,9,11,10,16,13,15,14,20,17,19,18,24,21,23,22,44,41,43,42,48,45,47,46,60,57,59,58,64,61,63,62,92,89,91,90,96,93,95,94,100,97,99,98,104,101,103,102,108,105,107,106,112,109,111,110,124,121,123,122,128,125,127,126,132,129,131,130,135,133,134,147,144,146,145,151,148,150,149)

```

 <br> 

Add the PC1 coordinate of each sample into the experimental design, to be used as factor in the model

```

BarretMyofCnt=RBarretTNFATGFBCntGMask[,samplesorder]
sampleTableMyof=cbind(sampleTableTNFATGFB[samplesorder,],pcabatchALL[samplesorder,1:2])
sampleTableMyof$PC1 <- scale(sampleTableMyof$PC1 , center = TRUE) 

#normalized and batch-corrected expression matrix 
#using the variance stabilizing transformation
BarretMYO_Batch_vsd=RBarretTNFATGFBCntGMaskBatch_vsd[,samplesorder]

write.csv(assay(BarretMYO_Batch_vsd),file="BarretMYO_Batch_vsd.csv")

```

 <br> 

Check gene name duplicates

```

Genenames = as.list(read.table("Gencode_33_Selected_Genename_GMask.txt",header=FALSE,as.is=TRUE))

#Find duplicate gene names
duplicated_genes <- Genenames$V1[duplicated(Genenames$V1)]
if(length(duplicated_genes) > 0){
  cat("Duplicate gene names found:", paste(duplicated_genes, collapse = ", "))
} else {
  cat("No duplicate gene names found.")
}

#Duplicate gene names found: TBCE, ATXN7, AHRR, MATR3, HSPA14, TMSB15B


#Find duplicate gene names at indices
duplicated_indices <- which(duplicated(Genenames$V1))
if(length(duplicated_indices) > 0){
  cat("Duplicate gene names found at indices:", paste(duplicated_indices, collapse = ", "))
} else {
  cat("No duplicate gene names found.")
}
#Duplicate gene names found at indices: 1530, 3024, 4133, 4576, 7638, 15459

```

```

#Correct the duplicated gene names 
#by appending "_1" to the end of each duplicate gene name
Genenames$V1[1530] <- "TBCE_1"
Genenames$V1[3024] <- "ATXN7_1"
Genenames$V1[4133] <- "AHRR_1"
Genenames$V1[4576] <- "MATR3_1"
Genenames$V1[7638] <- "HSPA14_1"
Genenames$V1[15459] <- "TMSB15B_1"

```

 <br> 
 
**First round of pairwise comparisons:**

All cell lines UNTREATED vs All cell lines TGFβ

All cell lines UNTREATED vs All cell lines TNFα

All cell lines UNTREATED vs All cell lines TNFα/TGFβ

We model the data correcting for PC1 and Line (patient-specific expression) and then test for the treatment effect

```

#The factor "Batch" is included as a covariate to account for potential batch effects, 
#and the factor "Treatment" is included as the variable of interest for differential expression analysis.

RBarretMYOFCntGMaskTreatment <- DESeqDataSetFromMatrix(BarretMyofCnt, colData= sampleTableMyof, design= ~ Batch + Treatment)
RBarretMYOFCntGMaskTreatment <- DESeq(RBarretMYOFCntGMaskTreatment)

```

<br>

The filterFun argument specifies a multiple testing correction method to apply to the results, in this case the independent hypothesis weighting (IHW) method.

```
#Calculate differential expression results for the pairwise comparisons of 
#the TG vs CC, TN vs CC, and TT vs CC treatment groups, respectively. 

#The resulting output for each comparison will contain a table of genes 
#with their corresponding LFCs, p-values, and adjusted p-values based on the specified multiple testing correction method.

RBarretMYOFCntGMaskTreatment_TG <- results(RBarretMYOFCntGMaskTreatment,contrast=c("Treatment", "TG", "CC"),filterFun=ihw)
RBarretMYOFCntGMaskTreatment_TN <- results(RBarretMYOFCntGMaskTreatment,contrast=c("Treatment", "TN", "CC"),filterFun=ihw)
RBarretMYOFCntGMaskTreatment_TT <- results(RBarretMYOFCntGMaskTreatment,contrast=c("Treatment", "TT", "CC"),filterFun=ihw)

```

<br>

Visually check the names of the most significant genes (very low adjusted p-value)

```

Genenames$V1[which(RBarretMYOFCntGMaskTreatment_TG$padj<0.000000000000001)] 
Genenames$V1[which(RBarretMYOFCntGMaskTreatment_TN$padj<0.000000000000001)] 
Genenames$V1[which(RBarretMYOFCntGMaskTreatment_TT$padj<0.000000000000001)] 

```

<br>

**Second round of pairwise comparisons:**

Untreated NON-FIBROTIC cell lines vs Untreated FIBROTIC cell lines

TNFα/TGFβ NON-FIBROTIC cell lines vs TNFα/TGFβ FIBROTIC cell lines

TGFβ NON-FIBROTIC cell lines vs TGFβ FIBROTIC cell lines

TNFα NON-FIBROTIC cell lines vs TNFα FIBROTIC cell lines

We model the data correcting for Line (patient-specific expression) and then test for the Factor (treatment+phenotype combination) effect.

```

RBarretMYOFCntGMaskFactor <- DESeqDataSetFromMatrix(BarretMyofCnt, colData= sampleTableMyof, design= ~ Line + Factor)
RBarretMYOFCntGMaskFactor <- DESeq(RBarretMYOFCntGMaskFactor)

```

```
#Compute differential expression analysis results for the four contrasts of interest "CC_F" vs "CC_N", "TG_F" vs "TG_N", "TN_F" vs "TN_N", and "TT_F" vs "TT_N" in the dataset.

RBarretMYOFCntGMaskFactor_CC_Pheno <- results(RBarretMYOFCntGMaskFactor,contrast=c("Factor", "CC_F", "CC_N"),filterFun=ihw)
RBarretMYOFCntGMaskFactor_TG_Pheno <- results(RBarretMYOFCntGMaskFactor,contrast=c("Factor", "TG_F", "TG_N"),filterFun=ihw)
RBarretMYOFCntGMaskFactor_TN_Pheno <- results(RBarretMYOFCntGMaskFactor,contrast=c("Factor", "TN_F", "TN_N"),filterFun=ihw)
RBarretMYOFCntGMaskFactor_TT_Pheno <- results(RBarretMYOFCntGMaskFactor,contrast=c("Factor", "TT_F", "TT_N"),filterFun=ihw)

```

```
Genenames$V1[which(RBarretMYOFCntGMaskFactor_CC_Pheno$padj<0.01)] 
Genenames$V1[which(RBarretMYOFCntGMaskFactor_TG_Pheno$padj<0.01)] 
Genenames$V1[which(RBarretMYOFCntGMaskFactor_TN_Pheno$padj<0.01)] 
Genenames$V1[which(RBarretMYOFCntGMaskFactor_TT_Pheno$padj<0.01)] 

```

Export and process results

```

write.csv(as.data.frame(RBarretMYOFCntGMaskTreatment_TG),file="RBarretMYOFCntGMaskTreatment_TG.txt")
write.csv(as.data.frame(RBarretMYOFCntGMaskTreatment_TN),file="RBarretMYOFCntGMaskTreatment_TN.txt")
write.csv(as.data.frame(RBarretMYOFCntGMaskTreatment_TT),file="RBarretMYOFCntGMaskTreatment_TT.txt")
write.csv(as.data.frame(RBarretMYOFCntGMaskFactor_CC_Pheno),file="RBarretMYOFCntGMaskFactor_CC_Pheno.txt")
write.csv(as.data.frame(RBarretMYOFCntGMaskFactor_TG_Pheno),file="RBarretMYOFCntGMaskFactor_TG_Pheno.txt")
write.csv(as.data.frame(RBarretMYOFCntGMaskFactor_TN_Pheno),file="RBarretMYOFCntGMaskFactor_TN_Pheno.txt")
write.csv(as.data.frame(RBarretMYOFCntGMaskFactor_TT_Pheno),file="RBarretMYOFCntGMaskFactor_TT_Pheno.txt")

```

<br>

#### Pairwise results shreadsheet

Paste Gencode\_33\_Selected\_Geneid\_GMask.txt, Gencode\_33\_Selected\_Genename\_GMask.txt, Gencode\_33\_Selected\_MappSS\_GMask.txt, RBarretMYOFCntGMaskFactor\_CC\_Pheno\_test.txt, RBarretMYOFCntGMaskFactor\_TG\_Pheno\_test.txt, RBarretMYOFCntGMaskFactor\_TN\_Pheno\_test.txt, RBarretMYOFCntGMaskFactor\_TT\_Pheno\_test.txt, 
RBarretMYOFCntGMaskTreatment\_TG\_test.txt, RBarretMYOFCntGMaskTreatment\_TN\_test.txt, and
RBarretMYOFCntGMaskTreatment\_TT\_test.txt to formulate a shreadsheet.


```
#In terminal

tail -n +2 RBarretMYOFCntGMaskFactor_CC_Pheno.txt | cut -d"," -f3,7 > RBarretMYOFCntGMaskFactor_CC_Pheno_test.txt

tail -n +2 RBarretMYOFCntGMaskFactor_TG_Pheno.txt | cut -d"," -f3,7 > RBarretMYOFCntGMaskFactor_TG_Pheno_test.txt

tail -n +2 RBarretMYOFCntGMaskFactor_TN_Pheno.txt | cut -d"," -f3,7 > RBarretMYOFCntGMaskFactor_TN_Pheno_test.txt

tail -n +2 RBarretMYOFCntGMaskFactor_TT_Pheno.txt | cut -d"," -f3,7 > RBarretMYOFCntGMaskFactor_TT_Pheno_test.txt

tail -n +2 RBarretMYOFCntGMaskTreatment_TG.txt | cut -d"," -f3,7 > RBarretMYOFCntGMaskTreatment_TG_test.txt

tail -n +2 RBarretMYOFCntGMaskTreatment_TN.txt | cut -d"," -f3,7 > RBarretMYOFCntGMaskTreatment_TN_test.txt

tail -n +2 RBarretMYOFCntGMaskTreatment_TT.txt | cut -d"," -f3,7 > RBarretMYOFCntGMaskTreatment_TT_test.txt

paste Gencode_33_Selected_Geneid_GMask.txt Gencode_33_Selected_Genename_GMask.txt Gencode_33_Selected_MappSS_GMask.txt  RBarretMYOFCntGMaskFactor_CC_Pheno_test.txt RBarretMYOFCntGMaskFactor_TG_Pheno_test.txt RBarretMYOFCntGMaskFactor_TN_Pheno_test.txt RBarretMYOFCntGMaskFactor_TT_Pheno_test.txt RBarretMYOFCntGMaskTreatment_TG_test.txt RBarretMYOFCntGMaskTreatment_TN_test.txt RBarretMYOFCntGMaskTreatment_TT_test.txt > Barret_Myofibroblast_TGFTNF_PAIRWISEResults.txt

sed 's/,/\t/g' Barret_Myofibroblast_TGFTNF_PAIRWISEResults.txt > Barret_Myofibroblast_TGFTNF_PAIRWISEResults.txt

```

In Figure 13, it shows a clear formatted sheet containing Gene ID, Gene name, Mappability, fold change, and adjusted p-value. Then, you can sort each adjusted p-value column to focus on genes you are interested in. You can find the pairwise results spreadsheet [here](/spreadsheet/Barret_Myofibroblast_TGFTNF_PAIRWISEResults.xlsx).

![](/Pics/Figure_13.png)

<br>

#### Pathway analysis

To run a first test on the reliability of the differential expression results above, we retrieve the top most significant genes (by adjusted p-value) and run functional enrichment analysis using [Metascape](https://metascape.org/). An example of the enrichment analyses that can be used for validation are below, using the top 500 most significant genes by their overall response to the TNFa-TGFb combination.

The following figure shows the top functional terms ranked by significance (hypergeomtric p-value):
![](/Pics/Figure_16.png)

<br>

Two observations are noted: these results recapitulate the expected response to TNFa, but on the other hand too many categories are related with rather synonymous functional terms. The following figure shows an alternative representation, an enrichment network of the same significant pathways and functional terms. The network is arranged by the similarity between the genes contributing to each term, and colored according to functional groups with high overlap. A highly connected network indicates that most of the top significant pathways are highly redundant, and only a few of them are needed for further consideration.
![](/Pics/Figure_14.png)

<br>

The following shows the same enrichment network colored by statistical significance (hypergeomtric p-value):
![](/Pics/Figure_15.png)

<br>

Finally, we can screen how our differential results capture previously known protein-protein interactions. The figures below aggregate treatment-responsive genes in interaction modules. Many of the interactions are expected (e.g. groups of inflammation-related genes, modules of genes known to interact in the control of cell-cyle or NFKB-mediated regulation), and facilitate the validation of the experimental results and the identification of the most regulated genes in specific pathways.
![](/Pics/Figure_17.png)

<br>

## Future works 

To further our understanding of these complex signaling pathways, my future activities in this project will include:

* The identification of patient-specific responses.
* The identification of phenotype-specific responses, to better understand if our iPSC models can shed light into the predisposition of some patients to develop fibrotic complications.
* The integration of this experiment with a matched dataset (same iPSC lines, same treatments) obtained from iPSC lines differentiated into epithelial organoids (currently being analyzed). My goal is to integrate data from both the myofibroblast experiment and the epithelial experiment. By combining these datasets, I aim to identify patient-specific and phenotype-specific epithelial-mesenchymal interactions under inflammatory (TNFa) and pro-fibrotic (TGFb) conditions. These epithelial-mesenchymal interactions are known to be involved in the recurrence of fibrotic complication in some IBD patients. One possible way to do this is to merge both datasets to analyze differential changes in the levels of ligand-receptor interactions between both cell types. If successful, leveraging this combined dataset to model responses from iPSCs lines could provide a powerful tool for studying personalized interventions for these diseases.

<br> 

## References 

1.	D'Alessio, S., et al., Revisiting fibrosis in inflammatory bowel disease: the gut thickens. Nat Rev Gastroenterol Hepatol, 2022. 19(3): p. 169-184.
2.	Brooks, I.R., et al., Functional genomics and the future of iPSCs in disease modeling. Stem Cell Reports, 2022. 17(5): p. 1033-1047.
3.	Workman, M.J., et al., Modeling Intestinal Epithelial Response to Interferon-gamma in Induced Pluripotent Stem Cell-Derived Human Intestinal Organoids. Int J Mol Sci, 2020. 22(1).
4.	Ihara, S., Y. Hirata, and K. Koike, TGF-beta in inflammatory bowel disease: a key regulator of immune cells, epithelium, and the intestinal microbiota. J Gastroenterol, 2017. 52(7): p. 777-787.
5.	Wang, Q., et al., Applications of human organoids in the personalized treatment for digestive diseases. Signal Transduct Target Ther, 2022. 7(1): p. 336.
6.	Carcamo-Orive, I., et al., Analysis of Transcriptional Variability in a Large Human iPSC Library Reveals Genetic and Non-genetic Determinants of Heterogeneity. Cell Stem Cell, 2017. 20(4): p. 518-532 e9.
7.	Anders, S. and W. Huber, Differential expression analysis for sequence count data. Genome Biol, 2010. 11(10): p. R106.



<br>
