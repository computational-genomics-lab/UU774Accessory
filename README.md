# UU774Accessory
commandlines and small scripts for analyzing UU774 data

RNAseq Data Normalization:

# Generate FeatureCount files using the following commands. Input types are .bam files (alignment files) and the annotation file in gtf format.

featureCounts -a Genbank.gtf -o genbankFeatureAll PATH/*.bam  PATH/74_hs_6h*.bam  PATH/*.bam  PATH/*.bam

# Now since multiple comparisons are not possible, one can create metadata files with pairs of conditions and also the feature count file should have matching number of columns. 
Example:  cut -f 1,13,14,15,16 genbankFeatureAll > nmData
# Read the featureCount file
>countData <-    read.table("/home/sutripa/Mastigocladus_laminosus_74_data/nmData",header=TRUE, sep="\t")

>head(countData)
        Geneid nitroNMd12A nitroNMd12B nitroNPd12A nitroNPd12B
1 BLD44_000005         134          94          77          60
2 BLD44_000010         182         214          43          55
3 BLD44_000015          12          30          14          14
4 BLD44_000020           8          16           6          17
5 BLD44_000025         109          94          71          82
6 BLD44_000030          30          29          20          27

# Create metadata file nmMetadata
id      dex     cell
nitroNMd12A     nmControl       nitroNMd12A
nitroNMd12B     nmControl       nitroNMd12B
nitroNPd12A     nmTreated       nitroNPd12A
nitroNPd12B     nmTreated       nitroNPd12B
>nmmetadata <- read.table("/home/sutripa/Mastigocladus_laminosus_74_data/nmMetadata", header=TRUE, sep="\t")
# Now read and Construct DESEQDataSet Object
>dds <- DESeqDataSetFromMatrix(countData=countData, colData=nmmetadata, design=~dex, tidy=TRUE)
# Now we’re ready to run DESEQ function
> dds <- DESeq(dds)
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing

#Creating an all metadata file for all the conditions:
One can also create an allmetadata file with the following contents: It is absolutely essential that the rows of metadata file are in the same order as the column of the count file!!
id      dex     cell
nitroNMd12A     nmControl       nitroNMd12A
nitroNMd12B     nmControl       nitroNMd12B
nitroNPd12A     nmTreated       nitroNPd12A
nitroNPd12B     nmTreated       nitroNPd12B
hsCtrlA hscontrol       hsCtrlA
hsCtrlB hscontrol       hsCtrlB
hs6HA   6hTreated       hs6HA
hs6HB   6hTreate        hs6HB
hs24HA  24hTreated      hs24HA
hs24HA  24hTreated      hs24HA
tpNMd0A tpcontrol       tpNMd0A
tpNMd0B tpcontrol       tpNMd0B
tpNMd12A        12dTreated      tpNMd12A
tpNMd12B        12dTreated      tpNMd12A
And also read the entire readCount file that has same number of columns as the number of rows + 1 of allmetadata file:
> head(countData)
        Geneid hs24HA hs24HA.1 hsCtrlA hsCtrlB hs6HA hs6HB nitroNMd12A
1 BLD44_000005   3563     3551     596    2169  2710  3570         134
2 BLD44_000010   3553     3549     546    1816  2361  3452         182
3 BLD44_000015   4011     4038     662    2269  2955  4146          12
4 BLD44_000020   3576     3702     570    2138  2603  3518           8
5 BLD44_000025   2226     2497     362    1395  1624  2254         109
6 BLD44_000030    530      536     144     220   429   652          30
  nitroNMd12B nitroNPd12A nitroNPd12B tpNMd0A tpNMd0B tpNMd12A tpNMd12B
1          94          77          60    1175    1134     1336      934
2         214          43          55    4986    5188     5651     5377
3          30          14          14    2552    2465     2159     2109
4          16           6          17     855     788      864      625
5          94          71          82   24048   23114    16188    11336
6          29          20          27    1897    2001     1359      804
And Create a DESEQ dataset as described above.
[dds <- DESeqDataSetFromMatrix(countData=countData, colData=allmetadata, design=~dex, tidy=TRUE)]
Normalization (Tutorial from : https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html )
To perform the median of ratios method of normalization, DESeq2 has a single estimateSizeFactors() function that will generate size factors for us. We will use the function in the example below, but in a typical RNA-seq analysis this step is automatically performed by the DESeq() function, which we will see later.
dds<- estimateSizeFactors(dds)
By assigning the results back to the dds object we are filling in the slots of the DESeqDataSet object with the appropriate information. We can take a look at the normalization factor applied to each sample using:
> sizeFactors(dds)
     hs24HA    hs24HA.1     hsCtrlA     hsCtrlB       hs6HA       hs6HB
 6.56384808  6.70724313  1.04587719  3.82713123  4.75508467  6.55832785
nitroNMd12A nitroNMd12B nitroNPd12A nitroNPd12B     tpNMd0A     tpNMd0B
 0.08165984  0.09611280  0.05112812  0.05310010  1.84543173  1.97906435
   tpNMd12A    tpNMd12B
 2.10469868  1.50487865

Now, to retrieve the normalized counts matrix from dds, we use the counts() function and add the argument normalized=TRUE.
normalized_counts <- counts(dds, normalized=TRUE)

We can save this normalized data matrix to file for later use:

write.table(normalized_counts, file="/home/sutripa//Mastigocladus_laminosus_74_data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)
# Now lets take a look at the results

Commands for Normalizing count  data files from count data with scaffold_20 removed
> countData <- read.table("/home/sutripa/Mastigocladus_laminosus_74_data/genbankFeatureCountForDeseq", header=TRUE, sep="\t")> metadata <- read.table("/home/sutripa/Mastigocladus_laminosus_74_data/allMetadata",header=TRUE, sep="\t")
> dds <- DESeqDataSetFromMatrix(countData=countData, colData=metadata, design=~dex, tidy=TRUE)
> dds <- DESeq(dds)
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
> dds<- estimateSizeFactors(dds)
> normalized_counts <- counts(dds, normalized=TRUE)
> write.table(normalized_counts, file="/home/sutripa//Mastigocladus_laminosus_74_data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

Run Differential Expression Analysis:
> dds <- DESeq(dds)
using pre-existing size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing

> res <- results(dds)
> head(results(dds, tidy=TRUE))
           row  baseMean log2FoldChange     lfcSE       stat       pvalue
1 BLD44_000005  788.7977    -0.05435005 0.1085615 -0.5006382 6.166257e-01
2 BLD44_000010 1500.2293    -0.23291719 0.1628016 -1.4306815 1.525215e-01
3 BLD44_000015  696.0551     0.11600279 0.1736073  0.6681906 5.040119e-01
4 BLD44_000020  405.2468     0.06059023 0.1235176  0.4905394 6.237522e-01
5 BLD44_000025 3377.6819     0.69869211 0.1046162  6.6786253 2.411940e-11
6 BLD44_000030  380.9431     0.78874713 0.3396410  2.3222969 2.021695e-02
          padj
1 7.598537e-01
2 2.864322e-01
3 6.689222e-01
4 7.654435e-01
5 4.426884e-10
6 6.133583e-02
> summary(res)

out of 5910 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 1194, 20%
LFC < 0 (down)     : 1031, 17%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

Sort summary list by p-value
> res <- res[order(res$padj),]
> head(res)
log2 fold change (MLE): dex tpcontrol vs 12dTreated
Wald test p-value: dex tpcontrol vs 12dTreated
DataFrame with 6 rows and 6 columns
                     baseMean    log2FoldChange             lfcSE
                    <numeric>         <numeric>         <numeric>
BLD44_013210 1540.36821026245 -3.63198667791499 0.126735378567501
BLD44_030075 515.606449920767 -4.49447163253932 0.161430410267406
BLD44_006070 753.790190093366  3.07822078897396 0.122215459660653
BLD44_013160 617.153440655524 -3.83503630129863 0.160166604037641
BLD44_028555 7525.72136751065 -6.40094585261329 0.269100866418636
BLD44_027395 359.216856196291 -3.32313049761085 0.161059226577685
                          stat                pvalue                  padj
                     <numeric>             <numeric>             <numeric>
BLD44_013210 -28.6580331314554 1.27308894707933e-180 7.52395567723881e-177
BLD44_030075 -27.8415425265557 1.36358075883668e-170 4.02938114236238e-167
BLD44_006070  25.1868364077756 5.58377946880986e-140 1.10000455535554e-136
BLD44_013160 -23.9440445425024  1.0659213189821e-126 1.57489874879606e-123
BLD44_028555 -23.7864185938942 4.61640402653229e-125 5.45658955936117e-122
BLD44_027395 -20.6329719086785  1.38848995709581e-94  1.36766260773937e-91


See References: https://lashlock.github.io/compbio/R_presentation.html for plotting and continuation
https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html Detailed vignettes
https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html for edgeR




Running CheckM


checkm  lineage_wf  checkmBin  checkmout  # CheckmBin is the directory that contains the genome fasta file ( A single file or multiple files can be accommodated). Output will be stored in checkmout directory.

Running Pyani for tetra nucleotide abundance and AAI and ANI comparisons.

average_nucleotide_identity.py -i ./separateScaffolds -o Mastigo -l result74.txt -m ANIm -g --gformat png,pdf --write_excel
average_nucleotide_identity.py -i ./separateScaffolds -o MastigoTetra -l result74.txt -m TETRA -g --gformat png,pdf --write_excel # Here the separate scaffolds directory contains scaffolds as separate files.


Genbank Parser File for 5'UTR and 3'UTR and gene length:

# This program uses bioperl to parse the genbank file into pieces
# Author: Sucheta

#!/usr/bin/perl
#use strict;
use Bio::SeqIO;


#Name of the input genbank file
my $gbkfile = shift;


if ($gbkfile eq '') {
  print "Please enter the input filename\nUsage: perl <prog_name> filename\n\n";
  exit(1);
}
my ($seqio,$seq_object) = ();

eval {

 $seqio = Bio::SeqIO->new(-file => "$gbkfile");
};

if ($@ ne '') {
  print "ERROR: in opening the genbak file".$@;
  exit(1);
}



print "Scaffold_id\tScaffold_len\tGeneId\tUTR5\tUTR3\tGeneLength\tAnnotation\n";
#Reading the seq_object from the genbank file
#we are interested in locusid, startpos, endpos, length_protein, protein_sequence, name,  product
# This loops over all the sequences(scaffolds in my case)
while($seq_object = $seqio -> next_seq){
        my @utr5;
        my @utr3;
        my @cdsLen;
        my $firstPass = 1;
        my @gene;
        my @annotation;
        my $end;
        my $len;
        my $start;

        # This returns feature objects for each sequence
        foreach my $feat_object($seq_object->get_SeqFeatures) {
        #the $feat_object->primary_tag returns the values stored in the
        # left hand most side for each feature

                if($feat_object->primary_tag eq 'CDS') {
                        foreach my $tag ($feat_object->get_all_tags) {
                        chomp($tag);
                        if($tag eq 'protein_id'){
                        push(@gene,$feat_object->get_tag_values("protein_id"));
                        push(@annotation, $feat_object->get_tag_values("product"));
                        if ($firstPass){
                        $start = $feat_object->location->start;
                        $end = $feat_object->location->end;

                        push (@utr5,$start);
                        $firstPass = 0;

                        $len = $end - $start;
                        push (@cdsLen, $len);


                        }
                        else{
                        $start = $feat_object->location->start;
                        push (@utr5, $start - $end);
                        push (@utr3, $start - $end);
                        $end = $feat_object->location->end;
                        $len = $end - $start;
                        push (@cdsLen, $len);

                        }
                        }
                        break;
                }



           }
        }
        push (@utr3, ($seq_object->length - $end));
my $i=0;
#print "length of annotation, protein,  utr5, utr3, length is", scalar(@annotation),"\t",  scalar(@gene),"\t", scalar(@utr5), "\t", scalar(@utr3),"\t", scalar(@cdsLen), "\n";
foreach my $gene_id(@gene){
        chomp($annotation[$i]);
        print $seq_object->id,"\t", length($seq_object->seq),"\t";
        print "$gene_id\t$utr5[$i]\t$utr3[$i]\t$cdsLen[$i]\t$annotation[$i]\n";
        $i++;
}

}


R scripts for Boxplot Plotting:

myDat <- read.table("/home/ada/TOOLBOX/perlScripts/genbank/74InterGenic", sep="\t", header=TRUE)

>png(file="/home/ada/TOOLBOX/perlScripts/genbank/74Boxplot.png")
> boxplot(X5UTR ~ GeneLength, data=myDat, xlab="Gene Length", ylab="5'UTR length", main="Boxplot for UU774 5ÚTR length per gene length")
> dev.off()


For running Synteny analysis using MCScanX:

From gff to file readable by MCScanx:

$awk ‘{print $1”\t”$8”\t”$3”\t”$4}’ inputFile > outputFile (The format is chromosome\tgene\tlocationsStart\tLocationend\t)

Run MCscanX:

./MCScanX dir/all ( output files will be all.collinearity and all.html directory)

Run Dual Synteny plot:
======================
java dual_synteny_plotter -g /home/sutripa/MCscanX/MCScanX/cyano/all.gff -s /home/sutripa/MCscanX/MCScanX/cyano/all.collinearity -c dualsynteny1.ctl -o dualall1.png

Where dualsynteny1.ctl file contains the following:

450     //plot width (in pixels)
1500    //plot height (in pixels)
MN14,MN15,MN16  //chromosomes in x axis
JH98    //chromosomes in y axis

Run Ka/Ks analysis using MCScanX:

perl add_ka_and_ks_to_collinearity.pl  -i ../cyano/all.collinearity -d ../cyano/all.cds -o all.collinearity.kaks

Adding Ka/ks ratio to the last column:


#!/usr/bin/perl -w

open FH, $ARGV[0] or die "can't open file for readin $!\n";

while(<FH>){

        chomp;
        my @arr = split(/\t/,$_);

        if(scalar(@arr) eq 6){
                if($arr[5] eq 0){
                        print $_, "\tNA", "\n";
                }
                else{
                        print $_, "\t", $arr[4]/$arr[5],"\n";
                }
        }
        else{
                print $_, "\n";
        }

}

close(FH);

