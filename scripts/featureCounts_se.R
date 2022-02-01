library(Rsubread)

counts <-featureCounts(files="outputs/gs_rnaseq_star/GCF_010509575.1-ERX4307280Aligned.sortedByCoord.out.bam",
                       
                       # annotation
                       annot.inbuilt=NULL,
                       annot.ext="outputs/gs_genomes_bakta/GCF_010509575.1/GCF_010509575.1.gtf",
                       isGTFAnnotationFile=T,
                       GTF.featureType="CDS",
                       GTF.attrType = "ID",
                       
                       # level of summarization
                       useMetaFeatures=FALSE,
                       
                       # overlap between reads and features
                       allowMultiOverlap=TRUE,
                       minOverlap=31,
                       largestOverlap=FALSE,
                       readExtension5=0,
                       readExtension3=0,
                       read2pos=NULL,
                       
                       # multi-mapping reads
                       countMultiMappingReads=T,
                       fraction=T,
                       
                       # read filtering
                       minMQS=0,
                       splitOnly=FALSE,
                       nonSplitOnly=FALSE,
                       primaryOnly=FALSE,
                       ignoreDup=FALSE,
                       
                       # strandness
                       strandSpecific=0,
                       # exon-exon junctions
                       juncCounts=TRUE,
                       genome=NULL,
                       
                       # parameters specific to paired end reads
                       isPairedEnd=FALSE,
                       requireBothEndsMapped=FALSE,
                       checkFragLength=FALSE,
                       minFragLength=20,
                       maxFragLength=600,
                       countChimericFragments=FALSE,
                       autosort=TRUE,
                       
                       # miscellaneous
                       nthreads=1,
                       maxMOp=10)

write.table(x=data.frame(counts$annotation, counts$counts, stringsAsFactors=FALSE),
            file=snakemake@output[["counts"]], row.names=F, quote = F)

write.table(x=counts$counts_junction, files=snakemake@output[['juncs']], row.names = F, quote = F)
