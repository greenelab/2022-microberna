# featureCounts is a function in the Rsubread package used to estimate the number of
# reads that overlap with a given feature in a genome from a BAM file of reads that
# have been mapped to that reference. The featureCounts() function is highly 
# parameterizable, which gives it a lot of flexibility for different sequencing 
# technologies, but it also makes it difficult to know whether the exact set of 
# correct parameters is always being used. 
# As invoked here, featureCounts:
#     1. requires a minimum overlap of 31 base pairs between a read and a feature.
#        This was somewhat arbitrary, but I selected it because it matches the 
#        default minimum k-mer overlap for a read and a transcript in salmon.
#     2. Counts at the CDS level.
#     3. Allows reads to overlap multiple features. This allows "spanning" reads 
#        that likely correspond to reads that overlap two genes, as occurs in 
#        polycistronic transcripts from multi-gene operons.
#     4. Counts reads fractionally if they overlap multiple features. 
#     5. Ignores strandedness
#     6. Counts as SE (reads were aligned to the reference as SE regardless of 
#        whether they were SE or PE.

library(Rsubread)

counts <-featureCounts(files=snakemake@input[['bam']],
                       
                       # annotation
                       annot.inbuilt=NULL,
                       annot.ext=snakemake@input[['gtf']],
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

write.table(x=counts$counts_junction, file=snakemake@output[['juncs']], row.names = F, quote = F)
