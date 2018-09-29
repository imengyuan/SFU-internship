
## script_notes.txt

```
7.18 Quantification
running in /home/meng/quant

meng@Meng[quant] ~/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/util/align_and_estimate_abundance.pl --transcripts ../Desktop/MengFiles/Transcriptome/WRCTranscriptome69K.fasta --seqType fq --samples_file sample.txt --est_method RSEM --aln_method bowtie2 --thread_count 8 --trinity_mode --prep_reference

ERR message

(ERR): bowtie2-align exited with value 141
Error, cmd: set -o pipefail && bowtie2 --no-mixed --no-discordant --gbar 1000 --end-to-end -k 200  -q -X 800 -x /home/meng/quant/../Desktop/MengFiles/Transcriptome/WRCTranscriptome69K.fasta.bowtie2 -1 /home/meng/Desktop/MengFiles/TrimmedReads/UC1-Low-F1/UC1-Low-F1_1P.fastq -2 /home/meng/Desktop/MengFiles/TrimmedReads/UC1-Low-F1/UC1-Low-F1_2P.fastq -p 8 | samtools view -F 4 -S -b | samtools sort -n -o bowtie2.bam  died with ret: 256 at /home/meng/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/util/align_and_estimate_abundance.pl line 790.


meng@Meng[quant] ~/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/util/align_and_estimate_abundance.pl --transcripts ../Desktop/MengFiles/Transcriptome/WRCTranscriptome69K.fasta --seqType fq --samples_file sample.txt --est_method RSEM --aln_method bowtie --thread_count 8 --trinity_mode --prep_reference

Error while writing string output; 2889082 characters in string, 0 written
Command: bowtie --wrapper basic-0 -q --all --best --strata -m 300 --chunkmbs 512 -X 800 -S -p 8 -1 /home/meng/Desktop/MengFiles/TrimmedReads/UC1-Low-F1/UC1-Low-F1_1P.fastq -2 /home/meng/Desktop/MengFiles/TrimmedReads/UC1-Low-F1/UC1-Low-F1_2P.fastq /home/meng/quant/../Desktop/MengFiles/Transcriptome/WRCTranscriptome69K.fasta.bowtie 
Error, cmd: set -o pipefail && bowtie -q --all --best --strata -m 300 --chunkmbs 512 -X 800 -S -p 8 /home/meng/quant/../Desktop/MengFiles/Transcriptome/WRCTranscriptome69K.fasta.bowtie -1 /home/meng/Desktop/MengFiles/TrimmedReads/UC1-Low-F1/UC1-Low-F1_1P.fastq -2 /home/meng/Desktop/MengFiles/TrimmedReads/UC1-Low-F1/UC1-Low-F1_2P.fastq | samtools view -F 4 -S -b | samtools sort -n -o bowtie.bam  died with ret: 256 at /home/meng/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/util/align_and_estimate_abundance.pl line 790.


7.19
Quantification ok(trouble shooting:imcompatible samtools version)
now move onto differential gene analysis


7.20 DE Analysis

# generate matrix
## not normalized(can be replaced by the one below)
/home/meng/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/util/abundance_estimates_to_matrix.pl --est_method RSEM --quant_files quant_files.txt --gene_trans_map none --out_prefix genes --name_sample_by_basedir --cross_sample_norm none

## default TMM normailized(include same output as above)
/home/meng/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/util/abundance_estimates_to_matrix.pl --est_method RSEM --quant_files quant_files.txt --gene_trans_map none --out_prefix genes_TMM --name_sample_by_basedir 

# run DE analysis
/home/meng/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix rsem.isoform.counts.matrix  --method edgeR --samples_file sample_DE.txt

# clustering
/home/meng/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../rsem_TMM.isoform.TMM.EXPR.matrix -P 1e-3 -C 2 --samples ../sample_DE.txt


/home/meng/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../rsem_TMM.isoform.TMM.EXPR.matrix  -C 4 --samples ../sample_DE.txt

## GO
/home/meng/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../rsem_TMM.isoform.TMM.EXPR.matrix --samples ../sample_DE.txt --examine_GO_enrichment --GO_annots WRC66KAnnotation.xls --gene_lengths 100

/home/meng/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl -RdiffExpr.P0.001_C2.matrix.RData --Ptree 60


# re run DE using transcripts
/home/meng/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../trans.isoform.TMM.EXPR.matrix --samples ../sample_DE.txt --examine_GO_enrichment --GO_annots ../Annotation/go_annotationsTranscripts.txt --gene_lengths

7.26
/home/meng/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/DifferentialExpression/PtR --matrix ../trans.isoform.TMM.EXPR.matrix --samples ../sample_DE.txt --compare_replicates 

/home/meng/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/DifferentialExpression/PtR --matrix ../trans.isoform.TMM.EXPR.matrix --samples ../sample_DE.txt --sample_cor_matrix


7.27

# requires the transdecoder.pep file
tmhmm --short < transdecoder.pep > tmhmm.out


# auto
/home/meng/Downloads/Trinotate-Trinotate-v3.1.1/auto/autoTrinotate.pl --Trinotate_sqlite F1ToF4WRCTrinotate.sqlite --transcripts Transcriptome/WRCTranscriptome69K.fasta --gene_to_trans_map WRCTranscriptome69K.fasta.gene_trans_map --conf /home/meng/Downloads/Trinotate-Trinotate-v3.1.1/auto/conf.txt --CPU 6   
```

## script.txt

```
#PtR
grep GO:0009820 trans.isoform.counts.matrix.F1_vs_F3.edgeR.DE_results.P0.001_C2.F1-UP.subset.GOseq.enriched > GO0009820F1UpC2F1vsF3.txt

cut -f 10 GO0009820F1UpC2F1vsF3.txt > GO0009820F1UpC2F1vsF3Transcripts.txt

sed 's/, /\n/g' GO0009820F1UpC2F1vsF3Transcripts.txt > GO0009820F1UpC2F1vsF3Transcriptslist.txt

grep -wFf GO0009820F1UpC2F1vsF3Transcriptslist.txt ../trans.isoform.TMM.EXPR.matrix > GO0009820F1UpC2F1vsF3.TMM.EXPR.matrix

gedit GO0009820F1UpC2F1vsF3.TMM.EXPR.matrix

/home/meng/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/DifferentialExpression/PtR --matrix ../../quant/edgeR.11730.dir/GO0009820F1UpC2F1vsF3.TMM.EXPR.matrix --samples ../sample_DE.txt --heatmap --log2 --center_rows --order_columns_by_samples_file

# Goeast
 sed 's!,! // !g' go_annotationsTranscripts.txt > go_annotationsTranscriptsGOeast.txt


#9.4
grep -c 
sort -u
 ~/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/DifferentialExpression/run_GOseq.pl --genes_single_factor ../venn/trans.F1_vs_F234.C2.F1-UP.list --lengths ../Annotation/WRCTranscriptome66K.fasta.seq_lengths.txt --background trans.all.C2.DE.sorted.list --GO_assignments ../Annotation/go_annotationsTranscripts.txt


#9.10

```

## heatmap_script.txt

```
grep PF00067 DE_results.P0.001_C2.F1-UP.subset.Annotation > PF00067.P0.001_C2.F1-UP.subset.Annotation

cut -f 2 PF00067.P0.001_C2.F1-UP.subset.Annotation> PF00067.P0.001_C2.F1-UP.transcript.list 

grep -wFf PF00067.P0.001_C2.F1-UP.transcript.list trans.isoform.TMM.EXPR.matrix > PF00067.P0.001_C2.F1-UP.TMM.EXPR.matrix


gedit PF00067.P0.001_C2.F1-UP.TMM.EXPR.matrix

~/Downloads/trinityrnaseq-Trinity-v2.7.0-PRERELEASE/Analysis/DifferentialExpression/PtR --matrix  PF00067.P0.001_C2.F1-UP.TMM.EXPR.matrix --samples  sample_DE.txt --heatmap --log2 --center_rows --order_columns_by_samples_file 
```
