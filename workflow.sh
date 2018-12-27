bedtools makewindows -g Bdi.chromSize -w 1000000 > Bdi_windows_1000000.bed
bedtools multicov -bams ChIP-Seq/bam/Input.bam ChIP-Seq/bam/A_H3K9me3.bam ChIP-Seq/bam/B_H3K9me3.bam ChIP-Seq/bam/C_H3K9me3.bam ChIP-Seq/bam/D_H3K9me3.bam ChIP-Seq/bam/E_H3K9me3.bam -bed Bdi_windows_1000000.bed > input_H3K9me3_cov_1000000.tsv &
bedtools multicov -bams ChIP-Seq/bam/Input.bam ChIP-Seq/bam/A_H3K4me3.bam ChIP-Seq/bam/B_H3K4me3.bam ChIP-Seq/bam/C_H3K4me3.bam ChIP-Seq/bam/D_H3K4me3.bam ChIP-Seq/bam/E_H3K4me3.bam -bed Bdi_windows_1000000.bed > input_H3K4me3_cov_1000000.tsv &
bedtools multicov -bams MeDIP/bam/Input.bam MeDIP/bam/A.bam MeDIP/bam/B.bam MeDIP/bam/C.bam MeDIP/bam/D.bam MeDIP/bam/E.bam -bed Bdi_windows_1000000.bed > input_MeDIP_cov_1000000.tsv &

$HOME/R/3.5.0/bin/Rscript circos_plot.R circos/input_H3K9me3_cov_1000000.tsv circos/input_H3K9me3_cov_1000000.pdf
$HOME/R/3.5.0/bin/Rscript circos_plot.R circos/input_H3K4me3_cov_1000000.tsv circos/input_H3K4me3_cov_1000000.pdf
$HOME/R/3.5.0/bin/Rscript circos_plot.R circos/input_MeDIP_cov_1000000.tsv circos/input_MeDIP_cov_1000000.pdf
