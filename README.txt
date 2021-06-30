emeraldBGC Repository
========================

------------------------


emeraldBGC - SMBGC detection tool -

How to install emeraldBGC?

  requires: 
    Linux OS
    Bioconda : https://bioconda.github.io/user/install.html

  installation:
    conda create -n emeraldbgc emeraldbgc -c santiagosanchezf

  basic use:
    conda activate emeraldbgc
    emeraldbgc <nucleotide fasta file>
    conda deactivate emerald

  output:
    GFF3 format file

    The fields in this header are as follows:

      seqname: SeqID for sequence, as in prodigal output.
      source: emeraldbgc version.
      feature: Feature type name, i.e. CLUSTER, CLUSTER_border, CDS.
      start: Start position of feature
      end: End position of feature
      score: empty
      strand: empty
      frame: empty
      attributes:
        ID: ordinal ID for the cluster, beginning with 1.
        nearest_MiBIG: MiBIG accession of the nearest BGC to the cluster in the MIBIG space, measured in Jaccard distance.
        nearest_MiBIG_class: BGC class of nearest_MiBIG.
        nearest_MiBIG_jaccardDistance: Jaccard distance between ID and nearest_MiBIG.
        partial: Indicates if a CLUSTER is at the edge of the contig. as in prodigal partial. "0" shows the cluster is not at the edge, whereas a "1" indicates is at that edge, (i.e. a partial cluster).

    sample:
     ##gff-version 3
     DS999642	EMERALDv0.9.0	CLUSTER	1	136970	.	.	.	ID=DS999642_emrld_1;nearest_MiBIG=BGC0001397;nearest_MiBIG_class=NRP Polyketide;nearest_MiBIG_jaccardDistance=0.561;partial=10
