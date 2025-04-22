## The pipelines ddradseq_pre.bash and ddradseq_dedup.bash were written to help preprocessing reads from ddRADseq experiments using sequencing of double digested of genomic DNA inserted sorrounded by inline barcode on the P5 adaptor read and DBR region on the P7 adaptor read, the procedure modified from Schweyen et al.
## The new solution gave finally up to 10x more Stacks with Coverage 2-6 times higher than when preprocessing reads with the original proposed approach with STACKS's process_radtags and clone_filter. Comparison results in ComparisonFinal.ods (https://github.com/rafalwoycicki/ddRADseq_reads/blob/main/ComparisonsFinal.ods) file and in the attached graphs.

Legend:
| cpa | after cut site filtering paired reads |
| rnp | not paired remained reads after cut site filetring |
| dpa | paired reads after deduplication based on DBR region |
| p5 | reads from p5 adaptor site (forward reads) |
| p7 | reads from p7 adaptor site (reverse reads) |


| Old dataset 130nt All   |         | p5_cpa | p7_cpa | p5_rnp | p7_rnp | p5_dpa | p7_dpa |
|-------------------------|---------|--------|--------|--------|--------|--------|--------|
| process_radtags eb2 ea2 | Average | 1,11%  | 1,11%  | 7,86%  | 0,23%  | 0,49%  | 0,49%  |
| ddradseq_reads eb1      | Average | 7,11%  | 7,11%  | 0,46%  | 0,22%  | 2,56%  | 2,56%  |

<p float="left">
  <img src="https://github.com/user-attachments/assets/c378ddb5-baf9-48f3-b73c-34e151a4e96f" width="500" />
  <img src="https://github.com/user-attachments/assets/3fda75cd-02e7-4701-8fba-a45811de448e" width="500" />
</p>


| Old dataset 130nt All   |         | Stacks     | Mean_cov |
|-------------------------|---------|------------|----------|
| process_radtags eb2 ea2 | Average | 120675,625 | 10,51    |
| ddradseq_reads eb1      | Average | 131449,75  | 62,66    |

<p float="left">
  <img src="https://github.com/user-attachments/assets/44daabc6-cb9c-4ffd-a245-0930cfe2e6c6" width="500" />
  <img src="https://github.com/user-attachments/assets/350d5cfe-671f-420c-add0-7e5ffa9e2486" width="500" />
</p>


| New P1 130 130nt All    |         | p5_cpa | p7_cpa | p5_rnp | p7_rnp | p5_dpa | p7_dpa |
|-------------------------|---------|--------|--------|--------|--------|--------|--------|
| process_radtags eb2 ea2 | Average | 0,06%  | 0,06%  | 0,56%  | 0,10%  | 0,04%  | 0,04%  |
| ddradseq_reads eb1      | Average | 0,44%  | 0,44%  | 0,04%  | 0,01%  | 0,33%  | 0,33%  |

<p float="left">
  <img src="https://github.com/user-attachments/assets/1a854d40-429f-45de-84f6-bd224324e1d4" width="500" />
  <img src="https://github.com/user-attachments/assets/75c08b55-aed1-46f5-9f01-22f50f90a8aa" width="500" />
</p>


| New P1 130 130nt All    |         | Stacks   | Mean_cov |
|-------------------------|---------|----------|----------|
| process_radtags eb2 ea2 | Average | 3886,01  | 6,06     |
| ddradseq_reads eb1      | Average | 58417,83 | 12,21    |

<p float="left">
  <img src="https://github.com/user-attachments/assets/6ffa0d82-ced1-4828-8623-61c06c59887d" width="500" />
  <img src="https://github.com/user-attachments/assets/33e3c8fc-b28a-4496-8980-732a400035d6" width="500" />
</p>

| New P2 130 130nt All    |         | p5_cpa | p7_cpa | p5_rnp | p7_rnp | p5_dpa | p7_dpa |
|-------------------------|---------|--------|--------|--------|--------|--------|--------|
| process_radtags eb2 ea2 | Average | 0,06%  | 0,06%  | 0,52%  | 0,09%  | 0,04%  | 0,04%  |
| ddradseq_reads eb1      | Average | 0,40%  | 0,40%  | 0,05%  | 0,14%  | 0,27%  | 0,27%  |

<p float="left">
  <img src="https://github.com/user-attachments/assets/2b2443bb-125b-40e1-a0ba-be87cbb27dd1" width="500" />
  <img src="https://github.com/user-attachments/assets/5737cae1-1bd5-4afa-ba4b-64376d067fac" width="500" />
</p>


| New P2 130 130nt All    |         | Stacks   | Mean_cov |
|-------------------------|---------|----------|----------|
| process_radtags eb2 ea2 | Average | 6994     | 6,07     |
| ddradseq_reads eb1      | Average | 55495,21 | 12,35    |

<p float="left">
  <img src="https://github.com/user-attachments/assets/b970e291-baf7-4878-9e9e-c1f42ecf311e" width="500" />
  <img src="https://github.com/user-attachments/assets/bb173278-40c1-4538-90b5-cf5557a69f52" width="500" />
</p>
