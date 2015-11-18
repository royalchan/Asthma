########################
# This script is to perform ID conversion with the biomaRt library.
########################
library(biomaRt)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
nimblegenid <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_RefSeqIDs.txt")
refseqmrna <- getBM(filters="refseq_mrna", attributes= c("refseq_mrna", "hgnc_symbol", "entrezgene"), values=nimblegenid, mart=mart, uniqueRows = TRUE)
refseqmrnapredicted <- getBM(filters="refseq_mrna_predicted", attributes= c("refseq_mrna_predicted", "hgnc_symbol", "entrezgene"), values=nimblegenid, mart=mart, uniqueRows = TRUE)
refseqncrna <- getBM(filters="refseq_ncrna", attributes= c("refseq_ncrna", "hgnc_symbol", "entrezgene"), values=nimblegenid, mart=mart, uniqueRows = TRUE)
refseqncrnapredicted <- getBM(filters="refseq_ncrna_predicted", attributes= c("refseq_ncrna_predicted", "hgnc_symbol", "entrezgene"), values=nimblegenid, mart=mart, uniqueRows = TRUE)
refseqpeptide <- getBM(filters="refseq_peptide", attributes= c("refseq_peptide", "hgnc_symbol", "entrezgene"), values=nimblegenid, mart=mart, uniqueRows = TRUE)
refseqpeptidepredicted <- getBM(filters="refseq_peptide_predicted", attributes= c("refseq_peptide_predicted", "hgnc_symbol", "entrezgene"), values=nimblegenid, mart=mart, uniqueRows = TRUE)
write.table(refseqmrna, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_RefSeqIDtoGeneSymbol_biomaRt.txt", sep = "\t", append = TRUE, row.names = FALSE)
write.table(refseqmrnapredicted, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_RefSeqIDtoGeneSymbol_biomaRt.txt", sep = "\t", append = TRUE, row.names = FALSE)
write.table(refseqncrna, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_RefSeqIDtoGeneSymbol_biomaRt.txt", sep = "\t", append = TRUE, row.names = FALSE)
write.table(refseqncrnapredicted, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_RefSeqIDtoGeneSymbol_biomaRt.txt", sep = "\t", append = TRUE, row.names = FALSE)
write.table(refseqpeptide, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_RefSeqIDtoGeneSymbol_biomaRt.txt", sep = "\t", append = TRUE, row.names = FALSE)
write.table(refseqpeptidepredicted, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/RESULTS/AA/AA_RNA_Seq_20140630/20140630-NewestRound/Analyses/Combined/AA37TxExp_RefSeqIDtoGeneSymbol_biomaRt.txt", sep = "\t", append = TRUE, row.names = FALSE)





########################
# APPENDIX.
########################
# EXAMPLE:
library(biomaRt)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
refseq <- c("NP_001165902", "NP_001166046", "NP_001166152")
getBM(filters="refseq_peptide", attributes="external_gene_id", values=refseq, mart=mart)
  external_gene_id
1            AMPD3
2            MCPH1
3           ZNF641

# ALL FILTERS:
> listFilters(mart)[grep(".",as.character(listFilters(mart)$description)),]
                                         name                                                                   description
1                             chromosome_name                                                               Chromosome name
2                                       start                                                               Gene Start (bp)
3                                         end                                                                 Gene End (bp)
4                                  band_start                                                                    Band Start
5                                    band_end                                                                      Band End
6                                marker_start                                                                  Marker Start
7                                  marker_end                                                                    Marker End
8                                        type                                                                          Type
9                               encode_region                                                                 Encode region
10                                     strand                                                                        Strand
11                         chromosomal_region                     Chromosome Regions (e.g 1:100:10000:-1,1:100000:200000:1)
12                       with_ox_arrayexpress                                                       with ArrayExpress ID(s)
13                 with_illumina_humanwg_6_v1                                        with Illumina HumanWG 6 v1 probe ID(s)
14                 with_illumina_humanwg_6_v2                                        with Illumina HumanWG 6 v2 probe ID(s)
15                 with_illumina_humanwg_6_v3                                        with Illumina HumanWG 6 v3 probe ID(s)
16                with_illumina_humanht_12_v3                                      with Illumina Human HT 12 v3 probe ID(s)
17                with_illumina_humanht_12_v4                                      with Illumina Human HT 12 v4 probe ID(s)
18                with_illumina_humanref_8_v3                                       with Illumina Human HT 8 v3 probe ID(s)
19                      with_phalanx_onearray                                             with Phalanx onearray probe ID(s)
20                     with_codelink_codelink                                                     with Codelink probe ID(s)
21                       with_agilent_cgh_44b                                              with Agilent CGH 44b probe ID(s)
22                              with_wikigene                                                           with WikiGene ID(s)
23                        with_affy_primeview                                    with Affymetrix Microarray primeview ID(s)
24     with_efg_agilent_sureprint_g3_ge_8x60k                                  with Efg agilent sureprint g3 ge 8x60k ID(s)
25  with_efg_agilent_sureprint_g3_ge_8x60k_v2                               with Efg agilent sureprint g3 ge 8x60k v2 ID(s)
26      with_efg_agilent_wholegenome_4x44k_v1                                   with Efg agilent wholegenome 4x44k v1 ID(s)
27      with_efg_agilent_wholegenome_4x44k_v2                                   with Efg agilent wholegenome 4x44k v2 ID(s)
28                                  with_hgnc                                                               with HGNC ID(s)
29              with_ox_clone_based_vega_gene                                              with clone based VEGA gene ID(s)
30           with_ox_clone_based_ensembl_gene                                           with clone based Ensembl gene ID(s)
31     with_ox_clone_based_ensembl_transcript                                     with clone based Ensembl transcript ID(s)
32        with_ox_clone_based_vega_transcript                                        with clone based VEGA transcript ID(s)
33                                  with_ottg                                                   with VEGA gene ID(s) (OTTG)
34                                  with_ottt                                             with VEGA transcript ID(s) (OTTT)
35          with_shares_cds_and_utr_with_ottt                         with HAVANA transcript (where ENST identical to OTTT)
36                  with_shares_cds_with_ottt                      with HAVANA transcript (where ENST shares CDS with OTTT)
37                               with_mirbase                                                            with miRBase ID(s)
38                                 with_rfamR                                                               with Rfam ID(s)
39                                  with_ccds                                                               with CCDS ID(s)
40                                with_chembl                                                             with ChEMBL ID(s)
41                                  with_embl                                                               with EMBL ID(s)
42                            with_entrezgene                                                         with EntrezGene ID(s)
43                                 with_go_id                                                     with GO Term Accession(s)
44                                with_merops                                                             with MEROPS ID(s)
45                              with_mim_gene                                                           with MIM gene ID(s)
46                            with_mim_morbid                                                        with MIM disease ID(s)
47                                  with_ucsc                                                               with UCSC ID(s)
48                                   with_pdb                                                                with PDB ID(s)
49                        with_refseq_peptide                                                     with RefSeq protein ID(s)
50              with_refseq_peptide_predicted                                           with RefSeq predicted protein ID(s)
51                            with_protein_id                                                   with protein(Genbank) ID(s)
52                       with_uniprotsptrembl                                            with UniProtKB/TrEMBL Accession(s)
53                               with_unigene                                                            with UniGene ID(s)
54                      with_uniprot_genename                                                   with UniProt Genename ID(s)
55      with_uniprot_genename_transcript_name                                      with Uniprot Genename Transcript Name(s)
56              with_exp_est_anatomicalsystem                                         with Expression est anatomical system
57                with_exp_est_associatedwith                                           with Expression est associated with
58                      with_exp_est_celltype                                                 with Expression est cell type
59              with_exp_est_developmentstage                                         with Expression est development stage
60         with_exp_est_experimentaltechnique                                    with Expression est experimental technique
61            with_exp_est_microarrayplatform                                       with Expression est microarray platform
62                     with_exp_est_pathology                                                 with Expression est pathology
63                       with_exp_est_pooling                                                   with Expression est pooling
64             with_exp_est_tissuepreparation                                        with Expression est tissue preparation
65                     with_exp_est_treatment                                                 with Expression est treatment
66                    with_exp_atlas_celltype                                                with GNF/Atlas cell type ID(s)
67                with_exp_atlas_organismpart                                            with GNF/Atlas organism part ID(s)
68                with_exp_atlas_diseasestate                                            with GNF/Atlas disease state ID(s)
69                                with_dbass3                                                             with DBASS3 ID(s)
70                                   with_hpa                                                with Human Protein Atlas ID(s)
71                          with_affy_hc_g110                             with Affymetrix Microarray hc g110 probeset ID(s)
72                         with_affy_hg_focus                            with Affymetrix Microarray hg Focus probeset ID(s)
73                         with_affy_hg_u133a                            with Affymetrix Microarray hg u133a probeset ID(s)
74                         with_affy_hg_u133b                            with Affymetrix Microarray hg u133b probeset ID(s)
75                       with_affy_hg_u133a_2                          with Affymetrix Microarray hg u133a 2 probeset ID(s)
76                   with_affy_hg_u133_plus_2                      with Affymetrix Microarray hg u133 plus 2 probeset ID(s)
77                          with_affy_hg_u95a                             with Affymetrix Microarray hg u95a probeset ID(s)
78                        with_affy_hg_u95av2                           with Affymetrix Microarray hg u95av2 ID(s) probeset
79                          with_affy_hg_u95b                             with Affymetrix Microarray hg u95b probeset ID(s)
80                          with_affy_hg_u95c                             with Affymetrix Microarray hg u95c probeset ID(s)
81                          with_affy_hg_u95d                             with Affymetrix Microarray hg u95d probeset ID(s)
82                          with_affy_hg_u95e                             with Affymetrix Microarray hg u95e probeset ID(s)
83                         with_affy_u133_x3p                            with Affymetrix Microarray u133 x3p probeset ID(s)
84                         with_affy_hugenefl                            with Affymetrix Microarray HuGeneFL probeset ID(s)
85                 with_affy_hugene_1_0_st_v1                    with Affymetrix Microarray hugene 1 0 st v1 probeset ID(s)
86                 with_affy_hugene_2_0_st_v1                             with Affymetrix Microarray hugene 2 0 st v1 ID(s)
87                   with_affy_huex_1_0_st_v2                      with Affymetrix Microarray huex 1 0 st v2 probeset ID(s)
88                         with_ox_goslim_goa                                                            with GOSlim GOA(s)
89                        with_ox_ens_hs_gene                                             with Ensembl to LRG link gene IDs
90                 with_ox_ens_hs_translation                                      with Ensembl to LRG link translation IDs
91                  with_ox_ens_hs_transcript                                       with Ensembl to LRG link transcript IDs
92                            with_ox_uniparc                                                            with UniParc ID(s)
93                       with_ox_ens_lrg_gene                                                   with Ensembl LRG gene ID(s)
94                 with_ox_ens_lrg_transcript                                             with Ensembl LRG transcript ID(s)
95               with_ox_hgnc_transcript_name                                                  with HGNC transcript name(s)
96               with_ox_rfam_transcript_name                                                  with Rfam transcript name(s)
97            with_ox_mirbase_transcript_name                                               with miRBase transcript name(s)
98                    with_ox_uniprotsptrembl                                                    with UniProtSPTREMBL ID(s)
99                   with_ox_uniprotswissprot                                                   with UniProtSWISSPROT ID(s)
100                                with_go_go                                                                 with GO ID(s)
101                       with_ox_refseq_mrna                                                        with RefSeq mRNA ID(s)
102             with_ox_refseq_mrna_predicted                                              with RefSeq mRNA predicted ID(s)
103                      with_ox_refseq_ncrna                                                       with RefSeq ncRNA ID(s)
104            with_ox_refseq_ncrna_predicted                                             with RefSeq ncRNA predicted ID(s)
105                           ensembl_gene_id                                     Ensembl Gene ID(s) [e.g. ENSG00000139618]
106                     ensembl_transcript_id                               Ensembl Transcript ID(s) [e.g. ENST00000380152]
107                        ensembl_peptide_id                                  Ensembl protein ID(s) [e.g. ENSP00000369497]
108                           ensembl_exon_id                                     Ensembl exon ID(s) [e.g. ENSE00001508081]
109                                   hgnc_id                                                       HGNC ID(s) [e.g. 43668]
110                               hgnc_symbol                                                     HGNC symbol(s) [e.g. ZFY]
111                      hgnc_transcript_name                                    HGNC transcript name(s) [e.g. QRSL1P2-001]
112                              arrayexpress                                     ArrayExpress ID(s) [e.g. ENSG00000241328]
113                                      ccds                                                   CCDS ID(s) [e.g. CCDS10187]
114                                    chembl                                       ChEMBL ID(s) ID(s) [e.g. CHEMBL1075092]
115             clone_based_ensembl_gene_name                            Clone based Ensembl gene name(s) [e.g. AL691479.1]
116       clone_based_ensembl_transcript_name                  Clone based Ensembl transcript name(s) [e.g. AL691479.1-201]
117                clone_based_vega_gene_name                             Clone based VEGA gene name(s) [e.g. RP11-815M8.1]
118          clone_based_vega_transcript_name                   Clone based VEGA transcript name(s) [e.g. RP11-815M8.1-001]
119                                  codelink                                          Codelink probe ID(s) [e.g. GE550734]
120                               dbass3_name                                                 DBASS3 Gene Name [e.g. PDE6B]
121                                      embl                                                      EMBL ID(s) [e.g. D87018]
122                               ens_hs_gene                           Ensembl to LRG link gene IDs [e.g. ENSG00000108821]
123                         ens_hs_transcript                     Ensembl to LRG link transcript IDs [e.g. ENST00000225964]
124                        ens_hs_translation                    Ensembl to LRG link translation IDs [e.g. ENSP00000376544]
125                              ens_lrg_gene                                     LRG to Ensembl link gene IDs [e.g. LRG_3]
126                        ens_lrg_transcript                           LRG to Ensembl link transcript IDs [e.g. LRG_226t1]
127                                entrezgene                                                EntrezGene ID(s) [e.g. 115286]
128                                     go_id                                        GO Term Accession(s) [e.g. GO:0005515]
129                      goslim_goa_accession                                    GOSlim GOA Accessions(s) [e.g. GO:0005623]
130                                       hpa                              Human Protein Atlas Antibody ID [e.g. HPA002549]
131                      shares_cds_with_ottt HAVANA transcript (where ENST shares CDS with OTTT) [e.g. OTTHUMT00000088063]
132              shares_cds_and_utr_with_ottt    HAVANA transcript (where ENST identical to OTTT) [e.g. OTTHUMT00000088055]
133                                    merops                                                   MEROPS ID(s) [e.g. C19.028]
134                        mim_gene_accession                                           MIM Gene Accession(s) [e.g. 611882]
135                      mim_morbid_accession                                         MIM Morbid Accession(s) [e.g. 540000]
136                                mirbase_id                                              miRBase ID(s) [e.g. hsa-mir-137]
137                         mirbase_accession                                         miRBase Accession(s) [e.g. MI0000454]
138                   mirbase_transcript_name                              miRBase transcript name [e.g. hsa-mir-877.2-201]
139                                       pdb                                                         PDB ID(s) [e.g. 1J47]
140                                protein_id                                       Protein (Genbank) ID(s) [e.g. BAA20017]
141                               refseq_mrna                                         Refseq mRNA ID(s) [e.g. NM_001195597]
142                     refseq_mrna_predicted                               Refseq Predicted mRNA ID(s) [e.g. XM_001125684]
143                              refseq_ncrna                                           Refseq ncRNA ID(s) [e.g. NR_002834]
144                    refseq_ncrna_predicted                                 Refseq Predicted ncRNA ID(s) [e.g. XR_111404]
145                            refseq_peptide                                      Refseq protein ID(s) [e.g. NP_001005353]
146                  refseq_peptide_predicted                            Refseq predicted protein ID(s) [e.g. XP_005263247]
147                                      rfam                                                     Rfam ID(s) [e.g. RF00432]
148                      rfam_transcript_name                                      Rfam transcript name(s) [e.g. U3.30-201]
149                                      ucsc                                                  UCSC ID(s) [e.g. uc010ajn.1]
150                          uniprot_sptrembl                                     UniProt/TrEMBL Accession(s) [e.g. A2MYD1]
151                         uniprot_swissprot                                    UniProt/Swissprot ID(s) [e.g. GAGE4_HUMAN]
152               uniprot_swissprot_accession                                  UniProt/Swissprot Accession(s) [e.g. Q13068]
153                                   unigene                                                UniGene ID(s) [e.g. Hs.602394]
154                          uniprot_genename                                            UniProt Genename ID(s) [e.g. V4-4]
155                                   uniparc                                            UniParc ID(s) [e.g. UPI0000000AA1]
156          uniprot_genename_transcript_name                       Uniprot Genename Transcript Name ID(s) [e.g. SEPT1-202]
157                                      ottg                              VEGA Gene ID(s) (OTTG) [e.g. OTTHUMG00000036159]
158                                      ottt                        VEGA Transcript ID(s) (OTTT) [e.g. OTTHUMT00000088063]
159                               wikigene_id                                                  WikiGene ID(s) [e.g. 115286]
160                             wikigene_name                                              WikiGene Name(s) [e.g. SLC25A26]
161                              affy_hc_g110                                   Affy hc g110 probeset ID(s) [e.g. 113_i_at]
162                             affy_hg_focus                                 Affy hg focus probeset ID(s) [e.g. 201612_at]
163                              affy_hg_u95a                                   Affy hg u95a probeset ID(s) [e.g. 32647_at]
164                            affy_hg_u95av2                                 Affy hg u95av2 probeset ID(s) [e.g. 32647_at]
165                              affy_hg_u95b                                   Affy hg u95b probeset ID(s) [e.g. 53925_at]
166                              affy_hg_u95c                                 Affy hg u95c probeset ID(s) [e.g. 61056_r_at]
167                              affy_hg_u95d                                   Affy hg u95d probeset ID(s) [e.g. 79632_at]
168                              affy_hg_u95e                                   Affy hg u95e probeset ID(s) [e.g. 79965_at]
169                           affy_hg_u133a_2                             Affy hg u133a 2 probeset ID(s) [e.g. 200874_s_at]
170                             affy_hg_u133a                               Affy hg u133a probeset ID(s) [e.g. 200874_s_at]
171                             affy_hg_u133b                                 Affy hg u133b probeset ID(s) [e.g. 227057_at]
172                       affy_hg_u133_plus_2                           Affy hg u133 plus 2 probeset ID(s) [e.g. 241843_at]
173                             affy_hugenefl                              Affy HuGene FL probeset ID(s) [e.g. M58525_s_at]
174                     affy_hugene_1_0_st_v1                           Affy HuGene 1_0 st v1 probeset ID(s) [e.g. 8016215]
175                     affy_hugene_2_0_st_v1                          Affy HuGene 2_0 st v1 probeset ID(s) [e.g. 16942487]
176                       affy_huex_1_0_st_v2                             Affy HuEx 1_0 st v2 probeset ID(s) [e.g. 4033465]
177                            affy_primeview                      Affymetrix Microarray Primeview ID(s) [e.g. 11763890_at]
178                             affy_u133_x3p                     Affy u133 x3p probeset ID(s) [e.g. Hs2.205326.1.A1_3p_at]
179                           agilent_cgh_44b                               Agilent CGH 44b probe ID(s) [e.g. A_14_P131077]
180         efg_agilent_sureprint_g3_ge_8x60k                Agilent Sureprint G3 GE 8x60k probe ID(s) [e.g. A_33_P3356022]
181      efg_agilent_sureprint_g3_ge_8x60k_v2              Agilent Sureprint G3 GE 8x60k v2 probe ID(s) [e.g. A_24_P368544]
182          efg_agilent_wholegenome_4x44k_v1                  Agilent WholeGenome 4x44k v1 probe ID(s) [e.g. A_32_P196615]
183          efg_agilent_wholegenome_4x44k_v2                 Agilent WholeGenome 4x44k v2 probe ID(s) [e.g. A_33_P3356022]
184                     illumina_humanwg_6_v1                           Illumina HumanWG 6 V1 probe ID(s) [e.g. 0000940471]
185                     illumina_humanwg_6_v2                         Illumina HumanWG 6 V2 probe ID(s) [e.g. ILMN_1748182]
186                     illumina_humanwg_6_v3                         Illumina HumanWG 6 v3 probe ID(s) [e.g. ILMN_2103362]
187                    illumina_humanref_8_v3                       Illumina Human Ref 8 v3 probe ID(s) [e.g. ILMN_1768251]
188                    illumina_humanht_12_v3                       Illumina Human HT 12 v3 probe ID(s) [e.g. ILMN_1672925]
189                    illumina_humanht_12_v4                       Illumina Human HT 12 v4 probe ID(s) [e.g. ILMN_1768251]
190                          phalanx_onearray                             Phalanx OneArray probe ID(s) [e.g. PH_hs_0031946]
191                          transcript_count                                                           Transcript count >=
192                                   biotype                                                                          Type
193                                    source                                                                 Source (gene)
194                         transcript_source                                                           Source (transcript)
195                                    status                                                                 Status (gene)
196                         transcript_status                                                           Status (transcript)
197                     phenotype_description                                                         Phenotype description
198                                event_type                                                                    Event Type
199                          go_evidence_code                                                              GO Evidence code
200                            go_parent_term                                                         Parent term accession
201                            go_parent_name                                                              Parent term name
202                    anatomical_system_term                                                             Anatomical system
203                    development_stage_term                                                             Development Stage
204                            cell_type_term                                                                     Cell type
205                            pathology_term                                                                     Pathology
209                         with_paralog_hsap                                                        Paralogous Human Genes
210                         with_homolog_vpac                                                      Orthologous Alpaca Genes
211                         with_homolog_acar                                                Orthologous Anole Lizard Genes
212                         with_homolog_dnov                                                   Orthologous Armadillo Genes
213                         with_homolog_gmor                                                Orthologous Atlantic Cod Genes
214                         with_homolog_ogar                                                    Orthologous Bushbaby Genes
215                         with_homolog_cele                                      Orthologous Caenorhabditis elegans Genes
216                         with_homolog_fcat                                                         Orthologous Cat Genes
217                         with_homolog_amex                                                   Orthologous Cave fish Genes
218                         with_homolog_ggal                                                     Orthologous Chicken Genes
219                         with_homolog_ptro                                                  Orthologous Chimpanzee Genes
220                         with_homolog_psin                                    Orthologous Chinese softshell turtle Genes
221                         with_homolog_cint                                          Orthologous Ciona intestinalis genes
222                         with_homolog_csav                                              Orthologous Ciona savignyi Genes
223                         with_homolog_lcha                                                  Orthologous Coelacanth Genes
224                         with_homolog_sara                                                Orthologous Common Shrew Genes
225                         with_homolog_btau                                                         Orthologous Cow Genes
226                         with_homolog_cfam                                                         Orthologous Dog Genes
227                         with_homolog_ttru                                                     Orthologous Dolphin Genes
228                         with_homolog_apla                                                        Orthologous Duck Genes
229                         with_homolog_dmel                                                  Orthologous Drosophila Genes
230                         with_homolog_lafr                                                    Orthologous Elephant Genes
231                         with_homolog_mfur                                                      Orthologous Ferret Genes
232                         with_homolog_falb                                                  Orthologous Flycatcher Genes
233                         with_homolog_trub                                                        Orthologous Fugu Genes
234                         with_homolog_nleu                                                      Orthologous Gibbon Genes
235                         with_homolog_ggor                                                     Orthologous Gorilla Genes
236                         with_homolog_cpor                                                  Orthologous Guinea Pig Genes
237                         with_homolog_eeur                                                    Orthologous Hedgehog Genes
238                         with_homolog_ecab                                                       Orthologous Horse Genes
239                         with_homolog_dord                                                Orthologous Kangaroo Rat Genes
240                         with_homolog_pmar                                                     Orthologous Lamprey Genes
241                         with_homolog_etel                                      Orthologous Lesser hedgehog tenrec Genes
242                         with_homolog_mmul                                                     Orthologous Macaque Genes
243                         with_homolog_cjac                                                    Orthologous Marmoset Genes
244                         with_homolog_olat                                                      Orthologous Medaka Genes
245                         with_homolog_pvam                                                     Orthologous Megabat Genes
246                         with_homolog_mluc                                                    Orthologous Microbat Genes
247                         with_homolog_mmus                                                       Orthologous Mouse Genes
248                         with_homolog_mmur                                                 Orthologous Mouse Lemur Genes
249                         with_homolog_onil                                                Orthologous Nile tilapia Genes
250                         with_homolog_mdom                                                     Orthologous Opossum Genes
251                         with_homolog_pabe                                                   Orthologous Orangutan Genes
252                         with_homolog_amel                                                       Orthologous Panda Genes
253                         with_homolog_sscr                                                         Orthologous Pig Genes
254                         with_homolog_opri                                                        Orthologous Pika Genes
255                         with_homolog_xmac                                                   Orthologous Platyfish Genes
256                         with_homolog_oana                                                    Orthologous Platypus Genes
257                         with_homolog_ocun                                                      Orthologous Rabbit Genes
258                         with_homolog_rnor                                                         Orthologous Rat Genes
259                         with_homolog_pcap                                                  Orthologous Rock Hyrax Genes
260                         with_homolog_oari                                                       Orthologous Sheep Genes
261                         with_homolog_chof                                                       Orthologous Sloth Genes
262                         with_homolog_locu                                                 Orthologous Spotted gar Genes
263                         with_homolog_itri                                                    Orthologous Squirrel Genes
264                         with_homolog_gacu                                                 Orthologous Stickleback Genes
265                         with_homolog_tsyr                                                     Orthologous Tarsier Genes
266                         with_homolog_shar                                             Orthologous Tasmanian Devil Genes
267                         with_homolog_tnig                                                   Orthologous Tetraodon Genes
268                         with_homolog_tbel                                                  Orthologous Tree Shrew Genes
269                         with_homolog_mgal                                                      Orthologous Turkey Genes
270                         with_homolog_meug                                                     Orthologous Wallaby Genes
271                         with_homolog_xtro                                                     Orthologous Xenopus Genes
272                         with_homolog_scer                                                       Orthologous Yeast Genes
273                         with_homolog_tgut                                                 Orthologous Zebra Finch Genes
274                         with_homolog_drer                                                   Orthologous Zebrafish Genes
275                              with_profile                                             with Protein feature pfscan ID(s)
276                                with_tmhmm                                              with Protein feature tmhmm ID(s)
277                              with_tigrfam                                            with Protein feature tigrfam ID(s)
278                          with_superfamily                                        with Protein feature superfamily ID(s)
279                                with_smart                                              with Protein feature smart ID(s)
280                              with_signalp                                            with Protein feature signalp ID(s)
281                       with_low_complexity                                                with Protein feature seg ID(s)
282                                with_pirsf                                              with Protein feature pirsf ID(s)
283                                 with_coil                                             with Protein feature ncoils ID(s)
284                             with_interpro                                                           with InterPro ID(s)
285                 with_protein_feature_pfam                                                               with PFAM ID(s)
286               with_protein_feature_prints                                                             with PRINTS ID(s)
287                                   tigrfam                                                TIGRfam ID(s) [e.g. TIGR00172]
288                               superfamily                                             Superfamily ID(s) [e.g. SSF47095]
289                                     smart                                                    SMART ID(s) [e.g. SM00398]
290                                     pirsf                                                PIRSF ID(s) [e.g. PIRSF037653]
291                                    family                       Ensembl Protein Family ID(s) [e.g. ENSFM00250000000002]
292                                      pfam                                                     PFAM ID(s) [e.g. PF00046]
293                                    prints                                                   PRINTS ID(s) [e.g. PR00194]
294                                   profile                                                  PROFILE ID(s) [e.g. PS50313]
295                                  interpro                                               Interpro ID(s) [e.g. IPR007087]
296                 with_transmembrane_domain                                                         Transmembrane domains
297                        with_signal_domain                                                                Signal domains
298                germ_line_variation_source                           limit to genes with germline variation data sources
299                  somatic_variation_source                            limit to genes with somatic variation data sources
300                        with_validated_snp                                                Associated with validated SNPs
301                            so_parent_name                                                              Parent term name

> listMarts()
                                 biomart                                                                                                version
1                                ensembl                                                                           ENSEMBL GENES 79 (SANGER UK)
2                                    snp                                                                       ENSEMBL VARIATION 79 (SANGER UK)
3                             regulation                                                                      ENSEMBL REGULATION 79 (SANGER UK)
4                                   vega                                                                                   VEGA 59  (SANGER UK)
5                          fungi_mart_26                                                                              ENSEMBL FUNGI 26 (EBI UK)
6                    fungi_variations_26                                                                    ENSEMBL FUNGI VARIATION 26 (EBI UK)
7                        metazoa_mart_26                                                                            ENSEMBL METAZOA 26 (EBI UK)
8                  metazoa_variations_26                                                                  ENSEMBL METAZOA VARIATION 26 (EBI UK)
9                         plants_mart_26                                                                             ENSEMBL PLANTS 26 (EBI UK)
10                  plants_variations_26                                                                   ENSEMBL PLANTS VARIATION 26 (EBI UK)
11                      protists_mart_25                                                                           ENSEMBL PROTISTS 25 (EBI UK)
12                protists_variations_26                                                                 ENSEMBL PROTISTS VARIATION 26 (EBI UK)
13                                   msd                                                                                           MSD (EBI UK)
14                            cg_mart_02                                                              PROTEOMICS (UNIVERSITY OF CAMBRIDGE - UK)
15                                  htgt                                                                WTSI MOUSE GENETICS PROJECT (SANGER UK)
16                                 WS220                                                                                 WORMBASE 220 (CSHL US)
17                         parasite_mart                                                                                 ParaSite Mart (EBI UK)
18                               biomart                                                                            MGI (JACKSON LABORATORY US)
19                               example                                                                   FANTOM5 phase1.1 (RIKEN CSLST Japan)
20                                 pride                                                                                         PRIDE (EBI UK)
21                      prod-intermart_1                                                                                      INTERPRO (EBI UK)
22                               unimart                                                                                       UNIPROT (EBI UK)
23                             biomartDB                                                                        PARAMECIUM GENOME (CNRS FRANCE)
24                              biblioDB                                                                  PARAMECIUM BIBLIOGRAPHY (CNRS FRANCE)
25                    Eurexpress Biomart                                                                          EUREXPRESS (MRC EDINBURGH UK)
26                        phytozome_mart                                                                                              Phytozome
27                         metazome_mart                                                                                               Metazome
28                          HapMap_rel27                                                                                    HAPMAP 27 (NCBI US)
29                           oncomodules                                                                                    INTOGEN ONCOMODULES
30                                  ikmc                                                                         IKMC GENES AND PRODUCTS (IKMC)
31                 EMAGE gene expression                                                                                  EMAGE GENE EXPRESSION
32                 EMAP anatomy ontology                                                                                  EMAP ANATOMY ONTOLOGY
33               EMAGE browse repository                                                                                EMAGE BROWSE REPOSITORY
34                            GermOnline                                                                                             GERMONLINE
35   Sigenae_Oligo_Annotation_Ensembl_61                                                                  SIGENAE OLIGO ANNOTATION (ENSEMBL 61)
36 Sigenae Oligo Annotation (Ensembl 59)                                                                  SIGENAE OLIGO ANNOTATION (ENSEMBL 59)
37 Sigenae Oligo Annotation (Ensembl 56)                                                                  SIGENAE OLIGO ANNOTATION (ENSEMBL 56)
38                        Breast_mart_69                                                           BCCTB Bioinformatics Portal (UK and Ireland)
39                          K562_Gm12878 Predictive models of gene regulation from processed high-throughput epigenomics data: K562 vs. Gm12878
40                             Hsmm_Hmec    Predictive models of gene regulation from processed high-throughput epigenomics data: Hsmm vs. Hmec
41                            Pancreas63                                             PANCREATIC EXPRESSION DATABASE (BARTS CANCER INSTITUTE UK)
42                    Public_OBIOMARTPUB        Multi-species: marker, QTL, SNP, gene, germplasm, phenotype, association, with Gene annotations
43                          Public_VITIS                               Grapevine 8x, stuctural annotation with Genetic maps (genetic markers..)
44                      Public_VITIS_12x             Grapevine 12x.0, stuctural and functional annotation with Genetic maps (genetic markers..)
45                            Prod_WHEAT                                      Wheat, stuctural annotation with Genetic maps (genetic markers..)
46                        Public_TAIRV10                                              Arabidopsis Thaliana TAIRV10, genes functional annotation
47                          Public_MAIZE                                                            Zea mays ZmB73, genes functional annotation
48                           Prod_TOMATO                                                            Tomato, stuctural and functional annotation
49                           Prod_POPLAR                                                       Populus trichocarpa, genes functional annotation
50                        Prod_POPLAR_V2                                                  Populus trichocarpa, genes functional annotation V2.0
51                     Prod_BOTRYTISEDIT                                                      Botrytis cinerea T4, genes functional annotation 
52                            Prod_BOFUB                                                   Botrytis cinerea B0510, genes functional annotation 
53                    Prod_LMACULANSEDIT                                                    Leptosphaeria maculans, genes functional annotation
54                     vb_gene_mart_1502                                                                                       VectorBase Genes
55                      vb_snp_mart_1502                                                                                   VectorBase Variation
56                            expression                                                                                  VectorBase Expression
57                    ENSEMBL_MART_PLANT                                                             GRAMENE 40 ENSEMBL GENES (CSHL/CORNELL US)
58                ENSEMBL_MART_PLANT_SNP                                                                 GRAMENE 40 VARIATION (CSHL/CORNELL US)

> listDatasets(mart = useMart("ensembl"))
                          dataset                                 description         version
1          oanatinus_gene_ensembl      Ornithorhynchus anatinus genes (OANA5)           OANA5
2         cporcellus_gene_ensembl             Cavia porcellus genes (cavPor3)         cavPor3
3         gaculeatus_gene_ensembl      Gasterosteus aculeatus genes (BROADS1)         BROADS1
4          lafricana_gene_ensembl          Loxodonta africana genes (loxAfr3)         loxAfr3
5  itridecemlineatus_gene_ensembl  Ictidomys tridecemlineatus genes (spetri2)         spetri2
6         choffmanni_gene_ensembl         Choloepus hoffmanni genes (choHof1)         choHof1
7          csavignyi_gene_ensembl              Ciona savignyi genes (CSAV2.0)         CSAV2.0
8             fcatus_gene_ensembl         Felis catus genes (Felis_catus_6.2) Felis_catus_6.2
9        rnorvegicus_gene_ensembl          Rattus norvegicus genes (Rnor_5.0)        Rnor_5.0
10         psinensis_gene_ensembl      Pelodiscus sinensis genes (PelSin_1.0)      PelSin_1.0
11          cjacchus_gene_ensembl   Callithrix jacchus genes (C_jacchus3.2.1)  C_jacchus3.2.1
12        ttruncatus_gene_ensembl          Tursiops truncatus genes (turTru1)         turTru1
13       scerevisiae_gene_ensembl    Saccharomyces cerevisiae genes (R64-1-1)         R64-1-1
14          celegans_gene_ensembl     Caenorhabditis elegans genes (WBcel235)        WBcel235
15          csabaeus_gene_ensembl       Chlorocebus sabaeus genes (ChlSab1.1)       ChlSab1.1
16        oniloticus_gene_ensembl     Oreochromis niloticus genes (Orenil1.0)       Orenil1.0
17         trubripes_gene_ensembl           Takifugu rubripes genes (FUGU4.0)         FUGU4.0
18        amexicanus_gene_ensembl        Astyanax mexicanus genes (AstMex102)       AstMex102
19          pmarinus_gene_ensembl     Petromyzon marinus genes (Pmarinus_7.0)    Pmarinus_7.0
20        eeuropaeus_gene_ensembl         Erinaceus europaeus genes (eriEur1)         eriEur1
21       falbicollis_gene_ensembl      Ficedula albicollis genes (FicAlb_1.4)      FicAlb_1.4
22      ptroglodytes_gene_ensembl          Pan troglodytes genes (CHIMP2.1.4)      CHIMP2.1.4
23         etelfairi_gene_ensembl            Echinops telfairi genes (TENREC)          TENREC
24     cintestinalis_gene_ensembl               Ciona intestinalis genes (KH)              KH
25       nleucogenys_gene_ensembl         Nomascus leucogenys genes (Nleu1.0)         Nleu1.0
26           sscrofa_gene_ensembl              Sus scrofa genes (Sscrofa10.2)     Sscrofa10.2
27        ocuniculus_gene_ensembl     Oryctolagus cuniculus genes (OryCun2.0)       OryCun2.0
28     dnovemcinctus_gene_ensembl      Dasypus novemcinctus genes (Dasnov3.0)       Dasnov3.0
29         pcapensis_gene_ensembl           Procavia capensis genes (proCap1)         proCap1
30          tguttata_gene_ensembl     Taeniopygia guttata genes (taeGut3.2.4)     taeGut3.2.4
31        mlucifugus_gene_ensembl            Myotis lucifugus genes (myoLuc2)         myoLuc2
32          hsapiens_gene_ensembl              Homo sapiens genes (GRCh38.p2)       GRCh38.p2
33          pformosa_gene_ensembl       Poecilia formosa genes (PoeFor_5.1.2)    PoeFor_5.1.2
34             mfuro_gene_ensembl  Mustela putorius furo genes (MusPutFur1.0)    MusPutFur1.0
35        tbelangeri_gene_ensembl            Tupaia belangeri genes (tupBel1)         tupBel1
36           ggallus_gene_ensembl               Gallus gallus genes (Galgal4)         Galgal4
37       xtropicalis_gene_ensembl           Xenopus tropicalis genes (JGI4.2)          JGI4.2
38         ecaballus_gene_ensembl              Equus caballus genes (EquCab2)         EquCab2
39           pabelii_gene_ensembl                  Pongo abelii genes (PPYG2)           PPYG2
40        xmaculatus_gene_ensembl   Xiphophorus maculatus genes (Xipmac4.4.2)     Xipmac4.4.2
41            drerio_gene_ensembl                     Danio rerio genes (Zv9)             Zv9
42        lchalumnae_gene_ensembl         Latimeria chalumnae genes (LatCha1)         LatCha1
43     tnigroviridis_gene_ensembl Tetraodon nigroviridis genes (TETRAODON8.0)    TETRAODON8.0
44      amelanoleuca_gene_ensembl      Ailuropoda melanoleuca genes (ailMel1)         ailMel1
45          mmulatta_gene_ensembl               Macaca mulatta genes (MMUL_1)          MMUL_1
46         pvampyrus_gene_ensembl           Pteropus vampyrus genes (pteVam1)         pteVam1
47           panubis_gene_ensembl              Papio anubis genes (PapAnu2.0)       PapAnu2.0
48        mdomestica_gene_ensembl       Monodelphis domestica genes (monDom5)         monDom5
49     acarolinensis_gene_ensembl       Anolis carolinensis genes (AnoCar2.0)       AnoCar2.0
50            vpacos_gene_ensembl               Vicugna pacos genes (vicPac1)         vicPac1
51         tsyrichta_gene_ensembl            Tarsius syrichta genes (tarSyr1)         tarSyr1
52        ogarnettii_gene_ensembl          Otolemur garnettii genes (OtoGar3)         OtoGar3
53     dmelanogaster_gene_ensembl       Drosophila melanogaster genes (BDGP6)           BDGP6
54          mmurinus_gene_ensembl          Microcebus murinus genes (micMur1)         micMur1
55         loculatus_gene_ensembl        Lepisosteus oculatus genes (LepOcu1)         LepOcu1
56          olatipes_gene_ensembl                Oryzias latipes genes (HdrR)            HdrR
57          ggorilla_gene_ensembl           Gorilla gorilla genes (gorGor3.1)       gorGor3.1
58         oprinceps_gene_ensembl         Ochotona princeps genes (OchPri2.0)       OchPri2.0
59            dordii_gene_ensembl             Dipodomys ordii genes (dipOrd1)         dipOrd1
60            oaries_gene_ensembl                 Ovis aries genes (Oar_v3.1)        Oar_v3.1
61         mmusculus_gene_ensembl              Mus musculus genes (GRCm38.p3)       GRCm38.p3
62        mgallopavo_gene_ensembl            Meleagris gallopavo genes (UMD2)            UMD2
63           gmorhua_gene_ensembl                Gadus morhua genes (gadMor1)         gadMor1
64    aplatyrhynchos_gene_ensembl     Anas platyrhynchos genes (BGI_duck_1.0)    BGI_duck_1.0
65          saraneus_gene_ensembl               Sorex araneus genes (sorAra1)         sorAra1
66         sharrisii_gene_ensembl       Sarcophilus harrisii genes (DEVIL7.0)        DEVIL7.0
67          meugenii_gene_ensembl           Macropus eugenii genes (Meug_1.0)        Meug_1.0
68           btaurus_gene_ensembl                   Bos taurus genes (UMD3.1)          UMD3.1
69       cfamiliaris_gene_ensembl          Canis familiaris genes (CanFam3.1)       CanFam3.1

> listAttributes(mart)
                                                              name                                              description
1                                                  ensembl_gene_id                                          Ensembl Gene ID
2                                            ensembl_transcript_id                                    Ensembl Transcript ID
3                                               ensembl_peptide_id                                       Ensembl Protein ID
4                                                  ensembl_exon_id                                          Ensembl Exon ID
5                                                      description                                              Description
6                                                  chromosome_name                                          Chromosome Name
7                                                   start_position                                          Gene Start (bp)
8                                                     end_position                                            Gene End (bp)
9                                                           strand                                                   Strand
10                                                            band                                                     Band
11                                                transcript_start                                    Transcript Start (bp)
12                                                  transcript_end                                      Transcript End (bp)
13                                                external_gene_id                                     Associated Gene Name
14                                          external_transcript_id                               Associated Transcript Name
15                                                external_gene_db                                       Associated Gene DB
16                                              transcript_db_name                                 Associated Transcript DB
17                                                transcript_count                                         Transcript count
18                                           percentage_gc_content                                             % GC content
19                                                    gene_biotype                                             Gene Biotype
20                                              transcript_biotype                                       Transcript Biotype
21                                                          source                                            Source (gene)
22                                               transcript_source                                      Source (transcript)
23                                                          status                                            Status (gene)
24                                               transcript_status                                      Status (transcript)
25                                           phenotype_description                                    Phenotype description
26                                                     source_name                                              Source name
27                                               study_external_id                                 Study External Reference
28                                                           go_id                                        GO Term Accession
29                                                       name_1006                                             GO Term Name
30                                                 definition_1006                                       GO Term Definition
31                                                 go_linkage_type                                    GO Term Evidence Code
32                                                  namespace_1003                                                GO domain
33                                            goslim_goa_accession                                  GOSlim GOA Accession(s)
34                                          goslim_goa_description                                   GOSlim GOA Description
35                                                    arrayexpress                                             ArrayExpress
36                                                          chembl                                             ChEMBL ID(s)
37                                   clone_based_ensembl_gene_name                            Clone based Ensembl gene name
38                             clone_based_ensembl_transcript_name                      Clone based Ensembl transcript name
39                                      clone_based_vega_gene_name                               Clone based VEGA gene name
40                                clone_based_vega_transcript_name                         Clone based VEGA transcript name
41                                                            ccds                                                  CCDS ID
42                                                       dbass3_id        Database of Aberrant 3' Splice Sites (DBASS3) IDs
43                                                     dbass3_name                                         DBASS3 Gene Name
44                                                            embl                                        EMBL (Genbank) ID
45                                                     ens_hs_gene                             Ensembl to LRG link gene IDs
46                                               ens_hs_transcript                       Ensembl to LRG link transcript IDs
47                                              ens_hs_translation                      Ensembl to LRG link translation IDs
48                                                    ens_lrg_gene                                 LRG to Ensembl link gene
49                                              ens_lrg_transcript                           LRG to Ensembl link transcript
50                                                      entrezgene                                            EntrezGene ID
51                                                             hpa                          Human Protein Atlas Antibody ID
52                                                            ottg                                   VEGA gene ID(s) (OTTG)
53                                                            ottt                             VEGA transcript ID(s) (OTTT)
54                                            shares_cds_with_ottt      HAVANA transcript (where ENST shares CDS with OTTT)
55                                    shares_cds_and_utr_with_ottt         HAVANA transcript (where ENST identical to OTTT)
56                                                         hgnc_id                                               HGNC ID(s)
57                                                     hgnc_symbol                                              HGNC symbol
58                                            hgnc_transcript_name                                     HGNC transcript name
59                                                          merops                                                MEROPS ID
60                                                             pdb                                                   PDB ID
61                                            mim_morbid_accession                                     MIM Morbid Accession
62                                          mim_morbid_description                                   MIM Morbid Description
63                                              mim_gene_accession                                       MIM Gene Accession
64                                            mim_gene_description                                     MIM Gene Description
65                                               mirbase_accession                                     miRBase Accession(s)
66                                                      mirbase_id                                            miRBase ID(s)
67                                         mirbase_transcript_name                                  miRBase transcript name
68                                                      protein_id                                     Protein (Genbank) ID
69                                                     refseq_mrna                          RefSeq mRNA [e.g. NM_001195597]
70                                           refseq_mrna_predicted                RefSeq mRNA predicted [e.g. XM_001125684]
71                                                    refseq_ncrna                            RefSeq ncRNA [e.g. NR_002834]
72                                          refseq_ncrna_predicted                  RefSeq ncRNA predicted [e.g. XR_108264]
73                                                  refseq_peptide                    RefSeq Protein ID [e.g. NP_001005353]
74                                        refseq_peptide_predicted          RefSeq Predicted Protein ID [e.g. XP_001720922]
75                                                            rfam                                                  Rfam ID
76                                            rfam_transcript_name                                     Rfam transcript name
77                                                            ucsc                                                  UCSC ID
78                                                         unigene                                               Unigene ID
79                                                uniprot_sptrembl                                 UniProt/TrEMBL Accession
80                                               uniprot_swissprot                                     UniProt/SwissProt ID
81                                     uniprot_swissprot_accession                              UniProt/SwissProt Accession
82                                                uniprot_genename                                        UniProt Gene Name
83                                uniprot_genename_transcript_name                         Uniprot Genename Transcript Name
84                                                         uniparc                                                  UniParc
85                                                   wikigene_name                                            WikiGene Name
86                                                     wikigene_id                                              WikiGene ID
87                                            wikigene_description                                     WikiGene Description
88                               efg_agilent_sureprint_g3_ge_8x60k                      Agilent SurePrint G3 GE 8x60k probe
89                            efg_agilent_sureprint_g3_ge_8x60k_v2                   Agilent SurePrint G3 GE 8x60k v2 probe
90                                efg_agilent_wholegenome_4x44k_v1                       Agilent WholeGenome 4x44k v1 probe
91                                efg_agilent_wholegenome_4x44k_v2                       Agilent WholeGenome 4x44k v2 probe
92                                                    affy_hc_g110                                    Affy HC G110 probeset
93                                                   affy_hg_focus                                   Affy HG FOCUS probeset
94                                             affy_hg_u133_plus_2                             Affy HG U133-PLUS-2 probeset
95                                                 affy_hg_u133a_2                                 Affy HG U133A_2 probeset
96                                                   affy_hg_u133a                                   Affy HG U133A probeset
97                                                   affy_hg_u133b                                   Affy HG U133B probeset
98                                                  affy_hg_u95av2                                  Affy HG U95AV2 probeset
99                                                    affy_hg_u95b                                    Affy HG U95B probeset
100                                                   affy_hg_u95c                                    Affy HG U95C probeset
101                                                   affy_hg_u95d                                    Affy HG U95D probeset
102                                                   affy_hg_u95e                                    Affy HG U95E probeset
103                                                   affy_hg_u95a                                    Affy HG U95A probeset
104                                                  affy_hugenefl                                  Affy HuGene FL probeset
105                                            affy_huex_1_0_st_v2                             Affy HuEx 1_0 st v2 probeset
106                                          affy_hugene_1_0_st_v1                           Affy HuGene 1_0 st v1 probeset
107                                          affy_hugene_2_0_st_v1                           Affy HuGene 2_0 st v1 probeset
108                                                 affy_primeview                                           Affy primeview
109                                                  affy_u133_x3p                                   Affy U133 X3P probeset
110                                                agilent_cgh_44b                                    Agilent CGH 44b probe
111                                                       codelink                                           Codelink probe
112                                          illumina_humanwg_6_v1                              Illumina HumanWG 6 v1 probe
113                                          illumina_humanwg_6_v2                              Illumina HumanWG 6 v2 probe
114                                          illumina_humanwg_6_v3                              Illumina HumanWG 6 v3 probe
115                                         illumina_humanht_12_v3                           Illumina Human HT 12 V3 probe 
116                                         illumina_humanht_12_v4                           Illumina Human HT 12 V4 probe 
117                                         illumina_humanref_8_v3                            Illumina Human Ref 8 V3 probe
118                                               phalanx_onearray                                   Phalanx OneArray probe
119                                              anatomical_system                            Anatomical System (egenetics)
120                                              development_stage                            Development Stage (egenetics)
121                                                      cell_type                                    Cell Type (egenetics)
122                                                      pathology                                    Pathology (egenetics)
123                                                 atlas_celltype                                      GNF/Atlas cell type
124                                             atlas_diseasestate                                  GNF/Atlas disease state
125                                             atlas_organismpart                                  GNF/Atlas organism part
126                                             family_description                               Ensembl Family Description
127                                                         family                             Ensembl Protein Family ID(s)
128                                                          pirsf                                     PIRSF SuperFamily ID
129                                                    superfamily                                           Superfamily ID
130                                                          smart                                                 SMART ID
131                                                        profile                                               PROFILE ID
132                                                         prints                                                PRINTS ID
133                                                           pfam                                                  PFAM ID
134                                                        tigrfam                                               TIGRFam ID
135                                                       interpro                                              Interpro ID
136                                     interpro_short_description                               Interpro Short Description
137                                           interpro_description                                     Interpro Description
138                                                 low_complexity                                           Low complexity
139                                           transmembrane_domain                                     Transmembrane domain
140                                                  signal_domain                                            Signal domain
141                                                         ncoils                                                   Ncoils
142                                                ensembl_gene_id                                          Ensembl Gene ID
143                                          ensembl_transcript_id                                    Ensembl Transcript ID
144                                             ensembl_peptide_id                                       Ensembl Protein ID
145                                                chromosome_name                                          Chromosome Name
146                                                 start_position                                          Gene Start (bp)
147                                                   end_position                                            Gene End (bp)
148                                               transcript_start                                    Transcript Start (bp)
149                                                 transcript_end                                      Transcript End (bp)
150                                                         strand                                                   Strand
151                                               external_gene_id                                     Associated Gene Name
152                                               external_gene_db                                       Associated Gene DB
153                                                    5_utr_start                                             5' UTR Start
154                                                      5_utr_end                                               5' UTR End
155                                                    3_utr_start                                             3' UTR Start
156                                                      3_utr_end                                               3' UTR End
157                                                     cds_length                                               CDS Length
158                                               transcript_count                                         Transcript count
159                                                    description                                              Description
160                                                   gene_biotype                                             Gene Biotype
161                                               exon_chrom_start                                      Exon Chr Start (bp)
162                                                 exon_chrom_end                                        Exon Chr End (bp)
163                                                is_constitutive                                        Constitutive Exon
164                                                           rank                                  Exon Rank in Transcript
165                                                          phase                                                    phase
166                                              cdna_coding_start                                        cDNA coding start
167                                                cdna_coding_end                                          cDNA coding end
168                                           genomic_coding_start                                     Genomic coding start
169                                             genomic_coding_end                                       Genomic coding end
170                                                ensembl_exon_id                                          Ensembl Exon ID
171                                                      cds_start                                                CDS Start
172                                                        cds_end                                                  CDS End
173                                                ensembl_gene_id                                          Ensembl Gene ID
174                                          ensembl_transcript_id                                    Ensembl Transcript ID
175                                             ensembl_peptide_id                                       Ensembl Protein ID
176                                                      name_1078                       Gene Name With Corresponding Event
177                                   splicing_event__dm_name_1059                                          Chromosome Name
178                                            splicing_event_type                                               Event Type
179                                                       name_106                                               Event Name
180                                          seq_region_start_1078                                    Seq Region Start (bp)
181                                            seq_region_end_1078                                      Seq Region End (bp)
182                                         seq_region_strand_1078                                                   Strand
183                                                ensembl_gene_id                                          Ensembl Gene ID
184                                          ensembl_transcript_id                                    Ensembl Transcript ID
185                                             ensembl_peptide_id                                       Ensembl Protein ID
186                                                chromosome_name                                          Chromosome Name
187                                                 start_position                                          Gene Start (bp)
188                                                   end_position                                            Gene End (bp)
189                                                         strand                                                   Strand
190                                                           band                                                     Band
191                                               external_gene_id                                     Associated Gene Name
192                                               external_gene_db                                       Associated Gene DB
193                                               transcript_count                                         Transcript count
194                                          percentage_gc_content                                             % GC content
195                                                    description                                              Description
196                                    vpacos_homolog_ensembl_gene                                   Alpaca Ensembl Gene ID
197                    vpacos_homolog_canomical_transcript_protein                       Canonical Protein or Transcript ID
198                                 vpacos_homolog_ensembl_peptide                                Alpaca Ensembl Protein ID
199                                      vpacos_homolog_chromosome                                   Alpaca Chromosome Name
200                                     vpacos_homolog_chrom_start                             Alpaca Chromosome Start (bp)
201                                       vpacos_homolog_chrom_end                               Alpaca Chromosome End (bp)
202                                  vpacos_homolog_orthology_type                                            Homology Type
203                                         vpacos_homolog_subtype                                                 Ancestor
204                            vpacos_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
205                                         vpacos_homolog_perc_id                    % Identity with respect to query gene
206                                      vpacos_homolog_perc_id_r1                   % Identity with respect to Alpaca gene
207                                              vpacos_homolog_dn                                                       dN
208                                              vpacos_homolog_ds                                                       dS
209                             acarolinensis_homolog_ensembl_gene                             Anole Lizard Ensembl Gene ID
210             acarolinensis_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
211                          acarolinensis_homolog_ensembl_peptide                          Anole Lizard Ensembl Protein ID
212                               acarolinensis_homolog_chromosome                             Anole Lizard Chromosome Name
213                              acarolinensis_homolog_chrom_start                       Anole Lizard Chromosome Start (bp)
214                                acarolinensis_homolog_chrom_end                         Anole Lizard Chromosome End (bp)
215                           acarolinensis_homolog_orthology_type                                            Homology Type
216                                  acarolinensis_homolog_subtype                                                 Ancestor
217                     acarolinensis_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
218                                  acarolinensis_homolog_perc_id                    % Identity with respect to query gene
219                               acarolinensis_homolog_perc_id_r1             % Identity with respect to Anole Lizard gene
220                                       acarolinensis_homolog_dn                                                       dN
221                                       acarolinensis_homolog_ds                                                       dS
222                             dnovemcinctus_homolog_ensembl_gene                                Armadillo Ensembl Gene ID
223             dnovemcinctus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
224                          dnovemcinctus_homolog_ensembl_peptide                             Armadillo Ensembl Protein ID
225                               dnovemcinctus_homolog_chromosome                                Armadillo Chromosome Name
226                              dnovemcinctus_homolog_chrom_start                          Armadillo Chromosome Start (bp)
227                                dnovemcinctus_homolog_chrom_end                            Armadillo Chromosome End (bp)
228                           dnovemcinctus_homolog_orthology_type                                            Homology Type
229                                  dnovemcinctus_homolog_subtype                                                 Ancestor
230                     dnovemcinctus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
231                                  dnovemcinctus_homolog_perc_id                    % Identity with respect to query gene
232                               dnovemcinctus_homolog_perc_id_r1                % Identity with respect to Armadillo gene
233                                       dnovemcinctus_homolog_dn                                                       dN
234                                       dnovemcinctus_homolog_ds                                                       dS
235                                   gmorhua_homolog_ensembl_gene                             Atlantic Cod Ensembl Gene ID
236                   gmorhua_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
237                                gmorhua_homolog_ensembl_protein                          Atlantic Cod Ensembl Protein ID
238                                     gmorhua_homolog_chromosome                             Atlantic Cod Chromosome Name
239                                    gmorhua_homolog_chrom_start                       Atlantic Cod Chromosome Start (bp)
240                                      gmorhua_homolog_chrom_end                         Atlantic Cod Chromosome End (bp)
241                                 gmorhua_homolog_orthology_type                                            Homology Type
242                                        gmorhua_homolog_subtype                                                 Ancestor
243                           gmorhua_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
244                                        gmorhua_homolog_perc_id                    % Identity with respect to query gene
245                                     gmorhua_homolog_perc_id_r1             % Identity with respect to Atlantic Cod gene
246                                ogarnettii_homolog_ensembl_gene                                 Bushbaby Ensembl Gene ID
247                ogarnettii_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
248                             ogarnettii_homolog_ensembl_peptide                              Bushbaby Ensembl Protein ID
249                                  ogarnettii_homolog_chromosome                                 Bushbaby Chromosome Name
250                                 ogarnettii_homolog_chrom_start                           Bushbaby Chromosome Start (bp)
251                                   ogarnettii_homolog_chrom_end                             Bushbaby Chromosome End (bp)
252                              ogarnettii_homolog_orthology_type                                            Homology Type
253                                     ogarnettii_homolog_subtype                                                 Ancestor
254                        ogarnettii_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
255                                     ogarnettii_homolog_perc_id                    % Identity with respect to query gene
256                                  ogarnettii_homolog_perc_id_r1                 % Identity with respect to Bushbaby gene
257                                          ogarnettii_homolog_dn                                                       dN
258                                          ogarnettii_homolog_ds                                                       dS
259                                   celegans_homolog_chrom_start             Caenorhabditis elegans Chromosome Start (bp)
260                                  celegans_homolog_ensembl_gene                   Caenorhabditis elegans Ensembl Gene ID
261                  celegans_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
262                               celegans_homolog_ensembl_peptide                Caenorhabditis elegans Ensembl Protein ID
263                                    celegans_homolog_chromosome                   Caenorhabditis elegans Chromosome Name
264                                     celegans_homolog_chrom_end               Caenorhabditis elegans Chromosome End (bp)
265                                celegans_homolog_orthology_type                                            Homology Type
266                                       celegans_homolog_subtype                                                 Ancestor
267                          celegans_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
268                                       celegans_homolog_perc_id                    % Identity with respect to query gene
269                                    celegans_homolog_perc_id_r1   % Identity with respect to Caenorhabditis elegans gene
270                                            celegans_homolog_dn                                                       dN
271                                            celegans_homolog_ds                                                       dS
272                                    fcatus_homolog_ensembl_gene                                      Cat Ensembl Gene ID
273                    fcatus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
274                                 fcatus_homolog_ensembl_peptide                                   Cat Ensembl Protein ID
275                                      fcatus_homolog_chromosome                                      Cat Chromosome Name
276                                     fcatus_homolog_chrom_start                                Cat Chromosome Start (bp)
277                                       fcatus_homolog_chrom_end                                  Cat Chromosome End (bp)
278                                  fcatus_homolog_orthology_type                                            Homology Type
279                                         fcatus_homolog_subtype                                                 Ancestor
280                            fcatus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
281                                         fcatus_homolog_perc_id                    % Identity with respect to query gene
282                                      fcatus_homolog_perc_id_r1                      % Identity with respect to Cat gene
283                                              fcatus_homolog_dn                                                       dN
284                                              fcatus_homolog_ds                                                       dS
285                                amexicanus_homolog_ensembl_gene                                Cave fish Ensembl Gene ID
286                amexicanus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
287                             amexicanus_homolog_ensembl_peptide                             Cave fish Ensembl Protein ID
288                                  amexicanus_homolog_chromosome                                Cave fish Chromosome Name
289                                 amexicanus_homolog_chrom_start                          Cave fish Chromosome Start (bp)
290                                   amexicanus_homolog_chrom_end                            Cave fish Chromosome End (bp)
291                              amexicanus_homolog_orthology_type                                            Homology Type
292                                     amexicanus_homolog_subtype                                                 Ancestor
293                        amexicanus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
294                                     amexicanus_homolog_perc_id                    % Identity with respect to query gene
295                                  amexicanus_homolog_perc_id_r1                % Identity with respect to Cave fish gene
296                                   ggallus_homolog_ensembl_gene                                  Chicken Ensembl Gene ID
297                   ggallus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
298                                ggallus_homolog_ensembl_peptide                               Chicken Ensembl Protein ID
299                                     ggallus_homolog_chromosome                                  Chicken Chromosome Name
300                                    ggallus_homolog_chrom_start                            Chicken Chromosome Start (bp)
301                                      ggallus_homolog_chrom_end                              Chicken Chromosome End (bp)
302                                 ggallus_homolog_orthology_type                                            Homology Type
303                                        ggallus_homolog_subtype                                                 Ancestor
304                           ggallus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
305                                        ggallus_homolog_perc_id                    % Identity with respect to query gene
306                                     ggallus_homolog_perc_id_r1                  % Identity with respect to Chicken gene
307                                             ggallus_homolog_dn                                                       dN
308                                             ggallus_homolog_ds                                                       dS
309                              ptroglodytes_homolog_ensembl_gene                               Chimpanzee Ensembl Gene ID
310              ptroglodytes_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
311                           ptroglodytes_homolog_ensembl_peptide                            Chimpanzee Ensembl Protein ID
312                                ptroglodytes_homolog_chromosome                               Chimpanzee Chromosome Name
313                               ptroglodytes_homolog_chrom_start                         Chimpanzee Chromosome Start (bp)
314                                 ptroglodytes_homolog_chrom_end                           Chimpanzee Chromosome End (bp)
315                            ptroglodytes_homolog_orthology_type                                            Homology Type
316                                   ptroglodytes_homolog_subtype                                                 Ancestor
317                      ptroglodytes_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
318                                   ptroglodytes_homolog_perc_id                    % Identity with respect to query gene
319                                ptroglodytes_homolog_perc_id_r1               % Identity with respect to Chimpanzee gene
320                                        ptroglodytes_homolog_dn                                                       dN
321                                        ptroglodytes_homolog_ds                                                       dS
322                                 psinensis_homolog_ensembl_gene                 Chinese softshell turtle Ensembl Gene ID
323                 psinensis_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
324                              psinensis_homolog_ensembl_peptide              Chinese softshell turtle Ensembl Protein ID
325                                   psinensis_homolog_chromosome                 Chinese softshell turtle Chromosome Name
326                                  psinensis_homolog_chrom_start           Chinese softshell turtle Chromosome Start (bp)
327                                    psinensis_homolog_chrom_end             Chinese softshell turtle Chromosome End (bp)
328                               psinensis_homolog_orthology_type                                            Homology Type
329                                      psinensis_homolog_subtype                                                 Ancestor
330                         psinensis_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
331                                      psinensis_homolog_perc_id                    % Identity with respect to query gene
332                                   psinensis_homolog_perc_id_r1 % Identity with respect to Chinese softshell turtle gene
333                                           psinensis_homolog_dn                                                       dN
334                                           psinensis_homolog_ds                                                       dS
335                             cintestinalis_homolog_ensembl_gene                       Ciona intestinalis Ensembl Gene ID
336             cintestinalis_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
337                          cintestinalis_homolog_ensembl_peptide                    Ciona intestinalis Ensembl Protein ID
338                               cintestinalis_homolog_chromosome                       Ciona intestinalis Chromosome Name
339                              cintestinalis_homolog_chrom_start                 Ciona intestinalis Chromosome Start (bp)
340                                cintestinalis_homolog_chrom_end                   Ciona intestinalis Chromosome End (bp)
341                           cintestinalis_homolog_orthology_type                                            Homology Type
342                                  cintestinalis_homolog_subtype                                                 Ancestor
343                     cintestinalis_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
344                                  cintestinalis_homolog_perc_id                    % Identity with respect to query gene
345                               cintestinalis_homolog_perc_id_r1       % Identity with respect to Ciona intestinalis gene
346                                       cintestinalis_homolog_dn                                                       dN
347                                       cintestinalis_homolog_ds                                                       dS
348                                 csavignyi_homolog_ensembl_gene                           Ciona savignyi Ensembl Gene ID
349                 csavignyi_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
350                              csavignyi_homolog_ensembl_peptide                        Ciona savignyi Ensembl Protein ID
351                                   csavignyi_homolog_chromosome                           Ciona savignyi Chromosome Name
352                                  csavignyi_homolog_chrom_start                     Ciona savignyi Chromosome Start (bp)
353                                    csavignyi_homolog_chrom_end                       Ciona savignyi Chromosome End (bp)
354                               csavignyi_homolog_orthology_type                                            Homology Type
355                                      csavignyi_homolog_subtype                                                 Ancestor
356                         csavignyi_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
357                                      csavignyi_homolog_perc_id                    % Identity with respect to query gene
358                                   csavignyi_homolog_perc_id_r1           % Identity with respect to Ciona savignyi gene
359                                           csavignyi_homolog_dn                                                       dN
360                                           csavignyi_homolog_ds                                                       dS
361                                lchalumnae_homolog_ensembl_gene                               Coelacanth Ensembl Gene ID
362                lchalumnae_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
363                             lchalumnae_homolog_ensembl_peptide                            Coelacanth Ensembl Protein ID
364                                  lchalumnae_homolog_chromosome                               Coelacanth Chromosome Name
365                                 lchalumnae_homolog_chrom_start                         Coelacanth Chromosome Start (bp)
366                                   lchalumnae_homolog_chrom_end                           Coelacanth Chromosome End (bp)
367                              lchalumnae_homolog_orthology_type                                            Homology Type
368                                     lchalumnae_homolog_subtype                                                 Ancestor
369                        lchalumnae_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
370                                     lchalumnae_homolog_perc_id                    % Identity with respect to query gene
371                                  lchalumnae_homolog_perc_id_r1               % Identity with respect to Coelacanth gene
372                                  saraneus_homolog_ensembl_gene                             Common Shrew Ensembl Gene ID
373                  saraneus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
374                               saraneus_homolog_ensembl_peptide                          Common Shrew Ensembl Protein ID
375                                    saraneus_homolog_chromosome                             Common Shrew Chromosome Name
376                                   saraneus_homolog_chrom_start                       Common Shrew Chromosome Start (bp)
377                                     saraneus_homolog_chrom_end                         Common Shrew Chromosome End (bp)
378                                saraneus_homolog_orthology_type                                            Homology Type
379                                       saraneus_homolog_subtype                                                 Ancestor
380                          saraneus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
381                                       saraneus_homolog_perc_id                    % Identity with respect to query gene
382                                    saraneus_homolog_perc_id_r1             % Identity with respect to Common Shrew gene
383                                            saraneus_homolog_dn                                                       dN
384                                            saraneus_homolog_ds                                                       dS
385                                   btaurus_homolog_ensembl_gene                                      Cow Ensembl Gene ID
386                   btaurus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
387                                btaurus_homolog_ensembl_peptide                                   Cow Ensembl Protein ID
388                                     btaurus_homolog_chromosome                                      Cow Chromosome Name
389                                    btaurus_homolog_chrom_start                                Cow Chromosome Start (bp)
390                                      btaurus_homolog_chrom_end                                  Cow Chromosome End (bp)
391                                 btaurus_homolog_orthology_type                                            Homology Type
392                                        btaurus_homolog_subtype                                                 Ancestor
393                           btaurus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
394                                        btaurus_homolog_perc_id                    % Identity with respect to query gene
395                                     btaurus_homolog_perc_id_r1                      % Identity with respect to Cow gene
396                                             btaurus_homolog_dn                                                       dN
397                                             btaurus_homolog_ds                                                       dS
398                               cfamiliaris_homolog_ensembl_gene                                      Dog Ensembl Gene ID
399               cfamiliaris_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
400                            cfamiliaris_homolog_ensembl_peptide                                   Dog Ensembl Protein ID
401                                 cfamiliaris_homolog_chromosome                                      Dog Chromosome Name
402                                cfamiliaris_homolog_chrom_start                                Dog Chromosome Start (bp)
403                                  cfamiliaris_homolog_chrom_end                                  Dog Chromosome End (bp)
404                             cfamiliaris_homolog_orthology_type                                            Homology Type
405                                    cfamiliaris_homolog_subtype                                                 Ancestor
406                       cfamiliaris_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
407                                    cfamiliaris_homolog_perc_id                    % Identity with respect to query gene
408                                 cfamiliaris_homolog_perc_id_r1                      % Identity with respect to Dog gene
409                                         cfamiliaris_homolog_dn                                                       dN
410                                         cfamiliaris_homolog_ds                                                       dS
411                                ttruncatus_homolog_ensembl_gene                                  Dolphin Ensembl Gene ID
412                ttruncatus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
413                             ttruncatus_homolog_ensembl_peptide                               Dolphin Ensembl Protein ID
414                                  ttruncatus_homolog_chromosome                                  Dolphin Chromosome Name
415                                 ttruncatus_homolog_chrom_start                            Dolphin Chromosome Start (bp)
416                                   ttruncatus_homolog_chrom_end                              Dolphin Chromosome End (bp)
417                              ttruncatus_homolog_orthology_type                                            Homology Type
418                                     ttruncatus_homolog_subtype                                                 Ancestor
419                        ttruncatus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
420                                     ttruncatus_homolog_perc_id                    % Identity with respect to query gene
421                                  ttruncatus_homolog_perc_id_r1                  % Identity with respect to Dolphin gene
422                                          ttruncatus_homolog_dn                                                       dN
423                                          ttruncatus_homolog_ds                                                       dS
424                             dmelanogaster_homolog_ensembl_gene                               Drosophila Ensembl Gene ID
425             dmelanogaster_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
426                          dmelanogaster_homolog_ensembl_peptide                            Drosophila Ensembl Protein ID
427                               dmelanogaster_homolog_chromosome                               Drosophila Chromosome Name
428                              dmelanogaster_homolog_chrom_start                         Drosophila Chromosome Start (bp)
429                                dmelanogaster_homolog_chrom_end                           Drosophila Chromosome End (bp)
430                           dmelanogaster_homolog_orthology_type                                            Homology Type
431                                  dmelanogaster_homolog_subtype                                                 Ancestor
432                     dmelanogaster_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
433                                  dmelanogaster_homolog_perc_id                    % Identity with respect to query gene
434                               dmelanogaster_homolog_perc_id_r1               % Identity with respect to Drosophila gene
435                                       dmelanogaster_homolog_dn                                                       dN
436                                       dmelanogaster_homolog_ds                                                       dS
437                            aplatyrhynchos_homolog_ensembl_gene                                     Duck Ensembl Gene ID
438            aplatyrhynchos_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
439                         aplatyrhynchos_homolog_ensembl_peptide                                  Duck Ensembl Protein ID
440                              aplatyrhynchos_homolog_chromosome                                     Duck Chromosome Name
441                             aplatyrhynchos_homolog_chrom_start                               Duck Chromosome Start (bp)
442                               aplatyrhynchos_homolog_chrom_end                                 Duck Chromosome End (bp)
443                          aplatyrhynchos_homolog_orthology_type                                            Homology Type
444                                 aplatyrhynchos_homolog_subtype                                                 Ancestor
445                    aplatyrhynchos_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
446                                 aplatyrhynchos_homolog_perc_id                    % Identity with respect to query gene
447                              aplatyrhynchos_homolog_perc_id_r1                     % Identity with respect to Duck gene
448                                      aplatyrhynchos_homolog_dn                                                       dN
449                                      aplatyrhynchos_homolog_ds                                                       dS
450                                 lafricana_homolog_ensembl_gene                                 Elephant Ensembl Gene ID
451                 lafricana_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
452                              lafricana_homolog_ensembl_peptide                              Elephant Ensembl Protein ID
453                                   lafricana_homolog_chromosome                                 Elephant Chromosome Name
454                                  lafricana_homolog_chrom_start                           Elephant Chromosome Start (bp)
455                                    lafricana_homolog_chrom_end                             Elephant Chromosome End (bp)
456                               lafricana_homolog_orthology_type                                            Homology Type
457                                      lafricana_homolog_subtype                                                 Ancestor
458                         lafricana_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
459                                      lafricana_homolog_perc_id                    % Identity with respect to query gene
460                                   lafricana_homolog_perc_id_r1                 % Identity with respect to Elephant gene
461                                           lafricana_homolog_dn                                                       dN
462                                           lafricana_homolog_ds                                                       dS
463                                     mfuro_homolog_ensembl_gene                                   Ferret Ensembl Gene ID
464                     mfuro_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
465                                  mfuro_homolog_ensembl_peptide                                Ferret Ensembl Protein ID
466                                       mfuro_homolog_chromosome                                   Ferret Chromosome Name
467                                      mfuro_homolog_chrom_start                             Ferret Chromosome Start (bp)
468                                        mfuro_homolog_chrom_end                               Ferret Chromosome End (bp)
469                                   mfuro_homolog_orthology_type                                            Homology Type
470                                          mfuro_homolog_subtype                                                 Ancestor
471                             mfuro_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
472                                          mfuro_homolog_perc_id                    % Identity with respect to query gene
473                                       mfuro_homolog_perc_id_r1                   % Identity with respect to Ferret gene
474                                               mfuro_homolog_dn                                                       dN
475                                               mfuro_homolog_ds                                                       dS
476                               falbicollis_homolog_ensembl_gene                               Flycatcher Ensembl Gene ID
477               falbicollis_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
478                            falbicollis_homolog_ensembl_peptide                            Flycatcher Ensembl Protein ID
479                                 falbicollis_homolog_chromosome                               Flycatcher Chromosome Name
480                                falbicollis_homolog_chrom_start                         Flycatcher Chromosome Start (bp)
481                                  falbicollis_homolog_chrom_end                           Flycatcher Chromosome End (bp)
482                             falbicollis_homolog_orthology_type                                            Homology Type
483                                    falbicollis_homolog_subtype                                                 Ancestor
484                       falbicollis_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
485                                    falbicollis_homolog_perc_id                    % Identity with respect to query gene
486                                 falbicollis_homolog_perc_id_r1               % Identity with respect to Flycatcher gene
487                                         falbicollis_homolog_dn                                                       dN
488                                         falbicollis_homolog_ds                                                       dS
489                                 trubripes_homolog_ensembl_gene                                     Fugu Ensembl Gene ID
490                 trubripes_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
491                              trubripes_homolog_ensembl_peptide                                  Fugu Ensembl Protein ID
492                                   trubripes_homolog_chromosome                                     Fugu Chromosome Name
493                                  trubripes_homolog_chrom_start                               Fugu Chromosome Start (bp)
494                                    trubripes_homolog_chrom_end                                 Fugu Chromosome End (bp)
495                               trubripes_homolog_orthology_type                                            Homology Type
496                                      trubripes_homolog_subtype                                                 Ancestor
497                         trubripes_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
498                                      trubripes_homolog_perc_id                    % Identity with respect to query gene
499                                   trubripes_homolog_perc_id_r1                     % Identity with respect to Fugu gene
500                                           trubripes_homolog_dn                                                       dN
501                                           trubripes_homolog_ds                                                       dS
502                               nleucogenys_homolog_ensembl_gene                                   Gibbon Ensembl Gene ID
503               nleucogenys_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
504                            nleucogenys_homolog_ensembl_peptide                                Gibbon Ensembl Protein ID
505                                 nleucogenys_homolog_chromosome                                   Gibbon Chromosome Name
506                                nleucogenys_homolog_chrom_start                             Gibbon Chromosome Start (bp)
507                                  nleucogenys_homolog_chrom_end                               Gibbon Chromosome End (bp)
508                             nleucogenys_homolog_orthology_type                                            Homology Type
509                                    nleucogenys_homolog_subtype                                                 Ancestor
510                       nleucogenys_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
511                                    nleucogenys_homolog_perc_id                    % Identity with respect to query gene
512                                 nleucogenys_homolog_perc_id_r1                   % Identity with respect to Gibbon gene
513                                         nleucogenys_homolog_dn                                                       dN
514                                         nleucogenys_homolog_ds                                                       dS
515                                  ggorilla_homolog_ensembl_gene                                  Gorilla Ensembl Gene ID
516                  ggorilla_homolog_canomical_transcript_protein                       Canonical Protein or Transcript ID
517                               ggorilla_homolog_ensembl_peptide                               Gorilla Ensembl Protein ID
518                                    ggorilla_homolog_chromosome                                  Gorilla Chromosome Name
519                                   ggorilla_homolog_chrom_start                            Gorilla Chromosome Start (bp)
520                                     ggorilla_homolog_chrom_end                              Gorilla Chromosome End (bp)
521                                ggorilla_homolog_orthology_type                                            Homology Type
522                                       ggorilla_homolog_subtype                                                 Ancestor
523                          ggorilla_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
524                                       ggorilla_homolog_perc_id                    % Identity with respect to query gene
525                                    ggorilla_homolog_perc_id_r1                  % Identity with respect to Gorilla gene
526                                            ggorilla_homolog_dn                                                       dN
527                                            ggorilla_homolog_ds                                                       dS
528                                cporcellus_homolog_ensembl_gene                               Guinea Pig Ensembl Gene ID
529                cporcellus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
530                             cporcellus_homolog_ensembl_peptide                            Guinea Pig Ensembl Protein ID
531                                  cporcellus_homolog_chromosome                               Guinea Pig Chromosome Name
532                                 cporcellus_homolog_chrom_start                         Guinea Pig Chromosome Start (bp)
533                                   cporcellus_homolog_chrom_end                           Guinea Pig Chromosome End (bp)
534                              cporcellus_homolog_orthology_type                                            Homology Type
535                                     cporcellus_homolog_subtype                                                 Ancestor
536                        cporcellus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
537                                     cporcellus_homolog_perc_id                    % Identity with respect to query gene
538                                  cporcellus_homolog_perc_id_r1               % Identity with respect to Guinea Pig gene
539                                          cporcellus_homolog_dn                                                       dN
540                                          cporcellus_homolog_ds                                                       dS
541                                eeuropaeus_homolog_ensembl_gene                                 Hedgehog Ensembl Gene ID
542                eeuropaeus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
543                             eeuropaeus_homolog_ensembl_peptide                              Hedgehog Ensembl Protein ID
544                                  eeuropaeus_homolog_chromosome                                 Hedgehog Chromosome Name
545                                 eeuropaeus_homolog_chrom_start                           Hedgehog Chromosome Start (bp)
546                                   eeuropaeus_homolog_chrom_end                             Hedgehog Chromosome End (bp)
547                              eeuropaeus_homolog_orthology_type                                            Homology Type
548                                     eeuropaeus_homolog_subtype                                                 Ancestor
549                        eeuropaeus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
550                                     eeuropaeus_homolog_perc_id                    % Identity with respect to query gene
551                                  eeuropaeus_homolog_perc_id_r1                 % Identity with respect to Hedgehog gene
552                                          eeuropaeus_homolog_dn                                                       dN
553                                          eeuropaeus_homolog_ds                                                       dS
554                                 ecaballus_homolog_ensembl_gene                                    Horse Ensembl Gene ID
555                 ecaballus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
556                              ecaballus_homolog_ensembl_peptide                                 Horse Ensembl Protein ID
557                                   ecaballus_homolog_chromosome                                    Horse Chromosome Name
558                                  ecaballus_homolog_chrom_start                              Horse Chromosome Start (bp)
559                                    ecaballus_homolog_chrom_end                                Horse Chromosome End (bp)
560                               ecaballus_homolog_orthology_type                                            Homology Type
561                                      ecaballus_homolog_subtype                                                 Ancestor
562                         ecaballus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
563                                      ecaballus_homolog_perc_id                    % Identity with respect to query gene
564                                   ecaballus_homolog_perc_id_r1                    % Identity with respect to Horse gene
565                                           ecaballus_homolog_dn                                                       dN
566                                           ecaballus_homolog_ds                                                       dS
567                                    dordii_homolog_ensembl_gene                             Kangaroo Rat Ensembl Gene ID
568                    dordii_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
569                                 dordii_homolog_ensembl_peptide                          Kangaroo Rat Ensembl Protein ID
570                                      dordii_homolog_chromosome                             Kangaroo Rat Chromosome Name
571                                     dordii_homolog_chrom_start                       Kangaroo Rat Chromosome Start (bp)
572                                       dordii_homolog_chrom_end                         Kangaroo Rat Chromosome End (bp)
573                                  dordii_homolog_orthology_type                                            Homology Type
574                                         dordii_homolog_subtype                                                 Ancestor
575                            dordii_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
576                                         dordii_homolog_perc_id                    % Identity with respect to query gene
577                                      dordii_homolog_perc_id_r1             % Identity with respect to Kangaroo Rat gene
578                                              dordii_homolog_dn                                                       dN
579                                              dordii_homolog_ds                                                       dS
580                                  pmarinus_homolog_ensembl_gene                                  Lamprey Ensembl Gene ID
581                  pmarinus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
582                               pmarinus_homolog_ensembl_peptide                               Lamprey Ensembl Protein ID
583                                    pmarinus_homolog_chromosome                                  Lamprey Chromosome Name
584                                   pmarinus_homolog_chrom_start                            Lamprey Chromosome Start (bp)
585                                     pmarinus_homolog_chrom_end                              Lamprey Chromosome End (bp)
586                                pmarinus_homolog_orthology_type                                            Homology Type
587                                       pmarinus_homolog_subtype                                                 Ancestor
588                          pmarinus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
589                                       pmarinus_homolog_perc_id                    % Identity with respect to query gene
590                                    pmarinus_homolog_perc_id_r1                  % Identity with respect to Lamprey gene
591                                 etelfairi_homolog_ensembl_gene                   Lesser hedgehog tenrec Ensembl Gene ID
592                 etelfairi_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
593                              etelfairi_homolog_ensembl_peptide                Lesser hedgehog tenrec Ensembl Protein ID
594                                   etelfairi_homolog_chromosome                   Lesser hedgehog tenrec Chromosome Name
595                                  etelfairi_homolog_chrom_start             Lesser hedgehog tenrec Chromosome Start (bp)
596                                    etelfairi_homolog_chrom_end               Lesser hedgehog tenrec Chromosome End (bp)
597                               etelfairi_homolog_orthology_type                                            Homology Type
598                                      etelfairi_homolog_subtype                                                 Ancestor
599                         etelfairi_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
600                                      etelfairi_homolog_perc_id                    % Identity with respect to query gene
601                                   etelfairi_homolog_perc_id_r1   % Identity with respect to Lesser hedgehog tenrec gene
602                                           etelfairi_homolog_dn                                                       dN
603                                           etelfairi_homolog_ds                                                       dS
604                                  mmulatta_homolog_ensembl_gene                                  Macaque Ensembl Gene ID
605                  mmulatta_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
606                               mmulatta_homolog_ensembl_peptide                               Macaque Ensembl Protein ID
607                                    mmulatta_homolog_chromosome                                  Macaque Chromosome Name
608                                   mmulatta_homolog_chrom_start                            Macaque Chromosome Start (bp)
609                                     mmulatta_homolog_chrom_end                              Macaque Chromosome End (bp)
610                                mmulatta_homolog_orthology_type                                            Homology Type
611                                       mmulatta_homolog_subtype                                                 Ancestor
612                          mmulatta_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
613                                       mmulatta_homolog_perc_id                    % Identity with respect to query gene
614                                    mmulatta_homolog_perc_id_r1                  % Identity with respect to Macaque gene
615                                            mmulatta_homolog_dn                                                       dN
616                                            mmulatta_homolog_ds                                                       dS
617                                  cjacchus_homolog_ensembl_gene                                 Marmoset Ensembl Gene ID
618                  cjacchus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
619                               cjacchus_homolog_ensembl_peptide                              Marmoset Ensembl Protein ID
620                                    cjacchus_homolog_chromosome                                 Marmoset Chromosome Name
621                                   cjacchus_homolog_chrom_start                           Marmoset Chromosome Start (bp)
622                                     cjacchus_homolog_chrom_end                             Marmoset Chromosome End (bp)
623                                cjacchus_homolog_orthology_type                                            Homology Type
624                                       cjacchus_homolog_subtype                                                 Ancestor
625                          cjacchus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
626                                       cjacchus_homolog_perc_id                    % Identity with respect to query gene
627                                    cjacchus_homolog_perc_id_r1                 % Identity with respect to Marmoset gene
628                                            cjacchus_homolog_dn                                                       dN
629                                            cjacchus_homolog_ds                                                       dS
630                                  olatipes_homolog_ensembl_gene                                   Medaka Ensembl Gene ID
631                  olatipes_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
632                               olatipes_homolog_ensembl_peptide                                Medaka Ensembl Protein ID
633                                    olatipes_homolog_chromosome                                   Medaka Chromosome Name
634                                   olatipes_homolog_chrom_start                             Medaka Chromosome Start (bp)
635                                     olatipes_homolog_chrom_end                               Medaka Chromosome End (bp)
636                                olatipes_homolog_orthology_type                                            Homology Type
637                                       olatipes_homolog_subtype                                                 Ancestor
638                          olatipes_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
639                                       olatipes_homolog_perc_id                    % Identity with respect to query gene
640                                    olatipes_homolog_perc_id_r1                   % Identity with respect to Medaka gene
641                                            olatipes_homolog_dn                                                       dN
642                                            olatipes_homolog_ds                                                       dS
643                                 pvampyrus_homolog_ensembl_gene                                  Megabat Ensembl Gene ID
644                 pvampyrus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
645                              pvampyrus_homolog_ensembl_peptide                               Megabat Ensembl Protein ID
646                                   pvampyrus_homolog_chromosome                                  Megabat Chromosome Name
647                                  pvampyrus_homolog_chrom_start                            Megabat Chromosome Start (bp)
648                                    pvampyrus_homolog_chrom_end                              Megabat Chromosome End (bp)
649                               pvampyrus_homolog_orthology_type                                            Homology Type
650                                      pvampyrus_homolog_subtype                                                 Ancestor
651                         pvampyrus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
652                                      pvampyrus_homolog_perc_id                    % Identity with respect to query gene
653                                   pvampyrus_homolog_perc_id_r1                  % Identity with respect to Megabat gene
654                                           pvampyrus_homolog_dn                                                       dN
655                                           pvampyrus_homolog_ds                                                       dS
656                                mlucifugus_homolog_ensembl_gene                                 Microbat Ensembl Gene ID
657                mlucifugus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
658                             mlucifugus_homolog_ensembl_peptide                              Microbat Ensembl Protein ID
659                                  mlucifugus_homolog_chromosome                                 Microbat Chromosome Name
660                                 mlucifugus_homolog_chrom_start                           Microbat Chromosome Start (bp)
661                                   mlucifugus_homolog_chrom_end                             Microbat Chromosome End (bp)
662                              mlucifugus_homolog_orthology_type                                            Homology Type
663                                     mlucifugus_homolog_subtype                                                 Ancestor
664                        mlucifugus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
665                                     mlucifugus_homolog_perc_id                    % Identity with respect to query gene
666                                  mlucifugus_homolog_perc_id_r1                 % Identity with respect to Microbat gene
667                                          mlucifugus_homolog_dn                                                       dN
668                                          mlucifugus_homolog_ds                                                       dS
669                                 mmusculus_homolog_ensembl_gene                                    Mouse Ensembl Gene ID
670                 mmusculus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
671                              mmusculus_homolog_ensembl_peptide                                 Mouse Ensembl Protein ID
672                                   mmusculus_homolog_chromosome                                    Mouse Chromosome Name
673                                  mmusculus_homolog_chrom_start                              Mouse Chromosome Start (bp)
674                                    mmusculus_homolog_chrom_end                                Mouse Chromosome End (bp)
675                               mmusculus_homolog_orthology_type                                            Homology Type
676                                      mmusculus_homolog_subtype                                                 Ancestor
677                         mmusculus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
678                                      mmusculus_homolog_perc_id                    % Identity with respect to query gene
679                                   mmusculus_homolog_perc_id_r1                    % Identity with respect to Mouse gene
680                                           mmusculus_homolog_dn                                                       dN
681                                           mmusculus_homolog_ds                                                       dS
682                                  mmurinus_homolog_ensembl_gene                              Mouse Lemur Ensembl Gene ID
683                  mmurinus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
684                               mmurinus_homolog_ensembl_peptide                           Mouse Lemur Ensembl Protein ID
685                                    mmurinus_homolog_chromosome                              Mouse Lemur Chromosome Name
686                                   mmurinus_homolog_chrom_start                        Mouse Lemur Chromosome Start (bp)
687                                     mmurinus_homolog_chrom_end                          Mouse Lemur Chromosome End (bp)
688                                mmurinus_homolog_orthology_type                                            Homology Type
689                                       mmurinus_homolog_subtype                                                 Ancestor
690                          mmurinus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
691                                       mmurinus_homolog_perc_id                    % Identity with respect to query gene
692                                    mmurinus_homolog_perc_id_r1              % Identity with respect to Mouse Lemur gene
693                                            mmurinus_homolog_dn                                                       dN
694                                            mmurinus_homolog_ds                                                       dS
695                                oniloticus_homolog_ensembl_gene                             Nile tilapia Ensembl Gene ID
696                oniloticus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
697                             oniloticus_homolog_ensembl_peptide                          Nile tilapia Ensembl Protein ID
698                                  oniloticus_homolog_chromosome                             Nile tilapia Chromosome Name
699                                 oniloticus_homolog_chrom_start                       Nile tilapia Chromosome Start (bp)
700                                   oniloticus_homolog_chrom_end                         Nile tilapia Chromosome End (bp)
701                              oniloticus_homolog_orthology_type                                            Homology Type
702                                     oniloticus_homolog_subtype                                                 Ancestor
703                        oniloticus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
704                                     oniloticus_homolog_perc_id                    % Identity with respect to query gene
705                                  oniloticus_homolog_perc_id_r1             % Identity with respect to Nile tilapia gene
706                                mdomestica_homolog_ensembl_gene                                  Opossum Ensembl Gene ID
707                mdomestica_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
708                             mdomestica_homolog_ensembl_peptide                               Opossum Ensembl Protein ID
709                                  mdomestica_homolog_chromosome                                  Opossum Chromosome Name
710                                 mdomestica_homolog_chrom_start                            Opossum Chromosome Start (bp)
711                                   mdomestica_homolog_chrom_end                              Opossum Chromosome End (bp)
712                              mdomestica_homolog_orthology_type                                            Homology Type
713                                     mdomestica_homolog_subtype                                                 Ancestor
714                        mdomestica_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
715                                     mdomestica_homolog_perc_id                    % Identity with respect to query gene
716                                  mdomestica_homolog_perc_id_r1                  % Identity with respect to Opossum gene
717                                          mdomestica_homolog_dn                                                       dN
718                                          mdomestica_homolog_ds                                                       dS
719                                   pabelii_homolog_ensembl_gene                                Orangutan Ensembl Gene ID
720                   pabelii_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
721                                pabelii_homolog_ensembl_peptide                             Orangutan Ensembl Protein ID
722                                     pabelii_homolog_chromosome                                Orangutan Chromosome Name
723                                    pabelii_homolog_chrom_start                          Orangutan Chromosome Start (bp)
724                                      pabelii_homolog_chrom_end                            Orangutan Chromosome End (bp)
725                                 pabelii_homolog_orthology_type                                            Homology Type
726                                        pabelii_homolog_subtype                                                 Ancestor
727                           pabelii_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
728                                        pabelii_homolog_perc_id                    % Identity with respect to query gene
729                                     pabelii_homolog_perc_id_r1                % Identity with respect to Orangutan gene
730                                             pabelii_homolog_dn                                                       dN
731                                             pabelii_homolog_ds                                                       dS
732                              amelanoleuca_homolog_ensembl_gene                                    Panda Ensembl Gene ID
733              amelanoleuca_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
734                           amelanoleuca_homolog_ensembl_peptide                                 Panda Ensembl Protein ID
735                                amelanoleuca_homolog_chromosome                                    Panda Chromosome Name
736                               amelanoleuca_homolog_chrom_start                              Panda Chromosome Start (bp)
737                                 amelanoleuca_homolog_chrom_end                                Panda Chromosome End (bp)
738                            amelanoleuca_homolog_orthology_type                                            Homology Type
739                                   amelanoleuca_homolog_subtype                                                 Ancestor
740                      amelanoleuca_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
741                                   amelanoleuca_homolog_perc_id                    % Identity with respect to query gene
742                                amelanoleuca_homolog_perc_id_r1                    % Identity with respect to Panda gene
743                                        amelanoleuca_homolog_dn                                                       dN
744                                        amelanoleuca_homolog_ds                                                       dS
745                                   sscrofa_homolog_ensembl_gene                                      Pig Ensembl Gene ID
746                   sscrofa_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
747                                sscrofa_homolog_ensembl_peptide                                   Pig Ensembl Protein ID
748                                     sscrofa_homolog_chromosome                                      Pig Chromosome Name
749                                    sscrofa_homolog_chrom_start                                Pig Chromosome Start (bp)
750                                      sscrofa_homolog_chrom_end                                  Pig Chromosome End (bp)
751                                 sscrofa_homolog_orthology_type                                            Homology Type
752                                        sscrofa_homolog_subtype                                                 Ancestor
753                           sscrofa_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
754                                        sscrofa_homolog_perc_id                    % Identity with respect to query gene
755                                     sscrofa_homolog_perc_id_r1                      % Identity with respect to Pig gene
756                                             sscrofa_homolog_dn                                                       dN
757                                             sscrofa_homolog_ds                                                       dS
758                                 oprinceps_homolog_ensembl_gene                                     Pika Ensembl Gene ID
759                 oprinceps_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
760                              oprinceps_homolog_ensembl_peptide                                  Pika Ensembl Protein ID
761                                   oprinceps_homolog_chromosome                                     Pika Chromosome Name
762                                  oprinceps_homolog_chrom_start                               Pika Chromosome Start (bp)
763                                    oprinceps_homolog_chrom_end                                 Pika Chromosome End (bp)
764                               oprinceps_homolog_orthology_type                                            Homology Type
765                                      oprinceps_homolog_subtype                                                 Ancestor
766                         oprinceps_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
767                                      oprinceps_homolog_perc_id                    % Identity with respect to query gene
768                                   oprinceps_homolog_perc_id_r1                     % Identity with respect to Pika gene
769                                           oprinceps_homolog_dn                                                       dN
770                                           oprinceps_homolog_ds                                                       dS
771                                xmaculatus_homolog_ensembl_gene                                Platyfish Ensembl Gene ID
772                xmaculatus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
773                             xmaculatus_homolog_ensembl_peptide                             Platyfish Ensembl Protein ID
774                                  xmaculatus_homolog_chromosome                                Platyfish Chromosome Name
775                                 xmaculatus_homolog_chrom_start                          Platyfish Chromosome Start (bp)
776                                   xmaculatus_homolog_chrom_end                            Platyfish Chromosome End (bp)
777                              xmaculatus_homolog_orthology_type                                            Homology Type
778                                     xmaculatus_homolog_subtype                                                 Ancestor
779                        xmaculatus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
780                                     xmaculatus_homolog_perc_id                    % Identity with respect to query gene
781                                  xmaculatus_homolog_perc_id_r1                % Identity with respect to Platyfish gene
782                                 oanatinus_homolog_ensembl_gene                                 Platypus Ensembl Gene ID
783                 oanatinus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
784                              oanatinus_homolog_ensembl_peptide                              Platypus Ensembl Protein ID
785                                   oanatinus_homolog_chromosome                                 Platypus Chromosome Name
786                                  oanatinus_homolog_chrom_start                           Platypus Chromosome Start (bp)
787                                    oanatinus_homolog_chrom_end                             Platypus Chromosome End (bp)
788                               oanatinus_homolog_orthology_type                                            Homology Type
789                                      oanatinus_homolog_subtype                                                 Ancestor
790                         oanatinus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
791                                      oanatinus_homolog_perc_id                    % Identity with respect to query gene
792                                   oanatinus_homolog_perc_id_r1                 % Identity with respect to Platypus gene
793                                           oanatinus_homolog_dn                                                       dN
794                                           oanatinus_homolog_ds                                                       dS
795                                ocuniculus_homolog_ensembl_gene                                   Rabbit Ensembl Gene ID
796                ocuniculus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
797                             ocuniculus_homolog_ensembl_peptide                                Rabbit Ensembl Protein ID
798                                  ocuniculus_homolog_chromosome                                   Rabbit Chromosome Name
799                                 ocuniculus_homolog_chrom_start                             Rabbit Chromosome Start (bp)
800                                   ocuniculus_homolog_chrom_end                               Rabbit Chromosome End (bp)
801                              ocuniculus_homolog_orthology_type                                            Homology Type
802                                     ocuniculus_homolog_subtype                                                 Ancestor
803                        ocuniculus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
804                                     ocuniculus_homolog_perc_id                    % Identity with respect to query gene
805                                  ocuniculus_homolog_perc_id_r1                   % Identity with respect to Rabbit gene
806                                          ocuniculus_homolog_dn                                                       dN
807                                          ocuniculus_homolog_ds                                                       dS
808                               rnorvegicus_homolog_ensembl_gene                                      Rat Ensembl Gene ID
809               rnorvegicus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
810                            rnorvegicus_homolog_ensembl_peptide                                   Rat Ensembl Protein ID
811                                 rnorvegicus_homolog_chromosome                                      Rat Chromosome Name
812                                rnorvegicus_homolog_chrom_start                                Rat Chromosome Start (bp)
813                                  rnorvegicus_homolog_chrom_end                                  Rat Chromosome End (bp)
814                             rnorvegicus_homolog_orthology_type                                            Homology Type
815                                    rnorvegicus_homolog_subtype                                                 Ancestor
816                       rnorvegicus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
817                                    rnorvegicus_homolog_perc_id                    % Identity with respect to query gene
818                                 rnorvegicus_homolog_perc_id_r1                      % Identity with respect to Rat gene
819                                         rnorvegicus_homolog_dn                                                       dN
820                                         rnorvegicus_homolog_ds                                                       dS
821                                 pcapensis_homolog_ensembl_gene                               Rock Hyrax Ensembl Gene ID
822                 pcapensis_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
823                              pcapensis_homolog_ensembl_peptide                            Rock Hyrax Ensembl Protein ID
824                                   pcapensis_homolog_chromosome                               Rock Hyrax Chromosome Name
825                                  pcapensis_homolog_chrom_start                         Rock Hyrax Chromosome Start (bp)
826                                    pcapensis_homolog_chrom_end                           Rock Hyrax Chromosome End (bp)
827                               pcapensis_homolog_orthology_type                                            Homology Type
828                                      pcapensis_homolog_subtype                                                 Ancestor
829                         pcapensis_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
830                                      pcapensis_homolog_perc_id                    % Identity with respect to query gene
831                                   pcapensis_homolog_perc_id_r1               % Identity with respect to Rock Hyrax gene
832                                           pcapensis_homolog_dn                                                       dN
833                                           pcapensis_homolog_ds                                                       dS
834                                    oaries_homolog_ensembl_gene                                    Sheep Ensembl Gene ID
835                    oaries_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
836                                 oaries_homolog_ensembl_peptide                                 Sheep Ensembl Protein ID
837                                      oaries_homolog_chromosome                                    Sheep Chromosome Name
838                                     oaries_homolog_chrom_start                              Sheep Chromosome Start (bp)
839                                       oaries_homolog_chrom_end                                Sheep Chromosome End (bp)
840                                  oaries_homolog_orthology_type                                            Homology Type
841                                         oaries_homolog_subtype                                                 Ancestor
842                            oaries_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
843                                         oaries_homolog_perc_id                    % Identity with respect to query gene
844                                      oaries_homolog_perc_id_r1                    % Identity with respect to Sheep gene
845                                              oaries_homolog_dn                                                       dN
846                                              oaries_homolog_ds                                                       dS
847                                choffmanni_homolog_ensembl_gene                                    Sloth Ensembl Gene ID
848                choffmanni_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
849                             choffmanni_homolog_ensembl_peptide                                 Sloth Ensembl Protein ID
850                                  choffmanni_homolog_chromosome                                    Sloth Chromosome Name
851                                 choffmanni_homolog_chrom_start                              Sloth Chromosome Start (bp)
852                                   choffmanni_homolog_chrom_end                                Sloth Chromosome End (bp)
853                              choffmanni_homolog_orthology_type                                            Homology Type
854                                     choffmanni_homolog_subtype                                                 Ancestor
855                        choffmanni_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
856                                     choffmanni_homolog_perc_id                    % Identity with respect to query gene
857                                  choffmanni_homolog_perc_id_r1                    % Identity with respect to Sloth gene
858                                          choffmanni_homolog_dn                                                       dN
859                                          choffmanni_homolog_ds                                                       dS
860                                 loculatus_homolog_ensembl_gene                              Spotted gar Ensembl Gene ID
861                 loculatus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
862                              loculatus_homolog_ensembl_peptide                           Spotted gar Ensembl Protein ID
863                                   loculatus_homolog_chromosome                              Spotted gar Chromosome Name
864                                  loculatus_homolog_chrom_start                        Spotted gar Chromosome Start (bp)
865                                    loculatus_homolog_chrom_end                          Spotted gar Chromosome End (bp)
866                               loculatus_homolog_orthology_type                                            Homology Type
867                                      loculatus_homolog_subtype                                                 Ancestor
868                         loculatus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
869                                      loculatus_homolog_perc_id                    % Identity with respect to query gene
870                                   loculatus_homolog_perc_id_r1              % Identity with respect to Spotted gar gene
871                         itridecemlineatus_homolog_ensembl_gene                                 Squirrel Ensembl Gene ID
872         itridecemlineatus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
873                      itridecemlineatus_homolog_ensembl_peptide                              Squirrel Ensembl Protein ID
874                           itridecemlineatus_homolog_chromosome                                 Squirrel Chromosome Name
875                          itridecemlineatus_homolog_chrom_start                           Squirrel Chromosome Start (bp)
876                            itridecemlineatus_homolog_chrom_end                             Squirrel Chromosome End (bp)
877                       itridecemlineatus_homolog_orthology_type                                            Homology Type
878                              itridecemlineatus_homolog_subtype                                                 Ancestor
879                 itridecemlineatus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
880                              itridecemlineatus_homolog_perc_id                    % Identity with respect to query gene
881                           itridecemlineatus_homolog_perc_id_r1                 % Identity with respect to Squirrel gene
882                                   itridecemlineatus_homolog_dn                                                       dN
883                                   itridecemlineatus_homolog_ds                                                       dS
884                                gaculeatus_homolog_ensembl_gene                              Stickleback Ensembl Gene ID
885                gaculeatus_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
886                             gaculeatus_homolog_ensembl_peptide                           Stickleback Ensembl Protein ID
887                                  gaculeatus_homolog_chromosome                              Stickleback Chromosome Name
888                                 gaculeatus_homolog_chrom_start                        Stickleback Chromosome Start (bp)
889                                   gaculeatus_homolog_chrom_end                          Stickleback Chromosome End (bp)
890                              gaculeatus_homolog_orthology_type                                            Homology Type
891                                     gaculeatus_homolog_subtype                                                 Ancestor
892                        gaculeatus_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
893                                     gaculeatus_homolog_perc_id                    % Identity with respect to query gene
894                                  gaculeatus_homolog_perc_id_r1              % Identity with respect to Stickleback gene
895                                          gaculeatus_homolog_dn                                                       dN
896                                          gaculeatus_homolog_ds                                                       dS
897                                 tsyrichta_homolog_ensembl_gene                                  Tarsier Ensembl Gene ID
898                 tsyrichta_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
899                              tsyrichta_homolog_ensembl_peptide                               Tarsier Ensembl Protein ID
900                                   tsyrichta_homolog_chromosome                                  Tarsier Chromosome Name
901                                  tsyrichta_homolog_chrom_start                            Tarsier Chromosome Start (bp)
902                                    tsyrichta_homolog_chrom_end                              Tarsier Chromosome End (bp)
903                               tsyrichta_homolog_orthology_type                                            Homology Type
904                                      tsyrichta_homolog_subtype                                                 Ancestor
905                         tsyrichta_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
906                                      tsyrichta_homolog_perc_id                    % Identity with respect to query gene
907                                   tsyrichta_homolog_perc_id_r1                  % Identity with respect to Tarsier gene
908                                           tsyrichta_homolog_dn                                                       dN
909                                           tsyrichta_homolog_ds                                                       dS
910                                 sharrisii_homolog_ensembl_gene                          Tasmanian Devil Ensembl Gene ID
911                 sharrisii_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
912                              sharrisii_homolog_ensembl_peptide                       Tasmanian Devil Ensembl Protein ID
913                                   sharrisii_homolog_chromosome                          Tasmanian Devil Chromosome Name
914                                  sharrisii_homolog_chrom_start                    Tasmanian Devil Chromosome Start (bp)
915                                    sharrisii_homolog_chrom_end                      Tasmanian Devil Chromosome End (bp)
916                               sharrisii_homolog_orthology_type                                            Homology Type
917                                      sharrisii_homolog_subtype                                                 Ancestor
918                         sharrisii_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
919                                      sharrisii_homolog_perc_id                    % Identity with respect to query gene
920                                   sharrisii_homolog_perc_id_r1          % Identity with respect to Tasmanian Devil gene
921                                           sharrisii_homolog_dn                                                       dN
922                                           sharrisii_homolog_ds                                                       dS
923                             tnigroviridis_homolog_ensembl_gene                                Tetraodon Ensembl Gene ID
924             tnigroviridis_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
925                          tnigroviridis_homolog_ensembl_peptide                             Tetraodon Ensembl Protein ID
926                               tnigroviridis_homolog_chromosome                                Tetraodon Chromosome Name
927                              tnigroviridis_homolog_chrom_start                          Tetraodon Chromosome Start (bp)
928                                tnigroviridis_homolog_chrom_end                            Tetraodon Chromosome End (bp)
929                           tnigroviridis_homolog_orthology_type                                            Homology Type
930                                  tnigroviridis_homolog_subtype                                                 Ancestor
931                     tnigroviridis_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
932                                  tnigroviridis_homolog_perc_id                    % Identity with respect to query gene
933                               tnigroviridis_homolog_perc_id_r1                % Identity with respect to Tetraodon gene
934                                       tnigroviridis_homolog_dn                                                       dN
935                                       tnigroviridis_homolog_ds                                                       dS
936                                tbelangeri_homolog_ensembl_gene                               Tree Shrew Ensembl Gene ID
937                tbelangeri_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
938                             tbelangeri_homolog_ensembl_peptide                            Tree Shrew Ensembl Protein ID
939                                  tbelangeri_homolog_chromosome                               Tree Shrew Chromosome Name
940                                 tbelangeri_homolog_chrom_start                         Tree Shrew Chromosome Start (bp)
941                                   tbelangeri_homolog_chrom_end                           Tree Shrew Chromosome End (bp)
942                              tbelangeri_homolog_orthology_type                                            Homology Type
943                                     tbelangeri_homolog_subtype                                                 Ancestor
944                        tbelangeri_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
945                                     tbelangeri_homolog_perc_id                    % Identity with respect to query gene
946                                  tbelangeri_homolog_perc_id_r1               % Identity with respect to Tree Shrew gene
947                                          tbelangeri_homolog_dn                                                       dN
948                                          tbelangeri_homolog_ds                                                       dS
949                                mgallopavo_homolog_ensembl_gene                                   Turkey Ensembl Gene ID
950                mgallopavo_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
951                             mgallopavo_homolog_ensembl_peptide                                Turkey Ensembl Protein ID
952                                  mgallopavo_homolog_chromosome                                   Turkey Chromosome Name
953                                 mgallopavo_homolog_chrom_start                             Turkey Chromosome Start (bp)
954                                   mgallopavo_homolog_chrom_end                               Turkey Chromosome End (bp)
955                              mgallopavo_homolog_orthology_type                                            Homology Type
956                                     mgallopavo_homolog_subtype                                                 Ancestor
957                        mgallopavo_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
958                                     mgallopavo_homolog_perc_id                    % Identity with respect to query gene
959                                  mgallopavo_homolog_perc_id_r1                   % Identity with respect to Turkey gene
960                                          mgallopavo_homolog_dn                                                       dN
961                                          mgallopavo_homolog_ds                                                       dS
962                                  meugenii_homolog_ensembl_gene                                  Wallaby Ensembl Gene ID
963                  meugenii_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
964                               meugenii_homolog_ensembl_peptide                               Wallaby Ensembl Protein ID
965                                    meugenii_homolog_chromosome                                  Wallaby Chromosome Name
966                                   meugenii_homolog_chrom_start                            Wallaby Chromosome Start (bp)
967                                     meugenii_homolog_chrom_end                              Wallaby Chromosome End (bp)
968                                meugenii_homolog_orthology_type                                            Homology Type
969                                       meugenii_homolog_subtype                                                 Ancestor
970                          meugenii_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
971                                       meugenii_homolog_perc_id                    % Identity with respect to query gene
972                                    meugenii_homolog_perc_id_r1                  % Identity with respect to Wallaby gene
973                                            meugenii_homolog_dn                                                       dN
974                                            meugenii_homolog_ds                                                       dS
975                               xtropicalis_homolog_ensembl_gene                                  Xenopus Ensembl Gene ID
976               xtropicalis_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
977                            xtropicalis_homolog_ensembl_peptide                               Xenopus Ensembl Protein ID
978                                 xtropicalis_homolog_chromosome                                  Xenopus Chromosome Name
979                                xtropicalis_homolog_chrom_start                            Xenopus Chromosome Start (bp)
980                                  xtropicalis_homolog_chrom_end                              Xenopus Chromosome End (bp)
981                             xtropicalis_homolog_orthology_type                                            Homology Type
982                                    xtropicalis_homolog_subtype                                                 Ancestor
983                       xtropicalis_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
984                                    xtropicalis_homolog_perc_id                    % Identity with respect to query gene
985                                 xtropicalis_homolog_perc_id_r1                  % Identity with respect to Xenopus gene
986                                         xtropicalis_homolog_dn                                                       dN
987                                         xtropicalis_homolog_ds                                                       dS
988                               scerevisiae_homolog_ensembl_gene                                    Yeast Ensembl Gene ID
989               scerevisiae_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
990                            scerevisiae_homolog_ensembl_peptide                                 Yeast Ensembl Protein ID
991                                 scerevisiae_homolog_chromosome                                    Yeast Chromosome Name
992                                scerevisiae_homolog_chrom_start                              Yeast Chromosome Start (bp)
993                                  scerevisiae_homolog_chrom_end                                Yeast Chromosome End (bp)
994                             scerevisiae_homolog_orthology_type                                            Homology Type
995                                    scerevisiae_homolog_subtype                                                 Ancestor
996                       scerevisiae_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
997                                    scerevisiae_homolog_perc_id                    % Identity with respect to query gene
998                                 scerevisiae_homolog_perc_id_r1                    % Identity with respect to Yeast gene
999                                         scerevisiae_homolog_dn                                                       dN
1000                                        scerevisiae_homolog_ds                                                       dS
1001                                 tguttata_homolog_ensembl_gene                              Zebra Finch Ensembl Gene ID
1002                 tguttata_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
1003                              tguttata_homolog_ensembl_peptide                           Zebra Finch Ensembl Protein ID
1004                                   tguttata_homolog_chromosome                              Zebra Finch Chromosome Name
1005                                  tguttata_homolog_chrom_start                        Zebra Finch Chromosome Start (bp)
1006                                    tguttata_homolog_chrom_end                          Zebra Finch Chromosome End (bp)
1007                               tguttata_homolog_orthology_type                                            Homology Type
1008                                      tguttata_homolog_subtype                                                 Ancestor
1009                         tguttata_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
1010                                      tguttata_homolog_perc_id                    % Identity with respect to query gene
1011                                   tguttata_homolog_perc_id_r1              % Identity with respect to Zebra Finch gene
1012                                           tguttata_homolog_dn                                                       dN
1013                                           tguttata_homolog_ds                                                       dS
1014                                   drerio_homolog_ensembl_gene                                Zebrafish Ensembl Gene ID
1015                   drerio_homolog_canonical_transcript_protein                       Canonical Protein or Transcript ID
1016                                drerio_homolog_ensembl_peptide                             Zebrafish Ensembl Protein ID
1017                                     drerio_homolog_chromosome                                Zebrafish Chromosome Name
1018                                    drerio_homolog_chrom_start                          Zebrafish Chromosome Start (bp)
1019                                      drerio_homolog_chrom_end                            Zebrafish Chromosome End (bp)
1020                                 drerio_homolog_orthology_type                                            Homology Type
1021                                        drerio_homolog_subtype                                                 Ancestor
1022                           drerio_homolog_orthology_confidence                     Orthology confidence [0 low, 1 high]
1023                                        drerio_homolog_perc_id                    % Identity with respect to query gene
1024                                     drerio_homolog_perc_id_r1                % Identity with respect to Zebrafish gene
1025                                             drerio_homolog_dn                                                       dN
1026                                             drerio_homolog_ds                                                       dS
1027                                 hsapiens_paralog_ensembl_gene                            Human Paralog Ensembl Gene ID
1028                 hsapiens_paralog_canonical_transcript_protein                       Canonical Protein or Transcript ID
1029                              hsapiens_paralog_ensembl_peptide                         Human Paralog Ensembl Protein ID
1030                                   hsapiens_paralog_chromosome                            Human Paralog Chromosome Name
1031                                  hsapiens_paralog_chrom_start                             Human Paralog Chr Start (bp)
1032                                    hsapiens_paralog_chrom_end                               Human Paralog Chr End (bp)
1033                               hsapiens_paralog_orthology_type                                            Homology Type
1034                                      hsapiens_paralog_subtype                                                 Ancestor
1035                          hsapiens_paralog_paralogy_confidence                      Paralogy confidence [0 low, 1 high]
1036                                      hsapiens_paralog_perc_id                    % Identity with respect to query gene
1037                                   hsapiens_paralog_perc_id_r1                    % Identity with respect to Human gene
1038                                           hsapiens_paralog_dn                                                       dN
1039                                           hsapiens_paralog_ds                                                       dS
1040                                               ensembl_gene_id                                          Ensembl Gene ID
1041                                         ensembl_transcript_id                                    Ensembl Transcript ID
1042                                            ensembl_peptide_id                                       Ensembl Protein ID
1043                                               chromosome_name                                          Chromosome Name
1044                                                start_position                                          Gene Start (bp)
1045                                                  end_position                                            Gene End (bp)
1046                                                        strand                                                   Strand
1047                                                          band                                                     Band
1048                                              external_gene_id                                     Associated Gene Name
1049                                              external_gene_db                                       Associated Gene DB
1050                                              transcript_count                                         Transcript count
1051                                         percentage_gc_content                                             % GC content
1052                                                   description                                              Description
1053                                                   source_name                                         Variation Source
1054                                            source_description                                       Source description
1055                                                   external_id                                             Reference ID
1056                                                        allele                                          Variant Alleles
1057                                                     validated                                          Evidence status
1058                                                     mapweight                                                Mapweight
1059                                                  minor_allele                                             Minor allele
1060                                             minor_allele_freq                                   Minor allele frequency
1061                                            minor_allele_count                                       Minor allele count
1062                                         clinical_significance                                    Clinical significance
1063                                           transcript_location                                 Transcript location (bp)
1064                                         snp_chromosome_strand                              Variation Chromosome Strand
1065                                              peptide_location                                    Protein location (aa)
1066                                           chromosome_location                                 Chromosome Location (bp)
1067                                      polyphen_prediction_2076                                      PolyPhen prediction
1068                                           polyphen_score_2076                                           PolyPhen score
1069                                          sift_prediction_2076                                          SIFT prediction
1070                                               sift_score_2076                                               SIFT score
1071                                   distance_to_transcript_2076                                   Distance to transcript
1072                                                cds_start_2076                                                CDS Start
1073                                                  cds_end_2076                                                  CDS End
1074                                                 peptide_shift                                           Protein Allele
1075                                             synonymous_status                  Consequence Type (Transcript Variation)
1076                                            allele_string_2076                              Consequence specific allele
1077                                           somatic_source_name                                           Variant Source
1078                                    somatic_source_description                                       Source description
1079                                          somatic_reference_id                                             Reference ID
1080                                                somatic_allele                                          Variant Alleles
1081                                             somatic_validated                                        Validation status
1082                                             somatic_mapweight                                                Mapweight
1083                                          somatic_minor_allele                                             Minor allele
1084                                     somatic_minor_allele_freq                                   Minor allele frequency
1085                                    somatic_minor_allele_count                                       Minor allele count
1086                                 somatic_clinical_significance                                    Clinical significance
1087                                   somatic_transcript_location                                 Transcript location (bp)
1088                                 somatic_snp_chromosome_strand                              Variation Chromosome Strand
1089                                      somatic_peptide_location                                    Protein location (aa)
1090                                   somatic_chromosome_location                                 Chromosome Location (bp)
1091        mart_transcript_variation_som__dm_sift_prediction_2076                                          SIFT prediction
1092             mart_transcript_variation_som__dm_sift_score_2076                                               SIFT score
1093    mart_transcript_variation_som__dm_polyphen_prediction_2076                                      PolyPhen prediction
1094         mart_transcript_variation_som__dm_polyphen_score_2076                                           PolyPhen score
1095 mart_transcript_variation_som__dm_distance_to_transcript_2076                                   Distance to transcript
1096                                         somatic_peptide_shift                                           Protein Allele
1097                                     somatic_synonymous_status                  Consequence Type (Transcript Variation)
1098                                        somatic_cds_start_2076                                                CDS Start
1099                                          somatic_cds_end_2076                                                  CDS End
1100          mart_transcript_variation_som__dm_allele_string_2076                              Consequence specific allele
1101                                        transcript_exon_intron                                   Unspliced (Transcript)
1102                                              gene_exon_intron                                         Unspliced (Gene)
1103                                              transcript_flank                                       Flank (Transcript)
1104                                                    gene_flank                                             Flank (Gene)
1105                                       coding_transcript_flank                         Flank-coding region (Transcript)
1106                                             coding_gene_flank                               Flank-coding region (Gene)
1107                                                          5utr                                                   5' UTR
1108                                                          3utr                                                   3' UTR
1109                                                     gene_exon                                           Exon sequences
1110                                                          cdna                                           cDNA sequences
1111                                                        coding                                          Coding sequence
1112                                                       peptide                                                  Protein
1113                                                upstream_flank                                           upstream_flank
1114                                              downstream_flank                                         downstream_flank
1115                                               ensembl_gene_id                                          Ensembl Gene ID
1116                                                   description                                              Description
1117                                              external_gene_id                                     Associated Gene Name
1118                                              external_gene_db                                       Associated Gene DB
1119                                               chromosome_name                                          Chromosome Name
1120                                                start_position                                          Gene Start (bp)
1121                                                  end_position                                            Gene End (bp)
1122                                                  gene_biotype                                             Gene Biotype
1123                                                        family                             Ensembl Protein Family ID(s)
1124                                             cdna_coding_start                                  CDS start (within cDNA)
1125                                               cdna_coding_end                                    CDS end (within cDNA)
1126                                                   5_utr_start                                             5' UTR Start
1127                                                     5_utr_end                                               5' UTR End
1128                                                   3_utr_start                                             3' UTR Start
1129                                                     3_utr_end                                               3' UTR End
1130                                         ensembl_transcript_id                                    Ensembl Transcript ID
1131                                            ensembl_peptide_id                                       Ensembl Protein ID
1132                                            transcript_biotype                                       Transcript Biotype
1133                                                        strand                                                   Strand
1134                                              transcript_start                                    Transcript Start (bp)
1135                                                transcript_end                                      Transcript End (bp)
1136                                                    cds_length                                               CDS Length
1137                                                     cds_start                                                CDS Start
1138                                                       cds_end                                                  CDS End
1139                                               ensembl_exon_id                                          Ensembl Exon ID
1140                                              exon_chrom_start                                      Exon Chr Start (bp)
1141                                                exon_chrom_end                                        Exon Chr End (bp)
1142                                                        strand                                                   Strand
1143                                                          rank                                  Exon Rank in Transcript
1144                                                         phase                                                    phase
1145                                             cdna_coding_start                                        cDNA coding start
1146                                               cdna_coding_end                                          cDNA coding end
1147                                          genomic_coding_start                                     Genomic coding start
1148                                            genomic_coding_end                                       Genomic coding end
1149                                               is_constitutive                                        Constitutive Exon
