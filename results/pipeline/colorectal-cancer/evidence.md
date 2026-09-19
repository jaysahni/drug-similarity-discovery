# Literature evidence - colorectal cancer

Pipeline task T3. Source: Europe PMC REST (`https://www.ebi.ac.uk/europepmc/webservices/rest/search`), searched 2026-09-19. Targets from `targets.json` (Open Targets MONDO_0005575, data 26.06).

**n = 25 targets searched, 73 PMIDs emitted, 73 verified resolving, 0 mismatches.** Every PMID below was re-fetched from Europe PMC by id in a second request and its title compared with the one written here; any that failed was removed and is listed in the verification section. Cross-source check against PubMed E-utilities: 73/73 titles matched PubMed esummary, 0 dropped.

Each quoted line is a **verbatim sentence from the downloaded abstract**, not a paraphrase - markup is rendered to text but no word is added, removed or reordered. Quotes are truncated at 400 characters. A quote that could not pair the target with an inhibitor/binder term is labelled underneath with its weaker claim tier.

## 1. MSH2 (P43246) - mutS homolog 2

Open Targets: overall **0.9239**, literature 0.934, known-drug (clinical) -, 0 drug/clinical candidates (0 approved).

Europe PMC hits: **36**, kept 3.

1. **PMID 42123360** - *Dysfunctional DNA Mismatch Repair Drives the Evolution of Gene Amplification in MTX-Resistant Human Colorectal Cancer Cells.* - International journal of molecular sciences, 2026. doi:10.3390/ijms27093774
   > MMR inhibition was achieved by depleting MSH2.

2. **PMID 42682916** - *Immune Checkpoint Inhibitor Resistance in Rectal Adenocarcinoma Containing a Mismatch Repair Deficiency Component: A Case Report.* - Cancer diagnosis & prognosis, 2026. doi:10.21873/cdp.10595
   > Immunohistochemically, MLH1 and MSH2 expression was lost in the high-grade adenocarcinoma component, but retained in the adenomatous and low-grade adenocarcinoma components.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

3. **PMID 41773741** - *Risk of Metachronous Colorectal Cancer after Segmental or Extended Resection in Patients with Lynch Syndrome.* - Journal of the American College of Surgeons, 2026. doi:10.1097/xcs.0000000000001892
   > Metachronous CRC was more common in patients with variants in high-risk ( MLH1 , MSH2 ) vs low-risk ( MSH6 , PMS2 ) genes: 80 of 326 patients (25%) vs 8 of 124 patients (6%) (p < 0.001).
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"MSH2" OR TITLE_ABS:"MSH-2" OR TITLE_ABS:"COCA1" OR TITLE_ABS:"FCC1" OR TITLE_ABS:"HNPCC1" OR TITLE_ABS:"LCFS2" OR TITLE_ABS:"LYNCH1") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 2. MSH6 (P52701) - mutS homolog 6

Open Targets: overall **0.9157**, literature 0.693, known-drug (clinical) -, 0 drug/clinical candidates (0 approved).

Europe PMC hits: **29**, kept 3.

1. **PMID 42466612** - *5-FU in combination with PARP inhibitor ABT-888 deregulates MGMT-dependent mismatch repair (MMR) pathway in MMR-proficient colorectal cancer stem cells by modulating MGMT/PARP1/MSH6 complex.* - Expert opinion on therapeutic targets, 2026. doi:10.1080/14728222.2026.2706483
   > Background: Previous study showed the PARP inhibitor ABT-888 potentiates the cytotoxicity of 5-fluorouracil (5-FU) by inhibiting PARP1-mediated mismatch repair (MMR) pathway via MSH6 deregulation in MMR-proficient colorectal cancer stem cells (CRC-CSCs).

2. **PMID 41773741** - *Risk of Metachronous Colorectal Cancer after Segmental or Extended Resection in Patients with Lynch Syndrome.* - Journal of the American College of Surgeons, 2026. doi:10.1097/xcs.0000000000001892
   > Metachronous CRC was more common in patients with variants in high-risk ( MLH1 , MSH2 ) vs low-risk ( MSH6 , PMS2 ) genes: 80 of 326 patients (25%) vs 8 of 124 patients (6%) (p < 0.001).
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

3. **PMID 40086819** - *Protein arginine methyltransferase 6 enhances immune checkpoint blockade efficacy via the STING pathway in MMR-proficient colorectal cancer.* - Journal for immunotherapy of cancer, 2025. doi:10.1136/jitc-2024-010639
   > MMR is the critical DNA repair pathway that maintains genomic integrity by correcting DNA mismatches, which is mediated by the MutSα or MutSβ complex consisting of MSH2 with MSH6 and MSH3, respectively.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"MSH6" OR TITLE_ABS:"MSH-6" OR TITLE_ABS:"HNPCC5" OR TITLE_ABS:"LYNCH5" OR TITLE_ABS:"MMRCS3" OR TITLE_ABS:"hMSH6" OR TITLE_ABS:"p160") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 3. MLH1 (P40692) - mutL homolog 1

Open Targets: overall **0.9142**, literature 0.986, known-drug (clinical) -, 0 drug/clinical candidates (0 approved).

Europe PMC hits: **70**, kept 3.

1. **PMID 42191394** - *MLH1 silencing reprograms PI3K/AKT and JAK/STAT signaling to enhance 5-Fluorouracil sensitivity in colorectal cancer cells.* - DNA repair, 2026. doi:10.1016/j.dnarep.2026.103941
   > Background/aim: This study aimed to evaluate the anticancer efficacy of 5-Fluorouracil (5-FU) in Caco-2 human colorectal adenocarcinoma cells treated with siRNA-mediated MLH1 gene inhibition, in terms of cell viability, apoptosis, and related gene/protein expression.

2. **PMID 41875719** - *Autophagy inhibition enhances PD-1 blockade efficacy in mismatch repair-proficient colorectal cancer.* - Cancer genetics, 2026. doi:10.1016/j.cancergen.2026.03.005
   > Methods: We utilized CRISPR/Cas9 to delete the mismatch repair gene Mlh1 in CT26 colorectal cancer cells, generating CT26-dMMR cells.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

3. **PMID 40458404** - *Fusobacterium nucleatum downregulated MLH1 expression in colorectal cancer by activating autophagy-lysosome pathway.* - Frontiers in immunology, 2025. doi:10.3389/fimmu.2025.1586146
   > RNA sequencing was used to evaluate the involved pathways, and non-targeted metabolomics was employed to analyze the metabolites regulating MLH1.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"MLH1" OR TITLE_ABS:"MLH-1" OR TITLE_ABS:"COCA2" OR TITLE_ABS:"FCC2" OR TITLE_ABS:"HNPCC2" OR TITLE_ABS:"LYNCH2" OR TITLE_ABS:"MMRCS1") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 4. PMS2 (P54278) - PMS1 homolog 2, mismatch repair system component

Open Targets: overall **0.8825**, literature 0.907, known-drug (clinical) -, 0 drug/clinical candidates (0 approved).

Europe PMC hits: **23**, kept 3.

1. **PMID 41773741** - *Risk of Metachronous Colorectal Cancer after Segmental or Extended Resection in Patients with Lynch Syndrome.* - Journal of the American College of Surgeons, 2026. doi:10.1097/xcs.0000000000001892
   > Metachronous CRC was more common in patients with variants in high-risk ( MLH1 , MSH2 ) vs low-risk ( MSH6 , PMS2 ) genes: 80 of 326 patients (25%) vs 8 of 124 patients (6%) (p < 0.001).
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

2. **PMID 40344511** - *PRMT5 Inhibitor Synergizes with Chemotherapy to Induce Resembling Mismatch Repair Deficiency and Enhance Anti-TIGIT Therapy in Microsatellite-Stable Colorectal Cancer.* - Advanced science (Weinheim, Baden-Wurttemberg, Germany), 2025. doi:10.1002/advs.202500271
   > Experiments in this study identify a DNA damage repair-related epigenetic gene, protein arginine methyltransferase 5 (PRMT5), whose inhibition enhances Irinotecan (CPT-11) sensitivity and synergistically induces a postmeiotic segregation increased 2 (PMS2)-deficient-like state, leading to the release of cytosolic double-stranded DNA.

3. **PMID 42405254** - *A retrospective study of differential prognostic factors in early-onset versus late-onset colorectal cancer: a comprehensive clinical and machine learning analysis.* - PeerJ, 2026. doi:10.7717/peerj.21484
   > In EO-CRC, distant metastasis, family history, Tumor, Node, and Metastasis (TNM) stage, PMS1 homolog 2, mismatch repair system component (PMS2), MutS Homolog 6 (MSH6), tumor size, concurrent polyps, and Ki-67 were major predictors.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"PMS2" OR TITLE_ABS:"PMS-2" OR TITLE_ABS:"HNPCC4" OR TITLE_ABS:"LYNCH4" OR TITLE_ABS:"MLH4" OR TITLE_ABS:"MMRCS4" OR TITLE_ABS:"PMSL2") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 5. APC (P25054) - APC regulator of Wnt signaling pathway

Open Targets: overall **0.8585**, literature 0.990, known-drug (clinical) -, 0 drug/clinical candidates (0 approved).

Europe PMC hits: **521**, kept 3.

1. **PMID 42495523** - *ACSS2 is required for colorectal cancer progression and a druggable target for colorectal cancer treatment.* - iScience, 2026. doi:10.1016/j.isci.2026.116860
   > Treatment with a small-molecule ACSS2 inhibitor markedly suppresses tumor growth in allograft/xenograft and Apc-mutant CRC models.

2. **PMID 41486293** - *Targeting PTPN13 with 11-amino-acid peptides of C-terminal APC prevents immune evasion of colorectal cancer.* - Cell research, 2026. doi:10.1038/s41422-025-01206-4
   > Colorectal cancer (CRC) remains largely refractory to immune-checkpoint blockade, with adenomatous polyposis coli (APC) mutations present in 80%-90% of cases.

3. **PMID 42283915** - *TELO2-interacting protein 1 (TTI1), a novel Wnt/β-catenin target gene, decreases chemo-sensitivity in colorectal cancer by modulating DNA damage responses.* - Molecular biomedicine, 2026. doi:10.1186/s43556-026-00475-8
   > Pharmacological suppression of the TTT complex by piperlongumine mimics TTI1 loss and enhances the anti-tumor activity of chemotherapy in cell lines, xenografts, Apc-mutant patient-derived organoids and Apcmin/+ colonic adenomas.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"APC" OR TITLE_ABS:"BTPS2" OR TITLE_ABS:"PPP1R46" OR TITLE_ABS:"DESMD") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 6. BRAF (P15056) - B-Raf proto-oncogene, serine/threonine kinase

Open Targets: overall **0.8385**, literature 0.996, known-drug (clinical) 0.952, 18 drug/clinical candidates (8 approved).

Europe PMC hits: **741**, kept 3.

1. **PMID 42101296** - *Reprogramming of Cellular Plasticity via ETS and MYC Core-Regulatory Circuits during Response to MAPK Inhibition in BRAF-Mutant Colorectal Cancer.* - Clinical cancer research : an official journal of the American Association for Cancer Research, 2026. doi:10.1158/1078-0432.ccr-25-4370
   > Building on this epigenetic vulnerability, bromodomain 2, a reader of H3K27ac-marked enhancers, was found to be synthetically lethal with BRAF + EGFR inhibition.

2. **PMID 42671662** - *BRAF V600E-Mutated Metastatic Colorectal Cancer: Practical Management in the Post-BREAKWATER Era.* - Journal of gastrointestinal cancer, 2026. doi:10.1007/s12029-026-01581-0
   > Therapeutic progress has been driven by the recognition that BRAF inhibition alone is inadequate in colorectal cancer because adaptive epidermal growth factor receptor (EGFR)-mediated feedback rapidly restores MAPK signaling.

3. **PMID 42487519** - *New Treatment Strategy and Future Research Direction for BRAF-Mutated Cancer.* - Cancer science, 2026. doi:10.1111/cas.70480
   > Treatment with BRAF inhibitor, and more recently, treatment with BRAF inhibitor plus MEK inhibitor was first developed in malignant melanoma.

<sub>query: `(TITLE_ABS:"BRAF" OR TITLE_ABS:"B-RAF1" OR TITLE_ABS:"BRAF-1" OR TITLE_ABS:"BRAF1" OR TITLE_ABS:"B-raf" OR TITLE_ABS:"RAFB1") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 7. TP53 (P04637) - tumor protein p53

Open Targets: overall **0.8347**, literature 0.997, known-drug (clinical) 0.106, 9 drug/clinical candidates (0 approved).

Europe PMC hits: **169**, kept 3.

1. **PMID 42409767** - *mTORC1 suppression by Trp53 mutation drives resistance to immune checkpoint blockade.* - Cell death & disease, 2026. doi:10.1038/s41419-026-09067-4
   > Mechanistically, Trp53 deletion downregulated mTORC1 inhibitor genes, leading to elevated mTORC1 signaling and diminished autophagy, which sensitized tumor cells to IFN-γ and TNF-α-induced apoptosis.
   <sub>matched via the alias `TRP53`, not the symbol TP53.</sub>

2. **PMID 42649998** - *Overcoming Resistance: Targeting Survivin-Driven Apoptotic Resistance Restores Irinotecan Sensitivity in TP53-Mutant Colorectal Cancer.* - Cancers, 2026. doi:10.3390/cancers18162687
   > Methods: Human CRC cell lines with varying TP53 status and murine tumor models were treated with irinotecan (or its active metabolite SN-38), the Survivin inhibitor YM-155, and CuET, alone or in combination.

3. **PMID 42297232** - *Regulation of EIF5A and hypusination by p53 determines colorectal cancer cell fitness.* - Cancer letters, 2026. doi:10.1016/j.canlet.2026.218682
   > Further, TP53 mutant CRC cells remained unresponsive to a DHPS inhibitor.

<sub>query: `(TITLE_ABS:"TP53" OR TITLE_ABS:"BCC7" OR TITLE_ABS:"BMFS5" OR TITLE_ABS:"LFS1" OR TITLE_ABS:"TRP53") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 8. PIK3CA (P42336) - phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha

Open Targets: overall **0.8178**, literature 0.990, known-drug (clinical) 0.186, 32 drug/clinical candidates (4 approved).

Europe PMC hits: **900**, kept 3.

1. **PMID 42357299** - *Synergistic Inhibition of Colorectal Cancer Growth by Combined PI3K and COX-2 Blockade in Cell Lines and Patient-Derived Organoids.* - Pharmaceutics, 2026. doi:10.3390/pharmaceutics18060683
   > Recent studies demonstrated a significant survival benefit from taking low-dose aspirin, a nonselective COX inhibitor, supporting further exploration of the synergistic effects of combined PI3Kα inhibitor (inavolisib) and COX-2 inhibitor (celecoxib) therapy.
   <sub>matched via the alias `PI3K`, not the symbol PIK3CA.</sub>

2. **PMID 42497063** - *DHODH is a synthetic lethal target in PIK3CA-mutant colorectal cancer.* - Cell reports, 2026. doi:10.1016/j.celrep.2026.117719
   > Here, we identify dihydroorotate dehydrogenase (DHODH) as a selective metabolic dependency in PIK3CA-mutant colorectal cancer (CRC).
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

3. **PMID 40635161** - *BCL-2 Family Inhibition Enhances mTORC1/2 Inhibition in PIK3CA-Mutant Colorectal Cancer.* - Molecular cancer therapeutics, 2025. doi:10.1158/1535-7163.mct-24-1096
   > However, resistance to single-agent PI3K pathway inhibitors has been observed across multiple clinical trials, necessitating the identification of combination therapies that overcome or prevent resistance to precision medicine strategies.
   <sub>matched via the alias `PI3K`, not the symbol PIK3CA.</sub>

<sub>query: `(TITLE_ABS:"PIK3CA" OR TITLE_ABS:"PI3K" OR TITLE_ABS:"CCM4" OR TITLE_ABS:"CWS5" OR TITLE_ABS:"PI3K-alpha" OR TITLE_ABS:"PI3Kalpha" OR TITLE_ABS:"p110-alpha") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 9. KRAS (P01116) - KRas proto-oncogene, GTPase

Open Targets: overall **0.8093**, literature 0.998, known-drug (clinical) 0.604, 3 drug/clinical candidates (2 approved).

Europe PMC hits: **1004**, kept 3.

1. **PMID 42209458** - *Efficacy of dual KRASG12D-EGFR blockade versus triple combinations in patient-derived models of KRASG12D-mutant colorectal cancer.* - Cell death & disease, 2026. doi:10.1038/s41419-026-08900-0
   > Likewise, dual therapy with trametinib + cetuximab was as effective as the triple regimen, suggesting functional redundancy between direct KRAS inhibition and downstream MEK blockade when EGFR is co-targeted.

2. **PMID 42276046** - *The cartography of KRAS inhibitor resistance in colorectal cancer.* - Cancer cell, 2026. doi:10.1016/j.ccell.2026.05.010
   > demonstrate that in patients with advanced colorectal cancer, there are regionally distinct genomic and transcriptomic adaptive responses to KRAS G12C inhibition, including epithelial cell state changes and pro-inflammatory pathway enrichment, and provide insights into key resistance mechanisms to KRAS G12C inhibition.

3. **PMID 42178060** - *Macrophage-hitchhiking nanomedicine codelivering gemcitabine and KRAS G12D inhibitor orchestrates chemo-immunotherapy for metastatic colorectal cancer.* - Journal of controlled release : official journal of the Controlled Release Society, 2026. doi:10.1016/j.jconrel.2026.115045
   > In this contribution, we report a macrophage-hitchhiking micellar nanomedicine (GemkiM) that exploits peritoneal macrophage chemotaxis to selectively deliver a KRAS G12D inhibitor and a gemcitabine prodrug to mCRC lung metastases.

<sub>query: `(TITLE_ABS:"KRAS" OR TITLE_ABS:"K-RAS2A" OR TITLE_ABS:"K-RAS2B" OR TITLE_ABS:"K-RAS4A" OR TITLE_ABS:"K-RAS4B" OR TITLE_ABS:"K-Ras4B" OR TITLE_ABS:"KRAS1") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 10. EPCAM (P16422) - epithelial cell adhesion molecule

Open Targets: overall **0.7899**, literature 0.950, known-drug (clinical) 0.780, 10 drug/clinical candidates (3 approved).

Europe PMC hits: **40**, kept 3.

1. **PMID 42032842** - *Simultaneous Establishment of Autologous Colorectal Cancer and Mesothelial Stromal Cell Lines from Malignant Ascites Reveals a Mesothelial-Stromal FGFR3 Axis as a Potential Vulnerability in Peritoneal Metastasis.* - Cancer medicine, 2026. doi:10.1002/cam4.71804
   > Lineage marker analysis using qPCR demonstrated that OMUCR-1 selectively expressed epithelial markers (EPCAM, KRT20), whereas CAmeso strongly expressed mesothelial-mesenchymal markers (ACTA2, MSLN) and lacked epithelial marker expression.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

2. **PMID 39630205** - *Plant-derived EpCAM-Fc fusion proteins induce in vivo immune response to produce IgGs inhibiting invasion and migration of colorectal cancer cells.* - Plant cell reports, 2024. doi:10.1007/s00299-024-03377-7
   > In the wound healing assay, EpCAM-FcKP IgG showed higher migration inhibition compared to EpCAM-FcP IgG in both cell types, with similar results to EpCAM-FcM IgG in SW620 cells.

3. **PMID 40316530** - *NAG-1/GDF15 as a tumor suppressor in colorectal cancer: inhibition of β-catenin and NF-κB pathways via interaction with EpCAM.* - Cell death & disease, 2025. doi:10.1038/s41419-025-07695-w
   > Mechanistically, the pro-NAG-1/GDF15 interacts with EpCAM, preventing its cleavage and nuclear translocation, thereby reducing β-catenin and NF-κB activity.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"EPCAM" OR TITLE_ABS:"17-1A" OR TITLE_ABS:"323/A3" OR TITLE_ABS:"Ep-CAM" OR TITLE_ABS:"Ber-Ep4" OR TITLE_ABS:"BerEp4" OR TITLE_ABS:"CD326") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 11. SMAD4 (Q13485) - SMAD family member 4

Open Targets: overall **0.7883**, literature 0.981, known-drug (clinical) -, 0 drug/clinical candidates (0 approved).

Europe PMC hits: **72**, kept 3.

1. **PMID 42014683** - *TGF-β/SMAD4/14-3-3σ/TFEB axis promotes mesenchymal-epithelial transition and inhibits autophagy in colorectal cancer.* - Cell death & disease, 2026. doi:10.1038/s41419-026-08733-x
   > Inhibition of autophagy and promotion of MET by SMAD4 was mediated by inhibition of TFEB via binding and sequestration of TFEB by 14-3-3σ.

2. **PMID 42366126** - *Recombinant bone morphogenetic protein-2 attenuates colorectal cancer progression by orchestrating Hippo pathway activation, Yes-associated protein inhibition, and epithelial-mesenchymal transition suppression.* - The Korean journal of physiology & pharmacology : official journal of the Korean Physiological Society and the Korean Society of Pharmacology, 2026. doi:10.4196/kjpp.25.369
   > This effect was associated with the increased expression of p53, p21, and Smad4, while the levels of cyclin D1, cyclin-dependent kinase 4 (CDK4), and CDK6 decreased.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

3. **PMID 42527030** - *IL-27 shapes NK cell heterogeneity and function in colorectal cancer.* - Journal for immunotherapy of cancer, 2026. doi:10.1136/jitc-2025-014667
   > Methods: We employed an orthotopic transplantation model of genetically engineered colorectal tumor organoids Apc-/-KrasG12D/+Trp53R172H/-Smad4-/- (AKPS) to investigate NK cell heterogeneity, maturation, and function during CRC progression.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"SMAD4" OR TITLE_ABS:"DPC4" OR TITLE_ABS:"MADH4" OR TITLE_ABS:"hSMAD4" OR TITLE_ABS:"MYHRS") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 12. FBXO11 (Q86XK2) - F-box protein 11

Open Targets: overall **0.7870**, literature 0.229, known-drug (clinical) -, 0 drug/clinical candidates (0 approved).

Europe PMC hits: **1**, kept 1.

1. **PMID 42302642** - *FBXO11 enhances colorectal cancer progression through TMBIM6 stabilization: Echinacoside as a natural inhibitor.* - Phytomedicine : international journal of phytotherapy and phytopharmacology, 2026. doi:10.1016/j.phymed.2026.158271
   > Structure-based virtual screening identified Echinacoside (ECH) as an FBXO11 inhibitor, with binding confirmed by CETSA and DARTS assays.

<sub>query: `(TITLE_ABS:"FBXO11" OR TITLE_ABS:"FBX11" OR TITLE_ABS:"PRMT9" OR TITLE_ABS:"UBR6" OR TITLE_ABS:"UG063H01" OR TITLE_ABS:"VIT-1" OR TITLE_ABS:"VIT1") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 13. ATM (Q13315) - ATM serine/threonine kinase

Open Targets: overall **0.7831**, literature 0.944, known-drug (clinical) -, 0 drug/clinical candidates (0 approved).

Europe PMC hits: **84**, kept 3.

1. **PMID 42484294** - *CRISPR Screening Identifies SMARCAL1 and MRN as Modulators of WRN Dependency in MSI-H Colorectal Cancer.* - Cancer research, 2026. doi:10.1158/0008-5472.can-26-1023
   > Acute disruption of the MRN complex conferred profound resistance to WRN inhibition, whereas ATM deficiency produced a more modest resistant phenotype.

2. **PMID 41661672** - *FASN Inhibition Enhances the Efficacy of Chemotherapy in Colorectal Cancer by Inhibiting the DNA Damage Response.* - Cancer research, 2026. doi:10.1158/0008-5472.can-25-1917
   > In addition, FASN inhibitor treatment blocked DDR by decreasing ATM expression and CHK2 phosphorylation.

3. **PMID 42283915** - *TELO2-interacting protein 1 (TTI1), a novel Wnt/β-catenin target gene, decreases chemo-sensitivity in colorectal cancer by modulating DNA damage responses.* - Molecular biomedicine, 2026. doi:10.1186/s43556-026-00475-8
   > Mechanistically, TTI1 maintains the integrity of the TELO2-TTI1-TTI2 complex and stabilizes the DNA damage response kinases ATM and ATR.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"ATM" OR TITLE_ABS:"TEL1" OR TITLE_ABS:"TELO1" OR TITLE_ABS:"ATDC") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 14. MUTYH (Q9UIF7) - mutY DNA glycosylase

Open Targets: overall **0.7826**, literature 0.478, known-drug (clinical) -, 0 drug/clinical candidates (0 approved).

Europe PMC hits: **6**, kept 3.

1. **PMID 31377904** - *Efficacy of immune checkpoint blockade in MUTYH-associated hereditary colorectal cancer.* - Investigational new drugs, 2020. doi:10.1007/s10637-019-00842-z
   > Given that these tumor features are associated with the response to immune checkpoint inhibitors, we administered nivolumab to a CRC patient who carried two inactive MUTYH alleles (p.Y179C and p.G396D) and previously experienced failure of chemotherapy.

2. **PMID 36862359** - *Genotype-Phenotype Correlations in Autosomal Dominant and Recessive APC Mutation-Negative Colorectal Adenomatous Polyposis.* - Digestive diseases and sciences, 2023. doi:10.1007/s10620-023-07890-9
   > The genetic predisposition to APC (-)/CAP has largely been associated with germline mutations in some susceptible genes, including the human mutY homologue (MUTYH) gene and the Nth-like DNA glycosylase 1 (NTHL1) gene, and DNA mismatch repair (MMR) can cause autosomal recessive APC (-)/CAP.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

3. **PMID 34721767** - *Molecular testing for colorectal cancer: Clinical applications.* - World journal of gastrointestinal oncology, 2021. doi:10.4251/wjgo.v13.i10.1288
   > The inactivation of DNA mismatch repair (MMR), or MUTYH gene, or DNA polymerase epsilon results in excessive tumor mutational burden; these CRCs are highly antigenic and therefore sensitive to immune checkpoint inhibitors.

<sub>query: `(TITLE_ABS:"MUTYH" OR TITLE_ABS:"hMYH") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 15. EGFR (P00533) - epidermal growth factor receptor

Open Targets: overall **0.7732**, literature 0.994, known-drug (clinical) 0.978, 82 drug/clinical candidates (28 approved).

Europe PMC hits: **1256**, kept 3.

1. **PMID 42209458** - *Efficacy of dual KRASG12D-EGFR blockade versus triple combinations in patient-derived models of KRASG12D-mutant colorectal cancer.* - Cell death & disease, 2026. doi:10.1038/s41419-026-08900-0
   > We systematically evaluated the biochemical, biological, and therapeutic activity of single, dual, and triple regimens combining the KRASG12D inhibitor MRTX1133 with cetuximab (EGFR inhibitor), alpelisib (PI3Kα inhibitor), or trametinib (MEK inhibitor) in a panel of patient-derived tumoroids and xenografts (PDXs) from metastatic CRC.

2. **PMID 42616657** - *RAF/MEK clamp inhibition with avutometinib and cetuximab in RAF1 S257L-mutated, anti-EGFR-refractory metastatic colorectal cancer in Gaucher's disease.* - The oncologist, 2026. doi:10.1093/oncolo/oyag335
   > After progression on fluoropyrimidine-, oxaliplatin-, and irinotecan-based regimens (the latter with the anti-EGFR antibody panitumumab) complicated by recurrent neutropenic sepsis, he received off-label avutometinib and defactinib (an oral RAF/MEK clamp inhibitor plus a FAK inhibitor) combined with cetuximab.

3. **PMID 42101296** - *Reprogramming of Cellular Plasticity via ETS and MYC Core-Regulatory Circuits during Response to MAPK Inhibition in BRAF-Mutant Colorectal Cancer.* - Clinical cancer research : an official journal of the American Association for Cancer Research, 2026. doi:10.1158/1078-0432.ccr-25-4370
   > Building on this epigenetic vulnerability, bromodomain 2, a reader of H3K27ac-marked enhancers, was found to be synthetically lethal with BRAF + EGFR inhibition.

<sub>query: `(TITLE_ABS:"EGFR" OR TITLE_ABS:"ERBB1" OR TITLE_ABS:"HER1" OR TITLE_ABS:"NISBD2" OR TITLE_ABS:"PIG61" OR TITLE_ABS:"ERBB" OR TITLE_ABS:"ERRP") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 16. FGFR2 (P21802) - fibroblast growth factor receptor 2

Open Targets: overall **0.7719**, literature 0.686, known-drug (clinical) 0.955, 27 drug/clinical candidates (11 approved).

Europe PMC hits: **11**, kept 3.

1. **PMID 41695355** - *Identification of small molecule inhibitors targeting FGFR through molecular docking-based screening.* - Frontiers in oncology, 2026. doi:10.3389/fonc.2026.1733391
   > In addition, ZINC000101867325 is also predicted to target FGFR2 mutations in colorectal cancer patients.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

2. **PMID 40871637** - *Synergistic Anticancer Effects of Fibroblast Growth Factor Receptor Inhibitor and Cannabidiol in Colorectal Cancer.* - Nutrients, 2025. doi:10.3390/nu17162609
   > Results: FGFR expression patterns were confirmed in various cancer cell lines, with NCI-H716 showing high FGFR2 expression.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

3. **PMID 40205008** - *FGFR as a Predictive Marker for Targeted Therapy in Gastrointestinal Malignancies: A Systematic Review.* - Journal of gastrointestinal cancer, 2025. doi:10.1007/s12029-025-01214-y
   > Alteration forms like FGFR2 fusion or rearrangement are associated with CC, while FGFR2 amplification and FGFR2b overexpression are associated with GC/OC.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"FGFR2" OR TITLE_ABS:"FGFR-2" OR TITLE_ABS:"BFR-1" OR TITLE_ABS:"CD332" OR TITLE_ABS:"CEK3" OR TITLE_ABS:"CFD1" OR TITLE_ABS:"ECT1") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 17. POLE (Q07864) - DNA polymerase epsilon, catalytic subunit

Open Targets: overall **0.7718**, literature 0.946, known-drug (clinical) 0.500, 7 drug/clinical candidates (5 approved).

Europe PMC hits: **43**, kept 3.

1. **PMID 40710181** - *Unveiling ctDNA Response: Immune Checkpoint Blockade Therapy in a Patient with POLE Mutation-Associated Early-Onset Colon Cancer.* - Current oncology (Toronto, Ont.), 2025. doi:10.3390/curroncol32070370
   > DNA polymerase epsilon (POLE) mutations occur at a higher rate than average-onset colorectal cancer (AOCRC).
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

2. **PMID 42299417** - *Severe renal toxicity following adjuvant envafolimab in a patient with ultra-hypermutated (POLE) stage II colorectal cancer: a case report.* - AME case reports, 2026. doi:10.21037/acr-2025-268
   > However, the safety profile of immune checkpoint inhibitors (ICIs) in early-stage POLE-mutated CRC remains unclear, especially the association between ultra-high TMB and immune-related adverse events (irAEs) lacks clinical evidence.

3. **PMID 41728249** - *Aggressive Right-Sided Colon Cancer in a Young Adult: Triple-Whammy Mutations (POLE, KRAS, BRCA1/2) Highlight Emerging Genetic Associations.* - ACG case reports journal, 2026. doi:10.14309/crj.0000000000002016
   > We report a case of early-onset CRC in a young male harboring pathogenic variants in BRCA1, BRCA2, and POLE, with no personal or familial cancer history.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"POLE" OR TITLE_ABS:"POLE1" OR TITLE_ABS:"CRCS12" OR TITLE_ABS:"FILS" OR TITLE_ABS:"IMAGEI") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 18. TCF7L2 (Q9NQB0) - transcription factor 7 like 2

Open Targets: overall **0.7620**, literature 0.893, known-drug (clinical) -, 0 drug/clinical candidates (0 approved).

Europe PMC hits: **103**, kept 3.

1. **PMID 42201413** - *Berberine impedes the DNA damage repair to inhibit colorectal cancer by regulating the SOX17/TCF4/PIM3 axis.* - Molecular genetics and genomics : MGG, 2026. doi:10.1007/s00438-026-02439-7
   > BBR inhibited the DDR of CRC cells by inactivating the β-catenin/TCF4 pathway through the regulation of SOX17.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term; matched via the alias `TCF4`, not the symbol TCF7L2.</sub>

2. **PMID 42689907** - *Synthesis of benzenesulfonamide derivatives reveals anti-CCSCs activity by regulating Wnt/β-catenin pathway.* - Future medicinal chemistry, 2026. doi:10.1080/17568919.2026.2726190
   > Methods: Based on the reported structural characteristics of LF3 which is a small molecule inhibitor of β-catenin/Transcription Factor 4 (TCF4) and the significant role of the naphthoquinone group in anti-tumor activity, we designed and synthesized a series of benzenesulfonamide derivatives.
   <sub>matched via the alias `TCF4`, not the symbol TCF7L2.</sub>

3. **PMID 41418870** - *Mechanistic insights into hypoxia-induced TCF7L2 upregulation and its oncogenic effects on colorectal cancer.* - Experimental cell research, 2026. doi:10.1016/j.yexcr.2025.114868
   > Additionally, Western blot and experiments employing the PI3K inhibitor LY294002 have demonstrated that TCF7L2 activates the PI3K/AKT signaling pathway, thereby facilitating the proliferation of CRC cells.

<sub>query: `(TITLE_ABS:"TCF7L2" OR TITLE_ABS:"TCF-4" OR TITLE_ABS:"TCF4" OR TITLE_ABS:"hTCF-4") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 19. FBXW7 (Q969H0) - F-box and WD repeat domain containing 7

Open Targets: overall **0.7585**, literature 0.950, known-drug (clinical) -, 0 drug/clinical candidates (0 approved).

Europe PMC hits: **29**, kept 3.

1. **PMID 41312772** - *Maritoclax Overcomes FBW7 Deficiency-Driven Irinotecan Resistance in Colorectal Cancer by Targeting MCL1.* - Cancer medicine, 2025. doi:10.1002/cam4.71419
   > While irinotecan (via its active metabolite SN38) is a first-line TOP1 inhibitor for advanced CRC, the mechanistic link between FBW7 dysfunction and irinotecan resistance remains elusive.
   <sub>matched via the alias `FBW7`, not the symbol FBXW7.</sub>

2. **PMID 42319868** - *Chemical Induction of MYC Protein Degradation via MYC-MAX Disruption and 20S Proteasome Activation.* - ACS chemical biology, 2026. doi:10.1021/acschembio.6c00258
   > MYC loss is proteasome-dependent and persists following knockdown of FBXW7, indicating a degradation mechanism distinct from canonical SCF-FBXW7-mediated turnover and consistent with direct 20S proteasomal degradation.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

3. **PMID 41952439** - *PVT1-104aa derived from the 8q24 gene desert promotes colorectal cancer tumorigenesis.* - Clinical and translational medicine, 2026. doi:10.1002/ctm2.70654
   > Mechanistically, PVT1-104aa enhanced c-Myc phosphorylation at Ser62, disrupted the c-Myc-FBW7 interaction, and thereby inhibited ubiquitin-mediated degradation of c-Myc, as shown by accelerated c-Myc turnover upon PVT1-104aa knockdown.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term; matched via the alias `FBW7`, not the symbol FBXW7.</sub>

<sub>query: `(TITLE_ABS:"FBXW7" OR TITLE_ABS:"FBX30" OR TITLE_ABS:"FBXW6" OR TITLE_ABS:"CDC4" OR TITLE_ABS:"FBW6" OR TITLE_ABS:"FBW7" OR TITLE_ABS:"FLJ11071") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 20. POLD1 (P28340) - DNA polymerase delta 1, catalytic subunit

Open Targets: overall **0.7561**, literature 0.831, known-drug (clinical) 0.500, 7 drug/clinical candidates (5 approved).

Europe PMC hits: **51**, kept 3.

1. **PMID 40981426** - *The Molecular and Functional Landscape of Resistance to FOLFIRI Chemotherapy in Metastatic Colorectal Cancer.* - Cancer discovery, 2026. doi:10.1158/2159-8290.cd-24-0556
   > FOLFIRI-resistant models showed transcriptional upregulation of innate immunity and mitochondrial metabolism genes, together with reduced expression of the DNA polymerase POLD1.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

2. **PMID 41213597** - *[Surgical Management of Hereditary Colorectal Cancer Syndromes].* - Zentralblatt fur Chirurgie, 2026. doi:10.1055/a-2724-3658
   > For rare polyposis syndromes such as NTHL1-, POLE-, or POLD1-associated syndromes, evidence-based recommendations are lacking, and treatment should follow FAP/aFAP protocols.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

3. **PMID 40652318** - *Predictive biomarkers for immune checkpoint inhibition for patients with colorectal cancer: a comprehensive review.* - Immunotherapy, 2025. doi:10.1080/1750743x.2025.2530853
   > Except for the loss of mismatch repair protein (MMR) function and POLE/POLD1 mutations, most of the biomarkers of response are largely investigational.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"POLD1" OR TITLE_ABS:"POLD" OR TITLE_ABS:"CDC2" OR TITLE_ABS:"CRCS10" OR TITLE_ABS:"IMD120" OR TITLE_ABS:"MDPL") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 21. FLT4 (P35916) - fms related receptor tyrosine kinase 4

Open Targets: overall **0.7524**, literature 0.515, known-drug (clinical) 0.971, 47 drug/clinical candidates (18 approved).

Europe PMC hits: **12**, kept 3.

1. **PMID 41459512** - *Fruquintinib saddles tumor immune tolerance by curbing pro-tumoral immature myeloid cell populations.* - Frontiers in immunology, 2025. doi:10.3389/fimmu.2025.1699980
   > In this context, we investigated the effects of fruquintinib, a selective oral VEGFR3/FLT4/CD310 inhibitor with high affinity for VEGFR3/FLT4/CD310, on tumor growth and size on colorectal (MC38, CT26) and breast (4T1, E0771) tumors.

2. **PMID 39352649** - *Clinical research progress of fruquintinib in the treatment of malignant tumors.* - Investigational new drugs, 2024. doi:10.1007/s10637-024-01476-6
   > It effectively curtails tumor growth by binding to and inhibiting VEGFR-1, VEGFR-2, and VEGFR-3.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term; matched via the alias `VEGFR-3`, not the symbol FLT4.</sub>

3. **PMID 35681695** - *SAR131675, a VEGRF3 Inhibitor, Modulates the Immune Response and Reduces the Growth of Colorectal Cancer Liver Metastasis.* - Cancers, 2022. doi:10.3390/cancers14112715
   > This study tested the effectiveness of SAR131675, a selective VEGFR-3 tyrosine kinase inhibitor, to inhibit CLM in a murine model.
   <sub>matched via the alias `VEGFR-3`, not the symbol FLT4.</sub>

<sub>query: `(TITLE_ABS:"FLT4" OR TITLE_ABS:"FLT-4" OR TITLE_ABS:"FLT41" OR TITLE_ABS:"CHTD7" OR TITLE_ABS:"LMPH1A" OR TITLE_ABS:"LMPHM1" OR TITLE_ABS:"VEGFR-3") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 22. ERBB2 (P04626) - erb-b2 receptor tyrosine kinase 2

Open Targets: overall **0.7495**, literature 0.990, known-drug (clinical) 0.841, 48 drug/clinical candidates (18 approved).

Europe PMC hits: **235**, kept 3.

1. **PMID 42663748** - *HER2-Targeted Therapy for Metastatic Colorectal Cancer: Current Evidence, Patient Selection, and Practical Treatment Sequencing.* - Journal of gastrointestinal cancer, 2026. doi:10.1007/s12029-026-01564-1
   > Dual HER2 blockade with trastuzumab-based combinations and antibody-drug conjugates, particularly trastuzumab deruxtecan, have demonstrated clinically meaningful activity in treatment-refractory disease.
   <sub>matched via the alias `HER2`, not the symbol ERBB2.</sub>

2. **PMID 42117066** - *Comprehensive genomic landscape of ERBB2 in Chinese GI tumors: mutation-centered landscapes and precision treatment opportunities.* - Therapeutic advances in medical oncology, 2026. doi:10.1177/17588359261445706
   > Compared with amplification, oncogenic ERBB2 mutations are preferentially associated with higher TMB/MSI-H and characteristic co-mutation signatures, supporting the clinical evaluation of mutation-selective HER2 inhibitors and rational combinations with immune checkpoint blockade.

3. **PMID 42375093** - *A Novel DRD2 Antagonist, SD2-2305, Exerts Anticancer Effects in Colorectal Cancer Cells through G1 Arrest and Caspase-Dependent Apoptosis.* - Biomolecules & therapeutics, 2026. doi:10.4062/biomolther.2026.111
   > While key survival pathways (JAK2/STAT3, PI3K/Akt, and MAPK) remained relatively unaffected, SD2-2305 modulated growth factor receptors post-transcriptionally, decreasing HER2/ErbB2 and increasing TGF-beta receptor 1 expression.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"ERBB2" OR TITLE_ABS:"CD340" OR TITLE_ABS:"HER-2" OR TITLE_ABS:"HER-2/neu" OR TITLE_ABS:"HER2" OR TITLE_ABS:"MLN-19" OR TITLE_ABS:"MLN19") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 23. KDR (P35968) - kinase insert domain receptor

Open Targets: overall **0.7495**, literature 0.857, known-drug (clinical) 0.971, 75 drug/clinical candidates (24 approved).

Europe PMC hits: **231**, kept 3.

1. **PMID 42331382** - *IL-10 secretion by CD8 T cells orchestrates early NK cell recruitment and antitumor immunity after anti-VEGFR-2/PD-1 therapy.* - Journal for immunotherapy of cancer, 2026. doi:10.1136/jitc-2026-015243
   > Methods: Using preclinical models of colorectal cancer, we examined the immune changes elicited by combined vascular endothelial growth factor receptor-2 (VEGFR-2) and programmed cell death protein-1 (PD-1) blockade.
   <sub>matched via the alias `VEGFR-2`, not the symbol KDR.</sub>

2. **PMID 42657050** - *1-Butylindoline-2,3-dione derivatives as potent VEGFR-2 inhibitors: synthesis, mechanistic anticancer studies, and molecular dynamics simulations.* - RSC advances, 2026. doi:10.1039/d6ra05194d
   > A novel series of N-butyl isatin (indolin-2-one) derivatives was designed and synthesized as potential VEGFR-2-targeted anticancer agents, drawing inspiration from the oxindole-based inhibitor sunitinib.
   <sub>matched via the alias `VEGFR-2`, not the symbol KDR.</sub>

3. **PMID 42324832** - *Identification of Sulfonamide-Based Thiazolidine-2,4-Dione Derivatives as VEGFR-2 and Tumor-Associated hCA IX Inhibitors With Hypoxia-Targeted Anti-Colorectal Cancer Activity.* - Drug development research, 2026. doi:10.1002/ddr.70333
   > Given their critical roles in tumor progression and hypoxia adaptation, concurrent inhibition of carbonic anhydrase IX (hCA IX) and vascular endothelial growth factor receptor-2 (VEGFR-2) is a viable therapeutic approach.
   <sub>matched via the alias `VEGFR-2`, not the symbol KDR.</sub>

<sub>query: `(TITLE_ABS:"KDR" OR TITLE_ABS:"CD309" OR TITLE_ABS:"FLK-1" OR TITLE_ABS:"FLK1" OR TITLE_ABS:"VEGFR-2" OR TITLE_ABS:"VEGFR2" OR TITLE_ABS:"VEGFR") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 24. FGFR3 (P22607) - fibroblast growth factor receptor 3

Open Targets: overall **0.7483**, literature 0.458, known-drug (clinical) 0.817, 25 drug/clinical candidates (11 approved).

Europe PMC hits: **19**, kept 3.

1. **PMID 41237766** - *Paneth-like transition drives resistance to dual targeting of KRAS and EGFR in colorectal cancer.* - Cancer cell, 2026. doi:10.1016/j.ccell.2025.10.010
   > Genetic or pharmacological inhibition of FGFR3 prevents the Paneth-like transition, restores drug sensitivity, and synergizes with KRAS-EGFR inhibition across multiple preclinical models.

2. **PMID 42032842** - *Simultaneous Establishment of Autologous Colorectal Cancer and Mesothelial Stromal Cell Lines from Malignant Ascites Reveals a Mesothelial-Stromal FGFR3 Axis as a Potential Vulnerability in Peritoneal Metastasis.* - Cancer medicine, 2026. doi:10.1002/cam4.71804
   > Treatment with the FGFR inhibitor BGJ398 reduced tumor growth and decreased stromal FGFR3-positive components, suggesting that stromal FGFR3 may represent a potential microenvironmental vulnerability in CRC with peritoneal dissemination.

3. **PMID 41106751** - *PDZK1, a direct target of Gli2, promotes FGFR3 trafficking and cell proliferation in colorectal cancer.* - Cellular signalling, 2025. doi:10.1016/j.cellsig.2025.112175
   > In colorectal cancer (CRC), dysregulated Hedgehog (Hh) and fibroblast growth factor receptor 3 (FGFR3)-related pathways drive tumor progression, but their molecular crosstalk remains poorly understood.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"FGFR3" OR TITLE_ABS:"FGFR-3" OR TITLE_ABS:"CD333" OR TITLE_ABS:"CEK2" OR TITLE_ABS:"HSFGFR3EX" OR TITLE_ABS:"JTK4") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## 25. BRCA2 (P51587) - BRCA2 DNA repair associated

Open Targets: overall **0.7477**, literature 0.845, known-drug (clinical) -, 0 drug/clinical candidates (0 approved).

Europe PMC hits: **25**, kept 3.

1. **PMID 41728249** - *Aggressive Right-Sided Colon Cancer in a Young Adult: Triple-Whammy Mutations (POLE, KRAS, BRCA1/2) Highlight Emerging Genetic Associations.* - ACG case reports journal, 2026. doi:10.14309/crj.0000000000002016
   > We report a case of early-onset CRC in a young male harboring pathogenic variants in BRCA1, BRCA2, and POLE, with no personal or familial cancer history.
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

2. **PMID 39602992** - *Discovery of a potent PARP1 PROTAC as a chemosensitizer for the treatment of colorectal cancer.* - European journal of medicinal chemistry, 2025. doi:10.1016/j.ejmech.2024.117062
   > Targeting PARP1 blockade exhibit specific toxicity towards tumor cells with BRCA1 or BRCA2 mutations through synthetic lethality.

3. **PMID 38067284** - *Epigenetically Downregulated Breast Cancer Gene 2 through Acetyltransferase Lysine Acetyltransferase 2B Increases the Sensitivity of Colorectal Cancer to Olaparib.* - Cancers, 2023. doi:10.3390/cancers15235580
   > We report that downregulating histone acetyltransferase KAT2B decreases BRCA2 expression by reducing the acetylation of the 27th amino acid in histone H3 (H3K27) at the promoter of the BRCA2 gene in colorectal cancer (CRC).
   <sub>claim tier `target-only`: this sentence is the strongest the abstract offered - it does not pair the target with an inhibitor/binder/antagonist term.</sub>

<sub>query: `(TITLE_ABS:"BRCA2" OR TITLE_ABS:"BRCC2" OR TITLE_ABS:"BROVCA2" OR TITLE_ABS:"FAD1" OR TITLE_ABS:"FANCD1" OR TITLE_ABS:"GLM3" OR TITLE_ABS:"PNCA2") AND TITLE_ABS:"colorectal cancer" AND (TITLE_ABS:"inhibitor" OR TITLE_ABS:"inhibition" OR TITLE_ABS:"binder" OR TITLE_ABS:"antagonist" OR TITLE_ABS:"blockade") AND (HAS_ABSTRACT:Y) AND (SRC:MED)`</sub>

## Verification

73 PMIDs emitted, 73 passed the Europe PMC re-fetch, 73 survived both checks, 0 rejected.

No mismatches: every emitted PMID resolved and returned the same title from Europe PMC, and matched PubMed.

Cross-source check (PubMed E-utilities esummary): 73/73 titles matched PubMed esummary, 0 dropped.

Claim tiers across the 73 quotes: `target+term` 41, `target-only` 32. A `target-only` or `first-sentence` quote is still verbatim, but the abstract did not put the target and an inhibitor term in one sentence - usually because the target has no inhibitor, which is itself the finding.

1 title(s) matched only after markup normalisation - same words, different rendering of superscripts or symbols by the two databases. These were kept, and both spellings are shown:

| source | pmid | as written here | as returned |
|---|---|---|---|
| PubMed | 42209458 | Efficacy of dual KRASG12D-EGFR blockade versus triple combinations in patient-derived models of KRASG12D-mutant colorectal cancer. | Efficacy of dual KRAS(G12D)-EGFR blockade versus triple combinations in patient-derived models of KRAS(G12D)-mutant colorectal cancer. |

## Not evaluated / why

Every target searched produced at least one verified PMID.
