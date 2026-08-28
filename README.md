# TReXomeDB

**Tandem Repeat eXome DataBase**

TReXomeDB is a web-accessible database for exploring human proteins containing tandem repeat domains and identifying exons that may be of interest as potential targets for antisense oligonucleotide (ASO)-mediated exon skipping.

Tandem repeat domains consist of consecutive repetitive protein units and can, in some proteins, provide a degree of structural and functional redundancy. This raises the possibility that removal of an individual repeat may be tolerated while preserving sufficient protein structure or function. TReXomeDB was developed to support the initial identification and prioritisation of such candidate exons.

**Website:** https://trexomedb.rnatherapy.nl/  
**GitHub:** https://github.com/DCRT-LUMC/Tandem-Repeat-Domain-Database

---

## About TReXomeDB

Exon skipping is an established therapeutic strategy for selected genetic diseases, but exons containing functional protein domains are often considered unsuitable targets because removing the encoded domain may disrupt protein function.

TReXomeDB focuses on an important potential exception: **protein tandem repeat domains**. In proteins containing repeated structural units, removal of an individual repeat may sometimes be compatible with preservation of partial or sufficient protein function.

The database therefore provides exon-level information that can be used as an initial filter when assessing potential exon-skipping targets.

Importantly, TReXomeDB provides **candidate prioritisation rather than proof of therapeutic suitability**. Exon-level criteria alone cannot determine whether removal of a repeat will preserve protein folding, interactions, or disease-relevant function. Structural modelling and experimental validation remain essential.

---

## What does TReXomeDB contain?

TReXomeDB integrates information on:

- human proteins containing annotated tandem repeat domains
- tandem repeat type and repeat boundaries
- UniProt identifiers and protein annotation
- Ensembl transcript and exon information
- exon coding status and reading frame status
- relationships between repeats and exon boundaries
- UniProt review status
- DCRT gene-level therapeutic eligibility information
- protein structures from AlphaFold, where available

The database can be explored at both the **protein level** and the **exon level**.

---

## Candidate exon prioritisation

TReXomeDB provides filters that can be combined to identify potential exon-skipping candidates.

These include:

- tandem repeat type
- UniProt review status
- gene-level therapeutic eligibility
- minimum repeat count
- in-frame exons
- fully coding exons
- repeats contained within a single exon
- exons containing a single repeat

A candidate exon that satisfies the default criteria is generally:

1. in-frame;
2. fully coding;
3. associated with a tandem repeat domain;
4. not affected by a repeat spanning multiple exons;
5. and, where applicable, contains a single complete repeat without additional unrelated protein-domain annotations.

These criteria are intended as an **initial screening step** and should be followed by assessment of the protein architecture, repeat-array position, linker regions, molecular interactions, and disease mechanism.

---

## Database construction

The database was developed by integrating repeat annotations with gene, transcript, and exon information.

The workflow described for the current resource is:

```text
UCSC Genome Browser
UniProt Repeat track (GRCh38 / hg38)
              |
              v
      UniProt protein data
              |
              v
     Ensembl transcript data
              |
              v
       Exon-level mapping
              |
              v
     TReXomeDB annotations
              |
              v
       Web interface + API
