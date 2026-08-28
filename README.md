# TReXomeDB

### Tandem Repeat eXome DataBase

TReXomeDB is a web-accessible resource for exploring human proteins containing tandem repeat domains and for identifying exons that may be of interest as potential targets for antisense oligonucleotide (ASO)-mediated exon skipping.

Tandem repeat domains consist of consecutive repetitive protein units and can, in some proteins, provide a degree of structural and functional redundancy. This raises the possibility that removal of individual repeat units may sometimes be tolerated while preserving sufficient protein structure or function. TReXomeDB was developed to support the initial identification and prioritisation of such candidate exons.

**Website:** https://trexomedb.rnatherapy.nl/

**GitHub repository:** https://github.com/DCRT-LUMC/Tandem-Repeat-Domain-Database

---

## About TReXomeDB

Antisense oligonucleotide-mediated exon skipping is an established therapeutic strategy for selected genetic diseases. However, exons encoding functional protein domains are often considered unsuitable targets because removal of the encoded domain may disrupt protein function.

TReXomeDB focuses on a potential exception to this principle: **protein tandem repeat domains**. Because these domains are composed of repetitive structural units, removal of an individual repeat may, in selected proteins, be compatible with preservation of partial or sufficient protein function.

TReXomeDB links protein tandem repeat annotations to their corresponding human transcripts and exons, allowing users to explore exon-level characteristics that may be relevant when prioritising potential exon-skipping targets.

The database is intended as an **initial candidate-prioritisation resource**. Identification of an exon in TReXomeDB does not demonstrate that the exon can be therapeutically skipped or that the resulting protein will retain sufficient function. Structural analysis and experimental validation are required before therapeutic development.

---

## What does TReXomeDB contain?

TReXomeDB integrates information at the protein, transcript, and exon levels, including:

- protein tandem repeat annotations
- tandem repeat type and repeat boundaries
- UniProt identifiers and protein information
- UniProt review status
- Ensembl transcript and exon information
- exon coding status
- exon reading-frame status
- repeat and exon boundaries
- information on whether repeats span multiple exons
- DCRT gene-level therapeutic eligibility information
- AlphaFold protein structures, where available

The resource allows users to explore tandem repeat-containing proteins and inspect the genomic and exon-level context of individual repeat domains.

---

## Candidate exon prioritisation

TReXomeDB provides filters that can be combined to identify potential exon-skipping candidates.

Available filters include:

- repeat type
- UniProt review status
- gene-level therapeutic eligibility
- repeat count
- exon frame
- non-spanning repeats
- fully coding exons
- single-repeat exons

The default exon-selection criteria are intended to identify exons that are more compatible with the proposed exon-skipping strategy.

In general, the preferred candidate exon is:

1. in-frame, so that exon removal does not disrupt the translational reading frame;
2. fully coding;
3. associated with a tandem repeat domain;
4. not part of a repeat that spans multiple exons;
5. ideally encoding a complete repeat domain without additional non-redundant protein-domain regions.

These criteria represent an initial screen. The position of the repeat within the overall protein and repeat array, linker regions, protein structure, molecular interactions, and disease mechanism must also be considered.

---

## Database construction

The database was developed by integrating protein tandem repeat annotations with gene, transcript, and exon information.

The database construction described in the associated manuscript can be summarised as:

```text
UCSC Genome Browser
UniProt Repeat track (hg38)
          |
          v
    UniProt information
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
   Web interface / API
```

Repeat annotations were obtained from the **UCSC Genome Browser hg38 UniProt Repeat track**. At the time of access (5 February 2025), the source contained **18,834 repeat annotations**.

These annotations were integrated with gene-level information from UniProt and mapped to transcripts and exons using Ensembl. Additional exon-level annotations were then generated to support the prioritisation of potential exon-skipping candidates.

A detailed description of the database construction and annotation strategy is provided in the associated manuscript.

---

## Using TReXomeDB

The web interface can be used to:

1. search for a specific gene or protein;
2. browse proteins containing tandem repeat domains;
3. filter proteins according to repeat and gene-level characteristics;
4. open individual protein detail pages;
5. inspect tandem repeat regions in available protein structures;
6. examine genomic and exon-level information;
7. identify exons fulfilling selected candidate criteria.

Protein detail pages also provide links to external resources, including Ensembl, GeneCards, and the Protein Data Bank where applicable.

Where an AlphaFold prediction is available, the protein viewer can display the predicted structure and highlight repeat regions. Repeat boundaries and exon boundaries can be displayed to help interpret the structural and genomic context of the repeat domains.

---

## API

TReXomeDB provides a public read-only API for programmatic access to the database.

The current API provides endpoints for accessing repeat-level and protein-level information, including filtering and pagination.

Main endpoints include:

```text
GET /api/repeats
GET /api/proteins
```

The API source code is located in:

```text
api/
```

The API documentation is available through the deployed application.

---

## Repository structure

The repository contains the TReXomeDB web application, API, data-processing scripts, source data, and example material.

```text
Tandem-Repeat-Domain-Database/
│
├── README.md
├── LICENSE
├── .gitignore
├── requirements.txt
│
├── api/
├── R/
├── col_data/
├── data/
├── docs/
├── scripts/
├── Examples/
└── Archive/
```

### `api/`

Contains the FastAPI application used to provide programmatic access to TReXomeDB.

### `R/`

Contains R-based development code related to the database processing workflow.

### `col_data/`

Contains data-processing scripts and examples used during the preparation and annotation of the database data.

### `data/`

Contains source and reference data used by the project.

### `docs/`

Contains the static web interface and associated website files.

### `scripts/`

Contains utility scripts used for processing repeat information and generating or analysing exon-removal sequences.

### `Examples/`

Contains example material associated with structural and exon-removal analyses.

### `Archive/`

Contains older development and prototype material retained for reference but not required for the main TReXomeDB application.

---

## Data provenance

TReXomeDB integrates information derived from several external resources, including:

- UCSC Genome Browser
- UniProt
- Ensembl
- AlphaFold

The version of these external resources is important when reproducing or extending the database because their annotations are updated over time.

The database described in the associated manuscript was constructed using the **GRCh38/hg38** genome build and source data available at the time of database development.

The manuscript reports the use of the UCSC Genome Browser UniProt Repeat track accessed on **5 February 2025**, containing 18,834 repeat annotations at that time.

---

## Limitations

TReXomeDB is intended to support **initial candidate prioritisation** and should not be interpreted as evidence that an exon is experimentally validated or clinically suitable for exon skipping.

In particular, exon-level criteria alone cannot determine whether removal of a repeat will preserve protein folding, protein interactions, or disease-relevant function.

Additional considerations include:

- the position of a repeat within a repeat array;
- whether the repeat is located at the edge or within the middle of an array;
- the consequences of removing linker regions;
- the overall architecture of the protein;
- effects on protein-protein interactions;
- the biological function of the affected repeat;
- the underlying disease mechanism;
- the potential consequences of multi-exon skipping.

The manuscript therefore proposes TReXomeDB as an initial filtering and prioritisation tool. Candidate exons require further structural investigation and ultimately experimental validation before therapeutic development.

---

## Example applications

The associated study used TReXomeDB to identify and prioritise candidate exons in disease-associated genes and subsequently modelled selected exon-removal events using AlphaFold.

Two illustrative examples were:

- **KIDINS220**, containing ankyrin repeats, where selected exon removals retained a visually similar overall repeat-array structure in the in silico models.
- **IFT172**, containing WD and TPR domains, where removal of a candidate WD-repeat-containing exon demonstrated why repeat position and higher-order repeat architecture must also be considered.

These examples illustrate both the potential of targeting tandem repeat domains and the limitations of relying on exon-level criteria alone.

---

## Publication

TReXomeDB was developed as part of the following study:

**Morgan J., Faber O., Aartsma-Rus A., van Roon-Mom W.M.C., Lauffer M.C.**

*Protein Tandem Repeat Domains as Therapeutic Targets: Can we expand the Scope of Antisense Oligonucleotide Therapeutics?*

The manuscript describes the rationale for targeting tandem repeat domains through exon skipping, the development of TReXomeDB, the database construction process, and illustrative case studies.

---

## Citation

When using TReXomeDB in research, please cite the associated publication.

Please also acknowledge the TReXomeDB resource:

**TReXomeDB – Tandem Repeat eXome DataBase**  
Dutch Center for RNA Therapeutics  
Leiden University Medical Center

---

## License

TReXomeDB is distributed under the **GNU Affero General Public License v3.0 (AGPL-3.0)**.

See [`LICENSE`](LICENSE) for the full license text.
