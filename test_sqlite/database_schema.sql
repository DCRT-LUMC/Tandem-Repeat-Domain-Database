-- Database schema for Tandem Repeat Domain Database

-- Genes table
CREATE TABLE genes (
    gene_id INTEGER PRIMARY KEY,
    gene_name VARCHAR(50) NOT NULL,
    chromosome VARCHAR(10),
    location VARCHAR(50),
    gene_type VARCHAR(100),
    ensembl_gene_id VARCHAR(50)
);

-- Gene aliases table
CREATE TABLE gene_aliases (
    alias_id INTEGER PRIMARY KEY,
    gene_id INTEGER NOT NULL,
    alias_name VARCHAR(100) NOT NULL,
    FOREIGN KEY (gene_id) REFERENCES genes(gene_id)
);

-- Proteins table 
CREATE TABLE proteins (
    protein_id VARCHAR(20) PRIMARY KEY,
    gene_id INTEGER,
    length INTEGER,
    description TEXT,
    status VARCHAR(100),
    FOREIGN KEY (gene_id) REFERENCES genes(gene_id)
);

-- Repeats table
CREATE TABLE repeats (
    repeat_id INTEGER PRIMARY KEY,
    protein_id VARCHAR(20),
    repeat_type VARCHAR(50),
    start_pos INTEGER,
    end_pos INTEGER,
    sequence TEXT,
    chrom VARCHAR(10),
    chrom_start INTEGER,
    chrom_end INTEGER,
    strand CHAR(1),
    repeat_length INTEGER,
    position TEXT,
    reserved TEXT,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id)
);

-- Transcripts table
CREATE TABLE transcripts (
    transcript_id VARCHAR(20) PRIMARY KEY,
    gene_id INTEGER,
    description TEXT,
    ensembl_transcript_id VARCHAR(50),
    FOREIGN KEY (gene_id) REFERENCES genes(gene_id)
);

-- Relationship between repeats and transcripts
CREATE TABLE repeat_transcripts (
    repeat_id INTEGER,
    transcript_id VARCHAR(20),
    genomic_start INTEGER,
    genomic_end INTEGER,
    exon_mapping TEXT,
    PRIMARY KEY (repeat_id, transcript_id),
    FOREIGN KEY (repeat_id) REFERENCES repeats(repeat_id),
    FOREIGN KEY (transcript_id) REFERENCES transcripts(transcript_id)
);

-- Exons table
CREATE TABLE exons (
    exon_id INTEGER PRIMARY KEY,
    repeat_id INTEGER,
    block_count INTEGER,
    block_sizes TEXT,
    block_starts TEXT,
    ensembl_exon_id VARCHAR(50),
    phase INTEGER,
    end_phase INTEGER,
    frame_status VARCHAR(50),
    ensembl_info TEXT,
    FOREIGN KEY (repeat_id) REFERENCES repeats(repeat_id)
);

-- Genomic coordinates table to store detailed location information
CREATE TABLE genomic_coordinates (
    coordinate_id INTEGER PRIMARY KEY,
    repeat_id INTEGER,
    chrom_start INTEGER,
    chrom_end INTEGER,
    strand CHAR(1),
    FOREIGN KEY (repeat_id) REFERENCES repeats(repeat_id)
);