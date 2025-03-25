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
    chrom VARCHAR(10),
    chrom_start INTEGER,
    chrom_end INTEGER,
    strand CHAR(1),
    repeat_length INTEGER,
    position TEXT,
    reserved TEXT,
    block_count INTEGER,
    block_sizes TEXT,
    block_starts TEXT,
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id)
);

-- Transcripts table
CREATE TABLE transcripts (
    transcript_id VARCHAR(20) PRIMARY KEY,
    gene_id INTEGER,
    description TEXT,
    ensembl_transcript_id VARCHAR(50),
    versioned_transcript_id VARCHAR(50),
    transcript_name VARCHAR(100),
    is_canonical BOOLEAN,
    biotype VARCHAR(50),
    exon_count INTEGER,
    FOREIGN KEY (gene_id) REFERENCES genes(gene_id)
);

-- Relationship between repeats and transcripts
CREATE TABLE repeat_transcripts (
    repeat_id INTEGER,
    transcript_id VARCHAR(20),
    genomic_start INTEGER,
    genomic_end INTEGER,
    exon_mapping TEXT,
    location VARCHAR(50),
    PRIMARY KEY (repeat_id, transcript_id),
    FOREIGN KEY (repeat_id) REFERENCES repeats(repeat_id),
    FOREIGN KEY (transcript_id) REFERENCES transcripts(transcript_id)
);

-- Exons table
CREATE TABLE exons (
    exon_id INTEGER PRIMARY KEY,
    ensembl_exon_id VARCHAR(50),
    phase INTEGER,
    end_phase INTEGER,
    frame_status VARCHAR(50),
    coding_status VARCHAR(50),
    utr_status VARCHAR(50),
    coding_percentage FLOAT
);

-- Transcript exons junction table
CREATE TABLE transcript_exons (
    transcript_id VARCHAR(20),
    exon_id INTEGER,
    exon_number INTEGER,
    overlap_bp INTEGER,
    exon_position_in_transcript VARCHAR(50),
    overlap_percentage FLOAT,
    PRIMARY KEY (transcript_id, exon_id),
    FOREIGN KEY (transcript_id) REFERENCES transcripts(transcript_id),
    FOREIGN KEY (exon_id) REFERENCES exons(exon_id)
);

-- Relationship between repeats and exons
CREATE TABLE repeat_exons (
    repeat_id INTEGER,
    exon_id INTEGER,
    PRIMARY KEY (repeat_id, exon_id),
    FOREIGN KEY (repeat_id) REFERENCES repeats(repeat_id),
    FOREIGN KEY (exon_id) REFERENCES exons(exon_id)
);

-- Create indexes for improved query performance
CREATE INDEX idx_exons_ensembl_id ON exons(ensembl_exon_id);
CREATE INDEX idx_repeats_protein_id ON repeats(protein_id);
CREATE INDEX idx_repeats_type ON repeats(repeat_type);
CREATE INDEX idx_transcript_exons_transcript ON transcript_exons(transcript_id);
CREATE INDEX idx_transcript_exons_exon ON transcript_exons(exon_id);
CREATE INDEX idx_repeat_exons_repeat ON repeat_exons(repeat_id);
CREATE INDEX idx_repeat_exons_exon ON repeat_exons(exon_id);