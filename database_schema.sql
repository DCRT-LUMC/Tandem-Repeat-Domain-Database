-- Database schema for Tandem Repeat Domain Database

-- Genes table
CREATE TABLE genes (
    gene_id INTEGER PRIMARY KEY,
    gene_name VARCHAR(50) NOT NULL,
    aliases TEXT,
    chromosome VARCHAR(10),
    location VARCHAR(50)
);

-- Proteins table 
CREATE TABLE proteins (
    protein_id VARCHAR(20) PRIMARY KEY,
    gene_id INTEGER,
    length INTEGER,
    description TEXT,
    status VARCHAR(20),
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
    FOREIGN KEY (protein_id) REFERENCES proteins(protein_id)
);

-- Transcripts table
CREATE TABLE transcripts (
    transcript_id VARCHAR(20) PRIMARY KEY,
    gene_id INTEGER,
    description TEXT,
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
