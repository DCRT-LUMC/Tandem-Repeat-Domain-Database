class DataManager {
    constructor() {
        this.geneEligibilityMap = new Map();
        this.proteinData = [];
        this.uniqueProteins = {};
    }

    /**
     * Load and parse CSV data for gene eligibility
     */
    async loadGeneEligibilityData() {
        try {
            const response = await fetch('DCRT_Gene_list_August_2024.csv');
            if (!response.ok) {
                throw new Error(`Failed to fetch CSV: ${response.status}`);
            }
            
            const csvText = await response.text();
            this.parseGeneEligibilityCSV(csvText);
            
            console.log(`Loaded eligibility data for ${this.geneEligibilityMap.size} genes`);
        } catch (error) {
            console.error('Error loading gene eligibility data:', error);
            throw error;
        }
    }

    /**
     * Parse CSV text and populate gene eligibility map
     */
    parseGeneEligibilityCSV(csvText) {
        const lines = csvText.split('\n');
        const headers = lines[0].split('\t');
        
        const entityNameIndex = headers.findIndex(h => h.trim() === 'Entity Name');
        const eligibilityIndex = headers.findIndex(h => h.trim() === 'Eligibility');
        
        if (entityNameIndex === -1 || eligibilityIndex === -1) {
            console.warn('Required columns not found in CSV');
            return;
        }
        
        for (let i = 1; i < lines.length; i++) {
            const columns = lines[i].split('\t');
            if (columns.length > Math.max(entityNameIndex, eligibilityIndex)) {
                const geneName = columns[entityNameIndex]?.trim();
                const eligibility = columns[eligibilityIndex]?.trim();
                
                if (geneName && eligibility) {
                    this.geneEligibilityMap.set(geneName, eligibility);
                }
            }
        }
    }

    /**
     * Get eligibility status for a gene
     */
    getGeneEligibility(geneName) {
        return this.geneEligibilityMap.get(geneName) || 'N/A';
    }

    /**
     * Load protein data from JSON
     */
    async loadProteinData() {
        try {
            const response = await fetch('./all_annotated_repeats.json');
            if (!response.ok) {
                throw new Error(`Network response was not ok: ${response.status}`);
            }
            
            const data = await response.json();
            console.log(`Data loaded successfully: ${data.length} records found`);
            
            this.proteinData = data;
            this.processProteinData();
            
            return this.getUniqueProteinsArray();
        } catch (error) {
            console.error('Error loading protein data:', error);
            throw error;
        }
    }

    /**
     * Process raw protein data into unique proteins with metadata
     */
    processProteinData() {
        const uniqueProteins = {};
        
        this.proteinData.forEach(item => {
            const uniqueKey = item.uniProtId;
            const aliases = this.extractAliases(item);
            const proteinMetadata = this.createProteinMetadata(item, aliases);
            
          if (!uniqueProteins[uniqueKey]) {
            const initialRepeatType = item.repeatType || 'Unknown';
            uniqueProteins[uniqueKey] = {
                ...proteinMetadata,
                repeats: [item],
                repeatCount: 1,
                repeatTypeCounts: {
                    [initialRepeatType]: 1
                },
                eligibility: this.getGeneEligibility(item.geneName || 'Unknown')
            };
        }else {
                this.mergeProteinData(uniqueProteins[uniqueKey], item, aliases);
            }
        });
        
        // Calculate the strict combined property for each unique protein
        Object.values(uniqueProteins).forEach(protein => {
            protein.hasOnlyNonSpanningInFrameRepeats = this.checkHasPerfectTargetExon(protein.repeats);
        });
        
        this.uniqueProteins = uniqueProteins;
    }

    /**
     * Extract aliases from protein item
     */
    extractAliases(item) {
        const aliases = new Set();
        
        if (item.aliases) {
            if (typeof item.aliases === 'string') {
                item.aliases.split(',').forEach(alias => {
                    const trimmed = alias.trim();
                    if (trimmed) aliases.add(trimmed);
                });
            } else if (Array.isArray(item.aliases)) {
                item.aliases.forEach(alias => {
                    if (alias && alias.trim()) aliases.add(alias.trim());
                });
            }
        }
        
        return Array.from(aliases).filter(a => a && a !== 'Unknown');
    }

    /**
     * Create protein metadata object
     */
    createProteinMetadata(item, aliases) {
        const hasInFrameExons = this.checkInFrameExons(item);
        const spansExons = item.blockCount > 1;
        const nonSpanningRepeat = item.blockCount === 1;
        const hasNonSpanningInFrameRepeats = this.checkNonSpanningInFrameRepeats(item);
        
        return {
            gene: item.geneName || 'Unknown',
            aliases: aliases,
            aliasString: aliases.length > 0 ? aliases.join(', ') : '',
            searchableText: this.createSearchableText(item, aliases),
            uniprotId: item.uniProtId || 'Unknown',
            repeatType: item.repeatType || 'Unknown',
            status: item.status || 'Unknown',
            chromosome: item.chrom || 'Unknown',
            hasInFrameExons: hasInFrameExons,
            hasSpanningExons: spansExons,
            hasNonSpanningRepeats: nonSpanningRepeat,
            hasNonSpanningInFrameRepeats: hasNonSpanningInFrameRepeats
        };
    }

    /**
     * Check if protein has in-frame exons
     */
    checkInFrameExons(item) {
        return item.ensembl_exon_info && 
               item.ensembl_exon_info.transcripts && 
               item.ensembl_exon_info.transcripts.some(t => 
                   t.containing_exons && 
                   t.containing_exons.some(e => 
                       e.frame_status === 'in_frame'
                   )
               );
    }

    /**
     * Check if protein has non-spanning repeats within in-frame AND fully coding exons
     */
    checkNonSpanningInFrameRepeats(item) {
        // First check if the repeat itself is non-spanning (blockCount = 1)
        if (!item.blockCount || item.blockCount !== 1) {
            return false;
        }
        
        // Then check if it has exon info and is within in-frame AND fully coding exons
        return item.ensembl_exon_info && 
               item.ensembl_exon_info.transcripts && 
               item.ensembl_exon_info.transcripts.some(t => 
                   t.containing_exons && 
                   t.containing_exons.some(e => 
                       e.frame_status === 'in_frame' && 
                       e.coding_status === 'fully_coding'
                   )
               );
    }

    /**
     * Check if protein has at least one "perfect target exon" that meets all 4 criteria:
     * 1. In-frame (reading frame aligns with CDS)
     * 2. Fully coding (entire exon is within CDS boundaries) 
     * 3. Contains at least one repeat with blockCount=1 (non-spanning repeat)
     * 4. Does NOT contain any repeat with blockCount>1 (no spanning repeats in that exon)
     */
    checkHasPerfectTargetExon(allRepeats) {
        // Build a map of exon_id -> repeats in that exon
        const exonToRepeats = new Map();
        
        allRepeats.forEach(repeat => {
            if (repeat.ensembl_exon_info && repeat.ensembl_exon_info.transcripts) {
                repeat.ensembl_exon_info.transcripts.forEach(transcript => {
                    if (transcript.containing_exons) {
                        transcript.containing_exons.forEach(exon => {
                            const exonKey = `${transcript.transcript_id}_${exon.exon_id}`;
                            if (!exonToRepeats.has(exonKey)) {
                                exonToRepeats.set(exonKey, {
                                    exonInfo: exon,
                                    repeats: []
                                });
                            }
                            exonToRepeats.get(exonKey).repeats.push(repeat);
                        });
                    }
                });
            }
        });
        
        // Check each exon to see if it meets all 4 criteria
        for (const [exonKey, exonData] of exonToRepeats) {
            const exon = exonData.exonInfo;
            const repeatsInExon = exonData.repeats;
            
            // Criteria 1 & 2: Must be in-frame AND fully coding
            if (exon.frame_status !== 'in_frame' || exon.coding_status !== 'fully_coding') {
                continue;
            }
            
            // Criteria 3: Must contain at least one repeat with blockCount=1
            const hasNonSpanningRepeat = repeatsInExon.some(r => r.blockCount === 1);
            if (!hasNonSpanningRepeat) {
                continue;
            }
            
            // Criteria 4: Must NOT contain any repeat with blockCount>1
            const hasSpanningRepeat = repeatsInExon.some(r => r.blockCount > 1);
            if (hasSpanningRepeat) {
                continue;
            }
            
            // If we get here, this exon meets all 4 criteria!
            return true;
        }
        
        return false;
    }

    /**
     * Create searchable text for protein
     */
    createSearchableText(item, aliases) {
        return [
            item.geneName || '',
            item.uniProtId || '',
            ...aliases
        ].join(' ').toLowerCase();
    }

    /**
     * Merge additional protein data into existing entry
     */
    mergeProteinData(existingProtein, newItem, newAliases) {
        const repeatType = newItem.repeatType || 'Unknown';

        existingProtein.repeats.push(newItem);
        existingProtein.repeatCount++;
        
        existingProtein.repeatTypeCounts = existingProtein.repeatTypeCounts || {};
        existingProtein.repeatTypeCounts[repeatType] = (existingProtein.repeatTypeCounts[repeatType] || 0) + 1;
        
        // Merge aliases
        const existingAliases = new Set(existingProtein.aliases);
        newAliases.forEach(alias => existingAliases.add(alias));
        existingProtein.aliases = Array.from(existingAliases);
        existingProtein.aliasString = Array.from(existingAliases).length > 0 ? 
            Array.from(existingAliases).join(', ') : '';
        existingProtein.searchableText = this.createSearchableText(
            { geneName: existingProtein.gene, uniProtId: existingProtein.uniprotId },
            Array.from(existingAliases)
        );
        
        // Update properties if found
        if (this.checkInFrameExons(newItem)) {
            existingProtein.hasInFrameExons = true;
        }
        if (newItem.blockCount > 1) {
            existingProtein.hasSpanningExons = true;
        }
        if (newItem.blockCount === 1) {
            existingProtein.hasNonSpanningRepeats = true;
        }
        if (this.checkNonSpanningInFrameRepeats(newItem)) {
            existingProtein.hasNonSpanningInFrameRepeats = true;
        }
    }

    /**
     * Get unique proteins as array
     */
    getUniqueProteinsArray() {
        return Object.values(this.uniqueProteins);
    }

    /**
     * Get unique values for filters
     */
    getFilterOptions() {
        const repeatTypes = new Set();
        const statuses = new Set();
        const eligibilityValues = new Set();
        
        this.proteinData.forEach(item => {
            repeatTypes.add(item.repeatType || 'Unknown');
            statuses.add(item.status || 'Unknown');
        });
        
        this.getUniqueProteinsArray().forEach(protein => {
            eligibilityValues.add(protein.eligibility);
        });
        
        return {
            repeatTypes: Array.from(repeatTypes).sort(),
            statuses: Array.from(statuses),
            eligibilityValues: Array.from(eligibilityValues).sort((a, b) => {
                if (a === 'N/A') return 1;
                if (b === 'N/A') return -1;
                return a.localeCompare(b);
            })
        };
    }
}
