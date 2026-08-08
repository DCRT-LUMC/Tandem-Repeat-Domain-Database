export class DataProcessor {
    constructor() {
        this.originalRepeatRegions = [];
        this.processedRepeatRegions = [];
        this.processedExonData = null;
    }

    async fetchProteinData(uniprotId, dataUrl) {
        try {
            const response = await fetch(dataUrl);
            const data = await response.json();
            
            const proteinRepeats = data.filter(item => item.uniProtId === uniprotId);
            
            if (!proteinRepeats.length) {
                throw new Error(`No data found for protein ${uniprotId}`);
            }
            
            this.originalRepeatRegions = proteinRepeats;
            return this.processProteinData(proteinRepeats);
        } catch (error) {
            console.error('Error loading data:', error);
            throw error;
        }
    }

    processProteinData(repeats) {
        try {
            // Sort repeats by position
            repeats.sort((a, b) => {
                const getStart = r => r.protein_start || 
                    (r.position && r.position.match(/amino acids (\d+)-(\d+)/) ? 
                     parseInt(r.position.match(/amino acids (\d+)-(\d+)/)[1]) : 0);
                
                return getStart(a) - getStart(b);
            });
            
            const repeatTypesPresent = [...new Set(repeats.map(r => r.repeatType).filter(Boolean))];

            const repeatTypeCounts = repeats.reduce((acc, repeat) => {
                const type = repeat.repeatType || 'Unknown';
                acc[type] = (acc[type] || 0) + 1;
                return acc;
            }, {});
            
            const repeatTypeSummary = repeatTypesPresent.length > 0
                ? repeatTypesPresent.join(' / ')
                : 'Unknown';
            
            const proteinInfo = {
                uniProtId: repeats[0].uniProtId,
                canonicalUniProtId: repeats[0].canonicalUniProtId || repeats[0].uniProtId || 'Unknown',
                uniProtDescription: repeats[0].uniProtDescription || 'Unknown',
                geneName: repeats[0].geneName || 'Unknown',
                repeatType: repeatTypeSummary,
                repeatTypesPresent: repeatTypesPresent,
                repeatTypeCounts: repeatTypeCounts,
                status: repeats[0].status || 'Unknown',
                chrom: repeats[0].chrom || 'Unknown',
                strand: repeats[0].strand || 'Unknown',
                aliases: repeats[0].aliases || []
            };
            // Extract genomic range
            if (repeats.length > 0) {
                const chromStarts = repeats.map(r => r.chromStart).filter(s => s !== undefined);
                const chromEnds = repeats.map(r => r.chromEnd).filter(e => e !== undefined);
                
                if (chromStarts.length > 0 && chromEnds.length > 0) {
                    const minStart = Math.min(...chromStarts);
                    const maxEnd = Math.max(...chromEnds);
                    proteinInfo.genomicRange = `${minStart}-${maxEnd}`;
                }
            }
            
            // Process repeat regions
            const repeatRegions = this.extractRepeatRegions(repeats);
            this.processedRepeatRegions = repeatRegions;
            
            // Process exon data
            const exonData = this.extractExonData(repeatRegions);
            if (exonData) {
                this.processedExonData = exonData;
            }
            
            return { proteinInfo, repeatRegions, exonData };
        } catch (error) {
            console.error("Error processing protein data:", error);
            throw error;
        }
    }

    extractRepeatRegions(repeats) {
    const repeatRegions = [];
    const colorOptions = ["#ff5f5f", "#5fba7d", "#5f87ff", "#ffaf5f", "#bf5fff", "#dc3545", "#fd7e14", "#ffc107", "#20c997", "#0dcaf0"];
    const typeCounters = {};

    repeats.forEach((repeat, index) => {
        const start = repeat.protein_start;
        const end = repeat.protein_end;
        const repeatType = repeat.repeatType || 'Unknown';

        typeCounters[repeatType] = (typeCounters[repeatType] || 0) + 1;
        const label = `${repeatType}${typeCounters[repeatType]}`;

        if (start && end) {
            const length = end - start + 1;
            const color = colorOptions[index % colorOptions.length];
            const colorHex = color.replace("#", "0x");

            repeatRegions.push({
                start,
                end,
                length,
                color,
                colorHex,
                position: repeat.position,
                ensembl_exon_info: repeat.ensembl_exon_info,
                chromStart: repeat.chromStart,
                chromEnd: repeat.chromEnd,
                repeatType: repeatType,
                blockCount: repeat.blockCount,
                label: label
            });
        } else {
            const posMatch = repeat.position && repeat.position.match(/amino acids (\d+)-(\d+)/);
            if (posMatch) {
                const start = parseInt(posMatch[1]);
                const end = parseInt(posMatch[2]);
                const length = end - start + 1;
                const color = colorOptions[index % colorOptions.length];
                const colorHex = color.replace("#", "0x");

                repeatRegions.push({
                    start,
                    end,
                    length,
                    color,
                    colorHex,
                    position: repeat.position,
                    ensembl_exon_info: repeat.ensembl_exon_info,
                    chromStart: repeat.chromStart,
                    chromEnd: repeat.chromEnd,
                    repeatType: repeatType,
                    blockCount: repeat.blockCount,
                    label: label
                });
            }
        }
    });

    return repeatRegions;
    }

    extractExonData(repeatRegions) {
        const exonData = { transcript: null, exons: [] };
        let hasExonData = false;
        
        // Check for exon data
        for (const repeat of repeatRegions) {
            if (repeat.ensembl_exon_info && repeat.ensembl_exon_info.transcripts && 
                repeat.ensembl_exon_info.transcripts.length > 0) {
                hasExonData = true;
                break;
            }
        }
        
        if (!hasExonData) {
            return null;
        }
        
        // Find canonical transcript
        for (const repeat of repeatRegions) {
            if (!repeat.ensembl_exon_info || !repeat.ensembl_exon_info.transcripts) continue;
            
            const canonical = repeat.ensembl_exon_info.transcripts.find(t => t.is_canonical);
            if (canonical) {
                exonData.transcript = canonical;
                break;
            }
        }
        
        // If no canonical found, use first transcript
        if (!exonData.transcript) {
            for (const repeat of repeatRegions) {
                if (repeat.ensembl_exon_info && repeat.ensembl_exon_info.transcripts && 
                    repeat.ensembl_exon_info.transcripts.length > 0) {
                    exonData.transcript = repeat.ensembl_exon_info.transcripts[0];
                    break;
                }
            }
        }
        
        if (!exonData.transcript) {
            return null;
        }
        
        // Collect exons
        const transcriptId = exonData.transcript.transcript_id;
        const exonMap = new Map();
        
        for (const repeat of repeatRegions) {
            if (!repeat.ensembl_exon_info || !repeat.ensembl_exon_info.transcripts) continue;
            
            const transcript = repeat.ensembl_exon_info.transcripts.find(t => t.transcript_id === transcriptId);
            if (!transcript || !transcript.containing_exons) continue;
            
            const repeatStart = repeat.protein_start || (repeat.position && repeat.position.match(/amino acids (\d+)-(\d+)/) ? parseInt(repeat.position.match(/amino acids (\d+)-(\d+)/)[1]) : null);
            const repeatEnd = repeat.protein_end || (repeat.position && repeat.position.match(/amino acids (\d+)-(\d+)/) ? parseInt(repeat.position.match(/amino acids (\d+)-(\d+)/)[2]) : null);
            
            if (!repeatStart || !repeatEnd) continue;
            
            const repeatLabel = repeat.label || `${repeat.repeatType || 'Unknown'}1`;
            const hasBlockCount1 = repeat.blockCount === 1;
            
            for (const exon of transcript.containing_exons) {
                const exonKey = exon.exon_id;
                
                if (!exonMap.has(exonKey)) {
                    exonMap.set(exonKey, {
                        exon_number: exon.exon_number,
                        exon_id: exon.exon_id,
                        position: exon.position,
                        coding_status: exon.coding_status || 'Unknown',
                        frame_status: exon.frame_status || 'Unknown',
                        start_phase: exon.phase !== undefined ? exon.phase : 'Unknown',
                        end_phase: exon.end_phase !== undefined ? exon.end_phase : 'Unknown',
                        repeat_regions: [{
                            label: repeatLabel,
                            repeatType: repeat.repeatType || 'Unknown',
                            start: repeatStart,
                            end: repeatEnd,
                            blockCount: repeat.blockCount
                        }],
                        repeat_start: exon.protein_start,
                        repeat_end: exon.protein_end,
                        overlap_percentage: exon.overlap_percentage,
                        exon_genomic_start: exon.exon_start,
                        exon_genomic_end: exon.exon_end,
                        has_blockcount_1: hasBlockCount1,
                        has_higher_blockcount: repeat.blockCount > 1
                    });
                } else {
                    const existingExon = exonMap.get(exonKey);
                    const hasRepeat = existingExon.repeat_regions.some(r => 
                        r.start === repeatStart && r.end === repeatEnd);
                    
                    if (!hasRepeat) {
                     existingExon.repeat_regions.push({
                            label: repeatLabel,
                            repeatType: repeat.repeatType || 'Unknown',
                            start: repeatStart,
                            end: repeatEnd,
                            blockCount: repeat.blockCount
                        });
                        
                        if (hasBlockCount1) {
                            existingExon.has_blockcount_1 = true;
                        }
                        
                        if (repeat.blockCount > 1) {
                            existingExon.has_higher_blockcount = true;
                        }
                        
                        exonMap.set(exonKey, existingExon);
                    }
                }
            }
        }
        
        exonData.exons = Array.from(exonMap.values()).sort((a, b) => a.exon_number - b.exon_number);
        
        return exonData.exons.length > 0 ? exonData : null;
    }

    getProcessedRepeatRegions() {
        return this.processedRepeatRegions;
    }

    getProcessedExonData() {
        return this.processedExonData;
    }

    getOriginalRepeatRegions() {
        return this.originalRepeatRegions;
    }
}
