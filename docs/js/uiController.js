export class UIController {
    constructor() {
        this.activeButton = null;
    }

    displayProteinInfo(info, repeats) {
        document.getElementById('proteinTitle').textContent = 
            `Focused view of ${info.repeatType} repeat domains in protein ${info.uniProtId} (${info.geneName})`;
            
        document.getElementById('uniprotLink').textContent = info.uniProtId;
        document.getElementById('uniprotLink').href = `https://www.uniprot.org/uniprotkb/${info.uniProtId}`;
        
        document.getElementById('geneName').textContent = info.geneName;
        document.getElementById('repeatType').textContent = info.repeatType;
        document.getElementById('status').textContent = info.status;
        
        // Add aliases to the protein info section
        document.getElementById('aliases').textContent = Array.isArray(info.aliases) ? 
            info.aliases.join(', ') : (info.aliases || 'None');
        
        // Genomic information (without aliases now)
        document.getElementById('chromosome').textContent = info.chrom;
        document.getElementById('strand').textContent = info.strand === '+' ? 'Forward (+)' : 'Reverse (-)';
        document.getElementById('genomicRange').textContent = info.genomicRange || 'Not available';
        
        // Update repeat summary
        document.getElementById('repeatSummary').textContent = 
            `Protein ${info.uniProtId} (${info.geneName}) contains ${repeats.length} ${info.repeatType} repeats in positions:`;
        
        // Create repeat list items
        const repeatList = document.getElementById('repeatList');
        repeatList.innerHTML = '';
        repeats.forEach(repeat => {
            const li = document.createElement('li');
            li.textContent = `${repeat.start}-${repeat.end} (${repeat.length} amino acids)`;
            repeatList.appendChild(li);
        });
        
        // Update external links
        this.updateExternalLinks(info);
        
        // Create legend for repeats
        this.createLegend(repeats);
    }

    updateExternalLinks(info) {
        const ensemblLink = document.getElementById('ensemblLink');
        ensemblLink.href = `https://www.ensembl.org/Homo_sapiens/Search/Results?q=${info.geneName}`;
        
        const genecardsLink = document.getElementById('genecardsLink');
        genecardsLink.href = `https://www.genecards.org/cgi-bin/carddisp.pl?gene=${info.geneName}`;
        
        const pdbLink = document.getElementById('pdbLink');
        pdbLink.href = `https://www.rcsb.org/search?query=${info.uniProtId}`;
    }

    createLegend(repeats) {
        const legend = document.getElementById('viewerLegend');
        
        // Clear existing items except the "Non-repeat regions"
        while (legend.childNodes.length > 1) {
            legend.removeChild(legend.firstChild);
        }
        
        // Add legend items for each repeat
        repeats.forEach((repeat, index) => {
            const item = document.createElement('div');
            item.className = 'legend-item';
            
            const colorBox = document.createElement('div');
            colorBox.className = 'color-box';
            colorBox.style.backgroundColor = repeat.color;
            
            const label = document.createElement('span');
            label.textContent = `${repeat.label} (${repeat.start}-${repeat.end})`;
            
            item.appendChild(colorBox);
            item.appendChild(label);
            
            // Insert at the top
            legend.insertBefore(item, legend.firstChild);
        });
    }

    updateExonDisplay(exonData) {
        const exonDataText = document.getElementById('exonDataText');

        if (exonData && exonData.exons && exonData.exons.length > 0) {
            const transcript = exonData.transcript;
            let exonHtml = `
                <p>Canonical transcript: <strong>${transcript?.transcript_id || 'Unknown'}</strong> (${transcript?.transcript_name || 'Unknown'})</p>
                <p>Exons overlapping repeats: <strong>${exonData.exons.length}</strong></p>
                
                <div class="card mb-3">
                    <div class="card-header">
                        <strong>Exon Highlighting Options</strong>
                    </div>
                    <div class="card-body">
                        <p class="mb-2">Select criteria for starring exons:</p>
                        <div class="form-check">
                            <input class="form-check-input exon-filter" type="checkbox" value="blockcount1" id="blockcount1" checked>
                            <label class="form-check-label" for="blockcount1">
                                Has blockCount = 1 only (non-spanning)
                            </label>
                        </div>
                        <div class="form-check">
                            <input class="form-check-input exon-filter" type="checkbox" value="blockcount1mixed" id="blockcount1mixed" checked>
                            <label class="form-check-label" for="blockcount1mixed">
                                Has blockCount = 1 (can include other blockCounts)
                            </label>
                        </div>
                        <div class="form-check">
                            <input class="form-check-input exon-filter" type="checkbox" value="inframe" id="inframe" checked>
                            <label class="form-check-label" for="inframe">
                                In-frame exons
                            </label>
                        </div>
                        <div class="form-check">
                            <input class="form-check-input exon-filter" type="checkbox" value="fullycoding" id="fullycoding" checked>
                            <label class="form-check-label" for="fullycoding">
                                Fully coding exons
                            </label>
                        </div>
                    </div>
                </div>
                
                <div class="mb-2">
                    <span class="legend-item exon-legend-item">
                        <i class="fas fa-star text-success"></i>&nbsp;&nbsp;Exons matching selected criteria
                    </span>
                </div>
                <table class="table table-sm table-striped">
                    <thead>
                        <tr>
                            <th>Exon #</th>
                            <th>Exon Boundaries</th>
                            <th>Coding Status</th>
                            <th>Start Phase</th>
                            <th>End Phase</th>
                            <th>Frame Status</th>
                            <th>Overlaps Repeats</th>
                            <th>Overlap %</th>
                        </tr>
                    </thead>
                    <tbody id="exonTableBody">
                        ${this.generateExonTableRows(exonData.exons)}
                    </tbody>
                </table>
            `;
            
            exonDataText.innerHTML = exonHtml;
            
            // Add event listeners to checkboxes
            document.querySelectorAll('.exon-filter').forEach(checkbox => {
                checkbox.addEventListener('change', () => {
                    this.updateExonTableStars(exonData.exons);
                });
            });
            
            // Apply stars immediately
            setTimeout(() => this.updateExonTableStars(exonData.exons), 0);
        } else {
            exonDataText.textContent = "No exon data available for this protein or the repeat regions.";
        }
    }

    generateExonTableRows(exons) {
        let rowsHtml = '';
        
        exons.forEach(exon => {
            const repeatLabels = exon.repeat_regions ? 
                exon.repeat_regions.map(r => r.label).join(", ") : 
                "Unknown";
            
            const frameStatus = this.formatFrameStatus(exon.frame_status);
            const startPhase = this.formatPhase(exon.start_phase);
            const endPhase = this.formatPhase(exon.end_phase);
            
            rowsHtml += `
                <tr data-exon-id="${exon.exon_id}" 
                    data-blockcount1="${exon.has_blockcount_1 === true && !exon.has_higher_blockcount}"
                    data-inframe="${(exon.frame_status === 'in_frame' || exon.frame_status === 'in-frame')}" 
                    data-fullycoding="${exon.coding_status === 'fully_coding'}">
                    <td class="exon-number">${exon.exon_number}</td>
                    <td>${exon.repeat_start}-${exon.repeat_end}</td>
                    <td>${exon.coding_status || 'Unknown'}</td>
                    <td>${startPhase}</td>
                    <td>${endPhase}</td>
                    <td>${frameStatus}</td>
                    <td>${repeatLabels}</td>
                    <td>${exon.overlap_percentage ? exon.overlap_percentage + '%' : 'Unknown'}</td>
                </tr>
            `;
        });
        
        return rowsHtml;
    }

    updateExonTableStars(exons) {
        const blockcount1Selected = document.getElementById('blockcount1').checked;
        const blockcount1MixedSelected = document.getElementById('blockcount1mixed').checked;
        const inframeSelected = document.getElementById('inframe').checked;
        const fullycodingSelected = document.getElementById('fullycoding').checked;
        
        const rows = document.querySelectorAll('#exonTableBody tr');
        rows.forEach(row => {
            const exonNumberCell = row.querySelector('.exon-number');
            
            const meetsBlockCountCriteria = !blockcount1Selected || row.dataset.blockcount1 === "true";
            const meetsBlockCountMixedCriteria = !blockcount1MixedSelected || row.dataset.blockcount1 === "true" || row.dataset.blockcount1 === "false";
            const meetsInFrameCriteria = !inframeSelected || row.dataset.inframe === "true";
            const meetsFullyCodingCriteria = !fullycodingSelected || row.dataset.fullycoding === "true";
            
            const isSpecialExon = meetsBlockCountCriteria && meetsBlockCountMixedCriteria && meetsInFrameCriteria && meetsFullyCodingCriteria;
            
            if (isSpecialExon) {
                const exonNumber = exonNumberCell.textContent.trim().split(' ')[0];
                exonNumberCell.innerHTML = `${exonNumber} <i class="fas fa-star text-success" title="Matches selected criteria"></i>`;
            } else {
                const exonNumber = exonNumberCell.textContent.trim().split(' ')[0];
                exonNumberCell.textContent = exonNumber;
            }
        });
    }

    formatPhase(phase) {
        return phase !== undefined ? phase : 'Unknown';
    }

    formatFrameStatus(status) {
        if (!status) return 'Unknown';
        if (status === 'in-frame') return 'In-frame';
        if (status === 'out-of-frame') return 'Out-of-frame';
        if (status === 'maintained') return 'Maintained';
        if (status === 'not-maintained') return 'Not maintained';
        return status;
    }

    setActiveButton(button) {
        try {
            document.querySelectorAll('.view-button').forEach(btn => {
                if (btn.id !== 'exonBoundariesBtn') {
                    btn.classList.remove('active');
                }
            });
            
            if (button.id !== 'exonBoundariesBtn') {
                button.classList.add('active');
            }
            this.activeButton = button;
        } catch (error) {
            console.error("Error in setActiveButton:", error);
        }
    }

    showError(message) {
        document.getElementById('loading').style.display = 'none';
        document.getElementById('errorContainer').style.display = 'block';
        document.getElementById('errorMessage').textContent = message;
    }

    showProteinDetails() {
        document.getElementById('loading').style.display = 'none';
        document.getElementById('proteinDetails').style.display = 'block';
    }

    toggleSection(sectionId) {
        const section = document.getElementById(sectionId);
        const button = event.currentTarget;
        
        if (section.classList.contains('collapsed')) {
            section.classList.remove('collapsed');
            button.innerHTML = button.innerHTML.replace('Show', 'Hide').replace('fa-eye', 'fa-eye-slash');
        } else {
            section.classList.add('collapsed');
            button.innerHTML = button.innerHTML.replace('Hide', 'Show').replace('fa-eye-slash', 'fa-eye');
        }
    }
}
