class FilterManager {
    constructor(table, tableData) {
        this.table = table;
        this.tableData = tableData;
        this.customFilters = new Map();
    }

    /**
     * Initialize all filter controls
     */
    initializeFilters() {
        this.setupRepeatTypeFilter();
        this.setupStatusFilter();
        this.setupEligibilityFilter();
        this.setupRepeatCountFilter();
        this.setupInFrameExonFilter();
        this.setupExonSpanningFilter();
        this.setupFullyCodingFilter();
        this.setupSingleRepeatFilter();
        this.setupResetButton();
    }

    /**
     * Setup repeat type filter
     */
    setupRepeatTypeFilter() {
        $('#repeatTypeFilter').on('change', (event) => {
            const value = event.target.value;
            this.table.column(4).search(value).draw();
        });
    }

    /**
     * Setup status filter
     */
    setupStatusFilter() {
        $('#statusFilter').on('change', (event) => {
            const value = event.target.value;
            
            this.removeCustomFilter('statusFilter');
            
            if (value) {
                const statusFilter = (settings, searchData, index, rowData, counter) => {
                    const status = rowData.status ? rowData.status.toLowerCase() : '';
                    
                    if (value === 'reviewed') {
                        return status.includes('reviewed') && !status.includes('unreviewed');
                    } else if (value === 'unreviewed') {
                        return status.includes('unreviewed');
                    }
                    
                    return true;
                };
                
                this.addCustomFilter('statusFilter', statusFilter);
            }
            
            this.table.draw();
        });
    }

    /**
     * Setup eligibility filter
     */
    setupEligibilityFilter() {
        $('#eligibilityFilter').on('change', (event) => {
            const value = event.target.value;
            this.table.column(7).search(value).draw();
        });
    }

    /**
     * Setup repeat count filter
     */
    setupRepeatCountFilter() {
        $('#repeatCountCheck, #repeatCountValue').on('change', () => {
            this.removeCustomFilter('repeatCountFilter');
            
            if ($('#repeatCountCheck').is(':checked')) {
                const minCount = parseInt($('#repeatCountValue').val());
                
                const repeatCountFilter = (settings, data, dataIndex, rowData) => {
                    const counts = Object.values(rowData.repeatTypeCounts || {});
                    const maxCount = counts.length > 0 ? Math.max(...counts) : 0;
                    return maxCount >= minCount;
                };
                
                this.addCustomFilter('repeatCountFilter', repeatCountFilter);
            }
            
            this.table.draw();
        });
    }

    /**
     * Setup in-frame exon filter
     */
    setupInFrameExonFilter() {
        $('#inFrameExonCheck').on('change', (event) => {
            const filterActive = $(event.target).is(':checked');
            
            this.removeCustomFilter('inFrameExonFilter');
            
            if (filterActive) {
                console.log("Filter active. Searching for proteins with in-frame exons.");
                const proteinsWithInFrameExons = this.tableData.filter(p => p.hasInFrameExons).length;
                console.log(`Found ${proteinsWithInFrameExons} proteins with in-frame exons out of ${this.tableData.length} total.`);
                
                const inFrameExonFilter = (settings, searchData, index, rowData, counter) => {
                    return rowData.hasInFrameExons === true;
                };
                
                this.addCustomFilter('inFrameExonFilter', inFrameExonFilter);
            }
            
            this.table.draw();
        });
    }

    /**
     * Setup exon spanning filter
     */
    setupExonSpanningFilter() {
        $('#exonSpanningCheck').on('change', (event) => {
            this.removeCustomFilter('exonSpanningFilter');
            
            if ($(event.target).is(':checked')) {
                console.log("Non-spanning filter active. Showing proteins with non-spanning repeats (blockCount = 1).");
                const proteinsWithNonSpanningRepeats = this.tableData.filter(p => p.hasNonSpanningRepeats).length;
                console.log(`Found ${proteinsWithNonSpanningRepeats} proteins with non-spanning repeats out of ${this.tableData.length} total.`);
                
                const exonSpanningFilter = (settings, searchData, index, rowData, counter) => {
                    return rowData.hasNonSpanningRepeats === true;
                };
                
                this.addCustomFilter('exonSpanningFilter', exonSpanningFilter);
            }
            
            this.table.draw();
        });
    }

         /**
         * Setup fully coding filter
         */
        setupFullyCodingFilter() {
            $('#fullyCodingCheck').on('change', (event) => {
                this.removeCustomFilter('fullyCodingFilter');
                
                if ($(event.target).is(':checked')) {
                    console.log("Fully coding filter active. Showing proteins with fully coding exons.");
                    const proteinsWithFullyCodingExons = this.tableData.filter(p => p.hasFullyCodingExons).length;
                    console.log(`Found ${proteinsWithFullyCodingExons} proteins with fully coding exons out of ${this.tableData.length} total.`);
                    
                    const fullyCodingFilter = (settings, searchData, index, rowData, counter) => {
                        return rowData.hasFullyCodingExons === true;
                    };
                    
                    this.addCustomFilter('fullyCodingFilter', fullyCodingFilter);
                }
                
                this.table.draw();
            });
        }

            /**
         * Setup single-repeat exon filter
         */
        setupSingleRepeatFilter() {
            $('#singleRepeatCheck').on('change', (event) => {
                this.removeCustomFilter('singleRepeatFilter');
                
                if ($(event.target).is(':checked')) {
                    console.log("Single-repeat filter active. Showing proteins with at least one exon containing exactly one tandem repeat.");
                    
                    const singleRepeatFilter = (settings, searchData, index, rowData, counter) => {
                        return rowData.hasSingleRepeatExon === true;
                    };
                    
                    this.addCustomFilter('singleRepeatFilter', singleRepeatFilter);
                }
                
                this.table.draw();
            });
        }
    /**
     * Setup reset filters button
     */
    setupResetButton() {
        $('#resetFilters').on('click', () => {
            this.resetAllFilters();
        });
    }

    /**
     * Add a custom filter
     */
    addCustomFilter(name, filterFunction) {
        filterFunction.filterName = name;
        this.customFilters.set(name, filterFunction);
        $.fn.dataTable.ext.search.push(filterFunction);
    }

    /**
     * Remove a custom filter
     */
    removeCustomFilter(name) {
        if (this.customFilters.has(name)) {
            const filterFunction = this.customFilters.get(name);
            $.fn.dataTable.ext.search = $.fn.dataTable.ext.search.filter(
                func => func.filterName !== name
            );
            this.customFilters.delete(name);
        }
    }

    /**
     * Reset all filters
     */
    resetAllFilters() {
        // Reset form controls
        $('#repeatTypeFilter').val('');
        $('#statusFilter').val('');
        $('#eligibilityFilter').val('');
        $('#repeatCountCheck').prop('checked', false);
        $('#repeatCountValue').val('1');
        $('#inFrameExonCheck').prop('checked', false);
        $('#exonSpanningCheck').prop('checked', false);
        $('#fullyCodingCheck').prop('checked', false);
        $('#singleRepeatCheck').prop('checked', false);
        
        // Clear DataTable search and remove custom search functions
        this.table.search('').columns().search('').draw();
        
        // Clear all custom filters
        this.customFilters.clear();
        $.fn.dataTable.ext.search = $.fn.dataTable.ext.search.filter(
            func => !func.filterName
        );
    }

    /**
     * Populate filter dropdowns
     */
    populateFilterDropdowns(filterOptions) {
        this.populateRepeatTypeFilter(filterOptions.repeatTypes);
        this.populateEligibilityFilter(filterOptions.eligibilityValues);
    }

    /**
     * Populate repeat type filter dropdown
     */
    populateRepeatTypeFilter(repeatTypes) {
        const repeatTypeFilter = document.getElementById('repeatTypeFilter');
        repeatTypeFilter.innerHTML = '<option value="">All Types</option>';
        
        repeatTypes.forEach(type => {
            const option = document.createElement('option');
            option.value = type;
            option.textContent = type;
            repeatTypeFilter.appendChild(option);
        });
    }

    /**
     * Populate eligibility filter dropdown
     */
    populateEligibilityFilter(eligibilityValues) {
        const eligibilityFilter = document.getElementById('eligibilityFilter');
        eligibilityFilter.innerHTML = '<option value="">All</option>';
        
        eligibilityValues.forEach(value => {
            const option = document.createElement('option');
            option.value = value;
            option.textContent = value;
            eligibilityFilter.appendChild(option);
        });
    }
}
