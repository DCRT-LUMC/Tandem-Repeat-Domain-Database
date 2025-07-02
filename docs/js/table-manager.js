class TableManager {
    constructor() {
        this.table = null;
        this.tableData = [];
    }

    /**
     * Initialize DataTable with configuration
     */
    initializeTable(tableData) {
        this.tableData = tableData;
        
        this.table = $('#proteinTable').DataTable({
            data: tableData,
            columns: this.getColumnConfiguration(),
            order: [[1, 'asc']],
            pageLength: 15,
            dom: 'lfrtip',
            responsive: true,
            columnDefs: this.getColumnDefinitions(),
            language: this.getLanguageConfiguration()
        });

        this.setupCustomSearch();
        return this.table;
    }

    /**
     * Get column configuration for DataTable
     */
    getColumnConfiguration() {
        return [
            { 
                data: 'chromosome',
                title: 'Chrom.'
            },
            { 
                data: 'gene',
                title: 'Gene Name'
            },
            { 
                data: 'aliasString',
                title: 'Aliases',
                render: this.renderAliases.bind(this)
            },
            { 
                data: 'uniprotId',
                title: 'UniProt ID'
            },
            { 
                data: 'repeatType',
                title: 'Repeat Type'
            },
            { 
                data: 'repeatCount',
                title: 'Repeat Count',
                render: this.renderRepeatCount
            },
            { 
                data: 'status',
                title: 'Status',
                render: this.renderStatus
            },
            { 
                data: 'eligibility',
                title: 'Eligibility',
                render: this.renderEligibility
            },
            {
                data: 'uniprotId',
                title: 'Actions',
                render: this.renderActions,
                orderable: false
            }
        ];
    }

    /**
     * Get column width definitions
     */
    getColumnDefinitions() {
        return [
            { width: "7%", targets: 0 },
            { width: "10%", targets: 1 },
            { width: "12%", targets: 2 },
            { width: "13%", targets: 3 },
            { width: "13%", targets: 4 },
            { width: "9%", targets: 5, className: "text-center" },
            { width: "10%", targets: 6, className: "text-center" },
            { width: "10%", targets: 7, className: "text-center" },
            { width: "16%", targets: 8, className: "text-center", orderable: false }
        ];
    }

    /**
     * Get language configuration for DataTable
     */
    getLanguageConfiguration() {
        return {
            search: "<i class='fas fa-search'></i> Search (genes, aliases, proteins):",
            paginate: {
                first: '<i class="fas fa-angle-double-left"></i>',
                previous: '<i class="fas fa-angle-left"></i>',
                next: '<i class="fas fa-angle-right"></i>',
                last: '<i class="fas fa-angle-double-right"></i>'
            }
        };
    }

    /**
     * Render aliases column
     */
    renderAliases(data, type, row) {
        if (type === 'display') {
            if (!data || data === '') return '<span class="text-muted">-</span>';
            if (data.length > 50) {
                const truncated = data.substring(0, 47) + '...';
                return `<span title="${data}">${truncated}</span>`;
            }
            return data;
        }
        return row.searchableText || data || '';
    }

    /**
     * Render repeat count with badge
     */
    renderRepeatCount(data) {
        return `<span class="badge bg-primary">${data}</span>`;
    }

    /**
     * Render status with appropriate badge
     */
    renderStatus(data) {
        if (!data) return '<span class="badge bg-secondary">Unknown</span>';
        
        let displayText = 'Unknown';
        let badgeClass = 'bg-info';
        
        const lowercaseStatus = data.toLowerCase();
        if (lowercaseStatus === 'reviewed' || 
            (lowercaseStatus.includes('reviewed') && !lowercaseStatus.includes('unreviewed'))) {
            displayText = 'Reviewed';
            badgeClass = 'bg-success';
        } else if (lowercaseStatus.includes('unreviewed')) {
            displayText = 'Unreviewed';
            badgeClass = 'bg-secondary';
        }
        
        return `<span class="badge ${badgeClass}" title="${data}">${displayText}</span>`;
    }

    /**
     * Render eligibility with colored badges
     */
    renderEligibility(data, type, row) {
        if (type === 'display') {
            let badgeClass = 'bg-secondary';
            let text = data || 'N/A';
            
            if (data === '1') {
                badgeClass = 'bg-success';
            } else if (data === '2') {
                badgeClass = 'bg-warning text-dark';
            } else if (data === 'N/A') {
                badgeClass = 'bg-secondary';
            }
            
            return `<span class="badge ${badgeClass}">${text}</span>`;
        }
        return data || 'N/A';
    }

    /**
     * Render actions column
     */
    renderActions(data) {
        return `<a href="detail.html?id=${data}" class="action-button">
            <i class="fas fa-eye"></i> View Details
        </a>`;
    }

    /**
     * Setup custom search functionality
     */
    setupCustomSearch() {
        $.fn.dataTable.ext.search.push((settings, searchData, index, rowData, counter) => {
            if (settings.nTable.id !== 'proteinTable') {
                return true;
            }
            
            const searchTerm = settings.oPreviousSearch.sSearch;
            if (!searchTerm) {
                return true;
            }
            
            const searchableContent = [
                rowData.gene || '',
                rowData.uniprotId || '',
                rowData.repeatType || '',
                rowData.chromosome || '',
                rowData.searchableText || ''
            ].join(' ').toLowerCase();
            
            return searchableContent.includes(searchTerm.toLowerCase());
        });
    }

    /**
     * Get table instance
     */
    getTable() {
        return this.table;
    }

    /**
     * Get table data
     */
    getTableData() {
        return this.tableData;
    }
}
