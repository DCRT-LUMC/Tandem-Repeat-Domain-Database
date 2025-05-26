class TrexomeApp {
    constructor() {
        this.dataManager = new DataManager();
        this.tableManager = new TableManager();
        this.filterManager = null;
    }

    /**
     * Initialize the application
     */
    async init() {
        try {
            console.log("Document ready, initializing application...");
            
            // Show loading indicator
            this.showLoading();
            
            // Load data
            await this.loadData();
            
            // Initialize UI components
            await this.initializeComponents();
            
            // Hide loading and show content
            this.showContent();
            
            console.log("Application initialized successfully");
        } catch (error) {
            console.error('Error initializing application:', error);
            this.showError(error);
        }
    }

    /**
     * Load all required data
     */
    async loadData() {
        console.log("Loading data...");
        
        // Load gene eligibility data
        await this.dataManager.loadGeneEligibilityData();
        
        // Load protein data
        const tableData = await this.dataManager.loadProteinData();
        console.log(`Processed ${tableData.length} unique proteins`);
        
        return tableData;
    }

    /**
     * Initialize UI components
     */
    async initializeComponents() {
        console.log("Initializing UI components...");
        
        // Get processed data
        const tableData = this.dataManager.getUniqueProteinsArray();
        const filterOptions = this.dataManager.getFilterOptions();
        
        // Initialize DataTable
        const table = this.tableManager.initializeTable(tableData);
        
        // Initialize filters
        this.filterManager = new FilterManager(table, tableData);
        this.filterManager.initializeFilters();
        this.filterManager.populateFilterDropdowns(filterOptions);
        
        console.log("UI components initialized");
    }

    /**
     * Show loading indicator
     */
    showLoading() {
        $('#loading').show();
        $('#tableContainer').hide();
    }

    /**
     * Show main content
     */
    showContent() {
        $('#loading').hide();
        $('#tableContainer').show();
    }

    /**
     * Show error message
     */
    showError(error) {
        $('#loading').html(`
            <div class="alert alert-danger p-4" role="alert">
                <h4 class="alert-heading"><i class="fas fa-exclamation-triangle me-2"></i> Error Loading Data</h4>
                <p>${error.message || 'Unable to load protein data. Please check your connection and try again.'}</p>
                <hr>
                <p class="mb-0">If the problem persists, please check the browser console for more details or contact support.</p>
            </div>
        `);
    }
}

// Initialize application when document is ready
$(document).ready(() => {
    const app = new TrexomeApp();
    app.init();
});
