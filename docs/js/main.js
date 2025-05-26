import { ProteinViewer } from './proteinViewer.js';
import { DataProcessor } from './dataProcessor.js';
import { UIController } from './uiController.js';
import { BoundaryController } from './boundaryController.js';

class DetailPageApp {
    constructor() {
        this.proteinViewer = new ProteinViewer();
        this.dataProcessor = new DataProcessor();
        this.uiController = new UIController();
        this.boundaryController = null; // Will be initialized after proteinViewer
        
        this.currentProteinData = null;
        this.currentRepeatRegions = null;
        this.currentExonData = null;
    }

    async initialize() {
        try {
            // Get the UniProt ID from URL
            const urlParams = new URLSearchParams(window.location.search);
            const uniprotId = urlParams.get('id');
            
            if (!uniprotId) {
                this.uiController.showError("No protein ID specified");
                return;
            }
            
            // Define data URL
            const dataUrl = './all_annotated_repeats.json';
            
            // Load and process the protein data
            const result = await this.dataProcessor.fetchProteinData(uniprotId, dataUrl);
            
            this.currentProteinData = result.proteinInfo;
            this.currentRepeatRegions = result.repeatRegions;
            this.currentExonData = result.exonData;
            
            // Display the protein information
            this.uiController.displayProteinInfo(result.proteinInfo, result.repeatRegions);
            
            // Update exon display if available
            if (result.exonData) {
                this.uiController.updateExonDisplay(result.exonData);
            }
            
            // Try to initialize the 3D viewer
            try {
                await this.proteinViewer.initialize(uniprotId, result.repeatRegions);
                
                // Initialize boundary controller after proteinViewer is ready
                this.boundaryController = new BoundaryController(this.proteinViewer);
                
                // 3D model loaded successfully - ensure the viewer section is expanded
                console.log("3D model loaded successfully");
                this.expand3DViewerSection();
                
            } catch (modelError) {
                // No 3D model available - collapse the viewer section and show info message
                console.log("No 3D model available, collapsing viewer section");
                this.collapse3DViewerSection();
                this.showNo3DModelMessage();
            }
            
            // Show the protein details
            this.uiController.showProteinDetails();
            
            // Set up event listeners
            this.setupEventListeners();
            
        } catch (error) {
            console.error('Error initializing application:', error);
            this.uiController.showError("Error loading protein data");
        }
    }

    collapse3DViewerSection() {
        const viewerSection = document.getElementById('viewerContainer');
        const toggleButton = document.querySelector('[onclick="toggleSection(\'viewerContainer\')"]');
        
        if (viewerSection && !viewerSection.classList.contains('collapsed')) {
            viewerSection.classList.add('collapsed');
            
            // Update the toggle button text and icon
            if (toggleButton) {
                toggleButton.innerHTML = '<i class="fas fa-eye me-1"></i> Show';
            }
        }
        
        console.log("3D Structure Visualization section collapsed due to missing model");
    }

    expand3DViewerSection() {
        const viewerSection = document.getElementById('viewerContainer');
        const toggleButton = document.querySelector('[onclick="toggleSection(\'viewerContainer\')"]');
        
        if (viewerSection && viewerSection.classList.contains('collapsed')) {
            viewerSection.classList.remove('collapsed');
            
            // Update the toggle button text and icon
            if (toggleButton) {
                toggleButton.innerHTML = '<i class="fas fa-eye-slash me-1"></i> Hide';
            }
        }
        
        console.log("3D Structure Visualization section expanded");
    }

    showNo3DModelMessage() {
        const viewerContainer = document.getElementById('proteinViewer');
        const loadingIndicator = viewerContainer.querySelector('.loading-indicator');
        
        if (loadingIndicator) {
            loadingIndicator.innerHTML = `
                <div class="alert alert-warning mb-0">
                    <i class="fas fa-exclamation-triangle me-2"></i>
                    <strong>No 3D Structure Available</strong><br>
                    <small>This protein does not have an available AlphaFold model for 3D visualization.</small>
                </div>
            `;
        }
    }

    hide3DViewerSection() {
        const viewerSection = document.getElementById('viewerContainer');
        const viewerHeader = viewerSection.previousElementSibling; // The section header
        
        if (viewerSection) {
            viewerSection.style.display = 'none';
        }
        
        if (viewerHeader && viewerHeader.classList.contains('section-header')) {
            viewerHeader.style.display = 'none';
        }
        
        console.log("3D Structure Visualization section hidden due to missing model");
    }

    show3DViewerSection() {
        const viewerSection = document.getElementById('viewerContainer');
        const viewerHeader = viewerSection.previousElementSibling; // The section header
        
        if (viewerSection) {
            viewerSection.style.display = 'block';
        }
        
        if (viewerHeader && viewerHeader.classList.contains('section-header')) {
            viewerHeader.style.display = 'flex';
        }
        
        console.log("3D Structure Visualization section shown");
    }

    setupEventListeners() {
        // Only set up 3D viewer event listeners if the viewer is available
        if (this.proteinViewer.viewerInitialized) {
            // 3D Viewer controls
            document.getElementById('highlightRepeatsBtn')?.addEventListener('click', (e) => {
                e.preventDefault();
                this.uiController.setActiveButton(e.target);
                this.proteinViewer.highlightRepeats(this.currentRepeatRegions);
            });
            
            document.getElementById('standardViewBtn')?.addEventListener('click', (e) => {
                e.preventDefault();
                this.uiController.setActiveButton(e.target);
                this.proteinViewer.showStandardView(this.currentRepeatRegions);
            });
            
            document.getElementById('zoomRepeatsBtn')?.addEventListener('click', (e) => {
                e.preventDefault();
                this.uiController.setActiveButton(e.target);
                this.proteinViewer.zoomToRepeats(this.currentRepeatRegions);
            });
            
            document.getElementById('highlightExonsBtn')?.addEventListener('click', (e) => {
                e.preventDefault();
                this.uiController.setActiveButton(e.target);
                this.highlightExons();
            });
            
            document.getElementById('exonBoundariesBtn')?.addEventListener('click', (e) => {
                e.preventDefault();
                e.target.classList.toggle('active');
                if (e.target.classList.contains('active')) {
                    e.target.innerHTML = '<i class="fas fa-dna me-1"></i> Hide Exon Boundaries';
                    this.proteinViewer.showExonBoundaries(this.currentExonData);
                    if (this.boundaryController) {
                        this.boundaryController.updateLegendWithExons(this.currentExonData);
                    }
                } else {
                    e.target.innerHTML = '<i class="fas fa-dna me-1"></i> Show Exon Boundaries';
                    this.proteinViewer.hideExonBoundaries();
                }
            });
            
            document.getElementById('repeatBoundariesBtn')?.addEventListener('click', (e) => {
                e.preventDefault();
                e.target.classList.toggle('active');
                if (e.target.classList.contains('active')) {
                    e.target.innerHTML = '<i class="fas fa-border-style me-1"></i> Hide Repeat Boundaries';
                    if (this.boundaryController) {
                        this.boundaryController.showRepeatBoundaries(this.currentRepeatRegions);
                    }
                } else {
                    e.target.innerHTML = '<i class="fas fa-border-style me-1"></i> Show Repeat Boundaries';
                    if (this.boundaryController) {
                        this.boundaryController.hideRepeatBoundaries();
                    }
                }
            });
        }

        // Section toggle functionality (always available)
        window.toggleSection = (sectionId) => {
            this.uiController.toggleSection(sectionId);
        };
    }

    highlightExons() {
        // Only proceed if viewer is available and initialized
        if (!this.proteinViewer.viewer || !this.proteinViewer.viewerInitialized) {
            console.error("Viewer not initialized yet");
            return;
        }
        
        if (!this.currentExonData?.exons?.length) {
            console.warn("No exon data available for this protein");
            return;
        }
        
        try {
            this.proteinViewer.viewer.removeAllSurfaces();
            this.proteinViewer.viewer.removeAllShapes();
            this.proteinViewer.viewer.removeAllLabels();
            
            this.proteinViewer.viewer.setStyle({}, {
                cartoon: {color: '0xcccccc', opacity: 0.5},
                surface: {opacity: 0.6, color: '0xdddddd'}
            });
            
            const exonColors = [
                "#4daf4a", "#377eb8", "#ff7f00", "#984ea3", "#e41a1c", 
                "#a65628", "#f781bf", "#999999", "#dede00", "#004529",
                "#8c510a", "#542788", "#a6cee3", "#b2df8a", "#6a3d9a"
            ];
            
            const highlightedExons = [];
            
            this.currentExonData.exons.forEach((exon, index) => {
                if (!exon.repeat_start || !exon.repeat_end) return;
                
                const start = exon.repeat_start;
                const end = exon.repeat_end;
                const exonNumber = exon.exon_number;
                
                const colorIndex = index % exonColors.length;
                const color = exonColors[colorIndex];
                const colorHex = color.replace("#", "0x");
                
                const selection = {resi: `${start}-${end}`};
                
                this.proteinViewer.viewer.setStyle(selection, {
                    cartoon: {color: colorHex, opacity: 1.0},
                    surface: {opacity: 0.8, color: colorHex}
                });
                
                const labelText = `Exon ${exonNumber}`;
                this.proteinViewer.viewer.addLabel(labelText, {
                    position: {resi: Math.floor((start + end) / 2)},
                    backgroundColor: 'black',
                    backgroundOpacity: 0.7,
                    fontColor: 'white',
                    fontSize: 12
                });
                
                highlightedExons.push({
                    number: exonNumber,
                    start: start,
                    end: end,
                    color: color
                });
            });
            
            this.updateLegendWithExonColors(highlightedExons);
            this.proteinViewer.viewer.render();
        } catch (error) {
            console.error("Error in highlightExons:", error);
        }
    }

    updateLegendWithExonColors(exons) {
        const legend = document.getElementById('viewerLegend');
        
        while (legend.firstChild) {
            legend.removeChild(legend.firstChild);
        }
        
        const nonRepeatItem = document.createElement('div');
        nonRepeatItem.className = 'legend-item';
        const nonRepeatColorBox = document.createElement('div');
        nonRepeatColorBox.className = 'color-box';
        nonRepeatColorBox.style.backgroundColor = '#cccccc';
        nonRepeatItem.appendChild(nonRepeatColorBox);
        nonRepeatItem.appendChild(document.createTextNode('Non-exon regions'));
        legend.appendChild(nonRepeatItem);
        
        exons.sort((a, b) => a.number - b.number);
        
        exons.forEach(exon => {
            const item = document.createElement('div');
            item.className = 'legend-item';
            
            const colorBox = document.createElement('div');
            colorBox.className = 'color-box';
            colorBox.style.backgroundColor = exon.color;
            
            const label = document.createElement('span');
            label.textContent = `Exon ${exon.number} (${exon.start}-${exon.end})`;
            
            item.appendChild(colorBox);
            item.appendChild(label);
            legend.appendChild(item);
        });
    }
}

// Initialize the application when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    // Hide the summary section by default
    document.getElementById('summaryContent').classList.add('collapsed');
    
    // Initialize the main application
    const app = new DetailPageApp();
    app.initialize();
});
