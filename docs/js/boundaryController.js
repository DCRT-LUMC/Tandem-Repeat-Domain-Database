export class BoundaryController {
    constructor(viewer) {
        this.viewer = viewer;
    }

    showRepeatBoundaries(repeatRegions) {
        if (!this.viewer?.viewer) {
            console.error("Viewer not initialized");
            return;
        }
        
        console.log("Showing repeat boundaries");
        
        if (!repeatRegions || repeatRegions.length === 0) {
            console.warn("No repeat regions available");
            return;
        }
        
        // Clear existing shapes and labels
        this.viewer.viewer.removeAllShapes();
        this.viewer.viewer.removeAllLabels();
        
        // Create a map to track all positions and their boundary types
        const positionMap = new Map();
        
        // First pass: Collect all boundary positions and identify overlaps
        repeatRegions.forEach(repeat => {
            if (!repeat.start || !repeat.end) return;
            
            // Add start position
            if (!positionMap.has(repeat.start)) {
                positionMap.set(repeat.start, { labels: [] });
            }
            positionMap.get(repeat.start).labels.push(`${repeat.label} start`);
            
            // Add end position
            if (!positionMap.has(repeat.end)) {
                positionMap.set(repeat.end, { labels: [] });
            }
            positionMap.get(repeat.end).labels.push(`${repeat.label} end`);
        });
        
        // Second pass: Add spheres and labels for each position
        for (const [position, info] of positionMap.entries()) {
            // Get atoms at this position
            const atoms = this.viewer.viewer.getModel().selectedAtoms({resi: position, atom: "CA"});
            
            if (atoms.length === 0) continue;
            
            // Get the first CA atom coordinates
            const atom = atoms[0];
            
            // Use different colors for spheres based on whether it's a start, end, or both
            let sphereColor = '0x3498DB'; // Default to blue (start)
            
            const hasStart = info.labels.some(l => l.includes('start'));
            const hasEnd = info.labels.some(l => l.includes('end'));
            
            if (hasStart && hasEnd) {
                // If position has both start and end, use orange
                sphereColor = '0xF39C12';
            } else if (hasEnd) {
                // If only end, use yellow
                sphereColor = '0xF1C40F';
            }
            
            // Add a wireframe sphere at this position
            this.viewer.viewer.addSphere({
                center: {x: atom.x, y: atom.y, z: atom.z},
                radius: 1.5,
                color: sphereColor,
                wireframe: true,
                linewidth: 3,
                alpha: 0.9
            });
            
            // Add labels with slight vertical offsets to avoid overlap
            info.labels.forEach((label, idx) => {
                // Offset each label vertically to prevent overlap
                const offset = (idx - (info.labels.length - 1) / 2) * 2;
                
                // Determine label color (blue for start, yellow for end)
                const bgColor = label.includes('start') ? '0x2980B9' : '0xF39C12';
                
                this.viewer.viewer.addLabel(label, {
                    position: {x: atom.x, y: atom.y + offset, z: atom.z},
                    backgroundColor: bgColor,
                    backgroundOpacity: 0.8,
                    fontColor: 'white',
                    fontSize: 12
                });
            });
        }
        
        // Update legend with repeat boundary info
        this.updateLegendWithRepeatBoundaries();
        
        this.viewer.viewer.render();
        console.log("Repeat boundaries display complete");
    }

    hideRepeatBoundaries() {
        if (!this.viewer?.viewer) return;
        this.viewer.viewer.removeAllShapes();
        this.viewer.viewer.removeAllLabels();
        this.viewer.viewer.render();
        console.log("Repeat boundaries hidden");
    }

    updateLegendWithRepeatBoundaries() {
        const legend = document.getElementById('viewerLegend');
        
        // Remove existing boundary legend items
        const existingBoundaryItems = document.querySelectorAll('.boundary-legend-item');
        existingBoundaryItems.forEach(item => item.remove());
        
        // Add new header for repeat boundary legend section
        const boundaryHeader = document.createElement('div');
        boundaryHeader.className = 'legend-item boundary-legend-item';
        boundaryHeader.innerHTML = '<strong>Repeat Boundaries:</strong>';
        legend.appendChild(boundaryHeader);
        
        // Add start marker
        const startItem = document.createElement('div');
        startItem.className = 'legend-item boundary-legend-item';
        const startBox = document.createElement('div');
        startBox.className = 'color-box';
        startBox.style.backgroundColor = '#3498DB'; // Blue
        startItem.appendChild(startBox);
        startItem.appendChild(document.createTextNode('Repeat start'));
        legend.appendChild(startItem);
        
        // Add end marker
        const endItem = document.createElement('div');
        endItem.className = 'legend-item boundary-legend-item';
        const endBox = document.createElement('div');
        endBox.className = 'color-box';
        endBox.style.backgroundColor = '#F1C40F'; // Yellow
        endItem.appendChild(endBox);
        endItem.appendChild(document.createTextNode('Repeat end'));
        legend.appendChild(endItem);
        
        // Add combined marker
        const bothItem = document.createElement('div');
        bothItem.className = 'legend-item boundary-legend-item';
        const bothBox = document.createElement('div');
        bothBox.className = 'color-box';
        bothBox.style.backgroundColor = '#F39C12'; // Orange
        bothItem.appendChild(bothBox);
        bothItem.appendChild(document.createTextNode('Repeat start & end (shared position)'));
        legend.appendChild(bothItem);
    }

    updateLegendWithExons(exonData) {
        const legend = document.getElementById('viewerLegend');
        
        // Remove existing boundary legend items
        const existingBoundaryItems = document.querySelectorAll('.boundary-legend-item, .exon-boundary-legend-item');
        existingBoundaryItems.forEach(item => item.remove());
        
        // Add exon boundary legend
        const exonHeader = document.createElement('div');
        exonHeader.className = 'legend-item exon-boundary-legend-item';
        exonHeader.innerHTML = '<strong>Exon Boundaries:</strong>';
        legend.appendChild(exonHeader);
        
        // Add start marker (green)
        const startItem = document.createElement('div');
        startItem.className = 'legend-item exon-boundary-legend-item';
        const startBox = document.createElement('div');
        startBox.className = 'color-box';
        startBox.style.backgroundColor = '#00FF00'; // Green
        startItem.appendChild(startBox);
        startItem.appendChild(document.createTextNode('Exon start'));
        legend.appendChild(startItem);
        
        // Add end marker (red)
        const endItem = document.createElement('div');
        endItem.className = 'legend-item exon-boundary-legend-item';
        const endBox = document.createElement('div');
        endBox.className = 'color-box';
        endBox.style.backgroundColor = '#FF0000'; // Red
        endItem.appendChild(endBox);
        endItem.appendChild(document.createTextNode('Exon end'));
        legend.appendChild(endItem);
        
        // Add combined marker (purple)
        const bothItem = document.createElement('div');
        bothItem.className = 'legend-item exon-boundary-legend-item';
        const bothBox = document.createElement('div');
        bothBox.className = 'color-box';
        bothBox.style.backgroundColor = '#FF00FF'; // Purple
        bothItem.appendChild(bothBox);
        bothItem.appendChild(document.createTextNode('Exon start & end (shared position)'));
        legend.appendChild(bothItem);
    }
}
