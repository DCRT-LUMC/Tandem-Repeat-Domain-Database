#!/usr/bin/env python3
import json
import os
import sys
import datetime
from pathlib import Path

def get_pdb_filename(uniprot_id):
    """
    Determine the PDB filename for a UniProt ID.
    First checks if a local PDB file exists, otherwise returns a path to alphafold.
    """
    # Check for local PDB file first
    local_path = f"data/AF-{uniprot_id}-F1-model_v4.pdb"
    if os.path.exists(local_path):
        return local_path
    
    # Otherwise return the URL to AlphaFold DB
    return f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"


def extract_protein_info(uniprot_id, repeats_data):
    """
    Extract information about a specific protein from the repeats data.
    Returns a dictionary with protein details and its repeat regions.
    """
    protein_info = {
        "uniprot_id": uniprot_id,
        "gene_name": "",
        "protein_name": "",
        "repeat_type": "",
        "status": "",
        "repeat_regions": []
    }
    
    # Find all repeats for this protein
    protein_repeats = []
    for repeat in repeats_data:
        if repeat.get("uniProtId") == uniprot_id:
            protein_repeats.append(repeat)
            
            # Update protein info from the first matching repeat
            if not protein_info["gene_name"]:
                protein_info["gene_name"] = repeat.get("geneName", "")
                protein_info["repeat_type"] = repeat.get("repeatType", "")
                protein_info["status"] = repeat.get("status", "")
    
    # Extract start and end positions from position field for each repeat
    for repeat in protein_repeats:
        position = repeat.get("position", "")
        start_pos = None
        end_pos = None
        
        if position and "amino acids" in position:
            try:
                pos_part = position.split("amino acids ")[1].split(" on")[0]
                if "-" in pos_part:
                    start_pos, end_pos = map(int, pos_part.split("-"))
            except (IndexError, ValueError):
                pass
        
        if start_pos is not None and end_pos is not None:
            # Create a color based on repeat type and position
            color_options = ["#ff5f5f", "#5fba7d", "#5f87ff", "#ffaf5f", "#bf5fff", "#dc3545", "#fd7e14", "#ffc107", "#20c997", "#0dcaf0"]
            color_index = len(protein_info["repeat_regions"]) % len(color_options)
            
            # Add the repeat region
            protein_info["repeat_regions"].append({
                "start": start_pos,
                "end": end_pos,
                "type": repeat.get("repeatType", ""),
                "color": color_options[color_index],
                "color_hex": color_options[color_index].replace("#", "0x"),
                "position": position,
                "chrom": repeat.get("chrom", ""),
                "chromStart": repeat.get("chromStart", ""),
                "chromEnd": repeat.get("chromEnd", ""),
                "ensembl_exon_info": repeat.get("ensembl_exon_info", {})
            })
    
    return protein_info


def generate_html_page(protein_info):
    """
    Generate an HTML visualization page for the specified protein.
    """
    uniprot_id = protein_info["uniprot_id"]
    gene_name = protein_info["gene_name"]
    repeat_type = protein_info["repeat_type"]
    status = protein_info["status"]
    repeat_regions = protein_info["repeat_regions"]
    
    # Skip if no repeat regions found
    if not repeat_regions:
        print(f"No repeat regions found for {uniprot_id}, skipping...")
        return None
    
    # Sort repeat regions by start position
    repeat_regions.sort(key=lambda x: x["start"])
    
    # Get PDB URL for this protein
    pdb_path = get_pdb_filename(uniprot_id)
    
    # Create JavaScript array of repeats
    repeats_js = []
    for i, repeat in enumerate(repeat_regions):
        repeats_js.append(f'''{{
            start: {repeat["start"]},
            end: {repeat["end"]},
            color: "{repeat["color_hex"]}",
            label: "{repeat_type}{i+1}"
        }}''')
    
    repeats_js_array = ",\n        ".join(repeats_js)
    
    # Create legend HTML
    legend_html = []
    for i, repeat in enumerate(repeat_regions):
        legend_html.append(f'''
    <div class="legend-item">
        <div class="color-box" style="background-color: {repeat["color"]};"></div>
        <span>{repeat_type}{i+1} ({repeat["start"]}-{repeat["end"]})</span>
    </div>''')
    
    legend_html_str = "\n".join(legend_html)
    
    html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Tandem Repeat Viewer - {gene_name} ({uniprot_id})</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css">
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    
    <style>
        body {{
            font-family: Arial, sans-serif;
            max-width: 900px;
            margin: 0 auto;
            padding: 20px;
        }}
        
        .viewer-title {{
            background-color: #f0f8ff;
            padding: 10px 15px;
            border-left: 4px solid #4682b4;
            margin-bottom: 20px;
        }}
        
        .protein-info {{
            display: flex;
            justify-content: space-between;
            background-color: #f8f9fa;
            padding: 10px;
            border-radius: 5px;
            margin-bottom: 15px;
        }}
        
        /* 3D viewer styling */
        #proteinViewer {{
            width: 100%;
            height: 400px;
            position: relative;
            margin: 20px 0;
            border: 1px solid #ddd;
            border-radius: 5px;
        }}
        
        .viewer-container {{
            margin-top: 30px;
        }}
        
        .viewer-controls {{
            display: flex;
            margin-bottom: 10px;
            gap: 10px;
        }}
        
        .view-button {{
            padding: 8px 16px;
            background-color: #f0f0f0;
            cursor: pointer;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }}
        
        .view-button.active {{
            background-color: #4682b4;
            color: white;
        }}
        
        .loading-indicator {{
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            background: rgba(255, 255, 255, 0.8);
            padding: 20px;
            border-radius: 5px;
            text-align: center;
            box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
            z-index: 100;
        }}
        
        .viewer-legend {{
            margin-top: 10px;
            padding: 10px;
            background-color: #f8f9fa;
            border-radius: 5px;
            display: flex;
            flex-wrap: wrap;
            gap: 15px;
        }}
        
        .legend-item {{
            display: flex;
            align-items: center;
        }}
        
        .color-box {{
            width: 20px;
            height: 20px;
            margin-right: 8px;
            border-radius: 3px;
        }}
        
        .gene-info {{
            display: flex;
            flex-wrap: wrap;
            gap: 20px;
            margin-bottom: 20px;
        }}
        
        .gene-info-item {{
            flex: 1;
            min-width: 250px;
        }}
        
        .key-fact {{
            background-color: #e9f7ef;
            padding: 8px 12px;
            border-radius: 4px;
            margin: 10px 0;
            display: inline-block;
        }}
        
        .back-button {{
            margin-bottom: 20px;
        }}
        
        .key-fact {{
            background-color: #e9f7ef;
            padding: 8px 12px;
            border-radius: 4px;
            margin: 10px 0;
            display: inline-block;
        }}
    </style>
</head>
<body>
    <div class="back-button">
        <a href="index.html" class="btn btn-outline-secondary btn-sm">← Back to Protein List</a>
    </div>
    
    <div class="viewer-title">
        <h2>Tandem Repeat Domain Viewer</h2>
        <p>Focused view of {repeat_type} repeat domains in protein {uniprot_id} ({gene_name})</p>
    </div>
    
    <div class="protein-info">
        <div>
            <strong>UniProt:</strong> {uniprot_id}
            <br>
            <strong>Gene:</strong> {gene_name}
        </div>
        <div>
            <strong>Repeat Type:</strong> {repeat_type}
            <br>
            <strong>Status:</strong> {status}
        </div>
    </div>
    
    <!-- Add 3D Viewer Container -->
    <div class="viewer-container">
        <h3>3D Structure Visualization</h3>
        <div class="viewer-controls">
            <button class="view-button active" id="highlightRepeatsBtn">Highlight Repeats</button>
            <button class="view-button" id="standardViewBtn">Standard View</button>
            <button class="view-button" id="zoomRepeatsBtn">Zoom to Repeats</button>
        </div>
        <div id="proteinViewer">
            <div class="loading-indicator">Loading protein structure...</div>
        </div>
        <div class="viewer-legend">{legend_html_str}
            <div class="legend-item">
                <div class="color-box" style="background-color: #cccccc;"></div>
                <span>Non-repeat regions</span>
            </div>
        </div>
    </div>
    
    <div class="summary-container mt-4">
        <h3>Summary of Tandem Repeats</h3>
        <p>Protein {uniprot_id} ({gene_name}) contains {len(repeat_regions)} {repeat_type} repeats in positions:</p>
        <ul>
            {"".join(f'<li>{r["start"]}-{r["end"]} ({r["end"] - r["start"] + 1} amino acids)</li>' for r in repeat_regions)}
        </ul>
        <p class="text-muted">Data source: Tandem Repeat Domain Database</p>
    </div>
    
    <script>
        document.addEventListener('DOMContentLoaded', function() {{
            // Initialize the 3D structure viewer using 3Dmol.js
            let viewer = null;
            const proteinViewerElement = document.getElementById('proteinViewer');
            const loadingIndicator = document.querySelector('.loading-indicator');
            
            // Define repeat regions
            const repeatRegions = [
        {repeats_js_array}
            ];
            
            // Path to PDB file
            const pdbUrl = "{pdb_path}";
            
            // Initialize the 3Dmol viewer and load the model
            function initViewer() {{
                // Create the 3Dmol viewer
                viewer = $3Dmol.createViewer(proteinViewerElement, {{
                    backgroundColor: "white",
                    antialias: true,
                    powerPreference: "high-performance"
                }});
                
                // Load the PDB model
                jQuery.ajax(pdbUrl, {{
                    success: function(data) {{
                        // Hide loading indicator once the model is loaded
                        loadingIndicator.style.display = 'none';
                        
                        // Add the model to the viewer
                        let model = viewer.addModel(data, "pdb");
                        
                        // Apply the repeats highlight by default
                        highlightRepeats();
                        
                        // Zoom directly to the repeat regions instead of the whole protein
                        const firstRepeat = repeatRegions[0].start;
                        const lastRepeat = repeatRegions[repeatRegions.length-1].end;
                        viewer.zoomTo({{resi: firstRepeat + "-" + lastRepeat}});
                        
                        // Rotate to a good angle for viewing repeats
                        viewer.rotate(30, {{x:1}});
                        viewer.rotate(20, {{y:1}});
                        
                        // Enable slabbing and render
                        viewer.enableSlabbing();
                        viewer.render();
                    }},
                    error: function(xhr, status, error) {{
                        loadingIndicator.textContent = "Error loading structure: " + error;
                        console.error("Failed to load PDB file:", error);
                    }}
                }});
            }}
            
            // Highlight the repeat domains specifically
            function highlightRepeats() {{
                if (!viewer) return;
                
                viewer.removeAllSurfaces();
                viewer.removeAllShapes();
                viewer.removeAllLabels();
                
                // Set whole protein to a light gray surface
                viewer.setStyle({{}}, {{
                    cartoon: {{color: '0xcccccc', opacity: 0.5}},
                    surface: {{opacity: 0.6, color: '0xdddddd'}}
                }});
                
                // Highlight each repeat with colored surface
                repeatRegions.forEach(repeat => {{
                    const selection = {{resi: repeat.start + "-" + repeat.end}};
                    
                    viewer.setStyle(selection, {{
                        cartoon: {{color: repeat.color, opacity: 1.0}},
                        surface: {{opacity: 0.8, color: repeat.color}}
                    }});
                    
                    // Add label
                    viewer.addLabel(repeat.label, {{
                        position: {{resi: Math.floor((repeat.start + repeat.end) / 2)}},
                        backgroundColor: 'black',
                        backgroundOpacity: 0.7,
                        fontColor: 'white',
                        fontSize: 12
                    }});
                }});
                
                viewer.zoomTo();
                viewer.render();
            }}
            
            // Show standard view (cartoon representation)
            function showStandardView() {{
                if (!viewer) return;
                
                viewer.removeAllSurfaces();
                viewer.removeAllShapes();
                viewer.removeAllLabels();
                
                // Color by chain
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: 'chainHetatm', thickness: 0.8}}}});
                
                // Highlight repeat regions with subtle effect
                repeatRegions.forEach(repeat => {{
                    viewer.addStyle({{resi: repeat.start + "-" + repeat.end}}, {{
                        cartoon: {{opacity: 1.0, thickness: 1.2}}
                    }});
                }});
                
                viewer.zoomTo();
                viewer.render();
            }}
            
            // Zoom specifically to the repeat regions
            function zoomToRepeats() {{
                if (!viewer) return;
                
                // First selection includes all repeats to zoom
                const firstRepeat = repeatRegions[0].start;
                const lastRepeat = repeatRegions[repeatRegions.length-1].end;
                viewer.zoomTo({{resi: firstRepeat + "-" + lastRepeat}});
                
                // Rotate to a good angle for viewing repeats
                viewer.rotate(30, {{x:1}});
                viewer.rotate(20, {{y:1}});
                
                viewer.render();
            }}
            
            // Add event listeners to buttons
            document.getElementById('highlightRepeatsBtn').addEventListener('click', function() {{
                setActiveButton(this);
                highlightRepeats();
            }});
            
            document.getElementById('standardViewBtn').addEventListener('click', function() {{
                setActiveButton(this);
                showStandardView();
            }});
            
            document.getElementById('zoomRepeatsBtn').addEventListener('click', function() {{
                zoomToRepeats();
            }});
            
            // Helper to set active button class
            function setActiveButton(button) {{
                document.querySelectorAll('.view-button').forEach(btn => {{
                    btn.classList.remove('active');
                }});
                button.classList.add('active');
            }}
            
            // Initialize the viewer when the document is ready
            initViewer();
        }});
    </script>
</body>
</html>
'''
    
    return html_content


def generate_index_page(proteins):
    """
    Generate an index HTML page that lists all available proteins.
    """
    # Group proteins by repeat type
    proteins_by_type = {}
    for protein in proteins:
        repeat_type = protein["repeat_type"]
        if repeat_type not in proteins_by_type:
            proteins_by_type[repeat_type] = []
        proteins_by_type[repeat_type].append(protein)
    
    # Create rows of protein cards for each type
    type_sections = []
    for repeat_type, type_proteins in proteins_by_type.items():
        # Create protein cards for this type
        protein_cards = []
        for protein in type_proteins:
            repeat_count = len(protein["repeat_regions"])
            filename = f"{protein['uniprot_id']}.html"
            
            protein_cards.append(f'''
            <div class="col-md-4 col-sm-6 mb-4">
                <div class="card h-100">
                    <div class="card-header">
                        <strong>{protein["gene_name"]}</strong> ({protein["uniprot_id"]})
                    </div>
                    <div class="card-body">
                        <p class="mb-2"><strong>Repeat Type:</strong> {protein["repeat_type"]}</p>
                        <p class="mb-2"><strong>Repeats:</strong> {repeat_count}</p>
                        <p class="mb-0"><strong>Status:</strong> {protein["status"]}</p>
                    </div>
                    <div class="card-footer bg-transparent">
                        <a href="{filename}" class="btn btn-primary btn-sm">View Structure</a>
                    </div>
                </div>
            </div>
            ''')
        
        # Create a section for this repeat type
        protein_cards_html = "\n".join(protein_cards)
        type_sections.append(f'''
        <div class="mb-5">
            <h2 class="mb-4">{repeat_type} Repeat Proteins</h2>
            <div class="row">
                {protein_cards_html}
            </div>
        </div>
        ''')
    
    type_sections_html = "\n".join(type_sections)
    current_date = datetime.datetime.now().strftime('%Y-%m-%d')
    
    index_html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Tandem Repeat Domain Database</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.3/dist/css/bootstrap.min.css">
    <style>
        .hero-section {{
            background-color: #f8f9fa;
            padding: 40px 0;
            margin-bottom: 40px;
            text-align: center;
        }}
        
        .protein-search {{
            max-width: 500px;
            margin: 0 auto;
        }}
        
        .card-header {{
            background-color: #f0f8ff;
        }}
    </style>
</head>
<body>
    <div class="hero-section">
        <div class="container">
            <h1>Tandem Repeat Domain Database</h1>
            <p class="lead">Interactive 3D visualization of protein tandem repeat domains</p>
            
            <div class="protein-search mt-4">
                <input type="text" class="form-control" id="proteinSearch" placeholder="Search for protein (UniProt ID or Gene Name)">
            </div>
        </div>
    </div>
    
    <div class="container">
        <div class="row mb-4">
            <div class="col">
                <h2>Protein Structure Visualization</h2>
                <p>This database provides 3D structure visualization for proteins with tandem repeat domains. Each protein visualization page includes:</p>
                <ul>
                    <li>3D structure viewer with highlighted repeat domains</li>
                    <li>Summary of repeat regions</li>
                    <li>Options to view the protein in different styles</li>
                </ul>
                <p>Select a protein from the listings below to view its structure.</p>
            </div>
        </div>
        
        {type_sections_html}
    </div>
    
    <footer class="bg-light mt-5 py-3">
        <div class="container text-center">
            <p class="mb-0">Tandem Repeat Domain Database - Generated on {current_date}</p>
        </div>
    </footer>
    
    <script>
    document.addEventListener('DOMContentLoaded', function() {{
        const searchInput = document.getElementById('proteinSearch');
        const cards = document.querySelectorAll('.card');
        
        searchInput.addEventListener('input', function() {{
            const query = this.value.toLowerCase().trim();
            
            cards.forEach(card => {{
                const cardText = card.textContent.toLowerCase();
                if (query === '' || cardText.includes(query)) {{
                    card.closest('.col-md-4').style.display = '';
                }} else {{
                    card.closest('.col-md-4').style.display = 'none';
                }}
            }});
        }});
    }});
    </script>
</body>
</html>
'''
    
    return index_html


def main():
    """
    Main function to generate protein visualization pages.
    """
    # Get path to JSON file with repeat data
    if len(sys.argv) > 1:
        json_file = sys.argv[1]
    else:
        json_file = 'output/1000_test_exons_hg38_repeats.json'
    
    # Create output directory
    output_dir = Path('data/protein_pages')
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read the JSON data
    try:
        with open(json_file, 'r') as f:
            repeats_data = json.load(f)
        print(f"Loaded {len(repeats_data)} items from {json_file}")
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        return
    
    # Extract unique proteins
    uniprot_ids = set()
    for repeat in repeats_data:
        uniprot_id = repeat.get("uniProtId")
        if uniprot_id:
            uniprot_ids.add(uniprot_id)
    
    print(f"Found {len(uniprot_ids)} unique proteins")
    
    # Generate a page for each protein
    processed_proteins = []
    for uniprot_id in uniprot_ids:
        print(f"Processing protein {uniprot_id}...")
        
        # Extract protein info
        protein_info = extract_protein_info(uniprot_id, repeats_data)
        
        # Skip proteins with no repeat regions
        if not protein_info["repeat_regions"]:
            continue
        
        # Generate HTML page
        html_content = generate_html_page(protein_info)
        
        if html_content:
            # Write to file
            output_file = output_dir / f"{uniprot_id}.html"
            with open(output_file, 'w') as f:
                f.write(html_content)
            print(f"Created visualization page for {uniprot_id}")
            
            processed_proteins.append(protein_info)
    
    # Generate index page
    index_html = generate_index_page(processed_proteins)
    with open(output_dir / "index.html", 'w') as f:
        f.write(index_html)
    print(f"Created index page with {len(processed_proteins)} proteins")


if __name__ == "__main__":
    main()