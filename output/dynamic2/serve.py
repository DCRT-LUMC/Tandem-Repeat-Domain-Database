#!/usr/bin/env python3
"""
Simple HTTP server to serve the Tandem Repeat Domain Database files.
This helps avoid CORS issues when accessing the JSON data.
"""

import http.server
import socketserver
import os
import sys
import json

# Default port
PORT = 8020

class CustomHandler(http.server.SimpleHTTPRequestHandler):
    """Custom HTTP request handler with JSON file info endpoint"""
    
    def do_GET(self):
        """Handle GET requests"""
        # Add an endpoint to check JSON file info
        if self.path == '/fileinfo':
            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            
            # Get the path to the hierarchical JSON file
            current_dir = os.path.dirname(os.path.abspath(__file__))
            json_path = os.path.join(current_dir, 'all_hierarchical.json')
            
            # Collect file info
            file_info = {
                'exists': os.path.exists(json_path),
                'size': os.path.getsize(json_path) if os.path.exists(json_path) else 0,
                'modified': os.path.getmtime(json_path) if os.path.exists(json_path) else 0,
            }
            
            # Try to read the first few bytes to check format
            if file_info['exists'] and file_info['size'] > 0:
                try:
                    with open(json_path, 'r') as f:
                        # Read first 1000 characters
                        file_start = f.read(1000)
                        file_info['preview'] = file_start
                        
                        # Try to parse the entire file to verify it's valid JSON
                        f.seek(0)  # Go back to start of file
                        data = json.load(f)
                        file_info['valid_json'] = True
                        file_info['top_level_keys'] = list(data.keys())
                        
                        # Count chromosomes, genes, etc.
                        chromosomes = [k for k in data.keys() if k.startswith('chr')]
                        file_info['chromosome_count'] = len(chromosomes)
                        
                        # Count genes for first chromosome
                        if chromosomes:
                            first_chrom = chromosomes[0]
                            gene_count = len(data[first_chrom])
                            file_info['first_chromosome'] = first_chrom
                            file_info['gene_count_in_first_chrom'] = gene_count
                        
                except Exception as e:
                    file_info['valid_json'] = False
                    file_info['error'] = str(e)
            
            self.wfile.write(json.dumps(file_info, indent=2).encode())
            return
            
        # For all other requests, use the default handler
        return http.server.SimpleHTTPRequestHandler.do_GET(self)

def run_server(port):
    """Run a simple HTTP server on the specified port"""
    
    # Set directory to script directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(current_dir)
    
    # Create and configure the server with our custom handler
    with socketserver.TCPServer(("", port), CustomHandler) as httpd:
        print(f"\n✅ Server started at http://localhost:{port}")
        print("📂 Serving files from:", os.getcwd())
        print("🧬 To view the Tandem Repeat Domain Database:")
        print(f"   - Main viewer: http://localhost:{port}/index.html")
        print(f"   - Detail pages: http://localhost:{port}/detail.html?id=[UNIPROT_ID]")
        print(f"   - JSON file info: http://localhost:{port}/fileinfo")
        print("\nPress Ctrl+C to stop the server.\n")
        
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nServer stopped.")

if __name__ == "__main__":
    # Get port from command line args or use default
    port = PORT
    if len(sys.argv) > 1:
        try:
            port = int(sys.argv[1])
        except ValueError:
            print(f"Invalid port: {sys.argv[1]}")
            print(f"Using default port {PORT}")
            port = PORT
    
    run_server(port)