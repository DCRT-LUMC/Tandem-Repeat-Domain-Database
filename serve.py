#!/usr/bin/env python3
"""
Simple HTTP server to serve the Tandem Repeat Domain Database files.
This helps avoid CORS issues when accessing the JSON data.
"""

import http.server
import socketserver
import os
import sys

# Default port
PORT = 8000

def run_server(port):
    """Run a simple HTTP server on the specified port"""
    handler = http.server.SimpleHTTPRequestHandler
    
    # Set directory to project root
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    # Create and configure the server
    with socketserver.TCPServer(("", port), handler) as httpd:
        print(f"\n✅ Server started at http://localhost:{port}")
        print("📂 Serving files from:", os.getcwd())
        print("🧬 To view the Tandem Repeat Domain Database:")
        print(f"   - Main viewer: http://localhost:{port}/output/dynamic_viewer/index.html")
        print(f"   - Detail pages: http://localhost:{port}/output/protein_pages/[UNIPROT_ID].html")
        print("\nPress Ctrl+C to stop the server.\n")
        
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nServer stopped.")

if __name__ == "__main__":
    # Use command line argument for port if provided
    if len(sys.argv) > 1:
        try:
            PORT = int(sys.argv[1])
        except ValueError:
            print(f"Invalid port: {sys.argv[1]}. Using default port {PORT}.")
    
    run_server(PORT)