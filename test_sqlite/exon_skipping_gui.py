import sqlite3
import tkinter as tk
from tkinter import ttk, messagebox, filedialog, simpledialog
import pandas as pd
import os
import json
from datetime import datetime

class ExonSkippingGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Exon Skipping Targets Search Tool")
        self.root.geometry("1100x800")
        self.conn = None
        self.results_data = []
        self.db_path = "test_sqlite/repeats.db"
        
        # Create frames
        self.create_frames()
        
        # Create menu
        self.create_menu()
        
        # Create search criteria widgets
        self.create_search_criteria()
        
        # Create results area
        self.create_results_area()
        
        # Status bar
        self.status_var = tk.StringVar()
        self.status_var.set("Ready")
        self.status_bar = tk.Label(self.root, textvariable=self.status_var, bd=1, relief=tk.SUNKEN, anchor=tk.W)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        
        # Connect to database
        self.connect_to_database()

    def create_frames(self):
        # Main frames
        self.criteria_frame = ttk.LabelFrame(self.root, text="Search Criteria")
        self.criteria_frame.pack(fill=tk.X, padx=10, pady=10)
        
        self.results_frame = ttk.LabelFrame(self.root, text="Results")
        self.results_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Sub-frames for criteria
        self.gene_frame = ttk.Frame(self.criteria_frame)
        self.gene_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.repeat_frame = ttk.Frame(self.criteria_frame)
        self.repeat_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.exon_frame = ttk.Frame(self.criteria_frame)
        self.exon_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Buttons frame
        self.button_frame = ttk.Frame(self.criteria_frame)
        self.button_frame.pack(fill=tk.X, padx=5, pady=10)

    def create_menu(self):
        menubar = tk.Menu(self.root)
        
        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Connect to Database", command=self.connect_dialog)
        file_menu.add_command(label="Export Results", command=self.export_results)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)
        menubar.add_cascade(label="File", menu=file_menu)
        
        # Help menu
        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="About", command=self.show_about)
        menubar.add_cascade(label="Help", menu=help_menu)
        
        self.root.config(menu=menubar)

    def create_search_criteria(self):
        # Gene criteria
        ttk.Label(self.gene_frame, text="Gene Filters:", font=("Arial", 10, "bold")).grid(row=0, column=0, sticky=tk.W)
        
        ttk.Label(self.gene_frame, text="Gene Name:").grid(row=1, column=0, sticky=tk.W)
        self.gene_name_var = tk.StringVar()
        self.gene_name_entry = ttk.Entry(self.gene_frame, textvariable=self.gene_name_var)
        self.gene_name_entry.grid(row=1, column=1, sticky=tk.W, padx=5)
        
        # Repeat criteria
        ttk.Label(self.repeat_frame, text="Repeat Filters:", font=("Arial", 10, "bold")).grid(row=0, column=0, sticky=tk.W)
        
        ttk.Label(self.repeat_frame, text="Min Repeat Count:").grid(row=1, column=0, sticky=tk.W)
        self.min_repeats_var = tk.IntVar(value=5)
        self.min_repeats_spin = ttk.Spinbox(self.repeat_frame, from_=1, to=100, textvariable=self.min_repeats_var, width=5)
        self.min_repeats_spin.grid(row=1, column=1, sticky=tk.W, padx=5)
        
        ttk.Label(self.repeat_frame, text="Repeat Type:").grid(row=1, column=2, sticky=tk.W, padx=(15,0))
        self.repeat_type_var = tk.StringVar()
        self.repeat_type_combo = ttk.Combobox(self.repeat_frame, textvariable=self.repeat_type_var, width=15)
        self.repeat_type_combo.grid(row=1, column=3, sticky=tk.W, padx=5)
        self.repeat_type_var.set("Any")
        
        # Block count criteria
        ttk.Label(self.repeat_frame, text="Block Count:").grid(row=2, column=0, sticky=tk.W)
        self.block_count_var = tk.StringVar(value="block_count = 1")
        self.block_count_combo = ttk.Combobox(self.repeat_frame, textvariable=self.block_count_var, width=20)
        self.block_count_combo.grid(row=2, column=1, sticky=tk.W, padx=5, columnspan=2)
        self.block_count_combo['values'] = [
            "Any", 
            "block_count = 1", 
            "block_count = repeat_length", 
            "Custom..."
        ]
        
        # Exon criteria
        ttk.Label(self.exon_frame, text="Exon Filters:", font=("Arial", 10, "bold")).grid(row=0, column=0, sticky=tk.W)
        
        ttk.Label(self.exon_frame, text="Min Overlap %:").grid(row=1, column=0, sticky=tk.W)
        self.min_overlap_var = tk.IntVar(value=70)
        self.min_overlap_spin = ttk.Spinbox(self.exon_frame, from_=0, to=100, textvariable=self.min_overlap_var, width=5)
        self.min_overlap_spin.grid(row=1, column=1, sticky=tk.W, padx=5)
        
        ttk.Label(self.exon_frame, text="Frame Status:").grid(row=1, column=2, sticky=tk.W, padx=(15,0))
        self.frame_status_var = tk.StringVar(value="in_frame")
        self.frame_status_combo = ttk.Combobox(self.exon_frame, textvariable=self.frame_status_var, width=15)
        self.frame_status_combo.grid(row=1, column=3, sticky=tk.W, padx=5)
        self.frame_status_combo['values'] = ["Any", "in_frame", "out_of_frame"]
        
        # Advanced SQL criteria
        ttk.Label(self.exon_frame, text="Additional WHERE clause (advanced):").grid(row=2, column=0, sticky=tk.W, columnspan=2)
        self.custom_sql_var = tk.StringVar()
        self.custom_sql_entry = ttk.Entry(self.exon_frame, textvariable=self.custom_sql_var, width=50)
        self.custom_sql_entry.grid(row=2, column=1, sticky=tk.W+tk.E, padx=5, columnspan=3)
        
        # Buttons
        ttk.Button(self.button_frame, text="Search", command=self.search).pack(side=tk.LEFT, padx=5)
        ttk.Button(self.button_frame, text="Clear", command=self.clear_criteria).pack(side=tk.LEFT, padx=5)
        ttk.Button(self.button_frame, text="Export Results", command=self.export_results).pack(side=tk.RIGHT, padx=5)

    def create_results_area(self):
        # Create notebook for results and queries
        self.results_notebook = ttk.Notebook(self.results_frame)
        self.results_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Results tab
        self.results_tab = ttk.Frame(self.results_notebook)
        self.results_notebook.add(self.results_tab, text="Results")
        
        # Create treeview for results
        self.tree_frame = ttk.Frame(self.results_tab)
        self.tree_frame.pack(fill=tk.BOTH, expand=True)
        
        self.tree_scroll_y = ttk.Scrollbar(self.tree_frame)
        self.tree_scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.tree_scroll_x = ttk.Scrollbar(self.tree_frame, orient=tk.HORIZONTAL)
        self.tree_scroll_x.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.tree = ttk.Treeview(self.tree_frame, 
                                 columns=("gene", "repeat_type", "count", "exon_id", 
                                          "frame", "overlap", "transcript"),
                                 show="headings",
                                 yscrollcommand=self.tree_scroll_y.set,
                                 xscrollcommand=self.tree_scroll_x.set)
        
        self.tree.column("gene", width=100)
        self.tree.column("repeat_type", width=100)
        self.tree.column("count", width=50)
        self.tree.column("exon_id", width=150)
        self.tree.column("frame", width=70)
        self.tree.column("overlap", width=70)
        self.tree.column("transcript", width=150)
        
        self.tree.heading("gene", text="Gene")
        self.tree.heading("repeat_type", text="Repeat Type")
        self.tree.heading("count", text="Count")
        self.tree.heading("exon_id", text="Exon ID")
        self.tree.heading("frame", text="Frame Status")
        self.tree.heading("overlap", text="Overlap %")
        self.tree.heading("transcript", text="Transcript")
        
        self.tree.pack(fill=tk.BOTH, expand=True)
        self.tree_scroll_y.config(command=self.tree.yview)
        self.tree_scroll_x.config(command=self.tree.xview)
        
        # Bind double click to show details
        self.tree.bind("<Double-1>", self.show_details)
        
        # SQL Query tab
        self.sql_tab = ttk.Frame(self.results_notebook)
        self.results_notebook.add(self.sql_tab, text="SQL Query")
        
        self.sql_text = tk.Text(self.sql_tab, wrap=tk.WORD)
        self.sql_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Details tab for selected exon
        self.details_tab = ttk.Frame(self.results_notebook)
        self.results_notebook.add(self.details_tab, text="Details")
        
        self.details_text = tk.Text(self.details_tab, wrap=tk.WORD)
        self.details_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def connect_to_database(self):
        """Connect to the SQLite database"""
        try:
            if os.path.exists(self.db_path):
                self.conn = sqlite3.connect(self.db_path)
                self.conn.row_factory = sqlite3.Row
                self.status_var.set(f"Connected to {self.db_path}")
                
                # Populate combo boxes with data from database
                self.populate_combo_boxes()
            else:
                self.status_var.set(f"Database file not found: {self.db_path}")
                messagebox.showerror("Error", f"Database file not found: {self.db_path}")
        except sqlite3.Error as e:
            self.status_var.set(f"Database connection error: {e}")
            messagebox.showerror("Error", f"Database connection error: {e}")

    def connect_dialog(self):
        """Show dialog to connect to a different database"""
        db_file = filedialog.askopenfilename(
            title="Select SQLite Database", 
            filetypes=(("SQLite files", "*.db *.sqlite"), ("All files", "*.*"))
        )
        if db_file:
            self.db_path = db_file
            self.connect_to_database()

    def populate_combo_boxes(self):
        """Populate combo boxes with data from database"""
        try:
            cursor = self.conn.cursor()
            
            # Get repeat types
            cursor.execute("SELECT DISTINCT repeat_type FROM repeats ORDER BY repeat_type")
            repeat_types = ["Any"] + [row[0] for row in cursor.fetchall()]
            self.repeat_type_combo['values'] = repeat_types
            
            # Get frame statuses
            cursor.execute("SELECT DISTINCT frame_status FROM exons WHERE frame_status IS NOT NULL ORDER BY frame_status")
            frame_statuses = ["Any"] + [row[0] for row in cursor.fetchall()]
            self.frame_status_combo['values'] = frame_statuses
            
        except sqlite3.Error as e:
            self.status_var.set(f"Error populating combo boxes: {e}")

    def build_query(self):
        """Build the SQL query based on search criteria"""
        # Base query with the gene_repeat_counts CTE
        query_parts = ["""
        -- Step 1: Find genes with multiple repeats of the same type
        WITH gene_repeat_counts AS (
            SELECT g.gene_id, g.gene_name, r.repeat_type, COUNT(DISTINCT r.repeat_id) as repeat_count
            FROM genes g
            JOIN proteins p ON g.gene_id = p.gene_id
            JOIN repeats r ON p.protein_id = r.protein_id
            WHERE 1=1
        """]
        
        params = []
        
        # Add where clauses to the CTE
        where_clauses_cte = []
        
        # Add block count condition (to CTE)
        block_count_selection = self.block_count_var.get()
        if block_count_selection and block_count_selection != "Any":
            if block_count_selection == "block_count = 1":
                where_clauses_cte.append("r.block_count = 1")
            elif block_count_selection == "block_count = repeat_length":
                where_clauses_cte.append("r.block_count = r.repeat_length")
            elif block_count_selection.startswith("Custom"):
                # Get a custom block count value via a simple dialog
                custom_val = simpledialog.askstring("Custom Block Count", "Enter a custom block count condition:", 
                                                   initialvalue="r.block_count > 1")
                if custom_val:
                    where_clauses_cte.append(custom_val)
        
        # Add repeat type filter (to CTE)
        repeat_type = self.repeat_type_var.get()
        if repeat_type and repeat_type != "Any":
            where_clauses_cte.append("r.repeat_type = ?")
            params.append(repeat_type)
        
        # Add gene name filter (to CTE)
        gene_name = self.gene_name_var.get().strip()
        if gene_name:
            where_clauses_cte.append("g.gene_name LIKE ?")
            params.append(f"%{gene_name}%")
        
        # Append where clauses to the CTE
        if where_clauses_cte:
            query_parts.append("AND " + " AND ".join(where_clauses_cte))
        
        # Group by and having clause
        query_parts.append(f"""
            GROUP BY g.gene_id, r.repeat_type
            HAVING COUNT(DISTINCT r.repeat_id) > ?
        )
        """)
        
        min_repeats = self.min_repeats_var.get()
        params.append(min_repeats)
        
        # Main query part
        query_parts.append("""
        -- Step 2: Find exons with significant overlap and in-frame status
        SELECT 
            grc.gene_name,
            grc.repeat_type,
            grc.repeat_count,
            e.ensembl_exon_id,
            e.frame_status,
            r.repeat_id,
            r.position,
            te.overlap_percentage,
            te.overlap_bp,
            te.exon_position_in_transcript,
            te.transcript_id,
            t.transcript_name,
            te.exon_number,
            r.block_count,
            r.repeat_length
        FROM gene_repeat_counts grc
        JOIN genes g ON grc.gene_id = g.gene_id
        JOIN proteins p ON g.gene_id = p.gene_id
        JOIN repeats r ON p.protein_id = r.protein_id AND r.repeat_type = grc.repeat_type
        JOIN repeat_exons re ON r.repeat_id = re.repeat_id
        JOIN exons e ON re.exon_id = e.exon_id
        JOIN transcript_exons te ON e.exon_id = te.exon_id
        JOIN transcripts t ON te.transcript_id = t.transcript_id
        WHERE 1=1
        """)
        
        # Add where clauses to the main query
        where_clauses = []
        
        # Add overlap percentage filter
        min_overlap = self.min_overlap_var.get()
        if min_overlap > 0:
            where_clauses.append("te.overlap_percentage >= ?")
            params.append(min_overlap)
        
        # Add frame status filter
        frame_status = self.frame_status_var.get()
        if frame_status and frame_status != "Any":
            where_clauses.append("e.frame_status = ?")
            params.append(frame_status)
        
        # Add block count condition (to main query as well)
        if block_count_selection and block_count_selection != "Any":
            if block_count_selection == "block_count = 1":
                where_clauses.append("r.block_count = 1")
            elif block_count_selection == "block_count = repeat_length":
                where_clauses.append("r.block_count = r.repeat_length")
        
        # Add custom WHERE clause
        custom_sql = self.custom_sql_var.get().strip()
        if custom_sql:
            where_clauses.append(f"({custom_sql})")
        
        # Append where clauses to main query
        if where_clauses:
            query_parts.append("AND " + " AND ".join(where_clauses))
        
        # Add ORDER BY
        query_parts.append("""
        ORDER BY grc.gene_name, grc.repeat_count DESC, te.overlap_percentage DESC
        """)
        
        # Combine all parts
        full_query = "\n".join(query_parts)
        
        return full_query, params
    
    def search(self):
        """Execute the search query"""
        if not self.conn:
            messagebox.showerror("Error", "Not connected to a database")
            return
        
        # Clear previous results
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.sql_text.delete(1.0, tk.END)
        self.details_text.delete(1.0, tk.END)
        self.results_data = []
        
        # Build the query
        query, params = self.build_query()
        
        # Display the query in the SQL tab
        self.sql_text.insert(tk.END, query)
        self.sql_text.insert(tk.END, "\n\nParameters: " + str(params))
        
        # Execute the query
        try:
            cursor = self.conn.cursor()
            cursor.execute(query, params)
            results = cursor.fetchall()
            
            # Store results for later use
            self.results_data = [dict(row) for row in results]
            
            # Group results by gene
            gene_results = {}
            for row in self.results_data:
                gene_name = row['gene_name']
                
                if gene_name not in gene_results:
                    # Insert main gene row
                    gene_id = self.tree.insert("", tk.END, text=gene_name, values=(
                        gene_name,
                        row['repeat_type'],
                        row['repeat_count'],
                        f"{row['repeat_count']} repeats",
                        "",
                        "",
                        ""
                    ))
                    
                    gene_results[gene_name] = {
                        'tree_id': gene_id,
                        'exons': []
                    }
                
                # Add exon as child of gene
                exon_id = self.tree.insert(gene_results[gene_name]['tree_id'], tk.END, text=row['ensembl_exon_id'], values=(
                    "",
                    "",
                    "",
                    row['ensembl_exon_id'],
                    row['frame_status'],
                    f"{row['overlap_percentage']}%",
                    row['transcript_name']
                ))
                gene_results[gene_name]['exons'].append(exon_id)
            
            self.status_var.set(f"Found {len(gene_results)} genes with suitable exons")
            
            # Switch to results tab
            self.results_notebook.select(0)
            
        except sqlite3.Error as e:
            self.status_var.set(f"Search error: {e}")
            messagebox.showerror("Error", f"Search error: {e}")

    def show_details(self, event):
        """Show details for the selected item"""
        self.details_text.delete(1.0, tk.END)
        
        selected_item = self.tree.focus()
        if not selected_item:
            return
            
        # Get parent to determine if this is a gene or exon row
        parent = self.tree.parent(selected_item)
        
        if not parent:  # This is a gene row
            gene_name = self.tree.item(selected_item)['values'][0]
            
            # Find all results for this gene
            gene_results = [r for r in self.results_data if r['gene_name'] == gene_name]
            if not gene_results:
                return
                
            # Display gene details
            self.details_text.insert(tk.END, f"Gene: {gene_name}\n")
            self.details_text.insert(tk.END, f"Repeat Type: {gene_results[0]['repeat_type']}\n")
            self.details_text.insert(tk.END, f"Repeat Count: {gene_results[0]['repeat_count']}\n")
            self.details_text.insert(tk.END, f"In-frame Exons with high overlap: {len(gene_results)}\n\n")
            
            self.details_text.insert(tk.END, "Exon Summary:\n")
            self.details_text.insert(tk.END, "-----------------\n")
            
            for idx, row in enumerate(gene_results[:10], 1):  # Show first 10 exons
                self.details_text.insert(tk.END, f"{idx}. ID: {row['ensembl_exon_id']}\n")
                transcript_info = f" (in transcript {row['transcript_id']})"
                self.details_text.insert(tk.END, f"   Position: {row['exon_position_in_transcript']}{transcript_info}\n")
                self.details_text.insert(tk.END, f"   Exon Number: {row['exon_number']} in {row['transcript_name']}\n")
                self.details_text.insert(tk.END, f"   Overlap: {row['overlap_percentage']}% ({row['overlap_bp']} bp)\n")
                self.details_text.insert(tk.END, f"   Repeat Position: {row['position']}\n")
                self.details_text.insert(tk.END, f"   Block Count: {row['block_count']}\n\n")
            
            if len(gene_results) > 10:
                self.details_text.insert(tk.END, f"... and {len(gene_results) - 10} more exons\n")
            
        else:  # This is an exon row
            # Find the corresponding exon in the results
            exon_id = self.tree.item(selected_item)['values'][3]
            
            # Find this specific exon in the results
            exon_results = [r for r in self.results_data if r['ensembl_exon_id'] == exon_id]
            if not exon_results:
                return
                
            row = exon_results[0]
            
            # Display detailed exon information
            self.details_text.insert(tk.END, f"Exon ID: {row['ensembl_exon_id']}\n")
            self.details_text.insert(tk.END, f"Gene: {row['gene_name']}\n")
            self.details_text.insert(tk.END, f"Repeat Type: {row['repeat_type']}\n")
            self.details_text.insert(tk.END, f"Frame Status: {row['frame_status']}\n\n")
            
            self.details_text.insert(tk.END, "Transcript Details:\n")
            self.details_text.insert(tk.END, "-----------------\n")
            self.details_text.insert(tk.END, f"Transcript ID: {row['transcript_id']}\n")
            self.details_text.insert(tk.END, f"Transcript Name: {row['transcript_name']}\n")
            self.details_text.insert(tk.END, f"Exon Number: {row['exon_number']}\n")
            self.details_text.insert(tk.END, f"Position in Transcript: {row['exon_position_in_transcript']}\n\n")
            
            self.details_text.insert(tk.END, "Repeat Details:\n")
            self.details_text.insert(tk.END, "-----------------\n")
            self.details_text.insert(tk.END, f"Repeat ID: {row['repeat_id']}\n")
            self.details_text.insert(tk.END, f"Position: {row['position']}\n")
            self.details_text.insert(tk.END, f"Overlap: {row['overlap_percentage']}% ({row['overlap_bp']} bp)\n")
            self.details_text.insert(tk.END, f"Block Count: {row['block_count']}\n")
            self.details_text.insert(tk.END, f"Repeat Length: {row['repeat_length']}\n")
        
        # Switch to details tab
        self.results_notebook.select(2)

    def clear_criteria(self):
        """Clear all search criteria"""
        self.gene_name_var.set("")
        self.min_repeats_var.set(5)
        self.repeat_type_var.set("Any")
        self.block_count_var.set("block_count = 1")
        self.min_overlap_var.set(70)
        self.frame_status_var.set("in_frame")
        self.custom_sql_var.set("")

    def export_results(self):
        """Export results to CSV file"""
        if not self.results_data:
            messagebox.showinfo("Export", "No results to export")
            return
            
        filename = filedialog.asksaveasfilename(
            title="Save Results",
            filetypes=(("CSV files", "*.csv"), ("All files", "*.*")),
            defaultextension=".csv"
        )
        
        if not filename:
            return
            
        try:
            df = pd.DataFrame(self.results_data)
            df.to_csv(filename, index=False)
            self.status_var.set(f"Results exported to {filename}")
        except Exception as e:
            messagebox.showerror("Export Error", str(e))

    def show_about(self):
        """Show about dialog"""
        about_text = """
        Exon Skipping Targets Search Tool
        
        This tool allows you to search for potential exon skipping targets
        based on customizable criteria.
        
        Part of the Tandem Repeat Domain Database project.
        """
        messagebox.showinfo("About", about_text)

if __name__ == "__main__":
    root = tk.Tk()
    app = ExonSkippingGUI(root)
    root.mainloop()