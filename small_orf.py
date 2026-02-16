import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import re

# Standard genetic code stop codons
STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"

class SmallORFFinder(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Small ORF Finder - Genome Analyzer")
        self.geometry("1100x700")

        self.orfs = []  # List of found ORFs: dicts with details
        self.input_seqs = {}  # Loaded sequences for reference

        self.create_widgets()

    def create_widgets(self):
        # Top controls
        top_frame = ttk.Frame(self, padding=10)
        top_frame.pack(fill=tk.X)

        ttk.Button(top_frame, text="Load FASTA/FNA File", 
                  command=self.load_file).pack(side=tk.LEFT, padx=5)

        ttk.Label(top_frame, text="Target AA length:").pack(side=tk.LEFT, padx=(20, 0))
        self.target_aa_var = tk.IntVar(value=50)
        ttk.Entry(top_frame, textvariable=self.target_aa_var, width=8).pack(side=tk.LEFT, padx=5)

        ttk.Label(top_frame, text="Tolerance % (±):").pack(side=tk.LEFT, padx=(10, 0))
        self.tol_var = tk.DoubleVar(value=10.0)
        ttk.Entry(top_frame, textvariable=self.tol_var, width=8).pack(side=tk.LEFT, padx=5)

        ttk.Button(top_frame, text="Find ORFs", 
                  command=self.find_orfs).pack(side=tk.LEFT, padx=10)

        ttk.Button(top_frame, text="Save ORFs (FASTA)", 
                  command=self.save_orfs).pack(side=tk.RIGHT)

        # Results table
        columns = ("seq_id", "frame", "strand", "start", "end", "aa_len", "nt_seq", "aa_seq")
        self.tree = ttk.Treeview(self, columns=columns, show="headings", height=20)
        
        for col in columns:
            self.tree.heading(col, text=col.replace("_", " ").title())
            self.tree.column(col, width=100 if col in ("start", "end") else 120)

        self.tree.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        # Status
        self.status_var = tk.StringVar(value="Load a genome file and set target AA length")
        ttk.Label(self, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W).pack(fill=tk.X)

    def load_file(self):
        filename = filedialog.askopenfilename(
            title="Select FASTA/FNA file",
            filetypes=[("FASTA/FNA files", "*.fasta *.fna *.fa *.gz"), ("All files", "*.*")]
        )
        if not filename:
            return

        self.input_seqs.clear()
        self.orfs.clear()
        for item in self.tree.get_children():
            self.tree.delete(item)

        try:
            if filename.endswith('.gz'):
                handle = gzip.open(filename, "rt")
            else:
                handle = open(filename, "r")

            for record in SeqIO.parse(handle, "fasta"):
                self.input_seqs[record.id] = str(record.seq).upper()

            handle.close()
            self.status_var.set(f"Loaded {len(self.input_seqs)} sequences from {filename}")
        except Exception as e:
            messagebox.showerror("Load error", f"Failed to load {filename}:\n{e}")

    def find_orfs(self):
        if not self.input_seqs:
            messagebox.showwarning("No file", "Please load a FASTA/FNA file first.")
            return

        target_aa = self.target_aa_var.get()
        tolerance = self.tol_var.get() / 100.0
        min_aa = int(target_aa * (1 - tolerance))
        max_aa = int(target_aa * (1 + tolerance))

        self.orfs.clear()
        for item in self.tree.get_children():
            self.tree.delete(item)

        self.status_var.set("Scanning for ORFs...")
        self.update_idletasks()

        total_orfs = 0
        for seq_id, seq in self.input_seqs.items():
            total_orfs += self._scan_sequence(seq_id, seq, min_aa, max_aa)

        self.status_var.set(f"Found {len(self.orfs)} small ORFs (target {target_aa}±{tolerance*100:.0f}% AA)")

        if self.orfs:
            self._populate_table()
        else:
            messagebox.showinfo("No matches", f"No ORFs found in range {min_aa}-{max_aa} amino acids.")

    def _scan_sequence(self, seq_id, seq, min_aa, max_aa):
        """Scan one sequence in all 6 frames for ORFs matching AA range. [web:54][web:55][web:57]"""
        found = 0
        seq_len = len(seq)

        for strand, nuc in [(+1, seq), (-1, seq[::-1])]:  # Forward and reverse complement
            for frame in range(3):
                frame_seq = nuc[frame:]
                # Length must be multiple of 3
                for i in range(0, len(frame_seq) - 2, 3):
                    codon = frame_seq[i:i+3]
                    if codon == START_CODON:
                        # Found potential start; scan for stop
                        for j in range(i + 3, len(frame_seq) - 2, 3):
                            stop_codon = frame_seq[j:j+3]
                            if stop_codon in STOP_CODONS:
                                nt_orf = frame_seq[i:j+3]
                                aa_len = (j + 3 - i) // 3
                                if min_aa <= aa_len <= max_aa:
                                    aa_orf = str(Seq(nt_orf).translate(to_stop=True))
                                    self.orfs.append({
                                        "seq_id": seq_id,
                                        "frame": frame + 1,
                                        "strand": "+" if strand == +1 else "-",
                                        "start": i + frame if strand == +1 else seq_len - (i + frame + len(nt_orf)),
                                        "end": i + frame + len(nt_orf) if strand == +1 else seq_len - (i + frame),
                                        "aa_len": aa_len,
                                        "nt_seq": nt_orf,
                                        "aa_seq": aa_orf
                                    })
                                    found += 1
                                break  # Stop after first in-frame stop
        return found

    def _populate_table(self):
        """Fill the results table."""
        for orf in self.orfs:
            self.tree.insert("", tk.END, values=(
                orf["seq_id"],
                f"{orf['strand']}{orf['frame']}",
                orf["start"],
                orf["end"],
                orf["aa_len"],
                orf["nt_seq"][:50] + "..." if len(orf["nt_seq"]) > 50 else orf["nt_seq"],
                orf["aa_seq"][:30] + "..." if len(orf["aa_seq"]) > 30 else orf["aa_seq"]
            ))

    def save_orfs(self):
        if not self.orfs:
            messagebox.showwarning("No data", "No ORFs to save.")
            return

        filename = filedialog.asksaveasfilename(
            title="Save ORFs as FASTA",
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")]
        )
        if not filename:
            return

        try:
            with open(filename, "w") as f:
                for i, orf in enumerate(self.orfs, 1):
                    f.write(f">{orf['seq_id']}_ORF{i}_frame{orf['strand']}{orf['frame']}_len{orf['aa_len']}aa_pos{orf['start']}-{orf['end']}\n")
                    f.write(f"{orf['nt_seq']}\n")
            self.status_var.set(f"Saved {len(self.orfs)} ORFs to {filename}")
        except Exception as e:
            messagebox.showerror("Save error", f"Failed to save:\n{e}")

if __name__ == "__main__":
    app = SmallORFFinder()
    app.mainloop()
