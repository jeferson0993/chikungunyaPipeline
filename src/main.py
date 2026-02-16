import os
import re
import sys
import subprocess
import warnings
import matplotlib

# Configuração CRÍTICA para Docker (Backend sem interface gráfica)
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Configurações Globais
warnings.filterwarnings("ignore")
sns.set_theme(style="whitegrid")

class ChikungunyaPipeline:
    """
    Pipeline Dockerizado para análise evolutiva do Vírus Chikungunya.
    Author: Jeferson F Silva
    Contact: jeferson0993@gmail.com
    """

    def __init__(self, email, work_dir="/app/data"):
        self.email = email
        Entrez.email = email
        # No Docker, work_dir será o volume montado
        self.work_dir = work_dir
        os.makedirs(work_dir, exist_ok=True)
        self.metadata = pd.DataFrame()
        print(f"🚀 [Docker] Pipeline Inicializado. Output: {self.work_dir}")

    def fetch_data_ncbi(self, query="Chikungunya virus[Organism] AND Brazil[Geo Location] AND 2015:2024[Date - Collection]", max_records=50):
        print(f"📥 Buscando dados no NCBI: '{query}'...")
        try:
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_records, idtype="acc")
            record = Entrez.read(handle)
            ids = record["IdList"]
            
            if not ids:
                print("⚠️ Nenhuma sequência encontrada.")
                return None

            handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
            records = list(SeqIO.parse(handle, "genbank"))
            
            raw_data = []
            for r in records:
                meta = {'accession': r.id, 'date': None, 'country': 'Unknown', 'length': len(r.seq)}
                for feature in r.features:
                    if feature.type == 'source':
                        meta['country'] = feature.qualifiers.get('country', ['Unknown'])[0]
                        date_str = feature.qualifiers.get('collection_date', ['Unknown'])[0]
                        try:
                            meta['date'] = int(re.search(r'\d{4}', date_str).group())
                        except:
                            meta['date'] = None

                if meta['length'] > 100 and meta['date'] is not None:
                    meta['seq_object'] = r
                    raw_data.append(meta)

            self.metadata = pd.DataFrame(raw_data)
            if self.metadata.empty: return None
            
            self.metadata.sort_values('date', inplace=True)
            self.fasta_raw = os.path.join(self.work_dir, "raw_sequences.fasta")
            SeqIO.write([x['seq_object'] for x in raw_data], self.fasta_raw, "fasta")
            
            print(f"✅ Download: {len(self.metadata)} sequências salvas.")
            return self.metadata
        except Exception as e:
            print(f"❌ Erro Download: {e}")
            return None

    def align_sequences(self):
        print("⚙️ Executando MAFFT (System Call)...")
        output_file = os.path.join(self.work_dir, "aligned.fasta")
        # MAFFT deve estar instalado no container via apt-get
        cmd = f"mafft --auto --quiet {self.fasta_raw} > {output_file}"
        os.system(cmd)
        self.fasta_aligned = output_file
        self.alignment = AlignIO.read(output_file, "fasta")
        print(f"✅ Alinhado: {self.alignment.get_alignment_length()} pb.")

    def analyze_molecular_evolution(self):
        """Pipeline integrado de tradução e dN/dS simplificado"""
        print("🧬 Analisando Evolução Molecular...")
        
        records = []
        for record in self.alignment:
            seq_str = str(record.seq).upper().replace("-", "")
            remainder = len(seq_str) % 3
            if remainder != 0: seq_str = seq_str[:-remainder]
            records.append(SeqRecord(Seq(seq_str), id=record.id))
            
        if not records: return

        ref_seq = str(records[0].seq)
        results = []
        
        for rec in records[1:]:
            query_seq = str(rec.seq)
            min_len = min(len(ref_seq), len(query_seq))
            min_len = min_len - (min_len % 3)
            
            S, N = 0, 0
            for i in range(0, min_len, 3):
                codon_ref = ref_seq[i:i+3]
                codon_qry = query_seq[i:i+3]
                if codon_ref != codon_qry and 'N' not in codon_ref:
                    if Seq(codon_ref).translate() != Seq(codon_qry).translate():
                        N += 1
                    else:
                        S += 1
            
            ratio = (N/S) if S > 0 else 0
            meta_row = self.metadata[self.metadata['accession'] == rec.id]
            date = meta_row['date'].values[0] if not meta_row.empty else 0
            
            results.append({
                'id': rec.id, 'date': date, 'synonymous': S, 
                'nonsynonymous': N, 'total_mutations': S+N, 'dn_ds': ratio
            })
            
        self.mutation_df = pd.DataFrame(results)
        self.mutation_df.to_csv(os.path.join(self.work_dir, "mutation_data.csv"), index=False)

    def generate_plots(self):
        print("📊 Gerando Figuras...")
        if self.mutation_df.empty: return

        # Plot 1: Regressão Temporal
        plt.figure(figsize=(10, 6))
        sns.regplot(data=self.mutation_df, x='date', y='total_mutations', 
                    line_kws={'color':'red'}, scatter_kws={'alpha':0.6})
        plt.title('Molecular Clock: Mutations vs Time')
        plt.savefig(os.path.join(self.work_dir, "plot_temporal.png"), dpi=300)
        plt.close()

        # Plot 2: Árvore
        try:
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(self.alignment)
            tree = DistanceTreeConstructor().nj(dm)
            tree.root_at_midpoint()
            Phylo.write(tree, os.path.join(self.work_dir, "tree.xml"), "phyloxml")
            
            plt.figure(figsize=(12, 12))
            Phylo.draw(tree, do_show=False)
            plt.savefig(os.path.join(self.work_dir, "plot_tree.png"), dpi=300)
            plt.close()
        except Exception as e:
            print(f"Erro na árvore: {e}")

def run():
    # Pega email do env var ou usa default
    email = os.environ.get('NCBI_EMAIL', 'jeferson0993@gmail.com')
    pipeline = ChikungunyaPipeline(email)
    
    # Executa passos
    if pipeline.fetch_data_ncbi(max_records=30):
        pipeline.align_sequences()
        pipeline.analyze_molecular_evolution()
        pipeline.generate_plots()
        print(f"\n✅ Concluído! Verifique a pasta './data' no seu computador.")

if __name__ == "__main__":
    run()
