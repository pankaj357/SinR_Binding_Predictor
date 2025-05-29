#!/usr/bin/env python3
import os
import sys
import tempfile
import multiprocessing
from flask import Flask, render_template, request, redirect, url_for, send_file, flash
from werkzeug.utils import secure_filename
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from functools import partial
import warnings
warnings.filterwarnings('ignore')

# Initialize Flask app
app = Flask(__name__, static_folder='static', static_url_path='/static')
app.secret_key = 'your_secret_key_here'
UPLOAD_FOLDER = tempfile.gettempdir()
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# ==================== EXACT COPY FROM YOUR STANDALONE SCRIPT ====================
class DNAShapePredictor:
    """Predicts DNA shape features"""
    def __init__(self):
        self.shape_params = {
            'AA': {'mgw': 4.0, 'prop_tw': -14.0, 'roll': 0.6},
            'AT': {'mgw': 4.2, 'prop_tw': -13.3, 'roll': 1.1},
            'AC': {'mgw': 4.1, 'prop_tw': -14.5, 'roll': 0.9},
            'AG': {'mgw': 4.0, 'prop_tw': -14.0, 'roll': 1.0},
            'TA': {'mgw': 4.7, 'prop_tw': -11.8, 'roll': 3.3},
            'TT': {'mgw': 4.0, 'prop_tw': -14.0, 'roll': 0.6},
            'TC': {'mgw': 4.1, 'prop_tw': -14.8, 'roll': 0.7},
            'TG': {'mgw': 4.2, 'prop_tw': -14.3, 'roll': 1.2},
            'CA': {'mgw': 4.1, 'prop_tw': -14.8, 'roll': 0.7},
            'CT': {'mgw': 4.1, 'prop_tw': -14.5, 'roll': 0.9},
            'CC': {'mgw': 4.2, 'prop_tw': -15.0, 'roll': 0.8},
            'CG': {'mgw': 4.3, 'prop_tw': -14.8, 'roll': 1.3},
            'GA': {'mgw': 4.0, 'prop_tw': -14.0, 'roll': 1.0},
            'GT': {'mgw': 4.2, 'prop_tw': -14.3, 'roll': 1.2},
            'GC': {'mgw': 4.3, 'prop_tw': -14.8, 'roll': 1.3},
            'GG': {'mgw': 4.2, 'prop_tw': -14.5, 'roll': 1.1}
        }

    def predict_shape(self, sequence):
        try:
            features = {'minor_groove_width': [], 'propeller_twist': [], 'roll': []}
            if len(sequence) < 2:
                return {k: [0] for k in features.keys()}

            for i in range(len(sequence)-1):
                dinuc = sequence[i:i+2].upper()
                if dinuc in self.shape_params:
                    params = self.shape_params[dinuc]
                    for key in features:
                        features[key].append(params[key.split('_')[0] if '_' in key else key])

            if not any(features.values()):
                return {k: [0] for k in features.keys()}
            return features
        except Exception as e:
            print(f"Error in predict_shape: {e}")
            return {k: [0] for k in features.keys()}

class SinRBindingPredictor:
    """Predicts SinR binding sites"""
    def __init__(self, known_motifs=None):
        if known_motifs is None:
            self.known_motifs = [
                "GTTCTCT",
                "AGAAGAC",
                "GTTNNNNNNNNAAC",
                "CACGAAAT",
                "TGAAAT"
            ]
        else:
            self.known_motifs = known_motifs
        self.shape_predictor = DNAShapePredictor()

    def _calculate_palindrome_score(self, sequence):
        try:
            if not sequence:
                return 0
            rev_comp = str(Seq(sequence).reverse_complement())
            matches = sum(a == b for a, b in zip(sequence, rev_comp))
            return matches / len(sequence)
        except Exception:
            return 0

    def _calculate_conservation_score(self, sequence):
        try:
            if not sequence:
                return 0
            scores = []
            for motif in self.known_motifs:
                if len(sequence) < len(motif):
                    continue
                if 'N' in motif:
                    non_n_positions = [i for i, c in enumerate(motif) if c != 'N']
                    motif_bases = [motif[i] for i in non_n_positions]
                    seq_bases = [sequence[i] for i in non_n_positions if i < len(sequence)]
                    if seq_bases:
                        score = sum(a == b for a, b in zip(motif_bases, seq_bases)) / len(seq_bases)
                        scores.append(score)
                else:
                    matches = sum(a == b for a, b in zip(sequence, motif))
                    scores.append(matches / len(motif))
            return max(scores) if scores else 0
        except Exception:
            return 0

    def calculate_shape_score(self, sequence):
        try:
            shape_features = self.shape_predictor.predict_shape(sequence)
            if not any(shape_features.values()):
                return 0

            avg_mgw = np.mean(shape_features['minor_groove_width'])
            avg_ptw = np.mean(shape_features['propeller_twist'])
            avg_roll = np.mean(shape_features['roll'])

            shape_score = (
                0.4 * (avg_mgw / 5.0) +
                0.3 * (abs(avg_ptw) / 15.0) +
                0.3 * (avg_roll / 3.5)
            )
            return shape_score
        except Exception:
            return 0

    def predict_binding_affinity(self, sequence):
        try:
            if not sequence or len(sequence) < 6:
                return 0

            at_content = (sequence.count('A') + sequence.count('T')) / len(sequence)
            palindrome_score = self._calculate_palindrome_score(sequence)
            conservation_score = self._calculate_conservation_score(sequence)
            shape_score = self.calculate_shape_score(sequence)

            score = (
                0.25 * at_content +
                0.25 * palindrome_score +
                0.25 * conservation_score +
                0.25 * shape_score
            )
            return score
        except Exception:
            return 0

def process_sequence_chunk(chunk, known_motifs):
    """Process a single sequence chunk - identical to standalone"""
    try:
        predictor = SinRBindingPredictor(known_motifs=known_motifs)
        sequence = chunk['sequence']
        window_size = 20
        step_size = 10
        results = []
        seen = set()

        for i in range(0, len(sequence) - window_size + 1, step_size):
            window = sequence[i:i + window_size]
            position = chunk['start'] + i
            key = (position, window)

            if key in seen:
                continue
            seen.add(key)

            score = predictor.predict_binding_affinity(window)
            results.append({
                'position': position,
                'sequence': window,
                'score': score,
                'contig': chunk['contig']
            })

        return results
    except Exception as e:
        print(f"Error processing chunk: {e}")
        return []
# ==================== END OF STANDALONE SCRIPT LOGIC ====================

def run_parallel_processing(chunks, known_motifs):
    """Parallel processing with error handling"""
    try:
        with multiprocessing.Pool(processes=min(4, multiprocessing.cpu_count())) as pool:
            results = []
            for result in pool.imap_unordered(partial(process_sequence_chunk, known_motifs=known_motifs), chunks):
                results.extend(result)
            return results
    except Exception as e:
        print(f"Parallel processing failed: {str(e)}")
        return []

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        # Validate file uploads
        if 'protein_fasta' not in request.files or 'genome_fasta' not in request.files or 'known_motifs' not in request.form:
            flash("Missing required inputs")
            return redirect(request.url)
            
        protein_file = request.files['protein_fasta']
        genome_file = request.files['genome_fasta']
        motifs_text = request.form.get('known_motifs', '').strip()
        
        if protein_file.filename == '' or genome_file.filename == '' or not motifs_text:
            flash("No selected files or empty motifs")
            return redirect(request.url)

        # Save files with verification
        os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
        protein_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(protein_file.filename))
        genome_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(genome_file.filename))
        
        protein_file.save(protein_path)
        genome_file.save(genome_path)

        # Read motifs from textarea input, splitting by lines, stripping each motif
        known_motifs = [line.strip() for line in motifs_text.splitlines() if line.strip()]

        # Process protein file
        try:
            protein_record = next(SeqIO.parse(protein_path, "fasta"))
            protein_name = protein_record.id
        except Exception:
            flash("Invalid protein FASTA file")
            return redirect(request.url)

        # Create chunks with identical parameters to standalone
        chunks = []
        chunk_size = 1000
        step_size = 500
        
        for record in SeqIO.parse(genome_path, "fasta"):
            seq = str(record.seq).upper()  # Consistent uppercase conversion
            for i in range(0, len(seq) - chunk_size + 1, step_size):
                chunks.append({
                    'contig': record.id,
                    'start': i,
                    'sequence': seq[i:i + chunk_size]
                })

        # Process chunks with identical logic
        results = run_parallel_processing(chunks, known_motifs)
        
        # Process results identically to standalone
        df = pd.DataFrame(results)
        df.drop_duplicates(subset=['position', 'sequence', 'contig'], inplace=True)
        df = df.sort_values('score', ascending=False)

        # Generate output files
        csv_path = os.path.join(app.config['UPLOAD_FOLDER'], "results.csv")
        df.to_csv(csv_path, index=False)

        # Visualization
        plt.figure(figsize=(10, 6))
        sns.histplot(df["score"], bins=50)
        plt.title('Binding Score Distribution')
        score_plot_path = os.path.join(app.config['UPLOAD_FOLDER'], "score_distribution.png")
        plt.savefig(score_plot_path)
        plt.close()

        plt.figure(figsize=(15, 6))
        sns.scatterplot(data=df, x='position', y='score')
        plt.title(f'{protein_name} Binding Sites Across Genome')
        position_plot_path = os.path.join(app.config['UPLOAD_FOLDER'], "genome_position.png")
        plt.savefig(position_plot_path)
        plt.close()

        # Generate full HTML report with top 10 and plots
        html_report_path = os.path.join(app.config['UPLOAD_FOLDER'], "results_report.html")
        with open(html_report_path, "w") as f:
            f.write(f"""
            <!DOCTYPE html>
            <html lang="en">
            <head>
                <meta charset="UTF-8" />
                <meta name="viewport" content="width=device-width, initial-scale=1" />
                <title>SinR Binding Site Full Report for {protein_name}</title>
                <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet" />
                <style>
                    body {{
                        font-family: Arial, sans-serif;
                        padding: 20px;
                        background-color: #f9f9f9;
                    }}
                    h1, h2 {{
                        color: #2a3f54;
                    }}
                    table {{
                        border-collapse: collapse;
                        width: 100%;
                        margin-bottom: 30px;
                        box-shadow: 0 0 10px rgba(0,0,0,0.1);
                    }}
                    th, td {{
                        border: 1px solid #ddd;
                        padding: 8px;
                        text-align: left;
                    }}
                    th {{
                        background-color: #2980b9;
                        color: white;
                    }}
                    tbody tr:nth-child(even) {{
                        background-color: #f2f6fc;
                    }}
                    .plot-img {{
                        max-width: 100%;
                        height: auto;
                        margin-bottom: 40px;
                        border-radius: 8px;
                        box-shadow: 0 6px 15px rgba(0,0,0,0.1);
                    }}
                </style>
            </head>
            <body>
                <h1>SinR Binding Site Full Report</h1>
                <h2>Top 10 Predicted Binding Sites for {protein_name}</h2>
                {df.head(10).to_html(index=False, classes='table table-striped table-bordered')}
                <h2>Score Distribution Plot</h2>
                <img src="score_distribution.png" alt="Score Distribution" class="plot-img" />
                <h2>Genome Position Plot</h2>
                <img src="genome_position.png" alt="Genome Position Plot" class="plot-img" />
            </body>
            </html>
            """)

        return render_template("results.html",
                            data={'protein_name': protein_name},
                            results=df.head(10).to_dict(orient='records'),
                            csv_path=url_for('download_csv'),
                            score_plot_path=url_for('download_plot'),
                            position_plot_path=url_for('download_position_plot'),
                            report_path=url_for('download_report'))

    return render_template("index.html")

@app.route('/download/csv')
def download_csv():
    path = os.path.join(app.config['UPLOAD_FOLDER'], "results.csv")
    return send_file(path, as_attachment=True, download_name="results.csv")

@app.route('/download/plot')
def download_plot():
    path = os.path.join(app.config['UPLOAD_FOLDER'], "score_distribution.png")
    return send_file(path, mimetype='image/png', as_attachment=True, download_name="score_distribution.png")

@app.route('/download/position_plot')
def download_position_plot():
    path = os.path.join(app.config['UPLOAD_FOLDER'], "genome_position.png")
    return send_file(path, mimetype='image/png', as_attachment=True, download_name="genome_position.png")

@app.route('/download/report')
def download_report():
    path = os.path.join(app.config['UPLOAD_FOLDER'], "results_report.html")
    return send_file(path, mimetype='text/html', as_attachment=True, download_name="results_report.html")

if __name__ == '__main__':
    os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
    multiprocessing.set_start_method('spawn')
    app.run(debug=True, port=5050)

