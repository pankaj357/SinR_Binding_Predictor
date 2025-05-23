import os
import sys
import tempfile
from flask import Flask, render_template, request, redirect, url_for, send_file, flash
from werkzeug.utils import secure_filename
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Seq import Seq
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# Initialize Flask app with static folder configuration
app = Flask(__name__, static_folder='static', static_url_path='/static')
app.secret_key = 'your_secret_key_here'
UPLOAD_FOLDER = tempfile.gettempdir()
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

class DNAShapePredictor:
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
    def __init__(self):
        self.known_motifs = [
            "GTTCTCT", "AGAAGAC", "GTTNNNNNNNNAAC", "CACGAAAT", "TGAAAT"
        ]
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

def process_sequence_chunk(chunk):
    try:
        predictor = SinRBindingPredictor()
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

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        protein_file = request.files.get('protein_fasta')
        genome_file = request.files.get('genome_fasta')
        known_motifs_text = request.form.get('known_motifs')

        if not protein_file or not genome_file or not known_motifs_text:
            flash("Please upload all required files and provide known motifs.")
            return redirect(request.url)

        protein_filename = secure_filename(protein_file.filename)
        genome_filename = secure_filename(genome_file.filename)
        protein_path = os.path.join(app.config['UPLOAD_FOLDER'], protein_filename)
        genome_path = os.path.join(app.config['UPLOAD_FOLDER'], genome_filename)
        protein_file.save(protein_path)
        genome_file.save(genome_path)

        try:
            protein_record = next(SeqIO.parse(protein_path, "fasta"))
        except Exception:
            flash("Invalid protein FASTA file.")
            return redirect(request.url)

        chunks = []
        chunk_size = 1000
        step_size = 500
        try:
            for record in SeqIO.parse(genome_path, "fasta"):
                seq = str(record.seq)
                for i in range(0, len(seq) - chunk_size + 1, step_size):
                    chunks.append({
                        'contig': record.id,
                        'start': i,
                        'sequence': seq[i:i + chunk_size]
                    })
        except Exception:
            flash("Invalid genome FASTA file.")
            return redirect(request.url)

        results = []
        for chunk in chunks:
            results.extend(process_sequence_chunk(chunk))

        df = pd.DataFrame(results)
        if df.empty:
            flash("No binding sites detected.")
            return redirect(request.url)

        csv_path = os.path.join(app.config['UPLOAD_FOLDER'], "results.csv")
        plot_path = os.path.join(app.config['UPLOAD_FOLDER'], "score_distribution.png")
        html_report_path = os.path.join(app.config['UPLOAD_FOLDER'], "results_report.html")
        df.to_csv(csv_path, index=False)

        plt.figure(figsize=(8, 4))
        sns.histplot(df["score"], bins=30, kde=True)
        plt.title("Predicted Binding Affinity Score Distribution")
        plt.xlabel("Score")
        plt.ylabel("Frequency")
        plt.tight_layout()
        plt.savefig(plot_path)
        plt.close()

        with open(html_report_path, "w") as f:
            f.write("<html><head><title>SinR Binding Site Report</title></head><body>")
            f.write("<h1>Top Predicted Binding Sites</h1>")
            f.write(df.head(10).to_html(index=False, escape=False))
            f.write("<h2>Score Distribution Plot</h2>")
            f.write(f'<img src="{url_for("download_plot")}" alt="Score Plot">')
            f.write("</body></html>")

        return render_template("results.html",
                               results=df.head(10).to_dict(orient='records'),
                               csv_path=url_for('download_csv'),
                               plot_path=url_for('download_plot'),
                               report_path=url_for('download_report'))

    return render_template("index.html")

@app.route('/download/csv')
def download_csv():
    path = os.path.join(app.config['UPLOAD_FOLDER'], "results.csv")
    return send_file(path, as_attachment=True)

@app.route('/download/plot')
def download_plot():
    path = os.path.join(app.config['UPLOAD_FOLDER'], "score_distribution.png")
    return send_file(path, mimetype='image/png')

@app.route('/download/report')
def download_report():
    path = os.path.join(app.config['UPLOAD_FOLDER'], "results_report.html")
    return send_file(path, mimetype='text/html')

if __name__ == '__main__':
    import multiprocessing
    multiprocessing.set_start_method('spawn')
    app.run(debug=True, port=5050)