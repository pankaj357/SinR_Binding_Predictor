@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        protein_file = request.files.get('protein_fasta')  # match HTML name
        genome_file = request.files.get('genome_fasta')    # match HTML name
        known_motifs = request.form.get('known_motifs')    # match HTML name

        if not protein_file or not genome_file or not known_motifs:
            flash("Please upload all required files and provide known motifs.")
            return redirect(request.url)

        # TEMP: Just return success for now
        return render_template("results.html", results={"status": "Success!"})

    return render_template("index.html")



        protein_filename = secure_filename(protein_file.filename)
        genome_filename = secure_filename(genome_file.filename)
        protein_path = os.path.join(app.config['UPLOAD_FOLDER'], protein_filename)
        genome_path = os.path.join(app.config['UPLOAD_FOLDER'], genome_filename)
        protein_file.save(protein_path)
        genome_file.save(genome_path)

        # Read protein info