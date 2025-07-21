import os
import numpy as np
from Bio import SeqIO
import shutil
import requests
import pandas as pd

from flask import Flask, render_template, request, send_file, flash, redirect
from werkzeug.utils import secure_filename

app = Flask(__name__)
app.secret_key = 'supersecretkey'

# Configure upload folder
UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'ab1'}
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
os.makedirs(UPLOAD_FOLDER, exist_ok=True)


def movingaverage(data_in, window_size):
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(data_in, window, 'same')


def process_ab1_file(ab1_file, window_size=10, qual_cutoff=30):
    low_cutoff = 15
    sample_seq = SeqIO.read(ab1_file, "abi")
    sample_seq.id = sample_seq.name
    fasta_file = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".fasta")
    qual_file  = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".qual")
    SeqIO.write(sample_seq, fasta_file, "fasta")
    SeqIO.write(sample_seq, qual_file, "qual")

    sample_qual       = SeqIO.read(qual_file, "qual")
    sample_qual_score = sample_qual.letter_annotations["phred_quality"]
    sample_qual_MA    = np.array(movingaverage(sample_qual_score, window_size))
    max_q             = np.max(sample_qual_MA)

    if max_q >= qual_cutoff:
        active_cutoff = qual_cutoff
    else:
        active_cutoff = low_cutoff

    print(f"Max MA quality = {max_q:.1f}, using cutoff = {active_cutoff}")

    if max_q > low_cutoff:
        above_idxs      = np.where(sample_qual_MA > active_cutoff)[0]
        start, end      = above_idxs.min(), above_idxs.max()
        trimmed_seq     = sample_seq[start:end]
        trimmed_qual    = sample_qual[start:end]

        trim_fasta_file = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".trim.fasta")
        trim_qual_file  = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".trim.qual")
        SeqIO.write(trimmed_seq,   trim_fasta_file, "fasta")
        SeqIO.write(trimmed_qual,  trim_qual_file,  "qual")

        fa_file = os.path.join(UPLOAD_FOLDER, sample_seq.id + ".fa")
        shutil.copy(trim_fasta_file, fa_file)

        return {
            "fasta": fasta_file,
            "qual": qual_file,
            "trim_fasta": trim_fasta_file,
            "trim_qual": trim_qual_file,
            "fa": fa_file
        }
    else:
        print(f"Max quality {max_q:.1f} below low cutoff {low_cutoff}")
        return None


def read_fasta(fastafile):
    sequences = []
    with open(fastafile, "r") as f:
        for line in f:
            sequences.append(line.rstrip("\n"))
    headers = [l for l in sequences if l.startswith(">")]
    indices = [sequences.index(h) for h in headers]
    seq_dic = {}
    for i, h in enumerate(headers):
        start = indices[i] + 1
        end   = indices[i+1] if i+1 < len(indices) else len(sequences)
        seq_dic[h] = "".join(sequences[start:end])
    return seq_dic

def swap_dna(dnastring):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = []
    end = len(dnastring) - (len(dnastring) % 3) - 1
    for i in range(0, end, 3):
        codon = dnastring[i:i+3]
        protein.append(table.get(codon, "N"))
    return "".join(protein)

def rev_seq(seq):
    comp = {'A':'T','C':'G','G':'C','T':'A'}
    return "".join(comp.get(b,b) for b in seq)[::-1]

def frame_id(seq):
    frames = {}
    seq_rev = rev_seq(seq)
    for offset, names in zip(range(3), [('+1','-1'),('+2','+3'),('right','left')]):
        fwd = swap_dna(seq[offset:])
        rev = swap_dna(seq_rev[offset:])
        frames[names[0]] = fwd
        frames[names[1]] = rev
    return frames

def gen_frames(seq_dict):
    return {hdr: frame_id(seq) for hdr, seq in seq_dict.items()}

def oframe(amino):
    out = []
    for i, aa in enumerate(amino):
        if aa == 'M':
            segment = amino[i:]
            stop = segment.find('_')
            out.append(segment[:stop+1] if stop != -1 else segment)
    return out

def find_prots(frames_dict):
    prots = {}
    for hdr, frames in frames_dict.items():
        candidates = []
        for seq in frames.values():
            candidates += oframe(seq)
        # pick longest
        best = max(candidates, key=len, default="")
        prots[hdr] = best
    return prots

class annotate:
    def __init__(self, aaseq, scheme):
        self.aaseq = aaseq
        self.scheme = scheme
    
    def __repr__(self):
        return "Annotation of VH or VL sequence using Kabat, Chothia, Contact, or IMGT scheme"
    
    def output(self, chain, lst, regionlst):
        self.chain = chain
        self.lst = lst
        self.regionlst = regionlst

        self.regiondict, self.numberdict = {}, {}
        for i in range(0, len(self.lst), 2):
            self.numberdict[self.lst[i]] = self.lst[i+1]
        
        if self.scheme == "kabat":
            print("Annotation in Kabat scheme:")
        elif self.scheme == "chothia":
            print("Annotation in Chothia scheme:")
        elif self.scheme == "contact":
            print("Annotation in Contact scheme:")
        else:
            print("Annotation in IMGT scheme:")
        
        if self.chain == "L":
            print("L-FR1:  ", self.regionlst[0])
            print("L-CDR1: ", self.regionlst[1])
            print("L-FR2:  ", self.regionlst[2])
            print("L-CDR2: ", self.regionlst[3])
            print("L-FR3:  ", self.regionlst[4])
            print("L-CDR3: ", self.regionlst[5])
            print("L-FR4:  ", self.regionlst[6])
            
            for region, seq in zip(["L-FR1", "L-CDR1", "L-FR2", "L-CDR2", "L-FR3", "L-CDR3", "L-FR4"], self.regionlst):
                self.regiondict[region] = seq
            
            return [self.regiondict, self.numberdict]
        else:
            print("H-FR1:  ", self.regionlst[0])
            print("H-CDR1: ", self.regionlst[1])
            print("H-FR2:  ", self.regionlst[2])
            print("H-CDR2: ", self.regionlst[3])
            print("H-FR3:  ", self.regionlst[4])
            print("H-CDR3: ", self.regionlst[5])
            print("H-FR4:  ", self.regionlst[6])
            
            for region, seq in zip(["H-FR1", "H-CDR1", "H-FR2", "H-CDR2", "H-FR3", "H-CDR3", "H-FR4"], self.regionlst):
                self.regiondict[region] = seq
            
            return [self.regiondict, self.numberdict]
    
    def analyze(self, chain, lst):
        self.chain = chain
        self.lst = lst
        if self.chain == "L":
            self.L_FR1, self.L_CDR1, self.L_FR2, self.L_CDR2, self.L_FR3, self.L_CDR3, self.L_FR4 = ["" for i in range(7)]
            try:
                if self.scheme in ["kabat", "chothia"]:
                    self.L_FR1 = "".join([self.lst[i+1] for i in range(0, self.lst.index("L24"), 2)])
                    self.L_CDR1 = "".join([self.lst[i+1] for i in range(self.lst.index("L24"), self.lst.index("L35"), 2)])
                    self.L_FR2 = "".join([self.lst[i+1] for i in range(self.lst.index("L35"), self.lst.index("L50"), 2)])
                    self.L_CDR2 = "".join([self.lst[i+1] for i in range(self.lst.index("L50"), self.lst.index("L57"), 2)])
                    self.L_FR3 = "".join([self.lst[i+1] for i in range(self.lst.index("L57"), self.lst.index("L89"), 2)])
                    self.L_CDR3 = "".join([self.lst[i+1] for i in range(self.lst.index("L89"), self.lst.index("L98"), 2)])
                    self.L_FR4 = "".join([self.lst[i+1] for i in range(self.lst.index("L98"), len(self.lst), 2)])
                return [self.L_FR1, self.L_CDR1, self.L_FR2, self.L_CDR2, self.L_FR3, self.L_CDR3, self.L_FR4]
            except ValueError:
                print("Unable to retrieve complete V region. Make sure the sequence has complete V region")
            except:
                print("An error occurred in the `analyze()` method")
        else:
            self.H_FR1, self.H_CDR1, self.H_FR2, self.H_CDR2, self.H_FR3, self.H_CDR3, self.H_FR4 = ["" for i in range(7)]
            try:
                if self.scheme == "kabat":
                    self.H_FR1 = "".join([self.lst[i+1] for i in range(0, self.lst.index("H31"), 2)])
                    self.H_CDR1 = "".join([self.lst[i+1] for i in range(self.lst.index("H31"), self.lst.index("H36"), 2)])
                    self.H_FR2 = "".join([self.lst[i+1] for i in range(self.lst.index("H36"), self.lst.index("H50"), 2)])
                    self.H_CDR2 = "".join([self.lst[i+1] for i in range(self.lst.index("H50"), self.lst.index("H66"), 2)])
                    self.H_FR3 = "".join([self.lst[i+1] for i in range(self.lst.index("H66"), self.lst.index("H95"), 2)])
                    self.H_CDR3 = "".join([self.lst[i+1] for i in range(self.lst.index("H95"), self.lst.index("H103"), 2)])
                    self.H_FR4 = "".join([self.lst[i+1] for i in range(self.lst.index("H103"), len(self.lst), 2)])
                elif self.scheme == "chothia":
                    self.H_FR1 = "".join([self.lst[i+1] for i in range(0, self.lst.index("H26"), 2)])
                    self.H_CDR1 = "".join([self.lst[i+1] for i in range(self.lst.index("H26"), self.lst.index("H33"), 2)])
                    self.H_FR2 = "".join([self.lst[i+1] for i in range(self.lst.index("H33"), self.lst.index("H52"), 2)])
                    self.H_CDR2 = "".join([self.lst[i+1] for i in range(self.lst.index("H52"), self.lst.index("H57"), 2)])
                    self.H_FR3 = "".join([self.lst[i+1] for i in range(self.lst.index("H57"), self.lst.index("H95"), 2)])
                    self.H_CDR3 = "".join([self.lst[i+1] for i in range(self.lst.index("H95"), self.lst.index("H103"), 2)])
                    self.H_FR4 = "".join([self.lst[i+1] for i in range(self.lst.index("H103"), len(self.lst), 2)])
                return [self.H_FR1, self.H_CDR1, self.H_FR2, self.H_CDR2, self.H_FR3, self.H_CDR3, self.H_FR4]
            except ValueError:
                print("Unable to retrieve complete V region. Make sure the sequence has complete V region")
            except:
                print("An error occurred in the `analyze()` method")
    
    def retrieve(self):
        self.url = "http://www.bioinf.org.uk/abs/abnum/abnum.cgi"
        try:
            if self.scheme not in ["kabat", "chothia", "contact", "imgt"]:
                raise Exception
        except ValueError:
            print("Incorrect scheme mode. Must be one of: kabat, chothia, contact, imgt (in lowercase)")
            return None
        else:
            if self.scheme == "kabat":
                self.sche = "-k"
            else:
                self.sche = "-c"
        
        try:
            # Prepare params with cleaned sequence
            clean_seq = self.aaseq.replace('_', 'X').replace('N', 'X')
            if len(clean_seq) < 50:
                print("Protein sequence too short for reliable annotation.")
                return None
            
            self.d = {"plain": 1, "scheme": self.sche, "aaseq": clean_seq}
            self.myPage = requests.get(self.url, params=self.d)
            self.text = self.myPage.text
            self.lst = self.text.split()
            print(f"Parameters sent: {self.d}")
            print(f"Response split list: {self.lst}")
            
            if len(self.lst) > 1 and not self.lst[0].startswith('#'):
                self.chain = self.lst[0][0]
                self.result = self.output(self.chain, self.lst, self.analyze(self.chain, self.lst))
                return self.result
            else:
                print("No annotation retrieved or error from server. Possibly incomplete VH/VL sequence or ambiguous residues.")
                return None
        except Exception as e:
            print(f"An error occurred in the `retrieve()` method: {e}")
            return None

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    analysis_results = None

    if request.method == 'POST':
        file = request.files.get('file')
        if not file or file.filename == '':
            flash('No file selected')
            return redirect(request.url)

        scheme = request.form.get('annotation_scheme', 'chothia')
        if not allowed_file(file.filename):
            flash('Only .ab1 files allowed')
            return redirect(request.url)

        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        try:
            base   = os.path.splitext(filename)[0]
            result = process_ab1_file(filepath)

            if result:
                seqs   = read_fasta(result['fa'])
                frames = gen_frames(seqs)
                prots  = find_prots(frames)
                best   = list(prots.values())[0]
                cleaned = best.replace('N','X').replace('_','X')

                annotator = annotate(cleaned, scheme)
                ann_res    = annotator.retrieve()  # may be None or [regions, mapping]

                # Output filenames
                out_fa   = f"{base}_output_fasta.fa"
                best_fa  = f"{base}_best_frame.fa"
                ann_txt  = f"{base}_annotation_result.txt"
                excel    = f"{base}_annotation_summary.xlsx"
                vhvl_fa  = f"{base}_vh_vl_full_length.fa"

                # 1) Write all sequences FASTA
                with open(os.path.join(UPLOAD_FOLDER, out_fa), 'w') as f:
                    for h, s in seqs.items():
                        f.write(f"{h}\n{s}\n")

                # 2) Write best protein frame
                with open(os.path.join(UPLOAD_FOLDER, best_fa), 'w') as f:
                    f.write(f">{base}_Best_Protein_Frame\n{cleaned}\n")

                # 3) Write annotation text file
                with open(os.path.join(UPLOAD_FOLDER, ann_txt), 'w') as f:
                    if ann_res:
                        f.write(f"Annotation Scheme: {scheme}\n\n")
                        f.write("Regions:\n" + str(ann_res[0]) + "\n\n")
                        f.write("Number mapping:\n" + str(ann_res[1]) + "\n")
                    else:
                        f.write("Annotation failed or incomplete.\n")

                # 4) Write Excel & VH/VL FASTA (always)
                if ann_res:
                    regions_map = ann_res[0]
                    num_map     = ann_res[1]
                    chain_type  = 'Heavy' if any(k.startswith('H-') for k in regions_map) else 'Light'
                    cdr1 = regions_map.get(f"{chain_type[0]}-CDR1", '')
                    cdr2 = regions_map.get(f"{chain_type[0]}-CDR2", '')
                    cdr3 = regions_map.get(f"{chain_type[0]}-CDR3", '')
                    full = ''.join(num_map.values())
                else:
                    chain_type, cdr1, cdr2, cdr3, full = ('N/A',)*5

                # Excel summary
                df = pd.DataFrame({
                    'Name': [filename],
                    'Chain': [chain_type],
                    'FullSequence': [full],
                    'CDR1': [cdr1],
                    'CDR2': [cdr2],
                    'CDR3': [cdr3]
                })
                df.to_excel(os.path.join(UPLOAD_FOLDER, excel), index=False)

                # VH/VL FASTA
                with open(os.path.join(UPLOAD_FOLDER, vhvl_fa), 'w') as f:
                    f.write(f">{base}_VHVL\n{full}\n")

                # Prepare template context
                analysis_results = {
                    'original_fasta': seqs,
                    'best_protein':   cleaned,
                    'annotation':     ann_res or ['N/A', 'N/A'],
                    'fasta_file':     out_fa,
                    'best_frame':     best_fa,
                    'annotation_result': ann_txt,
                    'excel_file':     excel,
                    'vh_vl_fasta':    vhvl_fa
                }
            else:
                flash('Failed to process AB1: low overall quality.')
        except Exception as e:
            flash(f'Error: {e}')

    return render_template('upload.html', analysis_results=analysis_results)


@app.route('/download/<filename>')
def download_file(filename):
    return send_file(os.path.join(app.config['UPLOAD_FOLDER'], filename), as_attachment=True)


if __name__ == '__main__':
    app.run(debug=True)
