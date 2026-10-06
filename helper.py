from aifc import Error
import os
from pathlib import Path
from shiny import render, ui, reactive
import gui
import gffutils
import pandas as pd
from faicons import icon_svg
from jbrowse_anywidget import LinearGenomeView
from shinywidgets import render_widget, output_widget
from pyfaidx import Fasta
from Bio.Seq import Seq
from pymsaviz import MsaViz
import tempfile
from pyfamsa import Aligner, Sequence
from Bio import SeqIO
from io import StringIO
import requests
import matplotlib.pyplot as plt
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import io
from Bio.Blast import NCBIWWW
import xml.etree.ElementTree as ET
import json
import subprocess

#BASE_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = Path(__file__).parent
WWW_DIR = BASE_DIR / "www"
BLAST_DIR = BASE_DIR / "BLAST_LIN"

ID_LENGTH = 14

EXECUTABLES = [
    "blastn",
    "blastp",
    "blastx",
    "tblastn",
    "tblastx"
]

COLORS = {
    "Exon": "#576ae6",
    "5' UTR": "#20e709b6",
    "3' UTR": "#ff4b4bee",
}

FASTA = {
    "genome": Fasta("www/data/genome.fa.gz"),
    "transcript": Fasta("www/data/transcript.fa.gz"),
    "cds": Fasta("www/data/cds.fa.gz"),
    "peptide": Fasta("www/data/proteome.fa.gz")
}

FEATURES = {
    "genome": ["five_prime_utr", "three_prime_utr", "exon", "intron"],
    "transcript": ["five_prime_utr", "three_prime_utr", "exon"],
    "cds": ["CDS"],
    "peptide": ["CDS"]
}

BLAST_TYPE = {
        "blastn": "genomic",
        "tblastn": "cds",
        "blastp": "protein",
        "blastx": "cds",
        "tblastx": "cds"
}

genome_ann = gffutils.FeatureDB("www/data/genome_ann.db")
proteome_ann = pd.read_excel("www/data/proteome.xlsx")
gene_ids = proteome_ann["Protein_id"].astype(str).tolist()

align_example = {
    "AT1G01010": "ATGGAGGATCAAGTTGGGTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTTGGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATCAGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTTCCAGTCAAAGTACAAATCGAGAGATGCTATGTGGTACTTCTTCTCTCGTAGAGAAAACAACAAAGGGAATCGACAGAGCAGGACAACGGTTTCTGGTAAATGGAAGCTTACCGGAGAATCTGTTGAGGTCAAGGACCAGTGGGGATTTTGTAGTGAGGGCTTTCGTGGTAAGATTGGTCATAAAAGGGTTTTGGTGTTCCTCGATGGAAGATACCCTGACAAAACCAAATCTGATTGGGTTATCCACGAGTTCCACTACGACCTCTTACCAGAACATCAGAGGACATATGTCATCTGCAGACTTGAGTACAAGGGTGATGATGCGGACATTCTATCTGCTTATGCAATAGATCCCACTCCCGCTTTTGTCCCCAATATGACTAGTAGTGCAGGTTCTGTGGTCAACCAATCACGTCAACGAAATTCAGGATCTTACAACACTTACTCTGAGTATGATTCAGCAAATCATGGCCAGCAGTTTAATGAAAACTCTAACATTATGCAGCAGCAACCACTTCAAGGATCATTCAACCCTCTCCTTGAGTATGATTTTGCAAATCACGGCGGTCAGTGGCTGAGTGACTATATCGACCTGCAACAGCAAGTTCCTTACTTGGCACCTTATGAAAATGAGTCGGAGATGATTTGGAAGCATGTGATTGAAGAAAATTTTGAGTTTTTGGTAGATGAAAGGACATCTATGCAACAGCATTACAGTGATCACCGGCCCAAAAAACCTGTGTCTGGGGTTTTGCCTGATGATAGCAGTGATACTGAAACTGGATCAATGATTTTCGAAGACACTTCGAGCTCCACTGATAGTGTTGGTAGTTCAGATGAACCGGGCCATACTCGTATAGATGATATTCCATCATTGAACATTATTGAGCCTTTGCACAATTATAAGGCACAAGAGCAACCAAAGCAGCAGAGCAAAGAAAAGGTGATAAGTTCGCAGAAAAGCGAATGCGAGTGGAAAATGGCTGAAGACTCGATCAAGATACCTCCATCCACCAACACGGTGAAGCAGAGCTGGATTGTTTTGGAGAATGCACAGTGGAACTATCTCAAGAACATGATCATTGGTGTCTTGTTGTTCATCTCCGTCATTAGTTGGATCATTCTTGTTGGTTAA",
    "AT3G18550": "CCTCTCTCTCCATATAGCTCTCCCTCTCTCTCTCTTTTTCTCTCTCTTTCTCTCCATTCTAGTTCATAGCACATAGACTCATTTAAGTACGTCGCCTGTGTTGTCTTTCTCAACCATTAGTCCTTTGCTTCGATAAAGGTATTCTCTCTCTACACATTTAGACACACACATAGATTCTCTAGAACAAAGACATATTGTACTATCTGACTAATCACTATTACTGCTCTTTTTTTTTCTTGATTTTGCAGGCTAAGGAAAAAAAGAAATCTCTCTTTCTAATTTTCTTGTTTTCATTTTGGTTGAGACTTGAGATCAATGTCTTCGTCATTGGAACAGTAACCTTTCTCTTTGATTCAAACTTTCTTTCATTTTCTTCTCTTTTTTCTTGTGTGTCTGAAGAAGATTATAGTACAAGTGTCATTCTCAAAATTTTGTCTTAGGGTTTTTTAATTTTGTTACTTCAAAAACCCCTAAAAAGGCATGAACAACAACATTTTCAGTACTACTACCACCATCAATGACGACTACATGTTATTCCCTTATAATGACCATTATTCCTCACAACCATTGCTCCCTTTTAGCCCTTCTTCTTCCATTAACGACATCTTGATTCACTCCACCTCTAACACATCAAACAATCATCTTGACCATCATCATCAATTCCAACAACCTTCTCCTTTTTCTCACTTCGAATTTGCCCCGGACTGCGCCCTCCTCACCTCTTTCCACCCAGAAAACAATGGCCATGATGATAACCAAACCATCCCAAACGACAATCATCATCCATCACTTCACTTTCCCTTGAACAACACCATTGTAGAACAACCCACTGAGCCCTCGGAAACTATAAACCTTATAGAAGATTCCCAGAGAATCTCAACTTCTCAAGACCCAAAAATGAAAAAAGCCAAGAAACCCAGCAGAACGGACAGGCACAGCAAGATCAAAACGGCCAAAGGGACACGAGATCGTAGGATGAGACTCTCGCTAGATGTCGCCAAAGAGTTGTTTGGCTTACAAGACATGCTTGGATTTGACAAAGCCAGCAAAACCGTTGAATGGTTGCTTACACAAGCAAAACCTGAGATCATAAAGATCGCGACAACCCTTTCTCACCATGGCTGCTTCAGCAGCGGCGATGAGTCTCATATCCGTAAGTATCCATTCTTTTCTAGTTTTAAGTCTATTCAAGTTTTCAACTTCTATGAATTTTTTTATATCATAATCTGAAAATTGTAAATTTCCTAATTAATTAAGGACCGGTGTTAGGATCCATGGACACATCTTCTGATCTATGTGAACTTGCATCCATGTGGACGGTCGACGATAGAGGCAGCAATACTAACACGACCGGTACGACCTTAACCGTTGCCTTTCTGATTTACAATTATATATAGAGTACTTTCTATCAAAATGTTACAAAAGATAAATTATCTTATACATTAGAAACAAGAGGAAACAAGGTCGATGGGAGATCGATGAGAGGGAAGAGAAAGAGGCCAGAACCGCGAACGCCCATTTTAAAGAAGTTGTCCAAGGAGGAGAGAGCGAAAGCTAGAGAAAGAGCAAAGGGTAGAACAATGGAGAAAATGATGATGAAGATGAAAGGAAGATCACAATTAGTGAAAGTTGTGGAAGAAGACGCTCATGATCATGGTGAGATAATAAAGAATAATAATAGAAGCCAAGTGAATCGGAGTTCTTTTGAGATGACACACTGCGAAGACAAGATCGAAGAACTTTGCAAGAACGATCGTTTTGCAGTTTGCAACGAATTTATCATGAATAAGAAAGATCACATTTCAAATGAATCTTATGACTTAGTCAACTACAAACCGAACTCATCATTCCCAGTGATTAACCACCATCGCAGCCAAGGAGCAGCTAATTCCATTGAGGTACCTTATTTATATCTACGAGTCACATTTATAATTTTTATGAGTTTTTTTCCAATGGTTTTATCTAAAGCTATACACAGTTAATTCAAGTTTTAGAGGTAGGGTATATCAAATTCAATGACACTAATGCACCACCATCGATTAATTAATGCTATATATAGTTTAATTAGAAGTTAGTTTTTGTGTTTTGTGGATTGGTACTTGGTTATTACTATCTAGTTCCACATCAGCGATTCAAAGGCGCCGGTTTTTTTTAACTGACTTGGTTATACTATTATGGTTTTTATTGTTCAAATTTTCATCGATTTTTTTTTTTTTGGGGCTATGTACCTTGCAGCAGCATCAGTTTACGGATCTTCATTACTCCTTCGGCGCGAAACCAAGAGACCTCATGCACAACTATCAAAACATGTATTGAAACTGATTTCATTTAACATTAATATCTCAAATTAATGATATATAGTTTGGTGAGAGAGATATGTTAATCCATCTTATGTTTTTTTTTCTTAACATCGTAATTAATCATGTTTATATGAAGGCTTTTTCGTCATTTTAGTGTCAAATATTAAACCTATGTGAAGAAATTAATGCAATTTGTATAATATTTTGTTTTCATAGTTCGAGTGGAT"
}

def get_exes():
    return EXECUTABLES

def get_ui():
    return gui.app_ui

def get_colours():
    return COLORS

def get_jbrowse():
    assembly = {
            "name": "Ppinnata",
            "uri": "www/data/genome.fa.gz",
        }
        
    view = LinearGenomeView(
        assembly=assembly,
        location="1:1..1000",
    )

    view.add_track(
        {
            "uri": "www/data/genome_ann.sorted.gff3.gz",
            "name": "Annotation",
        }
    )

    return view

def up_and_down_validate(input_up, input_down):
    try:
        if (input_up >= 0 and input_down >= 0) and not (input_up == 0 and input_down == 0):
            up = input_up
            down = input_down

        else:
            up = 0
            down = 0

        return (up, down)

    except Exception:
        return (0, 0)

def get_sequence(query, type="genome", up=0, down=0):

    fasta = FASTA[type]

    if type != "genome":    
        sequence = Seq(str(fasta[query.id]))

        return [str(sequence), "", ""]
    
    else:
        chrom = query.seqid
        start = query.start
        end = query.end

        sequence = Seq(str(fasta[chrom][start-1:end]))
        flank_up = Seq(str(fasta[chrom][start-1-up:start-1]))
        flank_down = Seq(str(fasta[chrom][end:end+down]))

        if query.strand == "-":
            sequence = sequence.reverse_complement()
            flank_up = flank_up.reverse_complement()
            flank_down = flank_down.reverse_complement()
    
        return [str(sequence), str(flank_up), str(flank_down)]

def get_header(isoform, type="genome", up=0, down=0):
    if type != "genome":
        if type == "transcript":
            return isoform.id
        
        elif type == "cds":
            return isoform.id + " CDS"
        
        elif type == "peptide":
            return isoform.id + ".p"
    
    else:
        try:
            if not up and not down:
                flanking = ""
            else:
                flanking = f"|upstream={up}|downstream={down}"
        
        except Exception:
            flanking = ""

        return f"P.pinnata Huanan_v1.0|{isoform.id.split(".")[0]}|{isoform.seqid}:{isoform.start}..{isoform.end} {"forward" if isoform.strand == "+" else "reverse"}{flanking}"

def get_all_annotation(query):
        annotation = {}
        chrom = query.seqid
        fasta = FASTA["genome"]
        is_reverse = query.strand != "+"

        for transcript in genome_ann.children(
            query,
            featuretype=["mRNA"],
            order_by="start"):

            annotation[transcript.id] = {}

            for type in ["genome", "transcript", "cds"]:

                if type == "genome":
                    intr_count = 0
                    temp_ann = {}

                    whole_seq = get_sequence(query)[0]

                    #temp_ann["intronxx"] = whole_seq

                    for feature in genome_ann.children(
                        transcript,
                        order_by="start",
                        reverse=is_reverse
                    ):

                        if feature.featuretype in FEATURES[type]:
                            seq = Seq(str(fasta[chrom][feature.start-1:feature.end]))

                            if is_reverse:
                                seq = seq.reverse_complement()
                                
                            seq = str(seq)

                            before = whole_seq.split(seq, 1)[0]

                            temp_ann[f"introx{intr_count}"] = before
                            intr_count += 1

                            whole_seq = whole_seq[len(before) + len(seq):]

                            temp_ann[feature.id.split(".")[-1]] = str(seq)

                else:
                    temp_ann = {}
                    
                    for feature in genome_ann.children(
                        transcript,
                        order_by="start"
                    ):

                        if feature.featuretype in FEATURES[type]:
                            seq = Seq(str(fasta[chrom][feature.start-1:feature.end]))

                            if is_reverse:
                                seq = seq.reverse_complement()
                            
                            temp_ann[feature.id.split(".")[-1]] = str(seq)
                
                    if is_reverse:
                        temp_ann = dict(reversed(list(temp_ann.items())))
                
                annotation[transcript.id][type] = temp_ann
            
            temp_ann = {}
            current_pos = 0

            if len("".join(annotation[transcript.id]["cds"].values())) % 3 == 0:

                protein = str(FASTA["peptide"][transcript.id])

                for cds, seq in annotation[transcript.id]["cds"].items():
                    cds_length = int(len(seq) / 3)
                    temp_ann[cds] = protein[current_pos:current_pos+cds_length]
                    current_pos += cds_length
            
            else:
                protein = str(Seq("".join(annotation[transcript.id]["cds"].values())).translate(to_stop=True))

                for cds, seq in annotation[transcript.id]["cds"].items():
                    remainder = len(seq) % 3
                    
                    if remainder == 1:
                        cds_length = int((len(seq) - 1) / 3)
                        temp_ann[cds] = protein[current_pos:current_pos+cds_length-1]
                        current_pos += cds_length - 1

                    elif remainder == 2:
                        cds_length = int((len(seq) + 1) / 3)
                        temp_ann[cds] = protein[current_pos:current_pos+cds_length+1]
                        current_pos += cds_length + 1

                    else:
                        cds_length = int(len(seq) / 3)
                        temp_ann[cds] = protein[current_pos:current_pos+cds_length]
                        current_pos += cds_length
                
            annotation[transcript.id]["peptide"] = temp_ann.copy()

        return annotation

def get_html(annotation):
    colors = {
        "exo": ["#abb6fe", "#d5dbff"],
        "CDS": ["#abb6fe", "#d5dbff"],
        "5UT": ["#9afb8fb6", "#c3f8bc"],
        "3UT": ["#f77a7aee", "#fdbfbfec"],
        "int": ["#ffffffd8", "#ffffff"] 
    }

    shade_index = {k: 0 for k in colors}
    html = []
    first = next(iter(annotation))[:3]
    current_shade = shade_index[first]
    shade_index[first] = 1 - current_shade

    for ann in annotation:

        text = annotation[ann]

        html.append(
                f'<span style="background:{colors[ann[:3]][current_shade]}">{text}</span>'
        )

        current_shade = shade_index[ann[:3]]
        shade_index[ann[:3]] = 1 - current_shade
        
    return "".join(html)

def get_homology_df(seq, type, tool, evalue=1):
    try:
        if tool == "BLAST+":
            db_path = os.path.join(
                BASE_DIR,
                "www",
                "data",
                "blast",
                "pinnata_blast_cdna_db" if type in ["blastn", "tblastn", "tblastx"] else "pinnata_blast_protein_db"
            )

            try:
                blast_exe = os.path.join(
                    BASE_DIR,
                    "BLAST_WIN",
                    type
                )

                cmd = [
                    blast_exe,
                    "-query", "-",
                    "-db", db_path,

                    "-outfmt",
                    "6 sseqid pident length mismatch gapopen "
                    "qstart qend sstart send evalue bitscore",

                    "-evalue", evalue
                ]

                result = subprocess.run(
                    cmd,
                    input=seq,
                    text=True,
                    capture_output=True,
                    check=True
                )

            except OSError:
                blast_exe = os.path.join(
                    BASE_DIR,
                    "BLAST_LIN",
                    type
                )

                cmd = [
                    blast_exe,
                    "-query", "-",
                    "-db", db_path,
                    "-outfmt",
                    "6 sseqid pident length mismatch gapopen "
                    "qstart qend sstart send evalue bitscore",

                    "-evalue", evalue
                ]

                result = subprocess.run(
                    cmd,
                    input=seq,
                    text=True,
                    capture_output=True,
                    check=True
                )
        
        elif tool == "DIAMOND":
            seq = f">Input\n{seq}"

            with open("temp/diamond_input.fasta", "w") as f:
                f.write(seq)
            
            try:
                cmd = [
                    'diamond_win', 
                    type, 
                    '--db', 
                    'www/data/blast/pinnata_diamond_protein_db.dmnd', 
                    '--outfmt', '6', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 
                    '--query', 'temp/diamond_input.fasta',

                    "--evalue", evalue
                    ]

                result = subprocess.run(
                    cmd, 
                    text=True,
                    capture_output=True,
                    check=True,
                    shell=False
                )
            
            except OSError:
                diamond_exe = os.path.join(
                    BASE_DIR,
                    "diamond_lin"
                )

                cmd = [
                    diamond_exe, 
                    type, 
                    '--db', 
                    'www/data/blast/pinnata_diamond_protein_db.dmnd', 
                    '--outfmt', '6', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 
                    '--query', 'temp/diamond_input.fasta',

                    "--evalue", evalue
                    ]

                result = subprocess.run(
                    cmd, 
                    text=True,
                    capture_output=True,
                    check=True,
                    shell=False
                )
            
        columns = [
            "Gene_ID",
            "Identity",
            "Length",
            "Mismatch",
            "Gap",
            "Q_Start",
            "Q_End",
            "S_Start",
            "S_End",
            "E_Value",
            "Bit_Score"
        ]

        if not result.stdout.strip():
            return pd.DataFrame(columns=columns)

        df = pd.read_csv(
            pd.io.common.StringIO(result.stdout),
            sep="\t",
            names=columns 
        )
            
        return df
    
    except Exception as e:
        print(e)
        return pd.DataFrame({"Status": ["Error: Invalid sequence."]})

def get_alignment(align_input):
    sequences = []

    for record in SeqIO.parse(StringIO(align_input), "fasta"):
        sequences.append(
            Sequence(
                str(record.id).encode('utf-8'), 
                str(record.seq).encode('utf-8'))
        )

    aligner = Aligner(guide_tree="upgma")
    msa = aligner.align(sequences)

    #for sequence in msa:print(sequence.id.decode().ljust(10), sequence.sequence.decode())

    return msa

def get_msa_to_fasta_string(input):
    fasta_string = ""
    
    for sequence in input:
        fasta_string += f">{str(sequence.id, "utf-8")}\n{str(sequence.sequence, "utf-8")}\n"

    return fasta_string

def render_alignment_image(align_input, consensus=True, color="Taylor"):
    msa = StringIO(align_input)

    mv = MsaViz(
        msa,
        color_scheme=color,
        wrap_length=70,
        show_count=True,
        show_consensus=consensus
    )

    # Temporary PNG
    output_file = tempfile.NamedTemporaryFile(
        suffix=".png",
        delete=False
    )
    output_file.close()

    mv.savefig(output_file.name)

    return {
        "src": output_file.name,
        "width": "100%",
        "height": "auto"
    } 

def render_alignment_phylo_tree(align_input):
    alignment_file = io.StringIO(align_input.strip())
    alignment = AlignIO.read(alignment_file, "fasta")

    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # 3. Construct the phylogenetic tree
    constructor = DistanceTreeConstructor()

    # Choose your method: 'nj' (Neighbor Joining) or 'upgma'
    tree = constructor.nj(distance_matrix) 

    # 4. Save and visualize the tree
    # Save to a standard Newick file
    Phylo.write(tree, "temp/temp_tree.nwk", "newick")

    n = len(tree.get_terminals())

    fig = plt.figure(
        figsize=(12, max(5, n * 0.3))
    )

    Phylo.draw(tree, axes=fig.gca(), do_show=False)

    return fig

def fetch_ensembl_gene(input, type="protein", isSingle=False):
    if not input:
        return
    
    if not isSingle:
        
        p = ui.Progress()
        p.set(40, message="Connecting to Ensembl", detail="Please wait..")

        try: 
            server = "https://rest.ensembl.org"
            ext = "/sequence/id"

            # Both content and accept headers must be JSON for POST requests
            headers = {
                "Content-Type": "application/json",
                "Accept": "application/json"
            }

            # Provide an array of IDs and enforce the translation type
            payload = {
                "ids": input,
                "type": type
            }

            response = requests.post(server + ext, headers=headers, json=payload)
            # Parse the JSON response
            batch_data = response.json()

            p.set(100, message="Sequence fetching successful!", detail="")
            p.close()

        except Exception as e:
            p.set(0, message="Error :(", detail="Please try again..")
            return

        result = {}

        for count in range(len(input)):
            result[input[count]] = batch_data[count]['seq']
    
    else:
        p = ui.Progress()

        p.set(20, message="ID submitted!", detail="Validating..")

        lookup_response = requests.get(
            f"https://rest.ensembl.org/lookup/id/{input}",
            headers={"Accept": "application/json"},
            timeout=30
        )

        if not lookup_response:
            p.set(80, message="Invalid ID!", detail="")
            return
        
        p.close()

        info = lookup_response.json()

        result = [
            info.get("display_name"),
            info.get("description"),
            fetch_ensembl_gene([input], type)[input]
        ]
        
    return result
