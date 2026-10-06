from shiny import App
from helper import *

for program in get_exes():
    path = os.path.join(BLAST_DIR, program)

    if os.path.exists(path):
        try:
            os.chmod(path, 0o755)
            print(f"Permission set: {program}")
        except Exception as e:
            print(f"Cannot chmod {program}: {e}")

os.chmod(os.path.join(BASE_DIR, "diamond_lin"), 0o755)

def server(input, output, session):
    query = ""

    current_gene = reactive.Value(None)
    index = reactive.Value(None)
    next_gene = reactive.Value(None)
    previous_gene = reactive.Value(None)

    current_page = reactive.Value("home")

    isoforms = reactive.Value(None)
    current_isoform = reactive.Value(None)
    sequence_type = reactive.Value("genome")
    fasta_header = reactive.Value(None)

    upstream = reactive.Value(0)
    downstream = reactive.Value(0)

    current_sequences = reactive.Value({})
    current_best_hit = reactive.Value({})

    current_homology_subject = reactive.Value("None")
    current_homology_dg = reactive.Value(None)
    current_homology_df = reactive.Value(None)
    current_homology_desc = reactive.Value("Valid gene description will appear here.")
    current_homology_query = reactive.Value(None)
    current_homology_parameters = reactive.Value(None)
    current_homology_selected = reactive.Value(None)
    current_homology_alignment = reactive.Value(None)

    current_align_input_type = reactive.Value(None)
    current_align_input = reactive.Value(None)
    current_align_input_string = reactive.Value(None)
    current_alignment = reactive.Value(None)

    @output
    @render.ui
    def ref_page():
        ref_texts = []

        with open("ref") as file:
            for line in file:
                ref_texts.append(ui.p(line))
        
        return ui.div(*ref_texts)

    @output
    @render.ui
    def about_page():
        about_texts = []

        with open("about") as file:
            for line in file:
                about_texts.append(
                    ui.h6(
                        line,
                        style="text-align: justify;"
                    )
                )
        
        return ui.div(*about_texts)

    @output
    @render.ui
    def logs_page():
        log_texts = []

        with open("logs") as file:
            for line in file:
                if line[0] == "v":
                    log_texts.append(ui.br())
                    log_texts.append(ui.p(ui.strong(line)))

                else:
                    log_texts.append(ui.p(line))
        
        return ui.div(*log_texts)
    
    @output
    @render.ui
    def homology_page_left_sidebar():
        return ui.card(
            ui.h4(ui.tags.u("Homology Analysis")),
            
            ui.input_radio_buttons(
                "homology_type_btn",
                "Choose an alignment type:",
                choices = {
                    "blastn": "blastn (BLAST+ only)",
                    "blastx": "blastx (DNA)",
                    "blastp": "blastp (amino acids)"
                }
            ),

            ui.panel_conditional(
                "input.homology_type_btn != 'blastn'",
                ui.input_radio_buttons(
                    "homology_tool_type_btn",
                    "Choose an alignment tool:",
                    choices = [
                        "BLAST+",
                        "DIAMOND"
                    ],
                ),
            ),

            ui.navset_card_tab(  
                ui.nav_panel(
                    "Sequence", 
                    ui.div(
                        ui.input_text_area(
                            "homology_text_input",
                            "Insert ONE sequence:",
                            rows=5,
                            width="100%",
                        ),

                        ui.div(
                            ui.tags.i(ui.h6("Example: ")),
                        
                            ui.tags.i(
                                ui.h6(
                                    ui.input_action_link(
                                    "example_alignment_link_1", 
                                    "AT1G01010"
                                    ),
                                ),
                            ),

                            ui.tags.i(ui.h6(", ")),

                            ui.tags.i(
                                ui.h6(
                                    ui.input_action_link(
                                    "example_alignment_link_2", 
                                    "AT3G18550"
                                    ),
                                ),
                            ),
                        style="display:flex; align-items:start;",
                        ),
                    ),
                    value="text",
                ),

                ui.nav_panel(
                    "Ensembl ID", 
                    ui.input_text(
                        "homology_ensembl_input",
                        "Insert an Ensembl Gene ID:",
                        value="AT1G01010"
                    ),

                    ui.div(
                        ui.tags.i(ui.h6("Example: ")),
                    
                        ui.tags.i(
                            ui.h6(
                                ui.input_action_link(
                                "example_alignment_link_3", 
                                "Glyma_12G040000"
                                ),
                            ),
                        ),

                        ui.tags.i(ui.h6(", ")),

                        ui.tags.i(
                            ui.h6(
                                ui.input_action_link(
                                "example_alignment_link_4", 
                                "LOC114162257"
                                ),
                            ),
                        ),
                    style="display:flex; align-items:start;",
                    ),

                    ui.output_ui("homology_ensembl_desc"),
                    value="ensembl"
                ),

                id="homology_method_tab"
            ),

            ui.card(
                ui.strong("Advanced Settings"),

                ui.input_numeric(
                    "homology_evalue",
                    "E-value",
                    value=1
                )
            ),

            ui.input_action_button(
                "homology_run_btn",
                "Submit"
            ),
        )
    
    @output
    @render.ui
    def homology_page_right_sidebar():
        return ui.card(
            ui.div(
                ui.h6("Showing results for: ", current_homology_subject()),
                
                style="display:space-between; align-items:start;",
            ),
            
            ui.output_ui("homology_output")
        )
    
    @output
    @render.ui
    def homology_page():
        return ui.div(
            ui.layout_columns(

                ui.output_ui(
                    "homology_page_left_sidebar",
                    height = "750px"),

                ui.output_ui("homology_page_right_sidebar"),
                
                col_widths=[4, 8],
            ),
            
            ui.output_ui("homology_align"),
        )
    
    @output
    @render.ui
    def homology_align():
        if not current_homology_dg():
            return 
        
        homology_align_type = {
            "blastn": "genome",
            "blastx": "cds",
            "blastp": "peptide"
        }

        df = current_homology_df()
        sel = homology_output_data_frame.cell_selection()

        out = [
            str(df.iloc[i, 0])
            for i in sel["rows"]
        ]

        hits_dict = {}
        type = current_homology_query()[2]

        for i in out:
            try:
                if type == "blastn":
                    query = genome_ann[i.split(".")[0]]
                
                else:
                    query = genome_ann[i]
            
            except Exception:
                return ui.card(ui.h6(ui.tags.i("Select a hits row to show alignment.")))

            hits_dict[i] = get_sequence(query, homology_align_type[type])[0]

        if not hits_dict:
            return ui.card(ui.h6(ui.tags.i("Select a hits row to show alignment.")))

        current_homology_selected.set({
            "hits": hits_dict,
            "type": type
        })

        ui.update_text(
            "gene_search",
            value=next(iter(hits_dict))
        )

        return ui.card(
            ui.output_ui("homology_render_display_align"),

            #ui.output_ui("input_align"),
             style="""
                border: 1px solid #ccc;
                border-radius: 5px;
            """
        )

    @output
    @render.ui
    def homology_render_display_align():
        input_string = []

        if current_homology_selected():
            hits = current_homology_selected()["hits"]

            for name, seq in hits.items():
                input_string.append(f">{name}\n{seq}")
        
        input_header = current_homology_subject()[2]
        
        input_string.append(f">{input_header}\n{current_homology_query()[1]}")

        align_input = "\n".join(input_string)

        try:
            msa = get_alignment(align_input)
            current_homology_alignment.set(get_msa_to_fasta_string(msa))

        except Exception as e:
            return e

        return ui.div(
                    ui.output_image("homology_render_align_image"),
                    style="""
                        width: 100%;
                        overflow-x: hidden;
                        overflow-y: auto;
                        max-height: 300px;
                        white-space: nowrap;
                    """
        )
    
    @output
    @render.image
    def homology_render_align_image():
        try:
            align_image = render_alignment_image(current_homology_alignment(), False, "Identity")

            return align_image

        except Exception as e:
            return e
        
    @reactive.effect
    @reactive.event(input.homology_run_btn)
    def _update_homology_query():
        current_homology_parameters.set(
            [
                str(input.homology_evalue())
            ]
        )

        current_homology_selected.set(None)

        if input.homology_method_tab() == "text":
            if not input.homology_text_input():
                return
            
            raw_input = input.homology_text_input()

            if ">" in raw_input:
                header = raw_input.split("\n")[0][1:]
                raw_input = "".join(raw_input.splitlines()[1:])
            
            else:
                raw_input = "".join(raw_input.splitlines())
                header = raw_input[:10] + "..." + f"[{len(raw_input)}]"
            
            current_homology_query.set([input.homology_tool_type_btn(), raw_input, input.homology_type_btn()])
            current_homology_subject.set([
                input.homology_tool_type_btn(), 
                " | ", 
                header,
                " | ", 
                input.homology_type_btn()]
                )
        
        elif input.homology_method_tab() == "ensembl":
            if not input.homology_ensembl_input():
                return

            try:
                query = input.homology_ensembl_input()
                result = fetch_ensembl_gene(query, BLAST_TYPE[input.homology_type_btn()], True)
                current_homology_query.set([input.homology_tool_type_btn(), result[2], input.homology_type_btn(), result[:2]])

                current_homology_subject.set([
                    input.homology_tool_type_btn(), 
                    " ", 
                    input.homology_ensembl_input(),
                    " ", 
                    input.homology_type_btn()]
                )

                current_homology_desc.set([
                    current_homology_query()[3][0],
                    current_homology_query()[3][1],
                ])
            
            except Exception:
                current_homology_query.set(None)
    
    @output
    @render.ui
    def homology_output():
        return ui.div(
            ui.output_data_frame("homology_output_data_frame"),
            #ui.download_button("download_homology_csv", "Download CSV", class_="btn-sm p-1"), style="align-items: end;"
        )
    
    @render.download(
        filename="homology_results.csv"
    )
    def download_homology_csv():
        if not current_homology_dg():
            yield current_homology_dg().to_csv(index=False),

    @output
    @render.data_frame
    def homology_output_data_frame():
        df_height = "750px"

        if not current_homology_query():
            return render.DataGrid(pd.DataFrame({"Status": ["Results will show here"]}), height=df_height)
        
        tool = current_homology_query()[0]   
        seq = current_homology_query()[1]
        type = current_homology_query()[2]

        if type == "blastn":
            tool = "BLAST+"
        
        blast_df = get_homology_df(
            seq, type, tool,
           current_homology_parameters()[0]
        )

        if blast_df.empty:
            blast_df = pd.DataFrame({"Status": ["No hits"]})
        
        blast_dg = render.DataGrid(blast_df, selection_mode="row", height=df_height)

        current_homology_df.set(blast_df)
        current_homology_dg.set(blast_dg)

        return blast_dg
        
    @output
    @render.ui
    def homology_ensembl_desc():
        desc = current_homology_desc()

        if len(desc) == 2:
            return ui.card(
                        ui.h6("Gene Name: ", desc[0]),
                        ui.h6("Description: ", desc[1])
            )
        
        return ui.card(ui.h6(desc))

    @reactive.effect
    @reactive.event(input.gene_search_btn)
    def gene_validate():
        try:
            check = genome_ann[input.gene_search()]

            if input.gene_search().count(".") != 0:
                check = genome_ann[input.gene_search().split(".")[0]]
            
            set_new_current_gene(check)
        
        except Exception as e:
            return ui.p(f"Gene not found ({e})")
    
    def set_new_current_gene(gene):
        upstream.set(0)
        downstream.set(0)

        if not gene:
            current_gene.set(None)
            current_isoform.set(None)
            isoforms.set(None)
            current_sequences.set(None)
            next_gene.set(None)
            previous_gene.set(None)
        
        else:
            current_gene.set(gene)
            sequence_type.set("genome")
            index.set(gene_ids.index(current_gene().id))

            if len(gene_ids[index() + 1]) == ID_LENGTH:
                next_gene.set(gene_ids[index() + 1])
            else:
                next_gene.set(None)
                
            if len(gene_ids[index() - 1]) == ID_LENGTH:
                previous_gene.set(gene_ids[index() - 1])
            else:
                previous_gene.set(None)

            current_sequences.set(get_all_annotation(current_gene()))
            isoforms.set(list(current_sequences().keys()))
            current_isoform.set(genome_ann[next(iter(current_sequences().keys()))])
            current_align_input_type.set(None)
            current_align_input.set(None)
            current_align_input_string.set(None)
            current_alignment.set(None)

            result = proteome_ann[
                proteome_ann["Protein_id"] == current_gene().id
            ]

            row = result.iloc[0]

            column = {
                "Soybean_best_hit": "Soybean",
                "Ath_best_hit": "Ath"
            }
        
            best_hits = {}

            for col, label in column.items():
                try:
                    value = row[col]

                    if value and len(value) > 4:
                        if label == "Soybean":
                            best_hits["Soybean"] = "_".join(["Glyma", row["Soybean_best_hit"].split(".")[1]])
                        
                        else:
                            best_hits[label] = row[col]
                    
                except Exception as e:
                    print(e)
            
            hits = fetch_ensembl_gene(list(best_hits.values()))
            current_best_hit.set(hits)

    @reactive.effect
    @reactive.event(input.example_A)
    def example_A():
        ui.update_text(
            "gene_search",
            value="Ponpi04G041600"
        )
        set_new_current_gene(genome_ann["Ponpi04G041600"])
    
    @reactive.effect
    @reactive.event(input.example_B)
    def example_B():
        ui.update_text(
            "gene_search",
            value="Ponpi01G000100"
        )
        set_new_current_gene(genome_ann["Ponpi01G000100"])
    
    def gene_nav():
        return ui.div(
            previous_button(),

            next_button(),

            style="""
                display: flex;
                justify-content: space-between;
                gap: 10px;
            """
        )
    
    @output
    @render.ui
    def gene_result():
        if not current_gene():
            return
        
        return ui.card(
            gene_nav(),

            ui.tags.u(ui.h5("Gene Information")), 
            gene_report(),
            
            ui.tags.u(ui.h5("Sequences Display")),
            ui.output_ui("gene_display_bar"),
            
            ui.tags.u(ui.h5("Protein Homologs")), 
            ui.output_ui("gene_align"),
        )
    
    def gene_report():
        query = current_gene()

        if not query:
            return None
        
        strand = "forward" if query.strand == "+" else "reverse"
        location = f"{query.seqid}:{query.start}..{query.end} {strand}"

        result = proteome_ann[
            proteome_ann["Protein_id"] == query.id
        ]

        row = result.iloc[0]

        return ui.div(
            gene_desc(query.id, location, isoforms(), row),

            gene_jbrowse(),
        )
                
    @output
    @render.ui
    def gene_align():
        return ui.div(
            ui.output_ui("render_display_align"),

            #ui.output_ui("input_align"),
             style="""
                border: 1px solid #ccc;
                border-radius: 5px;
            """
        )

    @output
    @render.ui
    def render_display_align():
        gene = genome_ann[current_gene()]
        input_string = []
         
        for isoform in isoforms():
            seq = "".join(get_all_annotation(gene)[isoform]["peptide"].values())
            input_string.append(f">{isoform}\n{seq}")

        current_align_input.set("\n".join(input_string))

        return ui.div(
            ui.output_ui("homolog_align")
        )
    
    @output
    @render.ui
    def homolog_align():
        homolog_string = []

        if current_best_hit():
            for id, seq in current_best_hit().items():
                homolog_string.append(f">{id}\n{seq}")
        
            align_input = current_align_input() + "\n" + "\n".join(homolog_string)
        
        else:
            align_input = current_align_input()

        #print(align_input)

        try:
            msa = get_alignment(align_input)
            current_alignment.set(get_msa_to_fasta_string(msa))

        except Exception as e:
            return e

        return ui.navset_card_tab(  
            ui.nav_panel(
                "Multiple Sequence Alignment", 
                ui.div(
                    ui.output_image("render_align_image"),
                    style="""
                        width: 100%;
                        overflow-x: hidden;
                        overflow-y: auto;
                        max-height: 300px;
                        white-space: nowrap;
                    """
                    ),
            ),

            ui.nav_panel("Phylogenetic Tree", ui.output_plot("render_phylo_tree")),

            id="align_output_tab"
        )
    
    @reactive.effect
    @reactive.event(input.example_alignment_link_1)
    def example_alignment_1():
        #seq1 = get_ensembl_gene("AT1G01010", "cds")
        #seq2 = get_ensembl_gene("AT3G18550", "genome")

        placeholder = align_example["AT1G01010"]

        ui.update_text_area(
            "homology_text_input",
            value=placeholder
        )

    @reactive.effect
    @reactive.event(input.example_alignment_link_2)
    def example_alignment_2():
        #seq1 = get_ensembl_gene("AT1G01010", "cds")
        #seq2 = get_ensembl_gene("AT3G18550", "genome")

        placeholder = align_example["AT3G18550"]
    
        ui.update_text_area(
            "homology_text_input",
            value=placeholder
        )
    
    @reactive.effect
    @reactive.event(input.example_alignment_link_3)
    def example_alignment_3():
        #seq1 = get_ensembl_gene("AT1G01010", "cds")
        #seq2 = get_ensembl_gene("AT3G18550", "genome")

        ui.update_text_area(
            "homology_ensembl_input",
            value="Glyma_12G040000"
        )
    
    @reactive.effect
    @reactive.event(input.example_alignment_link_4)
    def example_alignment_4():
        #seq1 = get_ensembl_gene("AT1G01010", "cds")
        #seq2 = get_ensembl_gene("AT3G18550", "genome")

        ui.update_text_area(
            "homology_ensembl_input",
            value="LOC114162257"
        )
            
    @reactive.effect
    @reactive.event(input.submit_align_btn)
    def _update_align_input_type():
        gene = genome_ann[input.isoform_align_btn().split(".")[0]]
        isoform = genome_ann[input.isoform_align_btn()]
        seq = "".join(get_all_annotation(gene)[isoform.id][input.display_align_btn()].values())
        seq_length = len(seq)

        current_align_input.set(f">{isoform.id} length={seq_length}\n{seq}\n")
        current_align_input_type.set(input.align_input_type())

        if current_align_input_type() == "text":
            align_input = input.query_align_text_input()
            
            if not align_input.strip():
                return
        
        elif current_align_input_type() == "file":
            file_input = input.query_align_file_input()

            if not file_input:
                return
            
            path = file_input[0]["datapath"]

            with open(path, "r", encoding="utf-8") as f:
                align_input = f.read()
        
        current_align_input_string.set(align_input)
    
    @output
    @render.ui
    def output_align():
        if not current_align_input_string():
            return
        
        align_input = current_align_input() + current_align_input_string()

        try:
            msa = get_alignment(align_input)
            current_alignment.set(get_msa_to_fasta_string(msa))

        except Exception as e:
            return e

        return ui.navset_card_tab(  
            ui.nav_panel(
                "Multiple Sequence Alignment", 
                ui.div(
                    ui.output_image("render_align_image"),
                    style="""
                        width: 100%;
                        overflow-x: hidden;
                        overflow-y: auto;
                        max-height: 300px;
                        white-space: nowrap;
                    """
                    ),
            ),

            ui.nav_panel("Phylogenetic Tree", ui.output_plot("render_phylo_tree")),

            id="align_output_tab"
        )
    
    @render.plot
    def render_phylo_tree():
        try:
            align_plot = render_alignment_phylo_tree(current_alignment())

            return align_plot

        except Exception as e:
            return e
    
    @output
    @render.image
    def render_align_image():
        try:
            align_image = render_alignment_image(current_alignment())

            return align_image

        except Exception as e:
            return e

    def gene_desc(id, location, isoforms, row):
        columns = {
            "Preferred_name": "Preferred Name: ",
            "Description": "Description: ",
            "Ath_best_hit": "A. thaliana Homolog: ",
            "Soybean_best_hit": "Soybean Homolog: ",
            "PFAMs": "Associated PlantFAMs: ",
        }

        column_ui = []

        for col, header in columns.items():
           
            text = row[col]

            try:
                if not text or len(text) <= 4:
                    text = "-"
                
                else:
                    if col == "Soybean_best_hit":
                        text = "_".join(text.split("."))
                    
                    if col == "PFAMs":
                        text = ", ".join(text.split(","))
                
            except Exception:
                text = "-"

            column_ui.append(
                ui.h6(
                    ui.tags.i(header), 
                    text
                ),
            )

        return ui.div(
            ui.h6(
                ui.tags.i("Gene Identifier: "), 
                id
            ),

            ui.h6(
                ui.tags.i("Location: "), 
                location, f" [len={len(get_sequence(current_gene())[0])}]"
            ),
            
            ui.h6(
                ui.tags.i("Transcript(s): "), 
                ", ".join(isoforms)
            ),

            ui.div(*column_ui),

            ui.br()
        )
    
    def gene_jbrowse():
        return ui.div(
            output_widget("render_jbrowse"),
            style="""
            width: 100%;
            height: 260px;
            overflow: hidden;
            position: relative;
            """
        ),
    
    @render_widget
    def render_jbrowse():
        if not current_gene():
            return
        
        gene = current_gene()

        loc = (
                f"{gene.seqid}:"
                f"{gene.start-500}.."
                f"{gene.end+500}"
        )

        view = get_jbrowse()

        view.location = loc
        
        return view
    
    @output
    @render.ui
    def gene_display_bar():
        transcripts = isoforms()
        colors = get_colours()

        return ui.div(
            ui.div(
                ui.input_select(
                    id="seq_type_select",
                    label=None,
                    choices={
                        "genome": "Genomic Sequence",
                        "transcript": "Transcript Sequence",
                        "cds": "CDS Sequence",
                        "peptide": "Peptide Sequence"
                    },
                    selected=sequence_type()
                ),

                ui.input_action_button(
                    id="seq_type_refresh",
                    label="Go",
                    class_="btn-sm",
                    width="90px",
                ).add_style(
                    "height: 36px;"
                ),
                style="display:flex; align-items:start;",
            ),

            ui.panel_conditional(
                    "input.seq_type_select !== 'genome'",
                    ui.input_select(
                        "isoform_select",
                        None,
                        choices=transcripts,
                        selected=current_isoform().id
                    )
            ),

            ui.panel_conditional(
                "input.seq_type_select === 'genome'",
                ui.div(
                    ui.div(
                        ui.h6(
                            ui.tags.em("Show flanking sequence: "),

                            ui.tags.em("Upstream - "),
                                
                            ui.input_numeric(
                                "upstream_box",
                                None,
                                value=upstream(),
                                width="80px",
                                min=0
                            ).add_style(
                                "height: 20px;"
                            ),

                            ui.tags.em("Downstream - "),

                            ui.input_numeric(
                                "downstream_box",
                                None,
                                value=downstream(),
                                width="80px",
                                min=0
                            ).add_style(
                                "height: 20px;"
                            ),
                            style="display: flex; align-items: center; gap: 10px;"
                        )
                    ),

                    ui.div(
                        *[
                            ui.span(
                                ui.span(
                                    style=f"""
                                    display:inline-block;
                                    width:11px;
                                    height:11px;
                                    background:{color};
                                    margin-right:5px;
                                    border:1px solid black;
                                    """
                                ),
                                name,
                                style="""
                                margin-right:15px;
                                align-items:flex-start;
                                """
                            )
                            for name, color in colors.items()
                        ],
                        style="""
                        padding:10px;
                        display:flex;
                        flex-wrap:wrap;
                        gap:10px;
                        """
                    ),

                    style="""
                        display: flex;
                        justify-content: space-between;
                        width: 100%;
                    """
                )
            ),
            ui.output_ui("render_sequence_display")
        )

    @reactive.effect
    @reactive.event(input.seq_type_refresh)
    def gene_sequence_display():
        sequence_type.set(input.seq_type_select())
        current_isoform.set(genome_ann[input.isoform_select()])
        current_sequences.set(get_all_annotation(current_gene()))
        
        validate = up_and_down_validate(input.upstream_box(), input.downstream_box())

        upstream.set(validate[0])
        downstream.set(validate[1])
    
    def update_sequence_display():
        if not current_gene() or not current_isoform():
            return

        if sequence_type() == "genome":
            query = current_gene()
        else:
             query = current_isoform()

        sequences = get_sequence(query, sequence_type(), upstream(), downstream())
        fasta_header.set(" ".join([get_header(query, sequence_type(), upstream(), downstream()), f"[len={len("".join(sequences))}]"]))
        annotation = current_sequences()[current_isoform().id][sequence_type()]
        html_ann = get_html(annotation)

        return [html_ann, sequences[1], sequences[2]]
    
    @output
    @render.ui
    def render_sequence_display(): 
        seq = update_sequence_display()
        
        if not seq:
            return

        highlight = seq[0]
        upstream = seq[1]
        downstream = seq[2]

        return ui.div(
                ui.HTML(
                    f"""
                    <div style="
                        font-family:monospace;
                        margin:0;
                        padding:0;
                        line-height:1;
                    ">
                        &gt;{fasta_header()}<br>
                        <span>{upstream}</span>{highlight}<span>{downstream}</span>
                    </div>
                    """
                ),
                style="height: 250px; overflow-y: auto;"
        )
    
    def previous_button():
        if previous_gene():
            return ui.input_action_button(
                "previous_gene_btn",
                label=[icon_svg("arrow-left"), " ", ui.tags.i(previous_gene())],
                class_="btn-sm",
                width="150px"
            ).add_style(
                "height: 36px;"
            ),
        
        else:
            return ui.div(
                style="visibility: hidden; width: 100px;"
            )
    
    def next_button():
        if next_gene():
            return ui.input_action_button(
                "next_gene_btn",
                label=[ui.tags.i(next_gene()), " ", icon_svg("arrow-right")],
                class_="btn-sm",
                width="150px"
            ).add_style(
                "height: 36px;"
            )
        
        else:
            return ui.div(
                style="visibility: hidden; width: 100px;"
            )
    
    @reactive.effect
    @reactive.event(input.next_gene_btn)
    def go_next_gene():
        if next_gene():
            gene = next_gene()
        else:
            return

        ui.update_text(
            "gene_search",
            value=gene
        )

        new_gene = genome_ann[gene]
        set_new_current_gene(new_gene)

    @reactive.effect
    @reactive.event(input.previous_gene_btn)
    def go_previous_gene():
        if previous_gene():
            gene = previous_gene()
        else:
            return 

        ui.update_text(
            "gene_search",
            value=gene
        )

        new_gene = genome_ann[gene]
        set_new_current_gene(new_gene)
    
    @reactive.effect
    def _update_page():
        return
        
        if current_page() != input.homepage():
            current_page.set(input.homepage())

            ui.update_text(
                "gene_search",
                value=""
            )
            set_new_current_gene(None)
        
        #print(current_page())

app = App(
    get_ui(), 
    server,
    static_assets={
        "/www": str(WWW_DIR.resolve())
    }
)
