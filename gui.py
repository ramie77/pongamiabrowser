from helper import *

app_ui = ui.page_fluid(

    ui.tags.style("""
        body {
            background-image: url("/www/bg.jpg");
            background-size: cover;
            background-position: center;
            background-repeat: no-repeat;
            background-attachment: fixed;
        }
        
        .page-overlay {
                background: rgba(255, 255, 255, 0.90);
                min-height: 100vh;
                padding: 30px;
            }
    """),

    ui.page_navbar(  
        ui.nav_panel(
            "Home",
            ui.div(
                ui.h2("Welcome to Pongamia Browser!"),
                ui.br(),
                ui.output_ui("logs_page"),
                class_="page-overlay"
            ),
            value="home"
        ),

        ui.nav_panel(
            "About",
            ui.div(
                ui.div(
                    ui.img(src='/www/flower.jpg', width='200px'), 

                    ui.h3(ui.tags.i(ui.strong("Pongamia pinnata"))),

                    style="""
                        display: flex;
                        justify-content: center;
                        align-items: center;
                        gap: 20px;
                    """,
                ),

                ui.br(),

                ui.output_ui("about_page"),

                class_="page-overlay"
            ),
            value="assembly"
        ),

        ui.nav_panel(
            "Gene Query",

            ui.div(

                ui.h6("Enter a Gene/Transcript ID:"),

                ui.div(
                    ui.input_text(  
                        "gene_search",  
                        "",   
                        value="",
                        width="300px"
                    ),

                    ui.input_action_button(
                        "gene_search_btn", 
                        "Search",
                        class_="btn-sm",
                        width="90px"
                    ).add_style(
                            "height: 36px;"
                    ),
                    style="display:flex; align-items:start;",
                ),

                ui.div(
                    ui.tags.i(ui.h6("Example: ")),
                    ui.tags.i(
                        ui.h6(
                            ui.input_action_link(
                                "example_A", 
                                "Ponpi04G041600"
                            ),
                        ),
                    ),
                    ui.tags.i(ui.h6(",")),
                    ui.tags.i(
                        ui.h6(
                            ui.input_action_link(
                                "example_B", 
                                "Ponpi01G000100"
                            ),
                        ),
                    ),
                    style="display:flex; align-items:start;",
                ),
                
                ui.output_ui("gene_result"), 

                class_="page-overlay"
            ),
            value="query"
        ),

        ui.nav_panel(
            "Homology Analysis",
            ui.div(
                ui.output_ui("homology_page"),
                class_="page-overlay",
            ),
            value="homology"
        ),

        ui.nav_panel(
            "References",
            ui.div(
                ui.output_ui("ref_page"),

                class_="page-overlay"
            ),
            value="ref"
        ),
        title="Pongamia Browser",  
        id="homepage"
    )
)