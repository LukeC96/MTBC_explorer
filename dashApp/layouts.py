from helpers import *
from dash import html, dcc, dash_table as dt
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto

# region network explorer page layout

hidden_content = html.Div(
    [     
        # label row
        dbc.Row([
            dbc.Col([
                html.H4('Expression over conditions')
            ],
            width = 12)
        ]),
        dbc.Row([
            # gene condition box plot
            dbc.Col([
                dcc.Loading(id="loading-boxplot", children = dcc.Graph(id='boxplot'))
            ],
            width = 12)
        ]),

        # label row
        dbc.Row([
            dbc.Col([
                html.H4('Cytoscape Network'),
                html.P('Configure the parameters below to adjust the network. Nodes can be coloured by either module colour or by level of expression of condition. Transcript types can be filtered and the edge weight threshold adjusted.')
            ],
            width = 12)
        ]),

        dbc.Row([
            dbc.Col(
                html.Label("Colour Mode:"),  # Label with spacing
                width="auto"
            ),
            dbc.Col(
                dbc.Switch(
                    id='colour-toggle-switch',
                    label="Module Colour / Condition Colour",
                    value=False,   # Default is Module Colour
                    style={'margin-top': '10px'}
                ),
                width="auto"
            ),
            dbc.Col(
                dcc.Dropdown(
                    id='condition-dropdown',
                    options=[{'label': cond, 'value': cond} for cond in dend_labels],
                    placeholder='Select condition'
                ),
                width="6"
            )
        ], align="center", justify="start"),  # Align items vertically centered and aligned to the left

        dbc.Row([
            # gene info 
            dbc.Col([
                dcc.Dropdown(
                id='transcript-type-filter',
                options=[{'label': t_type, 'value': t_type} for t_type in unique_transcript_types],
                multi=True,
                value = [],
                placeholder='Select transcript types'
                )
            ],
            width = 6),
            dbc.Col([
                dcc.Slider(
                    id='weight-threshold-slider',
                    min=0.1,
                    max=1,
                    step=0.1,
                    value=0.2,  # Initial value
                    marks=mark_labels,  # Slider marks from 0.1 to 1
                )
            ], width=6)
        ], style = {'margin-top': '20px'}),

        # cytoscape graph
        dbc.Row([ 
            dcc.Loading(
                id="loading-cytoscape-network",
                children = (
                    html.Div(style=cytoscape_styles['container'], children=[
                        html.Div(className='cy-container', style=cytoscape_styles['cy-container'], children=[
                            cyto.Cytoscape(
                                id='cytoscape-responsive-layout',
                                elements=[],
                                stylesheet=default_style_sheet,
                                style=cytoscape_styles['cytoscape'],
                                layout={
                                    'name': 'cose',
                                    'idealEdgeLength': 175,
                                    'nodeOverlap': 10,
                                    'refresh': 20,
                                    'fit': True,
                                    'padding': 30,
                                    'randomize': False,
                                    'componentSpacing': 100,
                                    'nodeRepulsion': 500000,
                                    'edgeElasticity': 100,
                                    'nestingFactor': 5,
                                    'gravity': 80,
                                    'numIter': 1000,
                                    'initialTemp': 200,
                                    'coolingFactor': 0.95,
                                    'minTemp': 1.0
                                },
                                responsive=True
                            )
                        ])
                    ])
                )
            )
        ], style = {'margin-top': '20px'}),

        #legends
        dbc.Row([
            dbc.Col(
                html.Div(id='transcript-type-counts')
                , width=6
            ),
            dbc.Col(
                html.Div(id='colour-scale'
                         , style={'display': 'none'}  # Hide by default
                )
                , width=6
                
            )
        ], style={'margin-top': '20px'}),
        html.Br(),
        
        # label row
        dbc.Row([
            dbc.Col([
                html.H4('Network gene table')
            ],
            width = 12)
        ]),

        # network transcript table
        dbc.Row([
            dcc.Loading(
                id="loading-network-transcript-table",
                children = dt.DataTable(
                    id = 'network-transcript-table',
                    columns = network_transcript_table_columns,
                    data = [],
                    page_size = 10,
                    page_current = 0,
                    sort_action= 'native'
                )
            )
        ])
    ]
)

networkExplorer = html.Div(
    [
        dcc.Store(id='gene-selected-store', data={'selected': False}),
        # gene table
        dbc.Row([
            dbc.Col([
                html.H4('Gene table')
            ],
            width = 6),
            dbc.Col([
                dbc.Input(type="search", placeholder="Search gene table...", id="search-input")
            ],
            width = 5),
            dbc.Col([
                dbc.Button("Search", color="dark", id="search-button"),
            ],
            width = 1)
        ], style = {'margin-top': '20px'}),
        html.Br(),
        dbc.Row([
            dcc.Loading(
                id="loading-gene-table",
                children = dt.DataTable(
                    id = 'gene-table',
                    columns = gene_table_columns,
                    data = complete_gene_info_df.to_dict('records'),
                    page_size = 10,
                    page_current = 0,
                    sort_action= 'native',
                    row_selectable= 'single'
                )
            )
        ]),
       
        html.Div(hidden_content, id='hidden-layout', style={'display': 'none'})
    ]   
)
#endregion


# region info page 

info = html.Div(
    [
        # gene table
        dbc.Row([
            dbc.Col([
                html.H4('Information')
            ])
        ], style = {'margin-top': '20px'}),
        dbc.Row([
            dbc.Col([
                html.P('This web application is designed to help users explore co-expression relationships in the Mycobacterium tuberculosis (Mtb) transcriptome. Utilisng Weighted Gene Co-expression Network Analysis (WGCNA), users can visualize and analyze networks of genes and non-coding RNAs (ncRNAs) to better understand their roles in the context of various biological conditions.')
            ])
        ], style = {'margin-top': '20px'}),
    ]   
)

#endregion