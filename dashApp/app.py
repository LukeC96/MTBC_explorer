import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
import plotly.graph_objs as go
from dash import Dash, html, dcc, Output, Input, State
from collections import Counter
from layouts import *
import re

external_stylesheets = [dbc.themes.BOOTSTRAP]
# start dash app
app = Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)
server = app.server

# region layout and navigation setup

navbar = dbc.NavbarSimple(
  [
    dbc.NavItem(
        dbc.NavLink(
            'Cytoscape network explorer', active=True, href=app.get_relative_path('/networkExplorer'), className="btn btn-light"
        ), style={"padding": "10px"}),
    dbc.NavItem(
        dbc.NavLink(
            'Info', active=True, href=app.get_relative_path('/info'), className="btn btn-light"
        ), style={"padding": "10px"})
  ], 
  brand="MTBC explorer",  
  brand_style={"color": "white"},
  color='dark', 
  dark=True,
  links_left=True
)

app.layout = html.Div([
  dcc.Location(id="page-url", refresh=False),
  navbar,
  dbc.Container(id="page-content"),
])

@app.callback(Output("page-content", "children"),
               [Input("page-url", "pathname")])
def display_layout(path_name):
    """
     Update the page layout based on the URL pathname

    Args:
        path_name (str): The pathname extracted from the URL indicating the current page being requested

    Returns:
        layout: The layout configured for the page in layouts.py
    """
    page_name = app.strip_relative_path(path_name)
    if page_name == 'networkExplorer':
        return networkExplorer
    if page_name == "info":
        return info

# endregion


# region callback functions

# following 2 callbacks ensure the rest of the layout below the gene table does not appear until a transcript is selected
@app.callback(
    Output('gene-selected-store', 'data'),
    [Input('gene-table', 'selected_rows')]  
)
def update_gene_selection(selected_rows):
    """
    Updates 'gene-selected-store' which is found in layouts.py

    Monitors if any row has been selected in the gene table, if a row is selected sets to true, else false

    Parameters:
    selected_rows (list[int]): The list of indices of selected rows in Gene table (can only select one though)

    Returns:
    dict: A dictionary with a key 'selected' boolean:  true if a row is selected, else false.
    """
    if selected_rows:
        return {'selected': True}
    return {'selected': False}

@app.callback(
    Output('hidden-layout', 'style'),
    [Input('gene-selected-store', 'data')]
)
def display_hidden_layout(selection_data):
    """
    Controls the visibility of a hidden layout section based on 'gene-selected-store'

    If value is true, displays rest of layout, else hides it.

    Parameters:
    selection_data (dict): A dictionary with a key 'selected' and boolean value

    Returns:
    dict: A style dictionary for the hidden layout: if a row is selected then display layout, else hide
    """
    if selection_data['selected']:
        return {'display': 'block',
            'marginBottom': '20px'}  # show rest of layout
    return {'display': 'none'}  # hide rest of layout


@app.callback(
    Output('gene-table', 'data'),
    [Input('search-button', 'n_clicks'),
     Input('search-input', 'n_submit')],
    [State('search-input', 'value')]
)
def filter_gene_table(n_clicks, n_submit, search_input):
    """
    Filter Gene table by user search input

    Args:
        n_clicks (int): Number of times the search button has been clicked
        n_submit (any): Allows enter to be pressed to trigger search
        search_input (str): User search input

    Returns:
        list[dict]: Filtered complete_gene_info_df
    """

    if search_input:
        filtered_data = complete_gene_info_df[
            complete_gene_info_df['gene_ID'].str.contains(search_input, case=False) |
            complete_gene_info_df['pred_name'].str.contains(search_input, case=False, na=False)
        ].to_dict('records')
        return filtered_data
    else:
        return complete_gene_info_df.to_dict('records')



# Callback to update box plot
@app.callback(
    Output('boxplot', 'figure'),
    [Input('gene-table', 'selected_rows'),
     Input('gene-table', 'data')]
)
def update_boxplot(selected_rows, filtered_data):
    """
    Update the condition boxplot by selected row in Gene table

    Args:
        selected_rows (list[int]): The list of indices of selected rows in Gene table (can only select one though)
        filtered_data (list[dict]): If Gene table has been filtered by user search then pass this in

    Returns:
        dict: Plotly figure object displaying the boxplots of normalised counts of expression for each condition for selected gene
                Returns empty dictionary if no gene selected
    """
    
    if selected_rows:
        if filtered_data:
        # Get gene_ID from selected row, if Gene table has been filtered use that input else use original dataframe
            selected_gene_id = filtered_data[selected_rows[0]]['gene_ID']
            selected_gene_name = filtered_data[selected_rows[0]]['pred_name']
        else:
            selected_gene_id = complete_gene_info_df.iloc[selected_rows[0]]['gene_ID']
            selected_gene_name = complete_gene_info_df.iloc[selected_rows[0]]['pred_name']

        # create dataframe of list of different conditions for labels
        txt_df = expression_conditions_df.loc[expression_conditions_df['gene_ID'].str.lower() == selected_gene_id.lower()].copy()
        txt_df['condition'] = pd.Categorical(txt_df['condition'], categories=dend_labels)

        # create plot object
        fig = px.box(txt_df, x='counts', y='condition', color='condition', title=selected_gene_name)
        fig.update_layout(xaxis_title='Normalised counts',yaxis_title='Condition', legend_title=None, height=700)
        fig.update_traces(showlegend=False)

        return fig
    else:
        return {}
    


# Callback to update cytoscape plot
@app.callback(
    Output('cytoscape-responsive-layout', 'elements'),
    [Input('gene-table', 'selected_rows'),
     Input('gene-table', 'data'),
     Input('transcript-type-filter', 'value'),
     Input('weight-threshold-slider', 'value')]
)
def update_cytoscape_plot(selected_rows, filtered_data, selected_transcript_types, selected_weight):
    """
    Update Cytopscape plot based on selected row in Gene table (use filtered date if applicable), transcript types selected and weight threshold set
    Creates nodes and edges that exceed the weight threshold using the adjacency matrix created in helpers.py from the R script

    Args:
        selected_rows (list[int]): The list of indices of selected rows in Gene table (can only select one though)
        filtered_data (list[dict]): If Gene table has been filtered by user search then pass this in
        selected_transcript_types (list[str]): List of selected transcript types for further filtering
        selected_weight (float): Weight threshold for filtering weights in the adjacency matrix

    Returns:
        list[dict]: A list of dictionaries representing the nodes and edges for the Cytoscape plot.
    """

    if selected_rows:
        if filtered_data:
            my_gene = filtered_data[selected_rows[0]]['gene_ID']
        else:
            my_gene = complete_gene_info_df.iloc[selected_rows[0]]['gene_ID']

        # get selected gene index and use to find all connected rows in adjacency matrix, filtered by selected weight 
        my_gene_index = complete_gene_info_df.index[complete_gene_info_df['gene_ID'] == my_gene].tolist()[0]
        weights_to_my_gene = numpy_adjacency_matrix[my_gene_index]
        connected_genes_indices = np.where(weights_to_my_gene > selected_weight)[0]
        connected_genes = complete_gene_info_df.iloc[connected_genes_indices]['gene_ID']
        connected_genes_list = connected_genes.tolist()
        
        nodes = []
        edges = []

        for gene_Id in connected_genes_list:

            # set node values and properties
            gene_info = complete_gene_info_df.loc[complete_gene_info_df['gene_ID'] == gene_Id]
            pred_name = gene_info['pred_name'].values[0]
            module_color = gene_info['moduleColor'].values[0]
            transcript_type = gene_info['transcript_type'].values[0]
            node_id = pred_name if pd.notnull(pred_name) else gene_Id

            # some module names contain numbers, remove these
            module_color = re.sub(r'\d+', '', module_color) 
            node_class = module_color

            # original selected gene has own class
            if gene_Id == my_gene:
                node_class += " selected"

            # can filter by transcript type so apply this
            if len(selected_transcript_types) == 0 or transcript_type in selected_transcript_types:
                nodes.append({'data': {'id': gene_Id, 'label': node_id, 'transcript_type': transcript_type, 'group': 'node'}, 'classes': node_class})

        node_degrees = {gene_Id: 0 for gene_Id in connected_genes_list}

        # Create edges from adjacency matrix
        for i in connected_genes_indices:
            for j in connected_genes_indices:
                if i != j and numpy_adjacency_matrix[i][j] > selected_weight:
                    source = complete_gene_info_df.iloc[i]['gene_ID']
                    target = complete_gene_info_df.iloc[j]['gene_ID']
                    if source in [node['data']['id'] for node in nodes] and target in [node['data']['id'] for node in nodes]:
                        weight = numpy_adjacency_matrix[i][j]
                        edges.append({'data': {'source': source, 'target': target, 'weight': weight, 'group': 'edge'}})

                        node_degrees[source] += 1
                        node_degrees[target] += 1

        max_degree = max(node_degrees.values(), default=1)
        
        # need to check if theres only one node
        if max_degree == 0:
            max_degree = 1
        min_size = 10
        max_size = 60

        for node in nodes:
            gene_id = node['data']['id']
            degree = node_degrees[gene_id]
    
            if max_degree == 1:  # if theres only one node or one edge
                normalised_size = (min_size + max_size) / 2
            else:
            # normalise degree based on min and max sizes
                normalised_size = min_size + ((degree / max_degree) * (max_size - min_size))
    
            # set normalised 
            node['data']['size'] = normalised_size

        weighted_elements = nodes + edges
       
        return weighted_elements
    else:
        return {}



@app.callback(
    [Output('network-transcript-table', 'data'),
     Output('transcript-type-counts', 'children')],
    [Input('gene-table', 'selected_rows'),
     Input('transcript-type-filter', 'value'),
     Input('cytoscape-responsive-layout', 'elements')]
)
def update_network_transcript_table(selected_rows, selected_transcript_types, elements):
    """
    Update the network transcript table which displays more information about the nodes in the network
    and display the counts of each transcript type within the generated network

    Args:
        selected_rows (list[int]): The list of indices of selected rows in Gene table (can only select one though)
        selected_transcript_types (list[str]): List of selected transcript types for further filtering
        elements (list[dict]): List of elements (Node and edges) in the Cytoscape plot (only need nodes but filter within function)

    Returns:
        tuple: A tuple containing:
            - list[dict]: List of dictionaries of further information for the transcripts within the network
            - str: Formatted string display counts of each transcript type in the network
    """

    if selected_rows:

        connected_genes_ids = [element['data']['id'] for element in elements if element['data'].get('group') == 'node']
        connected_genes_df = complete_gene_info_df.loc[complete_gene_info_df['gene_ID'].isin(connected_genes_ids)]

        # Filter by selected transcript types if provided
        if selected_transcript_types:
            connected_genes_df = connected_genes_df[connected_genes_df['transcript_type'].isin(selected_transcript_types)]

        # Count the number of each transcript type, formulate string to display counts
        transcript_type_counts = Counter(connected_genes_df['transcript_type'])

        # create legend item array
        legend_items = [
            {'shape': 'ellipse', 'label': 'Gene', 'count': transcript_type_counts.get('gene', 0)},
            {'shape': 'rectangle', 'label': 'Putative sRNA', 'count': transcript_type_counts.get('putative sRNA', 0)},
            {'shape': 'triangle', 'label': 'Putative UTR', 'count': transcript_type_counts.get('putative UTR', 0)}
        ]

        # configure legend dynamically based on counts of transcript types in current network
        legend_layout = html.Div([
        html.Div([
            html.Div(style={
                'display': 'inline-block',
                'width': '20px',
                'height': '20px',
                'backgroundColor': '#ccc',
                'marginRight': '10px',
                'border': '2px solid black',
                'borderRadius': '50%' if item['shape'] == 'ellipse' else '0',
                'clipPath': 'polygon(50% 0%, 0% 100%, 100% 100%)' if item['shape'] == 'triangle' else 'none',
            }),
            html.Span(f"{item['label']} ({item['count']})")
        ], style={'display': 'flex', 'alignItems': 'center', 'marginBottom': '10px'})
        for item in legend_items
        ], style={'border': '1px solid #ccc', 'padding': '10px', 'borderRadius': '5px'})

        return connected_genes_df.to_dict('records'), legend_layout
        
    else:
        return [], ""


@app.callback(
    Output('condition-dropdown', 'style'),
    Input('colour-toggle-switch', 'value')
)
def toggle_condition_dropdown(selected_value):
    if selected_value == True:
        return {'display': 'block'}
    return {'display': 'none'}

@app.callback(
    Output('cytoscape-responsive-layout', 'stylesheet'),
    [Input('colour-toggle-switch', 'value'),
     Input('condition-dropdown', 'value'),
     Input('cytoscape-responsive-layout', 'elements')]
)
def update_cytoscape_stylesheet(toggle_value, selected_condition, elements):
    """
    Update the Cytoscape stylesheet based on the selected colour input (if condition colour selected then also apply based on condition dropdown selection)
    Condition colour will scale from blue to red based on comparison of median transcript expression with overall maximum expression value found in the dataset

    Args:
        toggle_value (bool): Value of Colour Mode toggle switch, switches between Module Colour / Condition Colour
        selected_condition (str): The selected condition from the condition type dropdown
        elements (list[dict]): List of elements (Node and edges) in the Cytoscape plot (only need nodes but filter within function)

    Returns:
        list[dict]: A list of dictionaries representing the Cytoscape stylesheet with dynamic node colours based on either condition expression or module colour
    """


    # start with default style sheet
    dynamic_stylesheet = default_style_sheet.copy()
    node_elements = [element for element in elements if element['data'].get('group') == 'node']

    if toggle_value and selected_condition: 
        
        # get all condition counts for this condition
        condition_values = expression_conditions_df.loc[expression_conditions_df['condition'] == selected_condition].copy()
        
        # get the actual counts and find overall maximum
        count_column = condition_values.iloc[:, 2]
        max_value = count_column.max()

        for node in node_elements:
            gene_id = node['data']['id']
            gene_counts = condition_values.loc[condition_values['gene_ID'] == gene_id].iloc[:, 2].tolist()

            # use the median count for each node to calculate the colour relative to the max count value
            if gene_counts:
                gene_value = np.median(gene_counts)
                color_intensity = int((gene_value / max_value) * 255)
                color = f'rgb({color_intensity}, 0, {255 - color_intensity})'  # Example colour scale from blue to red
                dynamic_stylesheet.append({
                    'selector': f'node[id = "{gene_id}"]',
                    'style': {
                        'background-color': color,
                    }
                })

    return dynamic_stylesheet


@app.callback(
    Output('colour-scale', 'style'),
    Input('condition-dropdown', 'value')
)
def toggle_condition_selected(selected_value):
    if selected_value is not None:
        return {'display': 'block'}
    return {'display': 'none'}

@app.callback(
    Output('colour-scale', 'children'),
    [Input('condition-dropdown', 'value')]
)
def create_color_scale_legend(selected_condition):

    condition_values = expression_conditions_df.loc[expression_conditions_df['condition'] == selected_condition].copy()
    
    count_column = condition_values.iloc[:, 2]
    max_value = count_column.max()

    # define the color scale (from blue to red)
    colorscale = [[0, 'blue'], [1, 'red']]

    # dummy data array to generate the color gradient
    z = np.linspace(0, max_value, 100).reshape(1, 100)

    # heatmap for the color scale
    fig = go.Figure(data=go.Heatmap(
        z=z,
        colorscale=colorscale,
    ))

    # hide axes to focus on color scale only
    fig.update_layout(
        xaxis=dict(showticklabels=False, ticks=""),
        yaxis=dict(showticklabels=False, ticks=""),
        margin=dict(t=20, b=20, l=20, r=20),
        height=100,  
        width=500  
    )

    # return color scale figure as dcc.Graph component
    return dcc.Graph(figure=fig,
        config={
            'displayModeBar': False  # hide Plotly mode bar
        })


# endregion

if __name__ == '__main__':
    app.run(debug=True)

