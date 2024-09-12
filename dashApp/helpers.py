import pyreadr
import pandas as pd
import numpy as np
from rpy2.robjects import r



# region load data

# Load required RData and CSV files, extract as dataframes where necessary

# Main gene info data for Gene table and network
complete_gene_info = pyreadr.read_r("dashApp/R_data/complete_geneInfo.Rdata")
complete_gene_info_df = complete_gene_info[None]

# Extra info such as functional categories
cds_info = pd.read_csv('dashApp/R_data/cds_info.csv')

# Condition labels
dend_labels = pd.read_csv("dashApp/R_data/dend_labels.csv")
dend_labels = dend_labels['condition'].unique()

# Expression condition counts
expression_conditions = pyreadr.read_r("dashApp/R_data/expression_cond_df.RData")
expression_conditions_df = expression_conditions[None]

# endregion




# region data manipulation

# extract transcript type for column in Gene table
complete_gene_info_df['transcript_type'] = complete_gene_info_df['gene_ID'].str.split(':').str[0].str.replace('_', ' ')
# get unique transcript types
unique_transcript_types = complete_gene_info_df['transcript_type'].unique()


# define a function to retrieve the module membership value based on moduleColor
def get_module_membership_value(row):
    """Finds highest module membership value for each row 

    Args:
        row (pd.Series): Row from a Dataframe containing the 'moduleColor' on which module membership is based

    Returns:
        float: The highest module membership value to 2 decimal places
    """

    mm_column_name = f"MM{row['moduleColor']}"
    return round(row[mm_column_name], 2)

complete_gene_info_df['MM_value'] = complete_gene_info_df.apply(get_module_membership_value, axis=1)

# get unique functional categories
functional_categories = cds_info['Mycobrowser Functional Category'].unique()

# merge the dataframes based on the gene_ID column
merged_df = pd.merge(complete_gene_info_df, cds_info[['gene_ID', 'Mycobrowser Functional Category', 'Name']], on='gene_ID', how='left')

# update the 'pred_name' column where it is null
merged_df['pred_name'] = merged_df.apply(
    lambda row: row['Name'] if pd.isnull(row['pred_name']) else row['pred_name'],
    axis=1
)

# drop the 'Name' column as it's no longer needed
merged_df.drop(columns=['Name'], inplace=True)

# update with merged data
complete_gene_info_df = merged_df


# remove transcript types from gene id column as moved to separate column
complete_gene_info_df['gene_ID'] = complete_gene_info_df['gene_ID'].str.replace(r'.*:', '', regex=True)
expression_conditions_df['gene_ID'] = expression_conditions_df['gene_ID'].str.replace(r'.*:', '', regex=True)


# endregion


# region table column definitions

gene_table_column_mapping = {
    'gene_ID': 'Gene Id',
    'pred_name': 'Name/ Predicted Name',
    'transcript_type': 'Transcript Type',
    'moduleColor': 'Module Colour',
    'MM_value': 'Module Membership Value'
}

# define columns for datatable based on mapped names
gene_table_columns = [{"name": gene_table_column_mapping[key], "id": key} for key in gene_table_column_mapping.keys()]

# columns for network gene table
network_transcript_table_column_mapping = {
    'gene_ID': 'Gene Id',
    'pred_name': 'Name/ Predicted Name',
    'transcript_type': 'Transcript Type',
    'moduleColor': 'Module Colour',
    'Mycobrowser Functional Category': 'Mycobrowser Functional Category'
}

# define columns for datatable based on mapped names
network_transcript_table_columns = [{"name": network_transcript_table_column_mapping[key], "id": key} for key in network_transcript_table_column_mapping.keys()]

# endregion



# region load r script and data for network graph

# run r script to calculate adjacency matrix and return matrix
r.source("dashApp/wgcna_adjacency_matrix_script.R")
calculate_adjacency_matrix = r["calculate_adjacency_matrix"]
adjacency_matrix = calculate_adjacency_matrix()

# convert r matrix to numpy array
numpy_adjacency_matrix = np.array(adjacency_matrix)

# endregion


# region slider marks

# work out marks for slider (there is a bug where it doesnt display 1.0 correctly)

mark_values = np.linspace(0.1, 1, 10) # whatever computes the position of the marks
mark_labels = {}
for mark_val in mark_values:
    if abs(mark_val-round(mark_val)) < 1e-3: 
        mark_val = int(mark_val)
    mark_labels[mark_val] = {"label": str(round(mark_val, 2))}

# endregion




# region cytoscape network styles
cytoscape_styles = {
   'container': {
       #'position': 'fixed',
       'display': 'flex',
       'flex-direction': 'column',
       'height': '600px',
       'width': '100%'
   },
   'cy-container': {
       'flex': '1',
       'position': 'relative'
   },
   'cytoscape': {
       'position': 'absolute',
       'width': '100%',
       'height': '100%',
       'zIndex': 999,
       'label': 'data(label)',
       'border': '2px solid #000',  # container border
       'background-color': '#f0f0f0'  # container background colour
   }
}

default_style_sheet = [
    {
        'selector': 'node',
        'style': {
            'label': 'data(label)',
            'border-color': '#000',
            'border-width': '2px',
            "width": "mapData(size, 0, 100, 10, 60)",
            "height": "mapData(size, 0, 100, 10, 60)"
        }
    },
    {
        'selector': 'node.selected',
        'style': {
            'label': 'data(label)',
            'border-color': 'red',
            'border-width': '4px',
            'color': 'red'
        }
    },
    # change shape of node based on transcript type
    {
        'selector': 'node[transcript_type = "gene"]',
        'style': {
            'shape': 'ellipse' 
        }
    },
    {
        'selector': 'node[transcript_type = "putative sRNA"]',
        'style': {
            'shape': 'rectangle'  
        }
    },
    {
        'selector': 'node[transcript_type = "putative UTR"]',
        'style': {
            'shape': 'triangle' 
        }
    },
]
# get all unique colours from raw dataframe
unique_module_colors = complete_gene_info_df['moduleColor'].unique()

for color in unique_module_colors:
    default_style_sheet.append({
        'selector': f'node.{color}',
        'style': {
            'background-color': color
        }
    })

# default style for edges
default_style_sheet.append({
    'selector': 'edge',
    'style': {
        'width': 'mapData(weight, 0, 1, 1, 5)',  # edge width dependent on weight
        'line-color': '#ccc'
    }
})
# endregion