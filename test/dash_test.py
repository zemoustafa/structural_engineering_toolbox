import sys
sys.path.append(r'C:\_Github\structural_engineering_toolbox')
import dash
from dash import dcc, html, Input, Output, State, dash_table
import pandas as pd
from etabs_tools import etabs_api
from design_functions import as3600_column_design
from design_reports import design_reports
import dash_bootstrap_components as dbc

# Initialize the Dash app with a Bootstrap theme
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Initialize ETABS API
etabs_api_instance = etabs_api.etabs_api()

# Global variable to store designed piers DataFrame
designed_piers_df = None

# Define the layout of the app
app.layout = dbc.Container([
    dbc.Row([
        dbc.Col(html.H1("ETABS Piers Design Tool", className="text-center my-4"), width=12)
    ]),
    
    # Two-column layout
    dbc.Row([
        # Left column for input controls
        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    html.H4("Input Controls", className="card-title"),
                    html.Button("Connect to ETABS", id="connect-button", className="btn btn-primary mb-3"),
                    
                    # Dropdowns and controls (initially hidden)
                    html.Div(id="controls", style={"display": "none"}, children=[
                        dbc.Label("Unfactored EQ Envelope:"),
                        dcc.Dropdown(id="eq-env-1", className="mb-3"),
                        
                        dbc.Label("Factored EQ Envelope:"),
                        dcc.Dropdown(id="eq-env-2", className="mb-3"),
                        
                        dbc.Label("Wind Envelope:"),
                        dcc.Dropdown(id="wind-env", className="mb-3"),
                        
                        dbc.Label("Vertical Spacing (mm):"),
                        dcc.Dropdown(id="vertical-spacing", options=[150, 200, 250, 300], value=200, className="mb-3"),
                        
                        dbc.Label("Horizontal Spacing (mm):"),
                        dcc.Dropdown(id="horizontal-spacing", options=[150, 200, 250, 300], value=200, className="mb-3"),
                        
                        dbc.Label("Design Out of Plane Bending:"),
                        dbc.Checklist(id="out-of-plane", options=[{"label": "Yes", "value": True}], value=[], className="mb-3"),
                        
                        html.Button("Design Piers", id="design-button", className="btn btn-success mb-3"),
                    ]),
                ]),
            ]),
        ], width=4),  # Left column width
        
        # Right column for output
        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    html.Div(id="output-content", children=[
                        dcc.Markdown('''
                            **ETABS DESIGN PIERS AS COLUMNS**

                            **Introduction:**
                            - Tool is to design walls as columns to AS3600:2018 Section 10 for ULS load combinations
                            - ETABS model should be correct and run without errors
                            - Pier labels should be assigned to walls manually
                            - Piers assigned should be straight walls only, no T or L shaped walls are permitted
                            - Two walls meeting at a corner or junction should be assigned two separate pier labels
                            - Program will design the minimum bar size for the bar spacing selected in each direction
                            - This may not necessarily comply with code minimum requirements. This is to be checked manually
                            - X-X is defined as strong direction and Y-Y is to be the weak direction

                            **Load combinations:**
                            - ETABS load combinations must be set up as follows:
                            - eq_env_1 -> Envelope of load combinations with unfactored earthquake load for flexural design
                            - eq_env_2 -> Envelope of load combinations with factored earthquake load for shear design
                            - wind_env -> Envelope of all ULS wind combinations
                            - If designing non-ductile walls for earthquake actions, eq_env_1 and eq_env_2 will be the same
                        ''', style={"margin": "20px"}),
                    ]),
                ]),
            ]),
        ], width=8),  # Right column width
    ]),
    
    # Status label
    dbc.Row([
        dbc.Col(html.Div(id="status-label", className="my-3"), width=12)
    ]),
    
    # Export button (initially hidden)
    dbc.Row([
        dbc.Col(html.Button("Export to Excel", id="export-button", className="btn btn-warning", style={"display": "none"}), width=12)
    ]),
    
    # Table to display designed piers
    dbc.Row([
        dbc.Col(dash_table.DataTable(id="designed-piers-table", style_table={"margin": "20px"}), width=12)
    ]),
], fluid=True)

# Callback to handle the "Connect to ETABS" button
@app.callback(
    [Output("controls", "style"),
     Output("eq-env-1", "options"),
     Output("eq-env-2", "options"),
     Output("wind-env", "options")],
    Input("connect-button", "n_clicks")
)
def on_connect_button_clicked(n_clicks):
    if n_clicks == 0:
        return {"display": "none"}, [], [], []
    
    try:
        load_case_list = etabs_api_instance.get_load_cases_and_combinations()
        options = [{"label": case, "value": case} for case in load_case_list]
        return {"display": "block"}, options, options, options
    except Exception as e:
        return {"display": "none"}, [], [], []

# Callback to handle the "Design Piers" button
@app.callback(
    [Output("status-label", "children"),
     Output("export-button", "style"),
     Output("designed-piers-table", "data"),
     Output("designed-piers-table", "columns"),
     Output("output-content", "children")],
    [Input("design-button", "n_clicks")],
    [State("eq-env-1", "value"),
     State("eq-env-2", "value"),
     State("wind-env", "value"),
     State("vertical-spacing", "value"),
     State("horizontal-spacing", "value"),
     State("out-of-plane", "value")]
)
def on_design_button_clicked(n_clicks, eq_env_1, eq_env_2, wind_env, vertical_spacing, horizontal_spacing, out_of_plane):
    global designed_piers_df
    
    if n_clicks == 0:
        return "", {"display": "none"}, [], [], dash.no_update
    
    if not all([eq_env_1, eq_env_2, wind_env]):
        return "Error: Please select all load cases.", {"display": "none"}, [], [], dash.no_update
    
    try:
        piers = etabs_api_instance.get_piers(load_cases=[eq_env_1, eq_env_2, wind_env])
    except Exception as e:
        return f"Error fetching pier data: {e}", {"display": "none"}, [], [], dash.no_update
    
    try:
        designed_piers_df = as3600_column_design.design_all_piers(
            piers=piers,
            eq_env_1=eq_env_1,
            eq_env_2=eq_env_2,
            wind_env=wind_env,
            vertical_spacing=vertical_spacing,
            horizontal_spacing=horizontal_spacing
        )
        
        # Convert DataFrame to dictionary for Dash table
        data = designed_piers_df.to_dict("records")
        columns = [{"name": col, "id": col} for col in designed_piers_df.columns]
        
        return "Design complete", {"display": "block"}, data, columns, ""
    except Exception as e:
        return f"Error in design function: {e}", {"display": "none"}, [], [], dash.no_update

# Callback to handle the "Export to Excel" button
@app.callback(
    Output("status-label", "children", allow_duplicate=True),
    Input("export-button", "n_clicks"),
    prevent_initial_call=True
)
def on_export_button_clicked(n_clicks):
    global designed_piers_df
    
    if designed_piers_df is None:
        return "Error: No data to export. Please run the design first."
    
    try:
        design_reports.dataframe_to_xlsx(designed_piers_df)
        return "Export successful!"
    except Exception as e:
        return f"Error during export: {e}"

# Run the app
if __name__ == "__main__":
    app.run_server(debug=True)