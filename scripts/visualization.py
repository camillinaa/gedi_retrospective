# import libraries
from dash import Dash, dcc, html, Input, Output
import pandas as pd
import plotly.express as px


# initialize app
app = Dash(__name__)


# define the layout
app.layout = html.Div([
    html.H4('Interactive CNV annotation plot',
            style={'font-family': 'Arial, sans-serif', 'font-size': '24px', 'font-weight': 'bold'}),
    
    dcc.Graph(id="scatter-plot"),
    
    html.P("Filter by CNV overlap with hotspot:", style={'font-family': 'Arial, sans-serif'}),
    dcc.RangeSlider(
        id='range-slider',
        min=0, max=100, step=10,
        marks={0: {'label': '0', 'style': {'font-family': 'Arial', 'font-size': '14px'}},
               100: {'label': '100', 'style': {'font-family': 'Arial', 'font-size': '14px'}}},
        value=[0, 100]
    ),
    
    html.P("Select Y-axis:", style={'font-family': 'Arial, sans-serif'}),
    dcc.Dropdown(
        id="y-axis-dropdown",
        options=[
            {"label": "% High Signal Region", "value": "% High Signal Region"},
            {"label": "% Low Mappability Region", "value": "% Low Mappability Region"},
            {"label": "Gap Fraction", "value": "Gap Fraction"}
        ],
        value="% High Signal Region", 
        clearable=False
    ),
])

@app.callback(
    Output("scatter-plot", "figure"), 
    [Input("range-slider", "value"),
     Input("y-axis-dropdown", "value")]
)
def update_scatter_plot(slider_range, selected_y):
    try:
        df = pd.read_excel("results/genome_annotations.xlsx")

        if selected_y not in df.columns:
            return px.scatter(title="Error: Selected Y-axis not found in data")

        low, high = slider_range
        mask = (df['total_overlap'] > low) & (df['total_overlap'] < high)

        fig = px.scatter(
            df[mask], x="total_overlap", y=selected_y, 
            hover_data=['total_overlap'],
            labels={"total_overlap": "CNV Overlap", selected_y: selected_y},
            title=f"Scatter Plot: CNV Overlap vs {selected_y}"
        )

        return fig
    except Exception as e:
        print(f"Error in callback: {e}")
        return px.scatter(title="Error loading data")


# run app
if __name__ == "__main__":
    app.run_server(debug=True)
