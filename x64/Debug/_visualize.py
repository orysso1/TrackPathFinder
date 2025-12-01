import sys
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os

# 1. Überprüfung der Kommandozeilen-Argumente
if len(sys.argv) != 3:
    print("Verwendung: python _visualize.py <cones_file.csv> <path_file.csv>")
    sys.exit(1)

cones_file = sys.argv[1]
path_file = sys.argv[2]

# 2. Daten einlesen mit Pandas
try:
    # Annahme: cones.csv hat die Spalten: x_cone, y_cone, type
    cones_df = pd.read_csv(cones_file)
    # Annahme: mittelpfad.csv hat die Spalten: x, y
    path_df = pd.read_csv(path_file)
except FileNotFoundError as e:
    print(f"FEHLER: Datei nicht gefunden: {e.filename}")
    sys.exit(1)
except pd.errors.EmptyDataError:
    print("FEHLER: Eine der CSV-Dateien ist leer.")
    sys.exit(1)
except Exception as e:
    print(f"FEHLER beim Einlesen der Daten: {e}")
    sys.exit(1)


# 3. Plotly Visualisierungs-Logik

# Farbzuordnung für die Kegel
color_map = {0: 'blue', 1: 'yellow'} 

# Erster Plot: Scatter Plot der Kegel (Cones)
# 'type' muss eine Zeichenkette sein, damit Plotly die Farben korrekt zuordnet
cones_df['type_str'] = cones_df['type'].apply(lambda x: 'Blau (Links)' if x == 0 else 'Gelb (Rechts)')

fig = px.scatter(
    cones_df, 
    x='x_cone', 
    y='y_cone', 
    color='type_str',         
    color_discrete_map={'Blau (Links)': 'blue', 'Gelb (Rechts)': 'yellow'}, 
    title='Track Visualisierung: Kegel & Berechneter Pfad',
    # Verbessert die Darstellung, wenn über Datenpunkte gehovert wird
    hover_data={'type': False, 'type_str': False, 'x_cone': True, 'y_cone': True}, 
)

# Zweiter Plot: Linie des berechneten Pfads (Path) hinzufügen
# Wir verwenden go.Scatter für eine saubere Linien-Darstellung
fig.add_trace(
    go.Scatter(
        x=path_df['x_path'], 
        y=path_df['y_path'], 
        mode='lines', 
        name='Berechneter Pfad',
        line=dict(color='red', width=2)
    )
)

# 4. Layout anpassen (WICHTIG für Tracks)
fig.update_layout(
    xaxis_title='X-Koordinate (m)',
    yaxis_title='Y-Koordinate (m)',
    # Setzt das Seitenverhältnis auf 1:1, damit der Track nicht verzerrt aussieht
    yaxis=dict(scaleanchor="x", scaleratio=1), 
    showlegend=True
)

# 5. Interaktiven Plot im Browser anzeigen
print("Erzeuge interaktiven Plot...")
output_html_file = os.path.join(os.path.dirname(__file__), 'track_visualisierung.html')
fig.write_html(output_html_file)

print(f"Interaktive Visualisierung wurde gespeichert unter: {output_html_file}")

print("Visualisierung abgeschlossen.")
