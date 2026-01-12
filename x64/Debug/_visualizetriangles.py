import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd

# Daten laden
tri_data = pd.read_csv('triangles.csv')

fig, ax = plt.subplots()

# Über jedes Dreieck in der CSV iterieren
for index, row in tri_data.iterrows():
    # Punkte des Dreiecks definieren
    polygon = [
        (row['x1'], row['y1']),
        (row['x2'], row['y2']),
        (row['x3'], row['y3'])
    ]
    
    # Dreieck als Patch hinzufügen
    # facecolor='none' lässt es transparent, edgecolor='blue' zeichnet die Linien
    tri = patches.Polygon(polygon, closed=True, linewidth=0.5, edgecolor='blue', facecolor='lightblue', alpha=0.3)
    ax.add_patch(tri)

# Achsen automatisch anpassen
ax.autoscale()
plt.gca().set_aspect('equal', adjustable='box')
plt.show()