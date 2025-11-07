#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <CGAL/ch_graham_andrew.h>
#include <algorithm>

// --- CGAL ERFORDERLICHE HEADER ---
// Hinweis: Die genauen Header hängen von der verwendeten CGAL-Konfiguration ab.
// Diese sind Beispiele für eine einfache 2D-Triangulation.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

// Definiere den CGAL Kernel und die Triangulation
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_2<K> Triangulation;
typedef K::Point_2 Point;

// Struktur für die eingelesenen Kegeldaten
struct ConeData {
    double x;
    double y;
    int type; // 0 = Links (Blau), 1 = Rechts (Gelb)
};

// Struktur, die den Punkt und seinen Typ speichert (wird als Info an CGAL-Punkt gehängt)
struct ConePointInfo {
    int type;
};

struct Vector_2 {
    double x;
    double y;
};

// Berechnet das Skalarprodukt (Dot Product)
double dot_product(const Vector_2& v1, const Vector_2& v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

// Berechnet die quadratische Länge eines Vektors
double squared_length(const Vector_2& v) {
    return v.x * v.x + v.y * v.y;
}

double dot_product_coords(double x1, double y1, double x2, double y2) {
    return x1 * x2 + y1 * y2;
}

double squared_length_coords(double x, double y) {
    return x * x + y * y;
}

bool is_midpoint_inside_hull(const std::vector<CGAL::Point_2<K>>& hull_points,
    const CGAL::Point_2<K>& p)
{
    if (hull_points.size() < 3) return false;

    // Wir gehen davon aus, dass die Punkte gegen den Uhrzeigersinn (CCW) geordnet sind.
    // Für jeden Rand des Polygons muss der Testpunkt p "links" oder "kolinear" sein.
    // Wenn CGAL die Hülle im Uhrzeigersinn (CW) liefert, müssten wir nach "rechts" prüfen.
    // Standardmäßig liefert ch_graham_andrew CCW.

    for (size_t i = 0; i < hull_points.size(); ++i) {
        const CGAL::Point_2<K>& p_i = hull_points[i];
        const CGAL::Point_2<K>& p_next = hull_points[(i + 1) % hull_points.size()];

        // Prüfe die Orientierung: Orientierung von (p_i, p_next, p)
        CGAL::Orientation orientation = CGAL::orientation(p_i, p_next, p);

        // Wenn der Punkt rechts von der Kante liegt, ist er außerhalb der konvexen Hülle.
        if (orientation == CGAL::RIGHT) {
            return false;
        }
        // CGAL::LEFT und CGAL::COLLINEAR sind akzeptabel.
    }

    return true; // Der Punkt liegt innerhalb der Hülle oder auf dem Rand.
}

void smooth_path(std::vector<Point>& path, int window_size) {
    if (path.size() < 3 || window_size < 3) return;

    // Sicherstellen, dass die Fenstergröße ungerade ist
    if (window_size % 2 == 0) {
        window_size++;
    }

    int half_window = window_size / 2;
    std::vector<Point> smoothed_path = path; // Kopie des Originalpfads

    // Der Loop beginnt nach dem ersten 'half_window' und endet vor dem letzten 'half_window'
    for (size_t i = half_window; i < path.size() - half_window; ++i) {
        double sum_x = 0.0;
        double sum_y = 0.0;

        // Summiere die Koordinaten im Fenster
        for (int j = -half_window; j <= half_window; ++j) {
            sum_x += path[i + j].x();
            sum_y += path[i + j].y();
        }

        // Berechne den Durchschnitt und speichere ihn in der neuen Kopie
        smoothed_path[i] = Point(sum_x / window_size, sum_y / window_size);
    }

    // Ersetze den Originalpfad durch den geglätteten Pfad
    path = smoothed_path;
}

bool is_acceptable_turn(const Point& p_last, const Point& p_current, const Point& p_candidate) {

    // Vektor des letzten Segments (Pfadrichtung)
    Vector_2 v_current = { p_current.x() - p_last.x(), p_current.y() - p_last.y() };

    // Vektor zum Kandidatenpunkt
    Vector_2 v_next = { p_candidate.x() - p_current.x(), p_candidate.y() - p_current.y() };

    double dot_prod = dot_product(v_current, v_next);

    // Schwellenwert für den Winkel (z.B. cos(120 Grad) = -0.5)
    // Wenn das Skalarprodukt negativ wird, dreht sich der Pfad zurück (> 90 Grad).
    // Wenn es zu klein ist, ist es ein scharfer Knick.
    // Ein Wert von 0.25 erzwingt einen Winkel von weniger als ca. 75 Grad.
    const double COS_MAX_ANGLE_THRESHOLD = -0.50;

    // Normalisierte Prüfung: V1 * V2 / (|V1| * |V2|) > Schwellenwert
    // Wir verwenden die quadratischen Längen als Näherung, da die Wurzel teuer ist:
    double length_sq_current = squared_length(v_current);
    double length_sq_next = squared_length(v_next);

    // Vermeide Division durch Null für den unwahrscheinlichen Fall, dass die Vektoren Länge Null haben.
    if (length_sq_current < 1e-6 || length_sq_next < 1e-6) {
        return false;
    }

    // Wenn das normalisierte Skalarprodukt unter den Schwellenwert fällt, ist der Winkel zu scharf.
    // Skalarprodukt / (Länge * Länge) ist der Kosinus des Winkels.
    double cos_angle = dot_prod / std::sqrt(length_sq_current * length_sq_next);

    return cos_angle > COS_MAX_ANGLE_THRESHOLD;
}

int find_start_index(const std::vector<Point>& path_points) {
    if (path_points.empty()) return -1;

    int start_index = 0;
    // Verwenden Sie .x() statt .first
    double min_x = path_points[0].x();

    for (size_t i = 1; i < path_points.size(); ++i) {
        // Verwenden Sie .x() statt .first
        if (path_points[i].x() < min_x) {
            min_x = path_points[i].x();
            start_index = i;
        }
    }
    return start_index;
}

// --- Kernlogik für die Pfadbestimmung ---
void find_middle_path(const std::string& input_filepath, const std::string& output_filepath) {
    std::cout << "Starte Pfadfinder-Algorithmus...\n";
    std::vector<ConeData> cones;

    // 1. Daten einlesen (Cone-Positionen)
    std::ifstream inputFile(input_filepath);
    if (!inputFile.is_open()) {
        std::cerr << "Fehler: Eingabedatei " << input_filepath << " nicht gefunden.\n";
        return;
    }

    std::string line;
    int line_number = 1;
    if (std::getline(inputFile, line)) {
        line_number++;
    }
    else {
        std::cerr << "FEHLER: Cones-Datei ist leer.\n";
        return;
    }
    std::vector<std::pair<Point, ConePointInfo>> points_with_info;

    // Hilfsstruktur zur Speicherung der Cone-Infos
    std::map<Point, int> point_to_type;

    while (std::getline(inputFile, line)) {
        line_number++;

        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty() || line.find_first_not_of(" \t\n") == std::string::npos) {
            continue;
        }

        std::stringstream ss(line);
        std::string segment;
        ConeData cd;
        int value_count = 0;

        try {
            // Die Spalten in cones.csv sind: x_cone, y_cone, type
            if (std::getline(ss, segment, ',')) { cd.x = std::stod(segment); value_count++; }
            if (std::getline(ss, segment, ',')) { cd.y = std::stod(segment); value_count++; }
            if (std::getline(ss, segment, ',')) { cd.type = std::stoi(segment); value_count++; } // Hier muss stoi verwendet werden

            if (value_count == 3) {
                cones.push_back(cd); // HIER wird 'cones' befüllt!
            }
            else {
                std::cerr << "WARNUNG PathFinder: Nur " << value_count << " von 3 Spalten in Zeile " << line_number << " gefunden.\n";
            }
        }
        catch (const std::exception& e) {
            // Fangt alle Konvertierungsfehler (stod oder stoi) ab
            std::cerr << "FEHLER PathFinder: Konvertierung in Zeile " << line_number << " fehlgeschlagen: " << e.what() << "\n";
        }
        Point p(cd.x, cd.y);
        point_to_type[p] = cd.type;
    }
    inputFile.close();



    // Die Triangulation kann nur Punkte, keine Zusatzinformationen, direkt speichern.
    // Wir verwenden die map, um den Typ eines Punktes nach der Triangulation abzufragen.
    std::vector<Point> points_for_dt;
    for (const auto& pair : point_to_type) {
        points_for_dt.push_back(pair.first);
    }

    // 2. Delaunay-Triangulation durchführen
    Triangulation dt;
    dt.insert(points_for_dt.begin(), points_for_dt.end());
    std::cout << "Delaunay-Triangulation abgeschlossen. Dreiecke: " << dt.number_of_faces() << "\n";

    // 3. und 4. Filtern und Mittelpunkte extrahieren (Kern-Heuristik)
    std::vector<Point> raw_middle_points;

    // Iteriere über die Kanten der Triangulation (Edges)
    for (auto e_it = dt.finite_edges_begin(); e_it != dt.finite_edges_end(); ++e_it) {
        // Holen der Eckpunkte der Kante
        Triangulation::Vertex_handle v1 = e_it->first->vertex(dt.cw(e_it->second));
        Triangulation::Vertex_handle v2 = e_it->first->vertex(e_it->second);

        Point p1 = v1->point();
        Point p2 = v2->point();

        const double MAX_TRACK_WIDTH = 80.0; // Passen Sie diesen Wert an!
        double max_dist_sq = MAX_TRACK_WIDTH * MAX_TRACK_WIDTH;

        std::vector<CGAL::Point_2<K>> convex_hull_points;

        CGAL::ch_graham_andrew(points_for_dt.begin(), points_for_dt.end(),
            std::back_inserter(convex_hull_points));

        // 1. Distanzprüfung
        // CGAL::squared_distance(p1, p2) berechnet (x1-x2)^2 + (y1-y2)^2
        if (CGAL::squared_distance(p1, p2) >= max_dist_sq) {
            continue; // Kante ist zu lang (wilder Sprung), ignoriere sie
        }

        // Suche den Typ der beiden Punkte in unserer Map (Dies ist die Heuristik!)
        int type1 = point_to_type[p1];
        int type2 = point_to_type[p2];

        // Wir suchen nur Kanten, die einen blauen (0) und einen gelben (1) Kegel verbinden.
        if ((type1 == 0 && type2 == 1) || (type1 == 1 && type2 == 0)) {
            // Berechne den Mittelpunkt der Kante
            double mid_x = (p1.x() + p2.x()) / 2.0;
            double mid_y = (p1.y() + p2.y()) / 2.0;
            Point midpoint(mid_x, mid_y);

            if (is_midpoint_inside_hull(convex_hull_points, midpoint))
            {
                raw_middle_points.emplace_back(midpoint);
            }
        }
    }
    std::cout << "Gefundene potentielle Mittelpunkte: " << raw_middle_points.size() << "\n";

    // 5. Pfad-Folge bestimmen (Glätten und Sortieren)
    // HIER WIRD EINE ZUSÄTZLICHE HEURISTIK BENÖTIGT, 
    // um die raw_middle_points in die korrekte Reihenfolge zu bringen (z.B. mittels Nearest-Neighbor-Search oder Graph-Traversal).
    // Der Einfachheit halber verwenden wir hier die Rohdaten.
    if (raw_middle_points.empty()) {
        std::cerr << "WARNUNG: Es wurden keine Pfadpunkte generiert.\n";
        return;
    }

    if (raw_middle_points.size() < 1) {
        std::cerr << "WARNUNG: Es wurden keine Pfadpunkte fuer die Sortierung generiert.\n"; return;
    }

    std::vector<Point> sorted_path;
    std::vector<bool> visited(raw_middle_points.size(), false);

    // Startpunkt finden
    int current_index = find_start_index(raw_middle_points);

    if (current_index == -1) {
        // Dies sollte nicht passieren, wenn path_points.size() >= 1, aber zur Sicherheit
        std::cerr << "FEHLER: Startpunktindex konnte nicht gefunden werden.\n";
        return;
    }

    sorted_path.push_back(raw_middle_points[current_index]); // <-- HIER wird der erste Punkt gesetzt
    visited[current_index] = true;

    // if (current_index != -1) {
    //    sorted_path.push_back(raw_middle_points[current_index]);
    //    visited[current_index] = true;
    //}

    // Iteratives Sortieren
    const double LOCAL_SEARCH_RADIUS_SQ = 20.0 * 20.0;


    // Die robuste WHILE-Schleife: läuft, bis alle Punkte sortiert sind
    while (sorted_path.size() < raw_middle_points.size()) {

        int nearest_index = -1;
        double min_dist_sq = std::numeric_limits<double>::max();

        const Point& p_current = sorted_path.back();

        // p_last ist der vorletzte Punkt, nur wenn es mindestens 2 Punkte in sorted_path gibt
        const Point& p_last = (sorted_path.size() > 1) ? sorted_path[sorted_path.size() - 2] : p_current;

        for (size_t j = 0; j < raw_middle_points.size(); ++j) {
            if (!visited[j]) {

                double dist_sq = CGAL::squared_distance(p_current, raw_middle_points[j]);

                // 1. Lokale Suche (schneller Filter)
                if (dist_sq > LOCAL_SEARCH_RADIUS_SQ) {
                    continue;
                }

                // 2. Richtungsprüfung (Ab dem zweiten Punkt)
                double dx_current = p_current.x() - p_last.x();
                double dy_current = p_current.y() - p_last.y();
                double dx_next = raw_middle_points[j].x() - p_current.x();
                double dy_next = raw_middle_points[j].y() - p_current.y();

                

                // Der Kosinus des Winkels: 1.0 = perfekt gerade, -1.0 = perfekt rückwärts
                double length_sq_current = squared_length_coords(dx_current, dy_current);
                double length_sq_next = squared_length_coords(dx_next, dy_next);

                double cos_angle = 1.0;
                if (sorted_path.size() > 1 && length_sq_current > 1e-6 && length_sq_next > 1e-6) {
                    double dot_prod = dot_product_coords(dx_current, dy_current, dx_next, dy_next);
                    cos_angle = dot_prod / std::sqrt(length_sq_current * length_sq_next);

                }

                // Straffaktor: Bestraft starke Knicke und belohnt gerade Linien
                // (1 - cos_angle) ist 0 bei gerader Linie und 2 bei Rückwärtsbewegung
                double angle_penalty = 1.0 - cos_angle;

                // Berechne gewichtete Distanz: Distanz * (1 + Strafe)
                // Wenn Winkel schlecht ist, wird die Distanz künstlich vergrößert
                double weighted_dist_sq = dist_sq * (1.0 + 50.0 * angle_penalty); // 5.0 ist Gewichtungsfaktor!

                // 3. Nearest Neighbor Auswahl
                if (weighted_dist_sq < min_dist_sq) {
                    min_dist_sq = weighted_dist_sq;
                    nearest_index = j;
                }
            }
        }

        

        // FÜGE PUNKT HINZU ODER BRECHE AB
        if (nearest_index != -1) {
            sorted_path.push_back(raw_middle_points[nearest_index]);
            visited[nearest_index] = true;
        }
        else {
            // Abbruch, wenn kein geeigneter Nachbar mehr gefunden wird (Ende des Pfades)
            std::cerr << "WARNUNG: Nächster Nachbar konnte nicht gefunden werden. Sortierung abgebrochen bei " << sorted_path.size() << " Punkten.\n";
            break;
        }
    }

    if (!sorted_path.empty()) {
        const Point& last_p = sorted_path.back();
        std::cerr << "DEBUG: Letzter gefundener Punkt bei X=" << last_p.x()
            << ", Y=" << last_p.y() << ". Sortierung abgebrochen.\n";
    }

    if (sorted_path.empty()) {
        std::cerr << "FEHLER: Der sortierte Pfad ist leer. Keine Datei geschrieben.\n";
        return;
    }

    std::ofstream outputFile(output_filepath);

    if (!outputFile.is_open()) {
        std::cerr << "FEHLER: Konnte Ausgabedatei (" << output_filepath << ") nicht zum Schreiben öffnen.\n";
        return;
    }

    const int SMOOTHING_WINDOW = 5;

    smooth_path(sorted_path, SMOOTHING_WINDOW);
    std::cout << "INFO: Pfad geglättet mit Fenstergröße " << SMOOTHING_WINDOW << ".\n";

    // Header-Zeile schreiben
    outputFile << "x_path,y_path\n";

    // Schreibe alle Punkte des sortierten Pfades
    for (const auto& p : sorted_path) {
        // Da 'p' ein CGAL::Point_2 ist, verwenden wir .x() und .y()
        outputFile << p.x() << "," << p.y() << "\n";
    }

    outputFile.close();
    std::cout << "INFO: Pfad erfolgreich nach " << output_filepath << " geschrieben. Gesamtpunkte: " << sorted_path.size() << "\n";
}


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Verwendung: ./path_finder <input_cones_file> <output_path_file>\n";
        return 1;
    }
    find_middle_path(argv[1], argv[2]);
    return 0;
}
