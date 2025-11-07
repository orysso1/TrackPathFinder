#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <locale>

using namespace std;

// Struktur für die eingelesenen Track-Daten
struct TrackPoint {
    double x_m;
    double y_m;
    double w_tr_right_m;
    double w_tr_left_m;
};

// Struktur für die berechneten Cone-Positionen
struct ConeData {
    double x;
    double y;
    int type; // 0 = Links (Blau), 1 = Rechts (Gelb)
};

void convert_to_cones(const std::string& input_filepath, const std::string& output_filepath) {
    std::cout << "Starte Cone-Extraktion von " << input_filepath << "...\n";

    std::ifstream inputFile(input_filepath);
    if (!inputFile.is_open()) {
        std::cerr << "FEHLER: Eingabedatei " << input_filepath << " nicht gefunden.\n";
        return;
    }

    std::string line;
    std::vector<TrackPoint> track_points;
    int line_number = 0;

    // 1. Header überspringen
    if (std::getline(inputFile, line)) {
        line_number++;
        std::cout << "Header-Zeile (" << line_number << ") uebersprungen: " << line << "\n";
    }
    else {
        std::cerr << "FEHLER: Eingabedatei ist leer.\n";
        return;
    }

    // 2. Iteriere über die Datenzeilen und fange Fehler ab
    while (std::getline(inputFile, line)) {
        line_number++;

        // --- Robuste Zeilenprüfung ---

        // Entferne Carriage Returns und teste auf leere Zeilen
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        if (line.empty() || line.find_first_not_of(" \t\n") == std::string::npos) {
            continue;
        }

        // --- Konvertierung und Fehlerbehandlung ---

        std::stringstream ss(line);
        std::string segment;
        TrackPoint tp;
        int value_count = 0;

        try {
            // Lesen und Konvertieren der 4 Werte (Komma als Trennzeichen!)
            if (std::getline(ss, segment, ',')) { tp.x_m = std::stod(segment); value_count++; }
            if (std::getline(ss, segment, ',')) { tp.y_m = std::stod(segment); value_count++; }
            if (std::getline(ss, segment, ',')) { tp.w_tr_right_m = std::stod(segment); value_count++; }
            if (std::getline(ss, segment, ',')) { tp.w_tr_left_m = std::stod(segment); value_count++; }

            if (value_count == 4) {
                track_points.push_back(tp);
            }
            else {
                std::cerr << "WARNUNG in Zeile " << line_number << ": Nur " << value_count << " von 4 Spalten gefunden. Zeile ignoriert: " << line << "\n";
            }
        }
        catch (const std::invalid_argument& e) {
            std::cerr << "FEHLER in Zeile " << line_number << ": Konvertierung fehlgeschlagen (Pruefen Sie Dezimaltrennzeichen: Punkt vs. Komma). Zeile ignoriert: " << line << "\n";
        }
        catch (const std::out_of_range& e) {
            std::cerr << "FEHLER in Zeile " << line_number << ": Wert außerhalb des Bereichs. Zeile ignoriert.\n";
        }
    }
    inputFile.close();

    // --- Prüfung der Anzahl der Punkte ---

    if (track_points.size() < 2) {
        std::cerr << "FEHLER: Nicht genuegend gueltige Track-Punkte (" << track_points.size() << ") zum Berechnen der Cone-Positionen gefunden.\n";
        return;
    }

    std::cout << "Erfolgreich " << track_points.size() << " Track-Punkte eingelesen.\n";

    // --- Cone-Berechnungs-Logik (unveraendert) ---

    std::vector<ConeData> cones;

    for (size_t i = 0; i < track_points.size() - 1; ++i) {
        const auto& p1 = track_points[i];
        const auto& p2 = track_points[i + 1];

        // 1. Vektor des Segments
        double dx = p2.x_m - p1.x_m;
        double dy = p2.y_m - p1.y_m;

        // 2. Normalenvektor (90 Grad Drehung)
        // Rechts: (dy, -dx)
        // Links: (-dy, dx)
        double norm_x = dy;
        double norm_y = -dx;

        // 3. Normalisierung des Vektors
        double length = std::sqrt(norm_x * norm_x + norm_y * norm_y);
        if (length == 0.0) continue;

        double unit_norm_x = norm_x / length;
        double unit_norm_y = norm_y / length;

        // 4. Mittlere Breiten und Positionen berechnen (Mittelpunkt des Segments)
        double mid_x = (p1.x_m + p2.x_m) / 2.0;
        double mid_y = (p1.y_m + p2.y_m) / 2.0;
        double avg_w_right = (p1.w_tr_right_m + p2.w_tr_right_m) / 2.0;
        double avg_w_left = (p1.w_tr_left_m + p2.w_tr_left_m) / 2.0;

        // Cone Links (Blau)
        // Normalenvektor umgedreht
        cones.push_back({
            mid_x + unit_norm_x * avg_w_left, // x + (normal_x * width_left)
            mid_y + unit_norm_y * avg_w_left, // y + (normal_y * width_left)
            0 // Typ 0: Blau (Links)
            });

        // Cone Rechts (Gelb)
        // Normalenvektor * -1 (Richtung wechselt)
        cones.push_back({
            mid_x - unit_norm_x * avg_w_right, // x - (normal_x * width_right)
            mid_y - unit_norm_y * avg_w_right, // y - (normal_y * width_right)
            1 // Typ 1: Gelb (Rechts)
            });
    }

    // --- Ausgabe der Cone-Daten ---

    std::ofstream outputFile(output_filepath);
    outputFile << "x_cone,y_cone,type\n";
    for (const auto& c : cones) {
        outputFile << c.x << "," << c.y << "," << c.type << "\n";
    }
    outputFile.close();

    std::cout << "Erfolgreich " << cones.size() << " Cones nach " << output_filepath << " geschrieben.\n";
}

// Der main-Aufruf würde so aussehen:
int main(int argc, char* argv[]) {

    // NEU: Erzwinge die klassische C-Locale, um sicherzustellen, dass std::stod den PUNKT (.) als Dezimaltrennzeichen verwendet.
    // Dies löst Probleme mit Komma-Dezimaltrennzeichen in europäischen Regionen.
    std::locale::global(std::locale::classic());

    if (argc != 3) {
        std::cerr << "Verwendung: ./ConeExtractor.exe <input_tum_file> <output_cones_file>\n";
        return 1;
    }

    convert_to_cones(argv[1], argv[2]);

    return 0;
}
 