#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <algorithm>

#include <CGAL/ch_graham_andrew.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

// Define types for Delaunay (CGAL library)
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Delaunay_triangulation_2<Kernel> Delaunay;
typedef Delaunay::Point Point;

// Structure for ConeData with x,y coords and type (left or right)
struct ConeData {
    double x;
    double y;
    int type; // 0 = Left (Blue), 1 = Right (Yellow)
};

// Structure to safe type of point (extra info for CGAL point)
struct ConePointInfo {
    int type;
};

struct Vector_2 {
    double x;
    double y;
};

double find_min_squared_distance(const CGAL::Point_2<Kernel>& target,
    const std::vector<CGAL::Point_2<Kernel>>& cone_set) {

    double min_dist_sq = std::numeric_limits<double>::max();

    for (const auto& cone : cone_set) {
        double dist_sq = CGAL::squared_distance(target, cone);
        if (dist_sq < min_dist_sq) {
            min_dist_sq = dist_sq;
        }
    }
    return min_dist_sq;
}

// Calculate dot-product of two points
double dot_product(const Vector_2& v1, const Vector_2& v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

// Squared length of vector
double squared_length(const Vector_2& v) {
    return v.x * v.x + v.y * v.y;
}

// dot-product based from coords of points
double dot_product_coords(double x1, double y1, double x2, double y2) {
    return x1 * x2 + y1 * y2;
}

double squared_length_coords(double x, double y) {
    return x * x + y * y;
}

// Validate whether midpoint is inside hull/constraint based on orientation.
bool is_midpoint_inside_hull(const std::vector<CGAL::Point_2<Kernel>>& hull_points,
    const CGAL::Point_2<Kernel>& p)
{
    if (hull_points.size() < 3) return false;

    for (size_t i = 0; i < hull_points.size(); ++i) {
        const CGAL::Point_2<Kernel>& p_i = hull_points[i];
        const CGAL::Point_2<Kernel>& p_next = hull_points[(i + 1) % hull_points.size()];

        CGAL::Orientation orientation = CGAL::orientation(p_i, p_next, p);

        if (orientation == CGAL::RIGHT) {
            return false;
        }
    }
    return true;
}

// Smooths given path
void smooth_path(std::vector<Point>& path, int window_size) {
    if (path.size() < 3 || window_size < 3) return;

    // Make sure window size is odd
    if (window_size % 2 == 0) {
        window_size++;
    }

    int half_window = window_size / 2;
    std::vector<Point> smoothed_path = path;

    // Loop begins after first 'half_window' and ends before last 'half_window'
    for (size_t i = half_window; i < path.size() - half_window; ++i) {
        double sum_x = 0.0;
        double sum_y = 0.0;

        for (int j = -half_window; j <= half_window; ++j) {
            sum_x += path[i + j].x();
            sum_y += path[i + j].y();
        }

        smoothed_path[i] = Point(sum_x / window_size, sum_y / window_size);
    }

    path = smoothed_path;
}

// Validate candidate for next point (based on turn degree)
bool is_acceptable_turn(const Point& p_last, const Point& p_current, const Point& p_candidate) {

    Vector_2 v_current = { p_current.x() - p_last.x(), p_current.y() - p_last.y() };

    Vector_2 v_next = { p_candidate.x() - p_current.x(), p_candidate.y() - p_current.y() };

    double dot_prod = dot_product(v_current, v_next);

    // Magic number for maximum angle
    const double COS_MAX_ANGLE_THRESHOLD = -0.50;

    // Normalised check: V1 * V2 / (|V1| * |V2|) > Maximum (use of squares as approx to skip root calculation)
    double length_sq_current = squared_length(v_current);
    double length_sq_next = squared_length(v_next);

    // Avoid division by 0
    if (length_sq_current < 1e-6 || length_sq_next < 1e-6) {
        return false;
    }

    double cos_angle = dot_prod / std::sqrt(length_sq_current * length_sq_next);

    // If result is too small the angle is too sharp to be accepted.
    return cos_angle > COS_MAX_ANGLE_THRESHOLD;
}

int find_start_index(const std::vector<Point>& path_points) {
    if (path_points.empty()) return -1;

    int start_index = 0;
    double min_x = path_points[0].x();

    for (size_t i = 1; i < path_points.size(); ++i) {
        if (path_points[i].x() < min_x) {
            min_x = path_points[i].x();
            start_index = i;
        }
    }
    return start_index;
}

// Reads data from file (needs to be changed later to live in-out system) and calculates the middle path
void find_middle_path(const std::string& input_filepath, const std::string& output_filepath) {
    std::cout << "Begin pathfinding algorithm ...\n";
    std::vector<ConeData> cones;

    // Read cone positions as input (exchanged later for live-feed from Ros Node)
    std::ifstream inputFile(input_filepath);
    if (!inputFile.is_open()) {
        std::cerr << "ERROR: Input_file " << input_filepath << " not found.\n";
        return;
    }

    std::string line;
    int line_number = 1;
    if (std::getline(inputFile, line)) {
        line_number++;
    }
    else {
        std::cerr << "ERROR: cones-file empty.\n";
        return;
    }
    std::vector<std::pair<Point, ConePointInfo>> points_with_info;

    // Helping structure to get point (left, right) of point
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
            // Columns in cones.csv are: x_cone, y_cone, type
            if (std::getline(ss, segment, ',')) { cd.x = std::stod(segment); value_count++; }
            if (std::getline(ss, segment, ',')) { cd.y = std::stod(segment); value_count++; }
            if (std::getline(ss, segment, ',')) { cd.type = std::stoi(segment); value_count++; }

            // Add point if read correctly
            if (value_count == 3) {
                cones.push_back(cd);
            }
            else {
                std::cerr << "WARNING PathFinder: Only " << value_count << " of 3 columns in line " << line_number << " found.\n";
            }
        }
        catch (const std::exception& e) {
            std::cerr << "ERROR PathFinder: Conversion in line " << line_number << " failed: " << e.what() << "\n";
        }
        Point p(cd.x, cd.y);
        point_to_type[p] = cd.type;
    }
    inputFile.close();


    // Triangulation only safes coordinate data, not cone types
    // Use map to get cone type after triangulation
    std::vector<Point> points_for_dt;
    for (const auto& pair : point_to_type) {
        points_for_dt.push_back(pair.first);
    }

    std::vector<CGAL::Point_2<Kernel>> all_blue_cones;
    std::vector<CGAL::Point_2<Kernel>> all_yellow_cones;

    for (const auto& pair : point_to_type) {
        // pair.first is point (coords)
        // pair.second is type (0 for blue/left, 1 for yellow/right)

        if (pair.second == 0) { 
            all_blue_cones.push_back(pair.first);
        }
        else if (pair.second == 1) {
            all_yellow_cones.push_back(pair.first);
        }
    }

    std::set<CGAL::Point_2<Kernel>> raw_middle_points_set;

    // Delaunay-triangulation
    Delaunay dt;
    dt.insert(points_for_dt.begin(), points_for_dt.end());
    std::cout << "Delaunay-Triangulation finished. Triangles: " << dt.number_of_faces() << "\n";

	// Export edges for debugging
    double max_edge_length = 15.0; // Schwellenwert in Metern (anpassen!)
    double max_sq_dist = max_edge_length * max_edge_length;

    std::ofstream edge_file("filtered_edges.csv");
    edge_file << "x1,y1,x2,y2\n";

    for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
        // Hol dir die Segmente (Verbindungslinien) der Kante
        auto seg = dt.segment(eit);

        // CGAL nutzt oft quadrierte Distanzen (schneller, da kein Wurzelziehen)
        if (CGAL::squared_distance(seg.source(), seg.target()) <= max_sq_dist) {
            edge_file << seg.source().x() << "," << seg.source().y() << ","
                << seg.target().x() << "," << seg.target().y() << "\n";
        }
    }
    edge_file.close();

    // Magic numbers for configuration
    const double MAX_DISTANCE_RATIO_ERROR = 0.12;
    const double MAX_TRACK_WIDTH = 30.0; // Adjustable

    // Iterate over edges of triangulation
    for (auto e_it = dt.finite_edges_begin(); e_it != dt.finite_edges_end(); ++e_it) {
        Delaunay::Vertex_handle v1 = e_it->first->vertex(dt.cw(e_it->second));
        Delaunay::Vertex_handle v2 = e_it->first->vertex(e_it->second);

        Point p1 = v1->point();
        Point p2 = v2->point();
        
        double max_dist_sq = MAX_TRACK_WIDTH * MAX_TRACK_WIDTH;

        std::vector<CGAL::Point_2<Kernel>> convex_hull_points;

           
        // Check if distance is acceptable (avoid big jumps between points)
        if (CGAL::squared_distance(p1, p2) >= max_dist_sq) {
            continue;
        }

        int type1 = point_to_type[p1];
        int type2 = point_to_type[p2];

        // Only looking for vertex connecting left and right (blue and yellow) cone
        if ((type1 == 0 && type2 == 1) || (type1 == 1 && type2 == 0)) {
            double mid_x = (p1.x() + p2.x()) / 2.0;
            double mid_y = (p1.y() + p2.y()) / 2.0;
            Point midpoint(mid_x, mid_y);

            double dist_sq_to_blue = find_min_squared_distance(midpoint, all_blue_cones);

            double dist_sq_to_yellow = find_min_squared_distance(midpoint, all_yellow_cones);

            double average_sq = (dist_sq_to_blue + dist_sq_to_yellow) / 2.0;

            double diff_sq = std::abs(dist_sq_to_blue - dist_sq_to_yellow);

            // Check if distance between left and right cone is small enough
            if (average_sq > 0) {
                double error_ratio = diff_sq / average_sq;

                if (error_ratio <= MAX_DISTANCE_RATIO_ERROR  && dist_sq_to_blue < max_dist_sq && dist_sq_to_yellow < max_dist_sq) {
                    raw_middle_points_set.insert(midpoint);
                }
            }

        }
    }

    std::vector<CGAL::Point_2<Kernel>> raw_middle_points(raw_middle_points_set.begin(), raw_middle_points_set.end());
    std::cout << "Potential middle points found: " << raw_middle_points.size() << "\n";

    // Determine order of points (potentially has to be adjusted for live-feed)
    if (raw_middle_points.empty()) {
        std::cerr << "ERROR: No path points generated.\n";
        return;
    }

    std::vector<Point> sorted_path;
    std::vector<bool> visited(raw_middle_points.size(), false);

    int current_index = find_start_index(raw_middle_points);

    sorted_path.push_back(raw_middle_points[current_index]);
    visited[current_index] = true;

    const double LOCAL_SEARCH_RADIUS_SQ = 30.0 * 30.0;

    int error_rate = 0;
    while (sorted_path.size() < raw_middle_points.size()) {

        int nearest_index = -1;
        double min_dist_sq = std::numeric_limits<double>::max();

        const Point& p_current = sorted_path.back();

        // penultimate point
        const Point& p_last = (sorted_path.size() > 1) ? sorted_path[sorted_path.size() - 2] : p_current;

        for (size_t j = 0; j < raw_middle_points.size(); ++j) {
            if (!visited[j]) {

                double dist_sq = CGAL::squared_distance(p_current, raw_middle_points[j]);

                // Filter by local distance
                if (dist_sq > LOCAL_SEARCH_RADIUS_SQ) {
                    continue;
                }

                double dx_current = p_current.x() - p_last.x();
                double dy_current = p_current.y() - p_last.y();
                double dx_next = raw_middle_points[j].x() - p_current.x();
                double dy_next = raw_middle_points[j].y() - p_current.y();

                double length_sq_current = squared_length_coords(dx_current, dy_current);
                double length_sq_next = squared_length_coords(dx_next, dy_next);

                double cos_angle = 1.0;
                if (sorted_path.size() > 1 && length_sq_current > 1e-6 && length_sq_next > 1e-6) {
                    double dot_prod = dot_product_coords(dx_current, dy_current, dx_next, dy_next);
                    cos_angle = dot_prod / std::sqrt(length_sq_current * length_sq_next);
                }

                // Penalty for sharp angles based on angle's cosine
                double angle_penalty = 1.0 - cos_angle;

                // Magic number used to weight angle penalty
                const double weight_factor = 80.0;

                // Calculate distance weighted by angle's sharpness
                double weighted_dist_sq = dist_sq * (1.0 + weight_factor * angle_penalty);

                // Choice for best neighbouring edge
                if (weighted_dist_sq < min_dist_sq) {
                    min_dist_sq = weighted_dist_sq;
                    nearest_index = j;
                }
            }
        }
        
        // Backup behaviour if no fitting neighbouring edge was found
        if (nearest_index == -1) {
            error_rate = error_rate + 1;
            for (size_t j = 0; j < raw_middle_points.size(); ++j) {
                if (!visited[j]) {
                    double dist_sq = CGAL::squared_distance(p_current, raw_middle_points[j]);
                    if (dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                        nearest_index = static_cast<int>(j);
                    }
                }
            }
        }

        if (nearest_index != -1) {
            sorted_path.push_back(raw_middle_points[nearest_index]);
            visited[nearest_index] = true;
        }
        else {
            std::cerr << "WARNING: No valid next neighbour found. Exit sort with " << sorted_path.size() << " points.\n";
            break;
        }
    }

    if (!sorted_path.empty()) {
        const Point& last_p = sorted_path.back();
        std::cerr << "DEBUG: Last found point at X=" << last_p.x()
            << ", Y=" << last_p.y() << ". Sorting exited.\n";
    }

    if (sorted_path.empty()) {
        std::cerr << "ERROR: Sorting path empty. No file written.\n";
        return;
    }

    std::ofstream outputFile(output_filepath);

    if (!outputFile.is_open()) {
        std::cerr << "ERROR: Couldn't open output file (" << output_filepath << ") to write.\n";
        return;
    }

    const int SMOOTHING_WINDOW = 5;

    smooth_path(sorted_path, SMOOTHING_WINDOW);
    std::cout << "INFO: Path smoothend with window " << SMOOTHING_WINDOW << ".\n";

    // Header-Zeile schreiben
    outputFile << "x_path,y_path\n";

    // Schreibe alle Punkte des sortierten Pfades
    for (const auto& p : sorted_path) {
        // Da 'p' ein CGAL::Point_2 ist, verwenden wir .x() und .y()
        outputFile << p.x() << "," << p.y() << "\n";
    }

    outputFile.close();
    std::cout << "INFO: path successfully written to " << output_filepath << ". Total points: " << sorted_path.size() << " Error: " << error_rate << "\n";
}


int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: ./path_finder <input_cones_file> <output_path_file>\n";
        return 1;
    }
    find_middle_path(argv[1], argv[2]);
    return 0;
}
