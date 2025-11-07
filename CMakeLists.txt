# Die minimal erforderliche CMake-Version.
cmake_minimum_required(VERSION 3.16)

# Definiert den Projektnamen und die Sprache (C++).
project(PathFinder CXX)

# Setzt den C++-Standard auf C++11 (oder höher, falls gewünscht, z.B. 17).
set(CMAKE_CXX_STANDARD 11)

# **WICHTIG:** Sucht nach der CGAL-Bibliothek.
# Wenn CGAL auf Ihrem System korrekt installiert ist (einschließlich Boost, GMP/MPFR), 
# sollte dieser Befehl erfolgreich sein.
find_package(CGAL REQUIRED)

if (CGAL_FOUND)
    # Bindet die notwendigen CGAL-Einstellungen in das Projekt ein.
    include(${CGAL_USE_FILE})

    # Definiert die ausführbare Datei 'path_finder' aus der Quelldatei.
    add_executable(path_finder path_finder.cpp)

    # Verknüpft die CGAL-Bibliotheken mit der ausführbaren Datei.
    # Dies ist der Ersatz für die komplexen <CGAL_FLAGS> im g++-Befehl.
    target_link_libraries(path_finder PUBLIC ${CGAL_LIBRARIES})
else (CGAL_FOUND)
    message(FATAL_ERROR "CGAL wurde nicht gefunden! Bitte installieren Sie CGAL und seine Abhängigkeiten korrekt.")
endif (CGAL_FOUND)