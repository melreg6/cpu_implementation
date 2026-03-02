#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>

struct Vec3 { double x, y, z; };

std::vector<Vec3> read_pdb_coords(const std::string& fname) {
    std::vector<Vec3> coords;
    std::ifstream fin(fname);
    if (!fin) throw std::runtime_error("Cannot open " + fname);
    std::string line;
    while (std::getline(fin, line)) {
        if (line.rfind("ATOM", 0) == 0) {
            double x = std::stod(line.substr(30, 8));
            double y = std::stod(line.substr(38, 8));
            double z = std::stod(line.substr(46, 8));
            coords.push_back({x, y, z});

            //std::istringstream iss(line);
            //std::string token;
            //double x, y, z;
            // skip 6 tokens then read xyz
            //for (int i=0; i<6; ++i) iss >> token;
            //if (iss >> x >> y >> z)
                //coords.push_back({x, y, z});
        }
    }
    return coords;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: compare_pdb file1 file2\n";
        return 1;
    }
    auto a = read_pdb_coords(argv[1]);
    auto b = read_pdb_coords(argv[2]);
    if (a.size() != b.size()) {
        std::cerr << "Atom count mismatch: " << a.size() << " vs " << b.size() << "\n";
        return 1;
    }

    double rmsd = 0, rel_sum = 0, percenterror = 0;
    for (size_t i=0; i<a.size(); ++i) {
        double dx=a[i].x-b[i].x, dy=a[i].y-b[i].y, dz=a[i].z-b[i].z;
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
        rmsd += dist*dist;
        double mag = std::sqrt(a[i].x*a[i].x + a[i].y*a[i].y + a[i].z*a[i].z);
        if (mag > 1e-8) rel_sum += dist/mag;

        double pos1 =std::sqrt( (a[i].x*a[i].x + a[i].y*a[i].y + a[i].z*a[i].z) );
        double pos2 =std::sqrt( (b[i].x*b[i].x + b[i].y*b[i].y + b[i].z*b[i].z) );
        percenterror += (abs(pos1 - pos2)) / pos1;
    }
    rmsd = std::sqrt(rmsd / a.size());
    //double mean_rel = rel_sum / a.size();
    double mean_percenterror = percenterror/a.size();

    std::cout << "Atoms compared: " << a.size() << "\n"
              << "RMSD = " << rmsd << "\n"
              //<< "Mean relative deviation = " << mean_rel*100 << "%\n"
              << "Mean Percent Error = " << mean_percenterror*100 << "%\n";

    //for(int i = 0; i < a.size(); i++){
    //    std::cout << "checking reading: "<< a[i].x << " " << a[i].y << a[i].z << "\n";
    //}

}
