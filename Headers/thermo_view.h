#pragma once
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "inmost.h"

#define ERR(message) {std::cerr << "Error: " << message  << std::endl; exit(1);}

class ThermoMesh
{
public:
    ThermoMesh(int numCells_);

    void AssignNodeData(std::string varName,
                        std::vector<double> nodeData_);
    
    void AssignCellData(std::string varName,
                        std::vector<double> cellData_);

    void AssignCellThickness(std::vector<double> cellThickness);

    void SetWidth(double cellWidth_);

    void PlotSigma(std::string outputPath);

    void PlotVertical(std::string outputPath);
    
private:
    std::map<std::string, std::vector<double>> cellData;
    std::map<std::string, std::vector<double>> nodeData;

    std::vector<double> cellThickness;

    int numNodes;
    int numCells;
    double cellWidth = 0.1;

    void PlotMesh(std::string outputPath,
                   bool is_sigma);

    void AssignData(INMOST::Mesh* m, bool is_sigma);
};
