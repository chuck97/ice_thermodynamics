#pragma once
#include <vector>
#include <map>
#include <string>
#include <utility>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "defines.hpp"
#include "matvec.hpp"

namespace icethermo
{
    template <typename NumType>
    class Mesh
    {
    public:
        // Constructors
        Mesh(NumType thickness); // default constructor with 10 uniform layers with given thicknes
        Mesh(int n_uniform_layers, NumType thickness); // constructor with given number of uniform layers and given thickness
        Mesh(const std::vector<NumType>& unit_segment_decomposition, NumType thickness); // constructor with manual partition for sigma layer grid and thickness

        // cells and nodes number getters
        int GetCellsNum() const;
        int GetNodesNum() const;

        // Creators of cells and nodes data
        std::vector<NumType>& CreateCellData(const std::string& varname, bool visible = true);
        std::vector<NumType>& CreateNodeData(const std::string& varname, bool visible = true);

        // Deleters of cells and nodes data
        void DeleteCellData(const std::string& varname);
        void DeleteNodeData(const std::string& varname);

        // Getters of cells and nodes data
        std::vector<NumType>& GetCellData(const std::string& varname);
        std::vector<NumType>& GetNodeData(const std::string& varname);

        // Getter for cell thicknesses and total thickness
        std::vector<NumType>& GetCellsThickness();
        NumType GetTotalThickness() const;

        // Muters and Unmuters
        void MuteCellData(const std::string& varname);
        void MuteNodeData(const std::string& varname);
        
        void UnmuteCellData(const std::string& varname);
        void UnmuteNodeData(const std::string& varname);

        // Write mesh to File
        void SaveTXT(const std::string& filename) const;
        void SaveTXT(const std::string& filename, int postscript) const;

    private:
        // vector of cell thicknesses
        std::vector<NumType> cells_thickness;

        // container for cell data [name : {vector_of_values, is_visible}]
        std::map<std::string, std::pair<std::vector<NumType>, bool>> cells_data;

        // container for node data [name : {vector_of_values, is_visible}]
        std::map<std::string, std::pair<std::vector<NumType>, bool>> nodes_data;
    };
}