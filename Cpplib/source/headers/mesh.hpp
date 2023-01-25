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
#include <memory>

#ifdef USE_JSON_OUTPUT
#include <nlohmann/json.hpp>
using json = nlohmann::json;
#endif


#include "defines.hpp"
#include "matvec.hpp"

namespace icethermo
{
    template <typename NumType>
    class Mesh
    {
    public:
        // Constructors
        Mesh(); // default constructor for empty mesh
        Mesh(NumType thickness); // constructor with 10 uniform layers with given thicknes
        Mesh(int n_uniform_layers, NumType thickness); // constructor with given number of uniform layers and given thickness
        Mesh(const std::vector<NumType>& unit_segment_decomposition, NumType thickness); // constructor with manual partition for sigma layer grid and thickness
        Mesh(const Mesh<NumType>& other); // copy constructor
        ~Mesh(); // destructor

        // cells and nodes number getters
        int GetCellsNum() const;
        int GetNodesNum() const;

        // Creators of single, cells and nodes data
        std::shared_ptr<NumType> CreateSingleData(const std::string& varname, bool visible = true);
        std::shared_ptr<std::vector<NumType>> CreateCellsData(const std::string& varname, bool visible = true);
        std::shared_ptr<std::vector<NumType>> CreateNodesData(const std::string& varname, bool visible = true);

        // Deleters of single, cells and nodes data
        void DeleteSingleData(const std::string& varname);
        void DeleteCellsData(const std::string& varname);
        void DeleteNodesData(const std::string& varname);

        // Getters of single, cells and nodes data
        std::shared_ptr<NumType> GetSingleData(const std::string& varname);
        std::shared_ptr<std::vector<NumType>> GetCellsData(const std::string& varname);
        std::shared_ptr<std::vector<NumType>> GetNodesData(const std::string& varname);

        // Getter for cell thicknesses and total thickness
        std::shared_ptr<std::vector<NumType>> GetCellsThickness();
        NumType GetTotalThickness() const;

        // Muters and Unmuters
        void MuteSingleData(const std::string& varname);
        void MuteCellData(const std::string& varname);
        void MuteNodeData(const std::string& varname);
        
        void UnmuteSingleData(const std::string& varname);
        void UnmuteCellData(const std::string& varname);
        void UnmuteNodeData(const std::string& varname);

        // Write mesh to File
        void SaveTXT(const std::string& filename) const;
        void SaveTXT(const std::string& filename, int postscript) const;

#ifdef USE_JSON_OUTPUT

        void SaveJSON(const std::string& filename) const;
        void SaveJSON(const std::string& filename, int postscript) const;
#endif

        // Check existency of data
        bool CheckCellsDataExistency(const std::string& varname) const;
        bool CheckNodesDataExistency(const std::string& varname) const;
        bool CheckSingleDataExistency(const std::string& varname) const;

    private:
        // vector of cell thicknesses
        std::shared_ptr<std::vector<NumType>> cells_thickness;

        // container for stand-alone variables
        std::map<std::string, std::pair<std::shared_ptr<NumType>, bool>> single_data;

        // container for cell data [name : {vector_of_values, is_visible}]
        std::map<std::string, std::pair<std::shared_ptr<std::vector<NumType>>, bool>> cells_data;

        // container for node data [name : {vector_of_values, is_visible}]
        std::map<std::string, std::pair<std::shared_ptr<std::vector<NumType>>, bool>> nodes_data;
    };
}