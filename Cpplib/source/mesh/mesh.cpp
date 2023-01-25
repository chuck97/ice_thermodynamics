#include "mesh.hpp"

namespace icethermo
{
    template<typename NumType>
    Mesh<NumType>::Mesh()
    {
    };

    template<typename NumType>
    Mesh<NumType>::Mesh(int n_uniform_layers, NumType thickness)
    {
        if (n_uniform_layers <= 0)
        {
            THERMO_ERR("Number of layers in mesh should be greater than 1!");
        }

        std::vector<NumType> thicknesses(n_uniform_layers);
        
        for (int i =0; i < n_uniform_layers; ++i)
        {
            thicknesses[i] = thickness/n_uniform_layers;
        }
        cells_thickness = std::make_shared<std::vector<NumType>>(std::move(thicknesses));
    }

    template<typename NumType>
    Mesh<NumType>::Mesh(NumType thickness): Mesh(10, thickness)
    {
    }

    template<typename NumType>
    Mesh<NumType>::Mesh(const std::vector<NumType>& unit_segment_decomposition, NumType thickness)
    {   
        NumType sum_decomp = sum_vec(unit_segment_decomposition);
        
        if (std::abs(sum_decomp - 1.0) > 1e-5)
        {
            THERMO_ERR("Unit segment decomposition of length 1.0 should be given!");
        }

        std::vector<NumType> thicknesses = unit_segment_decomposition*thickness;
        cells_thickness = std::make_shared<std::vector<NumType>>(thicknesses);
    }

    template<typename NumType>
    Mesh<NumType>::Mesh(const Mesh<NumType>& other)
    {   
        this->cells_thickness = std::make_shared<std::vector<NumType>>(*(other.cells_thickness));

        for (auto item: other.single_data)
        {
            auto key = item.first;
            auto value = item.second;

            (this->single_data)[key] = 
            std::make_pair
            (
                std::make_shared<NumType>(*(value.first)),
                value.second
            );
        }

        for (auto item: other.cells_data)
        {
            auto key = item.first;
            auto value = item.second;

            (this->cells_data)[key] = 
            std::make_pair
            (
                std::make_shared<std::vector<NumType>>(*(value.first)),
                value.second
            );
        }

        for (auto item: other.nodes_data)
        {
            auto key = item.first;
            auto value = item.second;

            (this->nodes_data)[key] = 
            std::make_pair
            (
                std::make_shared<std::vector<NumType>>(*(value.first)),
                value.second
            );
        }
    }

    template<typename NumType>
    Mesh<NumType>::~Mesh()
    {   
        //delete cells_thickness.get();

        //std::cout << 111 << std::endl;
        
        //for (auto item: this->single_data)
        //{
        //    auto value = item.second;

        //    delete (value.first).get();
        //}

        //std::cout << 222 << std::endl;

        //for (auto item: this->cells_data)
        //{
        //    auto value = item.second;

        //    delete (value.first).get();
        //}

        //std::cout << 333 << std::endl;

        //for (auto item: this->nodes_data)
        //{
        //    auto value = item.second;

        //    delete (value.first).get();
        //}

        //std::cout << 444 << std::endl;
    }

    template<typename NumType>
    int Mesh<NumType>::GetCellsNum() const
    {
        return cells_thickness->size();
    }

    template<typename NumType>
    int Mesh<NumType>::GetNodesNum() const
    {
        return cells_thickness->size() + 1;
    }

    template<typename NumType>
    std::shared_ptr<NumType> Mesh<NumType>::CreateSingleData(const std::string& varname, bool visible)
    {
        if (single_data.count(varname) != 0)
        {
            THERMO_ERR("Variable \'" + varname + "\' already exists, could not create single variable!");
        }
        NumType zero_val = 0;
        single_data[varname] = {std::make_shared<NumType>(std::move(zero_val)), visible};
        return single_data[varname].first;
    }

    template<typename NumType>
    std::shared_ptr<std::vector<NumType>> Mesh<NumType>::CreateCellsData(const std::string& varname, bool visible)
    {
        if (cells_data.count(varname) != 0)
        {
            THERMO_ERR("Variable \'" + varname + "\' already exists, could not create cell variable!");
        }

        // make zero vector
        std::vector<NumType> zero_vec(cells_thickness->size());
        cells_data[varname] = {std::make_shared<std::vector<NumType>>(std::move(zero_vec)), visible};
        return cells_data[varname].first;
    }

    template<typename NumType>
    std::shared_ptr<std::vector<NumType>> Mesh<NumType>::CreateNodesData(const std::string& varname, bool visible)
    {
        if (nodes_data.count(varname) != 0)
        {
            THERMO_ERR("Variable \'" + varname + "\' already exists, could not create node variable!");
        }

        // make zero vector
        std::vector<NumType> zero_vec(cells_thickness->size() + 1);
        nodes_data[varname] = {std::make_shared<std::vector<NumType>>(std::move(zero_vec)), visible};
        return nodes_data[varname].first;
    }

    template<typename NumType>
    void Mesh<NumType>::DeleteSingleData(const std::string& varname)
    {
        if (single_data.count(varname) == 0)
        {
            return;
        }
        single_data[varname].first.reset();
        single_data.erase(varname);
    }

    template<typename NumType>
    void Mesh<NumType>::DeleteCellsData(const std::string& varname)
    {
        if (cells_data.count(varname) == 0)
        {
            return;
        }
        cells_data[varname].first.reset();
        cells_data.erase(varname);
    }

    template<typename NumType>
    void Mesh<NumType>::DeleteNodesData(const std::string& varname)
    {
        if (nodes_data.count(varname) == 0)
        {
            return;
        }
        nodes_data[varname].first.reset();
        nodes_data.erase(varname);
    }

    template<typename NumType>
    std::shared_ptr<NumType> Mesh<NumType>::GetSingleData(const std::string& varname)
    {
        if (single_data.count(varname) == 0)
        {
            THERMO_ERR("There is no single variable: \'" + varname + "\' - can't get!");
        }
        return single_data[varname].first;
    }

    template<typename NumType>
    std::shared_ptr<std::vector<NumType>> Mesh<NumType>::GetCellsData(const std::string& varname)
    {
        if (cells_data.count(varname) == 0)
        {
            THERMO_ERR("There is no cell variable: \'" + varname + "\' - can't get!");
        }
        return cells_data[varname].first;
    }

    template<typename NumType>
    std::shared_ptr<std::vector<NumType>> Mesh<NumType>::GetNodesData(const std::string& varname)
    {
        if (nodes_data.count(varname) == 0)
        {
            THERMO_ERR("There is no node variable: \'" + varname + "\' - can't get!");
        }
        return nodes_data[varname].first;
    }

    template<typename NumType>
    std::shared_ptr<std::vector<NumType>> Mesh<NumType>::GetCellsThickness()
    {
        return cells_thickness;
    }

    template<typename NumType>
    NumType Mesh<NumType>::GetTotalThickness() const
    {
        return sum_vec(*cells_thickness);
    }

    template<typename NumType>
    void Mesh<NumType>::MuteSingleData(const std::string& varname)
    {
        if (single_data.count(varname) == 0)
        {
            return;
        }
        single_data[varname].second = false;
    }

    template<typename NumType>
    void Mesh<NumType>::MuteCellData(const std::string& varname)
    {
        if (cells_data.count(varname) == 0)
        {
            return;
        }
        cells_data[varname].second = false;
    }

    template<typename NumType>
    void Mesh<NumType>::MuteNodeData(const std::string& varname)
    {
        if (nodes_data.count(varname) == 0)
        {
            return;
        }
        nodes_data[varname].second = false;
    }

    template<typename NumType>
    void Mesh<NumType>::UnmuteSingleData(const std::string& varname)
    {
        if (single_data.count(varname) == 0)
        {
            return;
        }
        single_data[varname].second = true;
    }

    template<typename NumType>
    void Mesh<NumType>::UnmuteCellData(const std::string& varname)
    {
        if (cells_data.count(varname) == 0)
        {
            return;
        }
        cells_data[varname].second = true;
    }

    template<typename NumType>
    void Mesh<NumType>::UnmuteNodeData(const std::string& varname)
    {
        if (nodes_data.count(varname) == 0)
        {
            return;
        }
        nodes_data[varname].second = true;
    }

    template<typename NumType>
    void Mesh<NumType>::SaveTXT(const std::string& filename) const
    {
        std::string filename_txt = filename + ".txt";
        std::fstream* ofs = new std::fstream;
        ofs->open(filename_txt, std::ios::out);

        if (!ofs->is_open())
        {
            THERMO_ERR("can't open file "+ filename_txt + " for logging!");
        }

        // cell thickness info
        *ofs << "### cells_thickness_array ###\n";
        for (int i = 0; i < cells_thickness->size(); i++)
        {
            if (i != cells_thickness->size() - 1)
            {
                *ofs << (*cells_thickness)[i] << " "; 
            }
            else
            {
                *ofs << (*cells_thickness)[i];
            }
        }
        *ofs << "\n";

        // single data
        *ofs << "#### Single data ###\n";
        for (auto item: single_data)
        {
            auto key = item.first;
            auto val = item.second;

            if (val.second)
            {
                *ofs << key + "\n";
                *ofs << *(val.first);
                *ofs << "\n";
            }
        }

        // cells data
        *ofs << "#### Cells data ###\n";
        for (auto item: cells_data)
        {
            auto key = item.first;
            auto val = item.second;

            if (val.second)
            {
                *ofs << key + "\n";
                for (int i = 0; i < val.first->size(); i++)
                {
                    if (i != val.first->size() - 1)
                    {
                        *ofs << (*val.first)[i] << " "; 
                    }
                    else
                    {
                        *ofs << (*val.first)[i];
                    }
                }
                *ofs << "\n";
            }
        }

        // nodes data
        *ofs << "#### Nodes data ###\n";
        for (auto item: nodes_data)
        {
            auto key = item.first;
            auto val = item.second;

            if (val.second)
            {
                *ofs << key + "\n";
                for (int i = 0; i < val.first->size(); i++)
                {
                    if (i != val.first->size() - 1)
                    {
                        *ofs << (*val.first)[i] << " "; 
                    }
                    else
                    {
                        *ofs << (*val.first)[i];
                    }
                }
                *ofs << "\n";
            }
        }

        ofs->close(); 
        delete ofs;

        std::cout << "Mesh saved to \'" + filename_txt + "\'\n";
    }

    template<typename NumType>
    void Mesh<NumType>::SaveTXT(const std::string& filename, int postscript) const
    {
        std::string file = filename;
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(5) << postscript;
        file += ss.str();
        SaveTXT(file);
    }

#ifdef USE_JSON_OUTPUT

    template<typename NumType>
    void Mesh<NumType>::SaveJSON(const std::string& filename) const
    {
        // create empty json object
        json j;

        // cells thickness
        j["cells_thickness_array"] = *cells_thickness;

        // single data
        for (auto item: single_data)
        {
            auto key = item.first;
            auto val = item.second;
            if (val.second)
            {
                j["single data"][key] = *(val.first);
            }
        }

        // cells data
        for (auto item: cells_data)
        {
            auto key = item.first;
            auto val = item.second;
            if (val.second)
            {
                j["cells data"][key] = *(val.first);
            }
        }

        // nodes data
        for (auto item: nodes_data)
        {
            auto key = item.first;
            auto val = item.second;

            if (val.second)
            {
                j["nodes data"][key] = *(val.first);
            }
        }

        // write json object to file
        std::string filename_json = filename + ".json";
        std::fstream* ofs = new std::fstream;
        ofs->open(filename_json, std::ios::out);

        if (!ofs->is_open())
        {
            THERMO_ERR("can't open file "+ filename_json + " for logging!");
        }

        *ofs << std::setw(4) << j;

        ofs->close(); 
        delete ofs;

        std::cout << "Mesh saved to \'" + filename_json + "\'\n";
    }

#endif

#ifdef USE_JSON_OUTPUT

    template<typename NumType>
    void Mesh<NumType>::SaveJSON(const std::string& filename, int postscript) const
    {
        std::string file = filename;
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(5) << postscript;
        file += ss.str();
        SaveJSON(file);
    }
#endif

    template<typename NumType>
    bool Mesh<NumType>::CheckCellsDataExistency(const std::string& varname) const
    {
        if (cells_data.count(varname) == 0)
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    template<typename NumType>
    bool Mesh<NumType>::CheckNodesDataExistency(const std::string& varname) const
    {
        if (nodes_data.count(varname) == 0)
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    template<typename NumType>
    bool Mesh<NumType>::CheckSingleDataExistency(const std::string& varname) const
    {
        if (single_data.count(varname) == 0)
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    // explicit instantiation of classes
    template class Mesh<float>;
    template class Mesh<double>;
}