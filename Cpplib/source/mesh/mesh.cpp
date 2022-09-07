#include "mesh.hpp"

namespace icethermo
{
    template<typename NumType>
    Mesh<NumType>::Mesh(int n_uniform_layers, NumType thickness)
    {
        if (n_uniform_layers <= 0)
        {
            THERMO_ERR("Number of layers in mesh should be greater than 1!");
        }

        if (thickness <= 0.0)
        {
            THERMO_ERR("mesh thickness should be greater than 0.0!");
        }

        std::vector<NumType> thicknesses(n_uniform_layers);
        
        for (int i =0; i < n_uniform_layers; ++i)
        {
            thicknesses[i] = thickness/n_uniform_layers;
        }
        cells_thickness = std::move(thicknesses);
    }

    template<typename NumType>
    Mesh<NumType>::Mesh(NumType thickness): Mesh(10, thickness)
    {
    }

    template<typename NumType>
    Mesh<NumType>::Mesh(const std::vector<NumType>& unit_segment_decomposition, NumType thickness)
    {
        if (thickness <= 0.0)
        {
            THERMO_ERR("mesh thickness should be greater than 0.0!");
        }
        
        NumType sum_decomp = std::accumulate(unit_segment_decomposition.begin(), unit_segment_decomposition.end(),
                                             NumType(0));
        
        if (std::abs(sum_decomp - 1.0) > 1e-5)
        {
            THERMO_ERR("Unit segment decomposition of length 1.0 should be given!");
        }

        std::vector<NumType> thicknesses = unit_segment_decomposition*thickness;
        cells_thickness = std::move(thicknesses);
    }

    template<typename NumType>
    int Mesh<NumType>::GetCellsNum() const
    {
        return cells_thickness.size();
    }

    template<typename NumType>
    int Mesh<NumType>::GetNodesNum() const
    {
        return cells_thickness.size() + 1;
    }

    template<typename NumType>
    std::vector<NumType>& Mesh<NumType>::CreateCellData(const std::string& varname, bool visible)
    {
        if (cells_data.count(varname) != 0)
        {
            THERMO_ERR("Variable \'" + varname + "\' already exists, could not create cell variable!");
        }

        // make zero vector
        std::vector<NumType> zero_vec(cells_thickness.size());
        cells_data[varname] = {zero_vec, visible};
        return cells_data[varname].first;
    }

    template<typename NumType>
    std::vector<NumType>& Mesh<NumType>::CreateNodeData(const std::string& varname, bool visible)
    {
        if (nodes_data.count(varname) != 0)
        {
            THERMO_ERR("Variable \'" + varname + "\' already exists, could not create node variable!");
        }

        // make zero vector
        std::vector<NumType> zero_vec(cells_thickness.size() + 1);
        nodes_data[varname] = {zero_vec, visible};
        return nodes_data[varname].first;
    }

    template<typename NumType>
    void Mesh<NumType>::DeleteCellData(const std::string& varname)
    {
        if (cells_data.count(varname) == 0)
        {
            return;
        }
        cells_data.erase(varname);
    }

    template<typename NumType>
    void Mesh<NumType>::DeleteNodeData(const std::string& varname)
    {
        if (nodes_data.count(varname) == 0)
        {
            return;
        }
        nodes_data.erase(varname);
    }

    template<typename NumType>
    std::vector<NumType>& Mesh<NumType>::GetCellData(const std::string& varname)
    {
        if (cells_data.count(varname) == 0)
        {
            THERMO_ERR("There is no cell variable: \'" + varname + "\' - can't get!");
        }
        return cells_data[varname].first;
    }

    template<typename NumType>
    std::vector<NumType>& Mesh<NumType>::GetNodeData(const std::string& varname)
    {
        if (nodes_data.count(varname) == 0)
        {
            THERMO_ERR("There is no node variable: \'" + varname + "\' - can't get!");
        }
        return nodes_data[varname].first;
    }

    template<typename NumType>
    std::vector<NumType>& Mesh<NumType>::GetCellsThickness()
    {
        return cells_thickness;
    }

    template<typename NumType>
    NumType Mesh<NumType>::GetTotalThickness() const
    {
        return std::accumulate(cells_thickness.begin(), cells_thickness.end(), NumType(0));
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
        *ofs << "### Cells thickness ###\n";
        for (int i = 0; i < cells_thickness.size(); i++)
        {
            if (i != cells_thickness.size() - 1)
            {
                *ofs << cells_thickness[i] << " "; 
            }
            else
            {
                *ofs << cells_thickness[i];
            }
        }
        *ofs << "\n";

        // cells data
        *ofs << "#### Cells data ###\n";
        for (auto [key, val]: cells_data)
        {
            if (val.second)
            {
                *ofs << key + "\n";
                for (int i = 0; i < val.first.size(); i++)
                {
                    if (i != val.first.size() - 1)
                    {
                        *ofs << val.first[i] << " "; 
                    }
                    else
                    {
                        *ofs << val.first[i];
                    }
                }
                *ofs << "\n";
            }
        }

        // nodes data
        *ofs << "#### Nodes data ###\n";
        for (auto [key, val]: nodes_data)
        {
            if (val.second)
            {
                *ofs << key + "\n";
                for (int i = 0; i < val.first.size(); i++)
                {
                    if (i != val.first.size() - 1)
                    {
                        *ofs << val.first[i] << " "; 
                    }
                    else
                    {
                        *ofs << val.first[i];
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


    // explicit instantiation of classes
    template class Mesh<float>;
    template class Mesh<double>;
}