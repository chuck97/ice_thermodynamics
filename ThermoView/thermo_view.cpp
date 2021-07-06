#include "thermo_view.h"

using namespace std;
using namespace INMOST;


ThermoMesh::ThermoMesh(int numCells_)
{
    numCells = numCells_;
    numNodes = numCells + 1;

    // create cell thickness vector
    vector<double> unit_thickness;
    for (int i = 0; i < numCells; ++i)
    {
        unit_thickness.push_back(1.0);
    }

    cellThickness = unit_thickness;
};

void ThermoMesh::AssignNodeData(std::string varName,
                                std::vector<double> nodeData_)
{
    nodeData[varName] = nodeData_;
};

void ThermoMesh::AssignCellData(std::string varName,
                             std::vector<double> cellData_)
{
    cellData[varName] = cellData_;
};

void ThermoMesh::AssignCellThickness(std::vector<double> cellThickness_)
{
    if (cellThickness_.size() != numCells)
    {
        ERR("thicknesses array size must be equal to number of cells");
    }

    cellThickness = cellThickness_;
};

void ThermoMesh::SetWidth(double cellWidth_)
{
    cellWidth = cellWidth_;
};

void ThermoMesh::AssignData(Mesh* m, bool is_sigma)
{
    // Assign Cell Data
    
    int counter = 0;

    for (auto& var : cellData)
    {
        Tag var_tag = m->CreateTag(var.first,
                                   DATA_REAL,
                                   CELL,
                                   NONE,
                                   1
                                   );
        
        for(Mesh::iteratorCell cellit = m->BeginCell();
                               cellit != m->EndCell();
                             ++cellit)
        {
            cellit->Real(var_tag) = var.second[counter];
            ++counter;
        }
        counter = 0;

        // plot cell data for nodes

        Tag node_var_tag = m->CreateTag(var.first,
                                        DATA_REAL,
                                        NODE,
                                        NONE,
                                        1
                                        );

        int node_counter = 0;
        int cell_counter = 0;

        for(Mesh::iteratorNode nodeit = m->BeginNode();
                               nodeit != m->EndNode();
                             ++nodeit)
        {
            if (cell_counter == 0)
            {
                nodeit->Real(node_var_tag) = var.second[0];
            }
            else if (cell_counter == numCells)
            {
                nodeit->Real(node_var_tag) = var.second[numCells-1];
            }
            else
            {
                nodeit->Real(node_var_tag) = 0.5*(var.second[cell_counter - 1] +
                                                  var.second[cell_counter]);
            }

            ++node_counter;

            if (node_counter%2 == 0)
            {
                ++cell_counter;
            }
        }
        counter = 0;

    }

    // Assign cell thickness 

    Tag z_tag = m->CreateTag("thickness",
                             DATA_REAL,
                             CELL,
                             NONE,
                             1
                             );

    for(Mesh::iteratorCell cellit = m->BeginCell();
                           cellit != m->EndCell();
                         ++cellit)
    {
        cellit->Real(z_tag) = cellThickness[counter];
        ++counter;
    }

    counter = 0;

    // node data

    for (auto& var : nodeData)
    {
        Tag var_tag = m->CreateTag(var.first,
                                   DATA_REAL,
                                   NODE,
                                   NONE,
                                   1
                                   );

        for(Mesh::iteratorNode nodeit = m->BeginNode();
                               nodeit != m->EndNode();
                             ++nodeit)
        {
            nodeit->Real(var_tag) = var.second[counter/2];
            ++counter;
        }
        counter = 0;
    }

};

void ThermoMesh::PlotMesh(std::string outputPath,
                           bool is_sigma)
{
    // create mesh
    Mesh *m;
    m = new Mesh();

    ElementArray<Node> nodes(m);
    nodes.reserve(numNodes*2);

    Storage::real coords[3];
    
    // bottom left node
    coords[0] = 0.0;
    coords[1] = 0.0;
    coords[2] = 0.0;
    Node node_bottom_left = m->CreateNode(coords);
    nodes.push_back(node_bottom_left);

    // bottom right node
    coords[0] = cellWidth;
    coords[1] = 0.0;
    coords[2] = 0.0;
    Node node_bottom_right = m->CreateNode(coords);
    nodes.push_back(node_bottom_right);

    double current_bottom = 0.0;

    for (int ncell = 0; ncell < numCells; ++ncell)
    {
        if (is_sigma)
        {
            current_bottom += 1.0/numCells;
        }
        else
        {
            current_bottom += cellThickness[ncell];
        }

        // left node
        coords[0] = 0.0;
        coords[1] = current_bottom;
        coords[2] = 0.0;
        Node node_left = m->CreateNode(coords);
        nodes.push_back(node_left);
        
        // right node
        coords[0] = cellWidth;
        coords[1] = current_bottom;
        coords[2] = 0.0;
        Node node_right = m->CreateNode(coords);
        nodes.push_back(node_right);

        // create rectangle cell
        ElementArray<Node> rect_verts(m);
        rect_verts.push_back(nodes[ncell*2]);
        rect_verts.push_back(nodes[ncell*2 + 1]);
        rect_verts.push_back(nodes[(ncell + 1)*2 + 1]);
        rect_verts.push_back(nodes[(ncell + 1)*2]);

        const INMOST_DATA_INTEGER_TYPE ne_face_nodes[8] = {0, 1, 1, 2, 2, 3, 3, 0};
        const INMOST_DATA_INTEGER_TYPE ne_num_nodes[4] = {2, 2, 2, 2};
        
        std::pair<Cell, bool> pair = m->CreateCell(rect_verts, ne_face_nodes, ne_num_nodes, 4);
    }

    AssignData(m, is_sigma);

    if (outputPath.substr(outputPath.size()-4, 4) != ".vtu")
    {
        ERR("Can't write mesh to this file. Filename should ended by .vtu")
    }
    
	m->Save(outputPath);
	
    cout << "Mesh saved to " << outputPath << std::endl;
   
    delete m;
};

void ThermoMesh::PlotSigma(std::string outputPath)
{
    PlotMesh(outputPath, true);
};

void ThermoMesh::PlotVertical(std::string outputPath)
{
    PlotMesh(outputPath, false);
};
