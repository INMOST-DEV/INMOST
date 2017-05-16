#include "octgrid.h"
#include <math.h>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
using namespace std;

double epsilon = 0.0001;

#define IS_GHOST(c) (c.GetStatus() == Element::Ghost)
#define IS_GHOST_P(c) (c->GetStatus() == Element::Ghost)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);

const int MAX_LEVEL = 10;
void default_transformation(double xyz[3]) { (void) xyz; }
int default_cell_should_unite(struct grid * g, int cell) { (void) g; (void) cell; return 0; }
int default_cell_should_split(struct grid * g, int cell) { (void) g; (void) cell; return 0; }
void mark_cell(struct grid* g, Cell cell);
void cellCoarse(struct grid* g, Cell cell);
void resolve_connections(grid* g);
void resolve_edges(grid* g, Cell* c = NULL);
inline bool equal(double a, double b);
void correct_brothers(struct grid* g, int size, int rank, int type);
int global = 0;

double time_u;
double time_p;
double time_r1;
double time_r2;
double ttt;


class CR_cells
{
    ElementArray<Cell>* d_cells[MAX_LEVEL];
    vector<bool>* d_valid[MAX_LEVEL];
    ElementArray<Cell> d_cr_cells; // Array for refine and coarse
    vector<bool> d_cr_valids;
    bool d_empty;
    grid* g;

    public: 
        bool active;
        CR_cells() 
        {
            d_empty = true; 
            active = false;
            for (int i = 0; i < MAX_LEVEL; i++) d_cells[i]  = NULL;
            for (int i = 0; i < MAX_LEVEL; i++) d_valid[i] = NULL;
        }
        void push(Cell c) 
        { 
            int l = c.Integer(g->c_tags.level);
            if (d_cells[l] == NULL) d_cells[l]  = new ElementArray<Cell>;
            if (d_valid[l] == NULL) d_valid[l] = new vector<bool>;

            d_cells[l]->push_back(c);
            d_valid[l]->push_back(true);
            c.Integer(g->c_tags.i) = d_cells[l]->size() - 1;

/*
            d_cr_cells.push_back(c); 
            d_cr_valids.push_back(true); 
            */
            d_empty = false; 
        }
        void set_grid(grid* _g) { g = _g; }
        
//        bool check(int i)        { return i >= d_cr_valids.size() ? false : d_cr_valids[i]; }
        bool check(int l, int i) { return d_valid[l] == NULL ? false : ( i >= d_valid[l]->size() ? false : (*d_valid[l])[i]); }
        
  //      void dis(int i)   { d_cr_valids[i] = false; }
        void dis(int l, int i)   { (*d_valid[l])[i] = false; }
        
    //    int  size()       {  return d_cr_cells.size(); }
        int  size(int l)  {  return d_cells[l] == NULL ? 0 : d_cells[l]->size(); }

      //  Cell operator[](int i) { return d_cr_cells[i]; }
        Cell get(int l, int i) 
        {
             if (d_cells[l] == NULL) TSNH;
             return (*(d_cells[l]))[i]; 
        }

        void clear()      
        {
             d_cr_cells.clear(); d_cr_valids.clear(); 
            for (int i = 0; i < MAX_LEVEL; i++)
            {
                if (d_cells[i] != NULL) delete d_cells[i];
                if (d_valid[i] != NULL) delete d_valid[i];
                d_cells[i] = NULL;
                d_valid[i] = NULL;
            }
            d_empty = true;
        }

      //  void set(Cell c, int i)   { d_cr_cells[i] = c; }
        void set1(Cell c, int i)   
        {
            int l = c.Integer(g->c_tags.level);
             (*d_cells[l])[i] = c; 
            c.Integer(g->c_tags.i) = i;
        }

        bool empty()      { return d_empty; }
       // bool check_i(int i )        { return (i >= 0 && i < d_cr_cells.size()); }
        bool check_i(int l, int i ) { return d_cells[l] == NULL ? false : (i >= 0 && i < d_cells[l]->size()); }

        void print_size()
        {
            for (int i = 0; i < MAX_LEVEL; i++)
                if (d_cells[i] != NULL) cout << i << ": " << d_cells[i]->size() << endl;
        }
};
CR_cells cr_cells;

/// Dump mesh to vtk file in folder "grids"
void dump_to_vtk(grid* g, const char* suffix)
{
	//thegrid.mesh->ResolveShared(); // Resolve duplicate nodes
	//thegrid.mesh->ExchangeGhost(2,NODE); // Construct Ghost cells in 2 layers connected via nodes
    int rank = g->mesh->GetProcessorRank(); // Get the rank of the current process
	int size = g->mesh->GetProcessorsNumber();

    std::stringstream filename;
    filename << "grids/grid_";
    filename << size;
    filename << suffix;
    if( size == 1 )
        filename << ".vtk";
    else
        filename << ".pvtk";
    g->mesh->Save(filename.str());
	cout << "Process " << rank << ": dumped mesh to file: " << filename.str() << endl;
}

void fill_proc_tag(grid* g)
{
    int rank = g->mesh->GetProcessorRank(); // Get the rank of the current process
    for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); it++)
    {
        it->Integer(g->c_tags.proc) = rank;
    }
}

int calc_sends(grid* g)
{
    int rank = g->mesh->GetProcessorRank(); // Get the rank of the current process
    int res = 0;
    for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); it++)
    {
        if (it->Integer(g->c_tags.proc) != rank) res++;
    }

    return res;
}

/// Redistribute grid by  partitioner
void redistribute(grid* g, int type)
{
    g->mesh->RemoveGhost();
    int rank = g->mesh->GetProcessorRank(); // Get the rank of the current process
	int size = g->mesh->GetProcessorsNumber();

	//LOG(2,"Process " << rank << ": redistribute. Cells: " << g->mesh->NumberOfCells())
    Partitioner * part = new Partitioner(g->mesh);
    
    // Specify the partitioner
    type = 1;
    if (type == 0) part->SetMethod(Partitioner::Parmetis, Partitioner::Partition);
    if (type == 1) part->SetMethod(Partitioner::Parmetis, Partitioner::Repartition);
    if (type == 2) part->SetMethod(Partitioner::Parmetis, Partitioner::Refine);
    
    try
    {
        part->Evaluate();
    }
    catch (INMOST::ErrorType er)
    {
        cout << "Exception: " << er << endl;
    }
    catch(...)
    {
    }
    delete part;

    g->mesh->RemoveGhost();
	correct_brothers(g,size,rank, 2);

    try
    {
        g->mesh->Redistribute(); 
    }
    catch (INMOST::ErrorType er)
    {
        cout << "Exception: " << er << endl;
    }
    catch(...)
    {
    }


    g->mesh->RemoveGhost();
    g->mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);
    //g->mesh->AssignGlobalID(CELL | EDGE | FACE | NODE);
//	LOG(2,"Process " << rank << ": redistribute completed")
}

/// Create tags
void init_mesh(struct grid* g)
{
    g->c_tags.busy     = g->mesh->CreateTag("busy"    ,DATA_INTEGER,CELL,false,1);
    g->c_tags.leaf     = g->mesh->CreateTag("leaf"    ,DATA_INTEGER,CELL,false,1);
    g->c_tags.level    = g->mesh->CreateTag("level"   ,DATA_INTEGER,CELL,false,1);
    g->c_tags.children = g->mesh->CreateTag("children",DATA_INTEGER,CELL,false,1<<DIM);
    g->c_tags.parent   = g->mesh->CreateTag("parent"  ,DATA_INTEGER,CELL,false,1);
    g->c_tags.is_valid = g->mesh->CreateTag("is_valid",DATA_INTEGER,CELL,false,1);
    g->c_tags.proc     = g->mesh->CreateTag("proc"    ,DATA_INTEGER,CELL,false,1);
    g->c_tags.i        = g->mesh->CreateTag("i"       ,DATA_INTEGER,CELL,false,1);

    // Cell
    g->c_tags.center   = g->mesh->CreateTag("center",DATA_REAL,CELL,false,3);
    g->c_tags.side     = g->mesh->CreateTag("side"  ,DATA_REAL,CELL,false,3);
    g->c_tags.floor    = g->mesh->CreateTag("floor" ,DATA_INTEGER,CELL,false,1);
    g->c_tags.chld_num = g->mesh->CreateTag("chld_num",DATA_INTEGER,CELL,false,1);
    g->c_tags.par_chld_nums = g->mesh->CreateTag("par_chld_nums",DATA_INTEGER,CELL,false,8);
    g->c_tags.to_split = g->mesh->CreateTag("to_split",DATA_INTEGER,CELL,false,1);
    g->c_tags.base_id  = g->mesh->CreateTag("base_id" ,DATA_INTEGER,CELL,false,1);
}   

int check_all_edges(grid* g, int n)
{
    int rank = g->mesh->GetProcessorRank(); // Get the rank of the current process
    int was = 0;
    for(Mesh::iteratorEdge edge = g->mesh->BeginEdge(); edge != g->mesh->EndEdge(); edge++)
    {
        Node n1 = edge->getBeg();
        Node n2 = edge->getEnd();
        int diff = 0;
        if (!equal(n1.RealArray(g->mesh->CoordsTag())[0], n2.RealArray(g->mesh->CoordsTag())[0])) diff++;
        if (!equal(n1.RealArray(g->mesh->CoordsTag())[1], n2.RealArray(g->mesh->CoordsTag())[1])) diff++;
        if (!equal(n1.RealArray(g->mesh->CoordsTag())[2], n2.RealArray(g->mesh->CoordsTag())[2])) diff++;
        if (diff > 1)
        {
            print_edge(g,edge->getAsEdge());
            was = 1;
        }
    }
    return was;
}

/// Debug
void print_node_center(struct grid* g, Node node)
{
	cout << node.RealArray(g->mesh->CoordsTag())[0] << " ";
	cout << node.RealArray(g->mesh->CoordsTag())[1] << " ";
	cout << node.RealArray(g->mesh->CoordsTag())[2] << endl;
}

/// Debug
void print_edge(struct grid* g, Edge edge)
{
    ElementArray<Node> nodes = edge.getNodes();
    for (ElementArray<Node>::iterator node = nodes.begin(); node != nodes.end(); node++)
        print_node_center(g, node->getAsNode());
    cout << endl;
}

/// Debug
void print_face_nodes(struct grid* g, Face face)
{
    ElementArray<Node> nodes = face.getNodes();
    for (ElementArray<Node>::iterator node = nodes.begin(); node != nodes.end(); node++)
        print_node_center(g, node->getAsNode());
    cout << endl;
}

/// Debug
void print_face_edges(struct grid* g, Face face)
{
    ElementArray<Edge> edges = face.getEdges();
    for (ElementArray<Edge>::iterator edge = edges.begin(); edge!= edges.end(); edge++)
        print_edge(g, edge->getAsEdge());
    cout << endl;
}

/// Debug
void print_cell_center(struct grid* g, Cell cell)
{
    double xyz[3];
    xyz[0] = cell.RealArrayDF(g->c_tags.center)[0];
    xyz[1] = cell.RealArrayDF(g->c_tags.center)[1];
    xyz[2] = cell.RealArrayDF(g->c_tags.center)[2];
    g->transformation(xyz);
    cout << xyz[0] << " " << xyz[1] << " " << xyz[2];
}

/// Makes correct order od edges in face
void reorder_edges(struct grid* g, ElementArray<Edge>* edges)
{
    ElementArray<Node> nodes_1;
    ElementArray<Node> nodes_2;
      
    for (int i = 0; i < edges->size(); i++)
        for (int j = i + 1; j < edges->size(); j++)
        {
           nodes_1 = (*edges)[i].getNodes();
           nodes_2 = (*edges)[j].getNodes();
           if (nodes_1[0] == nodes_2[0] || nodes_1[0] == nodes_2[1] ||
               nodes_1[1] == nodes_2[0] || nodes_1[1] == nodes_2[1])
           {
               if (i + 1 != j) 
               {
                   Edge tmp = (*edges)[j];
                   (*edges)[j] = (*edges)[i+1];
                   (*edges)[i+1] = tmp;
               }
               break;
           }
        }
}

inline bool equal(double a, double b)
{
    return (fabs(a - b) < epsilon);
}

int compare(double a, double b)
{
    if (equal(a,b)) return 0;
    if (a < b)      return -1;
    return 1;
}

// Compare nodes. Return -1 if node1 < node2, 1 if node1 > node2, else 0
///      (2)*-------*(6)  
///     z  /|      /|     
///     ^ /       / |     
///     |/  |    /  |     
///  (1)*-------*(5)|     
///     |   |   |   |     
///     |     y |   |     
///     |   |/  |   |     
///     |(3)*- -|- -*(7)  
///     |  /    |  /      
///     |       | /       
///  (0)*-------*(4) -> x 
int compare_nodes(struct grid* g, Node node1, Node node2, double y1, double y2)              
{                                                                                            
    double cx_1 = node1.RealArray(g->mesh->CoordsTag())[0];                                  
    double cy_1 = node1.RealArray(g->mesh->CoordsTag())[1];                                  
    double cz_1 = node1.RealArray(g->mesh->CoordsTag())[2];                                   
                                                                                             
    double cx_2 = node2.RealArray(g->mesh->CoordsTag())[0];                                  
    double cy_2 = node2.RealArray(g->mesh->CoordsTag())[1];                                  
    double cz_2 = node2.RealArray(g->mesh->CoordsTag())[2];                                  
                                                                                                                                                                                               
    if (fabs(cx_1 - cx_2) > epsilon)
    {                                                                                                  
        if (cx_1 < cx_2) return -1;                                                               
        if (cx_1 > cx_2) return  1;                                
    }

    if (fabs(cy_1 - cy_2) > epsilon)
    {
        if (cy_1 < cy_2) return -1;                                
        if (cy_1 > cy_2) return  1;                                
    }                                                          
    
    if (fabs(cz_1 - cz_2) > epsilon)
    {
        if (fabs(cy_1 - y2) > epsilon)
        {                                                          
            if (cz_1 < cz_2) return -1;                                 
            if (cz_1 > cz_2) return  1;                                 
        }
        else
        {
            if (cz_1 > cz_2) return -1;                                 
            if (cz_1 < cz_2) return  1;                                 
        }
    }

    return 0;
}

void reorder_1(struct grid* g, ElementArray<Node>* nodes)
{
    if (nodes->size() == 0) return;
    double y1 = nodes->begin()->RealArray(g->mesh->CoordsTag())[1];
    double y2 = y1;
   
    for (int i = 1; i < nodes->size(); i++)
    {
        y2 = (*nodes)[i].RealArray(g->mesh->CoordsTag())[1];
        if (fabs(y1 - y2) > epsilon) break;
    }
    if (y1 > y2) 
    {
        double t = y1; 
        y1 = y2; 
        y2 = t;
    }

    for (int i = 0; i < nodes->size(); i++)
        for (int j = i+1; j < nodes->size(); j++)
            if (compare_nodes(g,(*nodes)[i],(*nodes)[j],y1,y2 ) > 0)
            {
                Node tmp = (*nodes)[i];
                (*nodes)[i] = (*nodes)[j];
                (*nodes)[j] = tmp;
            }
}
                 
/// Create cube cell from verts
Cell CreateCubeElement(struct grid* g, ElementArray<Node> verts, int ghost) 
{
    /*      (5)*-------*(6)  */
	/*     z  /|      /|     */
	/*     ^ /       / |     */
	/*     |/  |    /  |     */
	/*  (4)*-------*(7)|     */
	/*     |   |   |   |     */
	/*     |     y |   |     */
    /*     |   |/  |   |     */
    /*     |(1)*- -|- -*(2)  */
    /*     |  /    |  /      */
    /*     |       | /       */
    /*     |/      |/        */
    /*  (0)*-------*(3) -> x */

    // Define six cube faces assuming verts are numerated in the way presented above
    const INMOST_DATA_INTEGER_TYPE face_nodes[24] = {0,4,5,1, 1,5,6,2, 2,6,7,3, 3,7,4,0, 0,1,2,3, 4,5,6,7};
    const INMOST_DATA_INTEGER_TYPE num_nodes[6]   = {4,       4,       4,       4,       4,       4};

    Cell res = g->mesh->CreateCell(verts,face_nodes,num_nodes,6).first;
    if (ghost) res.SetStatus(Element::Ghost);
    return res;
}

int l_n[3];
int get_num(int i, int j, int k)
{
    return i*(l_n[2]+1)*(l_n[1]+1) + j*(l_n[2]+1) + k;
}

void gridInit(struct grid * g, int n[3])
{
	int i,j,k,l,m, c[3], v, r, env[1<<DIM], flag;
    double dx,dy,dz;
    double sx,sy,sz;
    cr_cells.set_grid(g);

	g->mesh = new Mesh();
	g->mesh->SetDimensions(3);
    g->last_base_id = 0;

    g->mesh->SetCommunicator(INMOST_MPI_COMM_WORLD);
    MPI_Comm_set_errhandler(INMOST_MPI_COMM_WORLD,MPI_ERRORS_RETURN);
    int rank = g->mesh->GetProcessorRank(); // Get the rank of the current process
	int size = g->mesh->GetProcessorsNumber(); // Get the number of processors used in communicator comm
    
    int b = size / 2;
    int a = size - b;

    dx = 1.0/(double)n[0];
    dy = 1.0/(double)n[1];
    dz = 1.0/(double)n[2];

    if (size == 9 || size == 16)
    {
        int count_in = (size == 9 ? 3 : 4);
        l_n[2] = n[2];
        l_n[1] = n[1] / count_in;
        l_n[0] = n[0] / count_in;

        sy = (rank / count_in)*dy*l_n[1];
        sx = (rank % count_in)*dx*l_n[0];
        sz = 0;
    }
    else
    {
        l_n[2] = n[2];
        if (size > 1) l_n[1] = n[1]/2;
        else          l_n[1] = n[1];

        if (rank >= a) l_n[0] = n[0]/b;
        else           l_n[0] = n[0]/a;

        if (rank == size - 1) l_n[0] = n[0] - l_n[0]*(b-1);
        if (rank == a - 1)    l_n[0] = n[0] - l_n[0]*(rank);

        if (rank >= a) {
            sx = (rank - a)*(n[0]/b)*dx;
            sy = (n[1]/2)*dy;
        } else {
            sx = rank*(n[0]/a)*dx;
            sy = 0;
        }
        sz = 0;
    }

    init_mesh(g); 
    Storage::real xyz[3];
    ElementArray<Node> newverts(g->mesh);
    newverts.reserve(l_n[0]*l_n[1]*l_n[2]);
	for(i = 0; i <= l_n[0]; i++)
	{
		for(j = 0; j <= l_n[1]; j++)
		{
			for(k = 0; k <= l_n[2]; k++)
			{
                xyz[0] = sx + dx*i;
                xyz[1] = sy + dy*j;
                xyz[2] = sz + dz*k;
                g->transformation(xyz); 
                newverts.push_back(g->mesh->CreateNode(xyz)); 
            }
        }
    }
    
    /*      (5)*-------*(6)  */
	/*     z  /|      /|     */
	/*     ^ /       / |     */
	/*     |/  |    /  |     */
	/*  (4)*-------*(7)|     */
	/*     |   |   |   |     */
	/*     |     y |   |     */
    /*     |   |/  |   |     */
    /*     |(1)*- -|- -*(2)  */
    /*     |  /    |  /      */
    /*     |       | /       */
    /*     |/      |/        */
    /*  (0)*-------*(3) -> x */

    for(i = 0; i < l_n[0]; i++)
    {
        for(j = 0; j < l_n[1]; j++)
        {
            for(k = 0; k < l_n[2]; k++)
            {
                ElementArray<Node> e_nodes;
                e_nodes.push_back(newverts[get_num(i+0,j+0,k+0)]);
                e_nodes.push_back(newverts[get_num(i+0,j+1,k+0)]);
                e_nodes.push_back(newverts[get_num(i+1,j+1,k+0)]);
                e_nodes.push_back(newverts[get_num(i+1,j+0,k+0)]);

                e_nodes.push_back(newverts[get_num(i+0,j+0,k+1)]);
                e_nodes.push_back(newverts[get_num(i+0,j+1,k+1)]);
                e_nodes.push_back(newverts[get_num(i+1,j+1,k+1)]);
                e_nodes.push_back(newverts[get_num(i+1,j+0,k+1)]);
                Cell c = CreateCubeElement(g, e_nodes,0);

                c.RealArrayDF(g->c_tags.center)[0] = sx + dx*(i+0.5); 
                c.RealArrayDF(g->c_tags.center)[1] = sy + dy*(j+0.5); 
                c.RealArrayDF(g->c_tags.center)[2] = sz + dz*(k+0.5); 

                c.RealArrayDF(g->c_tags.side)[0] = dx;
                c.RealArrayDF(g->c_tags.side)[1] = dy;
                c.RealArrayDF(g->c_tags.side)[2] = dz;

                c.Integer(g->c_tags.floor) = k;
                c.Integer(g->c_tags.level) = 0;
                c.Integer(g->c_tags.chld_num) = -1;
                c.Integer(g->c_tags.base_id) = g->last_base_id++ + (rank * n[0]*n[1]*n[2]);;
                for (int i = 0; i < 8; i++)
                    c.IntegerArrayDF(g->c_tags.par_chld_nums)[i] = -1;

                c.Integer(g->c_tags.to_split) = 0;
            }
        }
    }

	g->mesh->ResolveShared(); // Resolve duplicate nodes
	//g->mesh->ExchangeGhost(1,NODE); // Construct Ghost cells in 2 layers connected via nodes
}

/// Return adjacent cell with number cell_num 
/// result sets to -1 if cell not found
/// 1 2     5 6
/// 0 3     4 7
Cell get_adjacent_cell(struct grid* g, Node node, int cell_num, int* result)
{   
    int num[8] = { 6,5,1,2,7,4,0,3 };

    ElementArray<Element> adjs = node->getAdjElements(CELL); 
    for (ElementArray<Element>::iterator adj = adjs.begin(); adj != adjs.end(); adj++)
    {
        ElementArray<Node> nodes = adj->getAsCell().getNodes();
        reorder_1(g, &nodes);
        if ((nodes[num[cell_num]].RealArray(g->mesh->CoordsTag())[0] == node.RealArray(g->mesh->CoordsTag())[0]) &&
            (nodes[num[cell_num]].RealArray(g->mesh->CoordsTag())[1] == node.RealArray(g->mesh->CoordsTag())[1]) &&
            (nodes[num[cell_num]].RealArray(g->mesh->CoordsTag())[2] == node.RealArray(g->mesh->CoordsTag())[2]))
            {
                *result = 0;
                return adj->getAsCell();
            }
   }
   *result = -1;
	return InvalidCell();
}

/// Returns true if edge is not a greatest in this cell
bool isSmallEdge(struct grid* g, Edge edge, Cell cell)
{
    ElementArray<Node> nodes = edge.getNodes(); 
    double c1[3];
    double c2[3];
    double s[3];
    c1[0] = nodes[0].RealArray(g->mesh->CoordsTag())[0];
    c1[1] = nodes[0].RealArray(g->mesh->CoordsTag())[1];
    c1[2] = nodes[0].RealArray(g->mesh->CoordsTag())[2];
    c2[0] = nodes[1].RealArray(g->mesh->CoordsTag())[0];
    c2[1] = nodes[1].RealArray(g->mesh->CoordsTag())[1];
    c2[2] = nodes[1].RealArray(g->mesh->CoordsTag())[2];
    g->rev_transformation(c1);
    g->rev_transformation(c2);
    
    if (nodes.size() != 2) TSNH
    double p = 0;
         if (!equal(nodes[0].RealArray(g->mesh->CoordsTag())[0], nodes[1].RealArray(g->mesh->CoordsTag())[0]))
        p = cell.RealArrayDF(g->c_tags.side)[0];
    else if (!equal(nodes[0].RealArray(g->mesh->CoordsTag())[1], nodes[1].RealArray(g->mesh->CoordsTag())[1]))
        p = cell.RealArrayDF(g->c_tags.side)[1];
    else if (!equal(nodes[0].RealArray(g->mesh->CoordsTag())[2], nodes[1].RealArray(g->mesh->CoordsTag())[2]))
        p = cell.RealArrayDF(g->c_tags.side)[2];
    else 
    {
        TSNH
    }

    s[0] = cell.RealArrayDF(g->c_tags.side)[0];
    s[1] = cell.RealArrayDF(g->c_tags.side)[1];
    s[2] = cell.RealArrayDF(g->c_tags.side)[2];

    double l = sqrt(((c1[0]-c2[0])*(c1[0]-c2[0]) + (c1[1]-c2[1])*(c1[1]-c2[1]) + (c1[2] - c2[2])* (c1[2] - c2[2]) ));
    bool res;
    if (equal(p,l)) res =  false;
    else res = true;

    return res;
}

/// Returns true if face is smallest in this cell
bool isSmallFace(struct grid* g, Face face, Cell cell)
{
    // If face is small, then face has 4 small edges
    ElementArray<Edge> edges = face.getEdges();
    if (edges.size() != 4) return false;

    for (ElementArray<Edge>::iterator edge = edges.begin(); edge != edges.end(); edge++)
    {
        if (!isSmallEdge(g,edge->getAsEdge(),cell)) return false;
    }

    return true;
}

/// Creates cell with faces based on cell. After deletes cell
Cell recreate_cell(struct grid* g, Cell cell, ElementArray<Face> faces)
{
    int ghost = IS_GHOST(cell);
    Cell res = g->mesh->CreateCell(faces).first;
    res.SetStatus(cell.GetStatus());

    res.RealArrayDF(g->c_tags.center)[0] = cell.RealArrayDF(g->c_tags.center)[0];
    res.RealArrayDF(g->c_tags.center)[1] = cell.RealArrayDF(g->c_tags.center)[1];
    res.RealArrayDF(g->c_tags.center)[2] = cell.RealArrayDF(g->c_tags.center)[2];
    
    res.RealArrayDF(g->c_tags.side)[0] = cell.RealArrayDF(g->c_tags.side)[0];
    res.RealArrayDF(g->c_tags.side)[1] = cell.RealArrayDF(g->c_tags.side)[1];
    res.RealArrayDF(g->c_tags.side)[2] = cell.RealArrayDF(g->c_tags.side)[2];

    res.Integer(g->c_tags.floor) = cell.Integer(g->c_tags.floor);
    res.Integer(g->c_tags.chld_num) = cell.Integer(g->c_tags.chld_num);
    int l = res.Integer(g->c_tags.level) = cell.Integer(g->c_tags.level);
    int cr_i = res.Integer(g->c_tags.i) = cell.Integer(g->c_tags.i);

/*
    if (cr_cells.check_i(l,cr_i))
        if (cr_cells.get(l.cr_i) == cell) 
            cr_cells.set1(res,cr_i);
*/
    if (cr_cells.active && l > 0)
        if (cr_cells.check_i(l,cr_i))
        {
            if (cr_cells.get(l,cr_i) == cell) 
                cr_cells.set1(res,cr_i);
            else
                TSNH
        }


    for (int i = 0; i < 8; i++)
        res.IntegerArrayDF(g->c_tags.par_chld_nums)[i] = cell.IntegerArrayDF(g->c_tags.par_chld_nums)[i];

    res.Integer(g->c_tags.to_split) = cell.Integer(g->c_tags.to_split);
    res.Integer(g->c_tags.base_id) = cell.Integer(g->c_tags.base_id);

    g->mesh->Destroy(cell);
    return res;
}

/// Split edge to 2
ElementArray<Edge> Split_Edge(struct grid* g, Edge edge, Cell* cell_for_split, Face* cor_face = NULL)
{
    int status = edge.GetStatus();
    ElementArray<Node> nodes = edge.getNodes();
    if (nodes.size() != 2) TSNH
    double c[3];
    c[0] = (nodes[0].RealArray(g->mesh->CoordsTag())[0] + nodes[1].RealArray(g->mesh->CoordsTag())[0])/2.0;
    c[1] = (nodes[0].RealArray(g->mesh->CoordsTag())[1] + nodes[1].RealArray(g->mesh->CoordsTag())[1])/2.0;
    c[2] = (nodes[0].RealArray(g->mesh->CoordsTag())[2] + nodes[1].RealArray(g->mesh->CoordsTag())[2])/2.0;

    Node center_node = g->mesh->CreateNode(c);
    center_node.SetStatus(edge.GetStatus());
    ElementArray<Node> ne1; ne1.push_back(nodes[0]); ne1.push_back(center_node);
    ElementArray<Node> ne2; ne2.push_back(nodes[1]); ne2.push_back(center_node);

    ElementArray<Edge> new_edges;
    Edge n_edge = g->mesh->CreateEdge(ne1).first;
    n_edge.SetStatus(status);
    new_edges.push_back(n_edge);

    n_edge = g->mesh->CreateEdge(ne2).first;
    n_edge.SetStatus(status);
    new_edges.push_back(n_edge);
    
    // For correct devide, neccesary recreate all adjacent faces and cells
    ElementArray<Face> faces = edge->getFaces();
    for (ElementArray<Face>::iterator face = faces.begin(); face != faces.end(); face++)
    {
        // Recreate face with new node
        ElementArray<Edge>     face_edges = face->getEdges();
        ElementArray<Edge> new_face_edges; 
        for (ElementArray<Edge>::iterator p = face_edges.begin(); p != face_edges.end(); p++)
            if (p->getAsEdge() != edge)
                new_face_edges.push_back(p->getAsEdge());

        new_face_edges.push_back(new_edges[0]); 
        new_face_edges.push_back(new_edges[1]); 
        
        reorder_edges(g,&new_face_edges);

        Face new_face = g->mesh->CreateFace(new_face_edges).first;
        new_face.SetStatus(face->GetStatus());
        if (cor_face && *cor_face == face->getAsFace()) *cor_face = new_face->getAsFace();
        
        // Face created. Now need recreate all cells connected with this face
        ElementArray<Cell> cells = face->getCells();
        for (ElementArray<Cell>::iterator cell = cells.begin(); cell != cells.end(); cell++)
        {
            ElementArray<Face>     cell_faces = cell->getFaces();
            ElementArray<Face> new_cell_faces;
            for (ElementArray<Face>::iterator p = cell_faces.begin(); p != cell_faces.end(); p++)
                if (p->getAsFace() != face->getAsFace())
                    new_cell_faces.push_back(p->getAsFace());
            new_cell_faces.push_back(new_face);

            if (cell->getAsCell() == *cell_for_split)
                *cell_for_split = recreate_cell(g, cell->getAsCell(), new_cell_faces);
            else
                recreate_cell(g, cell->getAsCell(), new_cell_faces);
        }

        g->mesh->Destroy(face->getAsFace());
    }

    g->mesh->Destroy(edge);
    return new_edges;
}

/// Reorder 8 nodes with this order
/// 2 4 7
/// 1   6
/// 0 3 5                                                                
void order_nodes_in_face(struct grid* g, ElementArray<Node>* nodes)
{
    bool swap = false;
    double n1[3], n2[3];
    for (int i = 0; i < nodes->size(); i++)
        for (int j = i+1; j < nodes->size(); j++)
        {
            swap = false;
            n1[2] = (*nodes)[i].RealArray(g->mesh->CoordsTag())[0]; n2[2] = (*nodes)[j].RealArray(g->mesh->CoordsTag())[0];
            n1[1] = (*nodes)[i].RealArray(g->mesh->CoordsTag())[1]; n2[1] = (*nodes)[j].RealArray(g->mesh->CoordsTag())[1];
            n1[0] = (*nodes)[i].RealArray(g->mesh->CoordsTag())[2]; n2[0] = (*nodes)[j].RealArray(g->mesh->CoordsTag())[2];

            if (compare(n2[2], n1[2])  < 0 ) swap = true;
            if (compare(n2[2], n1[2]) == 0 && compare(n2[1], n1[1])  < 0) swap = true;
            if (compare(n2[2], n1[2]) == 0 && compare(n2[1], n1[1]) == 0 && compare(n2[0], n1[0]) < 0) swap = true;
            
            if (swap)
            {
                Node tmp = (*nodes)[i];
                (*nodes)[i] = (*nodes)[j];
                (*nodes)[j] = tmp;
            }
        }
}

double dot(double* a, double* b, int n = 3)
{
    double res = 0;
    for (int i = 0; i < n; i++)
        res += a[i]*b[i];
    return res;
}

void order_nodes_in_face_b(struct grid* g, ElementArray<Node>* nodes)
{
    /*
    int rank = g->mesh->GetProcessorRank();
    for (int i = 0; i < nodes->size() - 1; i++)
    {
        double* nc = (*nodes)[i].RealArray(g->mesh->CoordsTag()).
        double min_dist = dist(nc,(*nodes)[i+1].RealArray(g->mesh->CoordsTag()));
        int min_j = i+1;
        for (int j = i+2; j < nodes->size(); j++)
        {
            double* ac = (*nodes)[j].RealArray(g->mesh->CoordsTag());
            if (dist(nc,ac) < min_dist)
            {
                min_dist = dist(nc,ac);
                min_j = j;
            }
        }
    }
    */
}
/// Reorder nodes with this order
/// 1 2 
/// 0 3                                                                 
void order_nodes_in_face_a(struct grid* g, ElementArray<Node>* nodes)
{

    int rank = g->mesh->GetProcessorRank();
    int I[2] = {0,2};

    if (nodes->size() == 0) return;
    double cen[3] = {0,0,0};
    double xyz[3] = {1,1,1};
    for (int i = 0; i < nodes->size(); i++)
    {
        if (i < nodes->size() - 1)
        {
            for (int j = 0;j < 3;j++)
                if (!equal((*nodes)[i].RealArray(g->mesh->CoordsTag())[j],(*nodes)[i+1].RealArray(g->mesh->CoordsTag())[j])) xyz[j] = 0;
        }
        cen[0] += (*nodes)[i].RealArray(g->mesh->CoordsTag())[0];
        cen[1] += (*nodes)[i].RealArray(g->mesh->CoordsTag())[1];
        cen[2] += (*nodes)[i].RealArray(g->mesh->CoordsTag())[2];
    }
    cen[0] /= nodes->size();
    cen[1] /= nodes->size();
    cen[2] /= nodes->size();
    
    if (xyz[0] == 0 && xyz[1] == 0) { I[0] = 0; I[1] = 1; }
    else if (xyz[1] == 0 && xyz[2] == 0) { I[0] = 1; I[1] = 2; }
    else if (xyz[0] == 0 && xyz[2] == 0) { I[0] = 0; I[1] = 2; }

    double a[3];
    a[0] = (*nodes)[0].RealArray(g->mesh->CoordsTag())[0] - cen[0];
    a[1] = (*nodes)[0].RealArray(g->mesh->CoordsTag())[1] - cen[1];
    a[2] = (*nodes)[0].RealArray(g->mesh->CoordsTag())[2] - cen[2];

    double* alphas = new double[nodes->size()];
    double t[3];
    double cos_v;
    for (int i = 0; i < nodes->size(); i++)
    {
        t[0] = (*nodes)[i].RealArray(g->mesh->CoordsTag())[0] - cen[0];
        t[1] = (*nodes)[i].RealArray(g->mesh->CoordsTag())[1] - cen[1];
        t[2] = (*nodes)[i].RealArray(g->mesh->CoordsTag())[2] - cen[2];

        cos_v = dot(t,a)/(sqrt(dot(t,t))*sqrt(dot(a,a))); 

        alphas[i] = acos(cos_v);
        double D = (t[I[0]]*a[I[1]]) - (t[I[1]]*a[I[0]]);
        if (D > 0) alphas[i] = 2*acos(-1) - alphas[i];
    }

    for (int i = 0; i < nodes->size(); i++)
    for (int j = nodes->size() - 1; j > i; j--)
    {
        if (alphas[j-1] > alphas[j])
        {
            double t = alphas[j-1];
            alphas[j-1] = alphas[j];
            alphas[j] = t;
           Node tmp = (*nodes)[j-1];
           (*nodes)[j-1] = (*nodes)[j];
           (*nodes)[j] = tmp;
        }
    }
    
    return;
}

/// Split face to 4. face don't delete
ElementArray<Face> split_face_to_4(struct grid* g, Face face, Node center_node, int ghost)
{
    // Reorder nodes in following order
    //   2 4 7
    //   1   6
    //   0 3 5
    ElementArray<Node> nodes = face.getNodes();
    
    if (nodes.size() != 8) { cout << nodes.size(); TSNH }
    
    order_nodes_in_face(g, &nodes);
  
    // Now create nodes
    int edges_nodes[4] = { 1,3,4,6 };
    ElementArray<Node> local_nodes;
    ElementArray<Edge> new_edges;
    for (int i = 0; i < 4; i++)
    {
        local_nodes.push_back(center_node);
        local_nodes.push_back(nodes[ edges_nodes[i] ]);
        Edge n_edge = g->mesh->CreateEdge(local_nodes).first;
        if (ghost) n_edge.SetStatus(Element::Ghost);
        new_edges.push_back(n_edge);

        local_nodes.clear();
    }

    // Create faces
    int faces_nodes[4][4] = { {1,2,4} , {4,7,6} , {6,5,3} , {3,0,1} };

    ElementArray<Face> result;
    for (int i = 0; i < 4; i++)
    {
        local_nodes.push_back(center_node);
        local_nodes.push_back(nodes[ faces_nodes[i][0] ]);
        local_nodes.push_back(nodes[ faces_nodes[i][1] ]);
        local_nodes.push_back(nodes[ faces_nodes[i][2] ]);
        result.push_back(g->mesh->CreateFace(local_nodes).first);
        local_nodes.clear();
    }

    return result;
}

/// Split cell to 8. Recreate adjacent cells
ElementArray<Cell> cellSplit_2(struct grid* g, Cell cell) 
{
    int rank = g->mesh->GetProcessorRank();
    int ghost = IS_GHOST(cell);

    ElementArray<Edge> edges = cell.getEdges();
    for (ElementArray<Edge>::iterator edge = edges.begin(); edge != edges.end(); edge++)
    {
        if (!isSmallEdge(g, edge->getAsEdge(),cell))
        {
            ElementArray<Edge> res = Split_Edge(g, edge->getAsEdge(), &cell);
        }
    }

    if (cell.getNodes().size() < 20) TSNH

    ElementArray<Face> faces = cell.getFaces();
    for (ElementArray<Face>::iterator face = faces.begin(); face != faces.end(); face++)
    {
        if (!isSmallFace(g,face->getAsFace(),cell)) 
        {
            double xyz[3] = {0,0,0};
            ElementArray<Node> face_nodes = face->getNodes(); 
            if (face_nodes.size() != 8) { cout << face_nodes.size(); TSNH}
            for (ElementArray<Node>::iterator node = face_nodes.begin(); node != face_nodes.end(); node++)
            {
                xyz[0] += node->RealArray(g->mesh->CoordsTag())[0];
                xyz[1] += node->RealArray(g->mesh->CoordsTag())[1];
                xyz[2] += node->RealArray(g->mesh->CoordsTag())[2];
            }
            xyz[0] /= face_nodes.size();
            xyz[1] /= face_nodes.size();
            xyz[2] /= face_nodes.size();

            Node new_node = g->mesh->CreateNode(xyz);

            int l_ghost = 1;

            ElementArray<Element> adjs = face->getAdjElements(CELL); 
            for (ElementArray<Element>::iterator adj = adjs.begin(); adj != adjs.end(); adj++)
            if (IS_GHOST_P(adj) == 0) l_ghost = 0;

            ElementArray<Face> new_faces = split_face_to_4(g,face->getAsFace(),new_node,l_ghost);
            if (new_faces.size() != 4) TSNH

            ElementArray<Cell> face_cells = face->getCells();
            if (face_cells.size() > 2) TSNH

            for (ElementArray<Cell>::iterator face_cell = face_cells.begin(); face_cell != face_cells.end(); face_cell++)
            {
                ElementArray<Face>     cell_faces = face_cell->getFaces();
                ElementArray<Face> new_cell_faces;
                for (ElementArray<Face>::iterator p = cell_faces.begin(); p != cell_faces.end(); p++)
                    if (p->getAsFace() != face->getAsFace())
                        new_cell_faces.push_back(p->getAsFace()); 
                new_cell_faces.push_back(new_faces[0]);
                new_cell_faces.push_back(new_faces[1]);
                new_cell_faces.push_back(new_faces[2]);
                new_cell_faces.push_back(new_faces[3]);

                if (face_cell->getAsCell() == cell)
                    cell = recreate_cell(g, face_cell->getAsCell() , new_cell_faces);
                else
                    recreate_cell(g, face_cell->getAsCell() , new_cell_faces);
            }
            
            g->mesh->Destroy(face->getAsFace());
        }
    }

    ElementArray<Node> nodes = cell.getNodes();
    if (cell.getFaces().size() != 24) TSNH  
    if (nodes.size()           != 26)
    {
         cout << nodes.size() << " "; 
         TSNH  
    }
        
    double xyz[3];
    xyz[0] = cell.RealArrayDF(g->c_tags.center)[0];
    xyz[1] = cell.RealArrayDF(g->c_tags.center)[1];
    xyz[2] = cell.RealArrayDF(g->c_tags.center)[2];
    g->transformation(xyz); 
    Node center_node = g->mesh->CreateNode(xyz);
    
    order_nodes_in_face(g, &nodes);

    nodes.push_back(center_node); 
    
    int cell_nodes[8][8] = { {0 ,3 ,12,9 ,1 ,4 ,26,10},
                             {3 ,6 ,14,12,4 ,7 ,15,26},
                             {12,14,23,20,26,15,24,21},
                             {9 ,12,20,17,10,26,21,18},
                             {1 ,4 ,26,10,2 ,5 ,13,11},
                             {4 ,7 ,15,26,5 ,8 ,16,13},
                             {26,15,24,21,13,16,25,22},
                             {10,26,21,18,11,13,22,19} };

    // Nodes order for function create cube 
    /*      (5)*-------*(6)  */
	/*     z  /|      /|     */
	/*     ^ /       / |     */
	/*     |/  |    /  |     */
	/*  (4)*-------*(7)|     */
	/*     |   |   |   |     */
	/*     |     y |   |     */
    /*     |   |/  |   |     */
    /*     |(1)*- -|- -*(2)  */
    /*     |  /    |  /      */
    /*     |       | /       */
    /*     |/      |/        */
    /*  (0)*-------*(3) -> x */

    int dirs[8][3] = { {-1,-1,-1} , {-1,1,-1} , {1,1,-1} , {1,-1,-1} , {-1,-1,1} , {-1,1,1} , {1,1,1} , {1,-1,1} };
    double cx = cell.RealArrayDF(g->c_tags.center)[0];
    double cy = cell.RealArrayDF(g->c_tags.center)[1];
    double cz = cell.RealArrayDF(g->c_tags.center)[2];
    double dx = cell.RealArrayDF(g->c_tags.side)[0];
    double dy = cell.RealArrayDF(g->c_tags.side)[1];
    double dz = cell.RealArrayDF(g->c_tags.side)[2];
    double ndx = dx/2.0;
    double ndy = dy/2.0;
    double ndz = dz/2.0;
    int level  = cell.Integer(g->c_tags.level);
    int floor  = cell.Integer(g->c_tags.floor);
    int base_id  = cell.Integer(g->c_tags.base_id);

    ElementArray<Node> local_nodes;
    ElementArray<Cell> result;
    for (int i = 0; i < 8; i++)
    {
        local_nodes.clear();
        for (int j = 0; j < 8; j++)
            local_nodes.push_back(nodes[cell_nodes[i][j]]);

        Cell c = CreateCubeElement(g, local_nodes,ghost);
        result.push_back(c);

        c.RealArrayDF(g->c_tags.center)[0] = cx + ndx*dirs[i][0]/2.0; 
        c.RealArrayDF(g->c_tags.center)[1] = cy + ndy*dirs[i][1]/2.0; 
        c.RealArrayDF(g->c_tags.center)[2] = cz + ndz*dirs[i][2]/2.0; 

        c.RealArrayDF(g->c_tags.side)[0] = ndx;
        c.RealArrayDF(g->c_tags.side)[1] = ndy;
        c.RealArrayDF(g->c_tags.side)[2] = ndz;

        c.Integer(g->c_tags.level) = level + 1;
        c.Integer(g->c_tags.floor) = floor;
        c.Integer(g->c_tags.chld_num) = i;
        c.Integer(g->c_tags.base_id) = base_id;
        bool ready = false;
        for (int k = 0; k < 8; k++)
        {
            c.IntegerArrayDF(g->c_tags.par_chld_nums)[k] = cell.IntegerArrayDF(g->c_tags.par_chld_nums)[k];
            if (!ready && c.IntegerArrayDF(g->c_tags.par_chld_nums)[k] == -1)
            {
                ready = true;
                c.IntegerArrayDF(g->c_tags.par_chld_nums)[k] = cell.Integer(g->c_tags.chld_num);
            }
        }
        c.Integer(g->c_tags.to_split) = 0;
    }

    ElementArray<Element> adjs = center_node.getAdjElements(EDGE); 

    g->mesh->Destroy(cell);
   
    return result;
}

int get_adj_count(Node node)
{
    ElementArray<Element> adjs = node.getAdjElements(CELL); 
    return adjs.size();
}

/// Unite array children to one cell
Cell cell_unite(struct grid* g, ElementArray<Cell> children, Node center)
{
    int rank = g->mesh->GetProcessorRank();
    double dx = children[0].RealArrayDF(g->c_tags.side)[0];
    double dy = children[0].RealArrayDF(g->c_tags.side)[1];
    double dz = children[0].RealArrayDF(g->c_tags.side)[2];
    double ndx = dx*2.0;
    double ndy = dy*2.0;
    double ndz = dz*2.0;
    int level  = children[0].Integer(g->c_tags.level);
    int floor  = children[0].Integer(g->c_tags.floor);
    int base_id  = children[0].Integer(g->c_tags.base_id);
    double xyz[3];
    xyz[0] = center->RealArray(g->mesh->CoordsTag())[0]; 
    xyz[1] = center->RealArray(g->mesh->CoordsTag())[1]; 
    xyz[2] = center->RealArray(g->mesh->CoordsTag())[2]; 
    int par_chld_nums[8];
    for(int i = 0; i < 8; i++)
        par_chld_nums[i] = children[0].IntegerArrayDF(g->c_tags.par_chld_nums)[i];

    /*      (2)*-------*(6)  */
	/*     z  /|      /|     */
	/*     ^ /       / |     */
	/*     |/  |    /  |     */
	/*  (1)*-------*(5)|     */
	/*     |   |   |   |     */
	/*     |     y |   |     */
    /*     |   |/  |   |     */
    /*     |(3)*- -|- -*(7)  */
    /*     |  /    |  /      */
    /*     |       | /       */
    /*     |/      |/        */
    /*  (0)*-------*(4) -> x */

    int verts[8][5] = { {0,1,2,3,-1} , {2,3, 6,7 ,-1} , {4,5,6 ,7 ,-1} , {0,1 ,5 ,4 ,3 } 
                       ,{1,2,6,5,-1} , {2,6,-1,-1,-1} , {6,5,-1,-1,-1} , {5,-1,-1,-1,-1} };

    ElementArray<Node> def_nodes;
    for (ElementArray<Cell>::iterator p = children.begin(); p != children.end(); p++)
    {
        ElementArray<Node> nodes = p->getNodes();
        if (nodes.size() != 8) { cout << rank << "|" << nodes.size(); TSNH }
        reorder_1(g,&nodes);
        int j = 0;
        while (j < 5 && verts[p->Integer(g->c_tags.chld_num)][j] != -1)
            def_nodes.push_back(nodes[verts[p->Integer(g->c_tags.chld_num)][j++]]->getAsNode());
    }
    if (def_nodes.size() != 26) TSNH
    
    order_nodes_in_face(g, &def_nodes);

    // Destroy cells 
    for (ElementArray<Cell>::iterator p = children.begin(); p != children.end(); p++)
    {
        int cr_i = p->Integer(g->c_tags.i);
        int l = p->Integer(g->c_tags.level);
        //cr_cells.dis(cr_i);
        cr_cells.dis(l,cr_i);
        g->mesh->Destroy(p->getAsCell());
    }
    g->mesh->Destroy(center);

    int x_nodes[6] = { 4,15,21,10,12,13 };
    int x_node_face_nodes[6][8] = {{0,1,2,5,8,7,6,3},{6,7,8,16,25,24,23,14},{17,18,19,22,25,24,23,20},{0,1,2,11,19,18,17,9},{0,3,6,14,23,20,17,9},{2,5,8,16,25,22,19,11}}; 
    ElementArray<Face> new_faces;
    for (int i = 0; i < 6; i++)
    {
        int adjs_count = get_adj_count(def_nodes[ x_nodes[i] ]);

        // This case means that adjacent cell splitted to 4. It is not neccessary to create
        if (adjs_count > 1)
        {
            ElementArray<Face> faces = def_nodes[ x_nodes[i] ].getFaces();
            int m = new_faces.size();
            for (ElementArray<Face>::iterator q = faces.begin(); q != faces.end(); q++)
                if (q->getCells().size() == 1)
                    new_faces.push_back(q->getAsFace());
            if ((new_faces.size() - m) != 4 ) TSNH
        }
        
        // This case means that adjacent cell in one cell. 
        if (adjs_count == 1)
        {
            ElementArray<Node> face_nodes;
            for (int j = 0; j < 8; j++)
                face_nodes.push_back(def_nodes[ x_node_face_nodes[i][j]] );
            Face new_face = g->mesh->CreateFace(face_nodes).first;

            ElementArray<Face> faces = def_nodes[ x_nodes[i] ].getFaces();
            if (faces.size() != 4) { cout << faces.size(); TSNH }
            if (faces[0].getCells().size() != 1) { cout << faces[0].getCells().size(); TSNH }

            if (def_nodes[x_nodes[i]].getCells().size() != 1) TSNH
            Cell neight_cell = (def_nodes[x_nodes[i]].getCells())[0];
            
            ElementArray<Face>     neight_faces = neight_cell.getFaces();
            ElementArray<Face> new_neight_faces;
            for (ElementArray<Face>::iterator q = neight_faces.begin(); q != neight_faces.end(); q++)
            {
                if (q->getAsFace() != faces[0] && 
                    q->getAsFace() != faces[1] && 
                    q->getAsFace() != faces[2] && 
                    q->getAsFace() != faces[3]) new_neight_faces.push_back(q->getAsFace());
            }
            new_neight_faces.push_back(new_face);
            if (neight_faces.size() - new_neight_faces.size() != 3) { TSNH }
            recreate_cell(g, neight_cell, new_neight_faces);

            new_faces.push_back(new_face);

            g->mesh->Destroy(faces[0]);
            g->mesh->Destroy(faces[1]);
            g->mesh->Destroy(faces[2]);
            g->mesh->Destroy(faces[3]);
            g->mesh->Destroy(def_nodes[x_nodes[i]]);
        }

        // This case means that cell has't adjacent cell
        if (adjs_count == 0)
        {
            ElementArray<Node> face_nodes;
            for (int j = 0; j < 8; j++)
                face_nodes.push_back(def_nodes[ x_node_face_nodes[i][j]] );
            Face new_face = g->mesh->CreateFace(face_nodes).first;

            new_faces.push_back(new_face);

            ElementArray<Face> faces = def_nodes[ x_nodes[i] ].getFaces();
            if (faces.size() != 4) TSNH

            g->mesh->Destroy(faces[0]);
            g->mesh->Destroy(faces[1]);
            g->mesh->Destroy(faces[2]);
            g->mesh->Destroy(faces[3]);
            g->mesh->Destroy(def_nodes[x_nodes[i]]);
        }
    }

    Cell c = g->mesh->CreateCell(new_faces).first;
    g->rev_transformation(xyz);
    c.RealArrayDF(g->c_tags.center)[0] = xyz[0]; 
    c.RealArrayDF(g->c_tags.center)[1] = xyz[1]; 
    c.RealArrayDF(g->c_tags.center)[2] = xyz[2]; 

    c.RealArrayDF(g->c_tags.side)[0] = ndx;
    c.RealArrayDF(g->c_tags.side)[1] = ndy;
    c.RealArrayDF(g->c_tags.side)[2] = ndz;

    c.Integer(g->c_tags.floor) = floor;
    c.Integer(g->c_tags.level) = level - 1;
    c.Integer(g->c_tags.base_id) = base_id;

    for (int i = 0; i < 8; i++)
    {
        c.IntegerArrayDF(g->c_tags.par_chld_nums)[i] = par_chld_nums[i]; 
    }

    for (int i = 0; i < 8; i++)
    {
        if (c.IntegerArrayDF(g->c_tags.par_chld_nums)[i] == -1)
        {
            if (i == 0) c.Integer(g->c_tags.chld_num) = -1;
            else 
            {
                c.Integer(g->c_tags.chld_num) = c.IntegerArrayDF(g->c_tags.par_chld_nums)[i-1]; 
                c.IntegerArrayDF(g->c_tags.par_chld_nums)[i-1] = -1;
            }
            break;
        }
    }

/*
    int cr_i = cr_cells.size();
    c.Integer(g->c_tags.i) = cr_i;
    cr_cells.push(c);
*/
    
    if (level - 1 > 0)
    {
        int cr_i = cr_cells.size(level - 1);
        c.Integer(g->c_tags.i) = cr_i;
        cr_cells.push(c);
    }

    
    return c;
}

/// Check that cell should split
bool cell_check_mark(struct grid* g, Cell cell)
{
    if (cell->Integer(g->c_tags.to_split) == 1) return false;
    if( g->cell_should_split(g,cell) )
    {
        mark_cell(g,cell);
        return true;
    }
    return false;
}

/// Check that cell should unite
bool cell_check_coarse(struct grid* g, Cell cell)
{
    if (cell->Integer(g->c_tags.level) == 0) return false;
    if( g->cell_should_unite(g,cell) )
    {
        cellCoarse(g,cell);
        return true;
    }
    return false;
}

/// Marks cell for split
void mark_cell(struct grid* g, Cell cell)
{
    // Check adjacent cell. Levels must be coordinated
    ElementArray<Node> nodes = cell.getNodes();
    int neight = false;
    for(ElementArray<Node>::iterator node = nodes.begin(); node != nodes.end(); node++)
    {
        ElementArray<Element> adjs = node->getAdjElements(CELL); 
        for (ElementArray<Element>::iterator adj = adjs.begin(); adj != adjs.end(); adj++)
        {
            if (adj->getAsCell().Integer(g->c_tags.floor) != cell.Integer(g->c_tags.floor))
                continue;

            if (adj->getAsCell().Integer(g->c_tags.level) < cell.Integer(g->c_tags.level))
            {
                neight = true; 
                mark_cell(g, adj->getAsCell());
            }
        }
    }

    if (!neight)
        cell.Integer(g->c_tags.to_split) = 1;
        
    return;
}
    
/// Return node with number num
Node get_corner_node(struct grid* g, Cell cell, int num)
{
    ElementArray<Node> nodes = cell.getNodes();
    Node result = nodes[0];

    int act[8][3] = { {-1,-1,-1} ,
                      {-1,-1, 1} ,
                      {-1, 1, 1} ,
                      {-1, 1,-1} ,
                      { 1,-1,-1} ,
                      { 1,-1, 1} ,
                      { 1, 1, 1} ,
                      { 1, 1,-1} };
    
    double n1[3], n2[3];
    bool swap;
    n1[0] = result.RealArray(g->mesh->CoordsTag())[0];
    n1[1] = result.RealArray(g->mesh->CoordsTag())[1];
    n1[2] = result.RealArray(g->mesh->CoordsTag())[2];
    
    for (int i = 0; i < nodes.size(); i++)
    {
        swap = false;
        n2[0] = nodes[i].RealArray(g->mesh->CoordsTag())[0];
        n2[1] = nodes[i].RealArray(g->mesh->CoordsTag())[1];
        n2[2] = nodes[i].RealArray(g->mesh->CoordsTag())[2];

             if (compare(n2[0],n1[0]) == act[num][0]) swap = true;
        else if (compare(n2[1],n1[1]) == act[num][1]) swap = true;
        else if (compare(n2[2],n1[2]) == act[num][2]) swap = true;

        if (swap)
        {
            result = nodes[i];
            n1[0] = n2[0];
            n1[1] = n2[1];
            n1[2] = n2[2];
        }
    }

    return result;
}

/// Unite cell with brothers
void cellCoarse(struct grid* g, Cell cell)
{
    if (cell.Integer(g->c_tags.level) == 0) return;
    int level = cell.Integer(g->c_tags.level);
    int chld_num = cell.Integer(g->c_tags.chld_num);

    /*      (2)*-------*(6)  */
	/*     z  /|      /|     */
	/*     ^ /       / |     */
	/*     |/  |    /  |     */
	/*  (1)*-------*(5)|     */
	/*     |   |   |   |     */
	/*     |     y |   |     */
    /*     |   |/  |   |     */
    /*     |(3)*- -|- -*(7)  */
    /*     |  /    |  /      */
    /*     |       | /       */
    /*     |/      |/        */
    /*  (0)*-------*(4) -> x */

    int center_node_num[8] = {6,5,1,2,7,4,0,3};                                              
    Node center_node = get_corner_node(g,cell,center_node_num[ chld_num ]);

    int i = 0;
    ElementArray<Cell> children;

    ElementArray<Element> adjs = center_node.getAdjElements(CELL); 
    if (adjs.size() != 8) { cout << adjs.size(); TSNH }
    for (ElementArray<Element>::iterator adj = adjs.begin(); adj != adjs.end(); adj++)
    {
        ElementArray<Node> local_nodes = adj->getAsCell().getNodes();
        for (ElementArray<Node>::iterator p = local_nodes.begin(); p != local_nodes.end(); p++)
        {
            ElementArray<Element> local_adjs = p->getAdjElements(CELL); 
            for (ElementArray<Element>::iterator q = local_adjs.begin(); q != local_adjs.end(); q++)
                if (q->getAsCell().Integer(g->c_tags.level) > adj->getAsCell().Integer(g->c_tags.level))
                {
                    int cr_i = cell.Integer(g->c_tags.i);
                    int l = cell.Integer(g->c_tags.level);
                    cellCoarse(g,q->getAsCell());
                    //if (cr_cells.check(i)) 
                    //{
                    //    cr_cells.push(cr_cells[i]);
                    //    cr_cells.dis(i);
                   // }
                    if (cr_cells.check(l,i)) 
                    {
                        cr_cells.push(cr_cells.get(l,i));
                        cr_cells.dis(l,i);
                    }
                    return;
                }
        }

        children.push_back(adj->getAsCell());
    }

    time_p += Timer() - ttt;
    ttt = Timer();
    Cell c = cell_unite(g,children, center_node);
    time_u += Timer() - ttt;
    ttt = Timer();
    resolve_edges(g,&c);
    time_r1 += Timer() - ttt;
    ttt = Timer();
   
    return;
}

void gridRefine(struct grid * g)
{
    bool all_ready = false;
    bool split_ready = false;
    while (1)
    {
        all_ready = true;
        /// Mark cells
        for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); it++)
        {   
     //       if (it->GetStatus() == Element::Ghost) continue; 
           if (cell_check_mark(g,it->getAsCell())) all_ready = false;
        }

        /// All cells was splitted
        if (all_ready) break; 
        
        /// Split
        while (1)
        {
            split_ready = true;
            for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); it++)
            {
       //         if (it->GetStatus() == Element::Ghost) continue; 
                if (it->getAsCell().Integer(g->c_tags.to_split) == 1)
                {
                    split_ready = false;
                    it->getAsCell().Integer(g->c_tags.to_split) = 0;
                    cellSplit_2(g,it->getAsCell());
                }
            }

            if (split_ready) break;
        }
    }
}

/// Unite 2 edges to 1
void unite_2_edges(struct grid* g, Node node)
{
    ElementArray<Edge> edges = node.getEdges();
    if (edges.size() != 2) TSNH
    ElementArray<Node> new_nodes;

    if (edges[0].getBeg() != node) new_nodes.push_back(edges[0].getBeg());
    if (edges[0].getEnd() != node) new_nodes.push_back(edges[0].getEnd());
    if (edges[1].getBeg() != node) new_nodes.push_back(edges[1].getBeg());
    if (edges[1].getEnd() != node) new_nodes.push_back(edges[1].getEnd());
    if (new_nodes.size() != 2) TSNH

    Edge new_edge = g->mesh->CreateEdge(new_nodes).first;

    ElementArray<Face> faces = edges[0].getFaces();
    for (ElementArray<Face>::iterator face = faces.begin(); face != faces.end(); face++)
    {
        ElementArray<Edge>     face_edges = face->getEdges();
        ElementArray<Edge> new_face_edges;
        for (ElementArray<Edge>::iterator face_edge = face_edges.begin(); face_edge != face_edges.end(); face_edge++)
            if (face_edge->getAsEdge() != edges[0] && face_edge->getAsEdge() != edges[1])
                new_face_edges.push_back(face_edge->getAsEdge());
        
        new_face_edges.push_back(new_edge);
        reorder_edges(g,&new_face_edges);
        Face new_face = g->mesh->CreateFace(new_face_edges).first;

        ElementArray<Cell>     face_cells = face->getCells();
        for (ElementArray<Cell>::iterator face_cell = face_cells.begin(); face_cell != face_cells.end(); face_cell++)
        {
            ElementArray<Face>     cell_faces = face_cell->getFaces();
            ElementArray<Face> new_cell_faces;
            for (ElementArray<Face>::iterator p = cell_faces.begin(); p != cell_faces.end(); p++)
                if (p->getAsFace() != face->getAsFace())
                    new_cell_faces.push_back(p->getAsFace());
            new_cell_faces.push_back(new_face);
            
            recreate_cell(g, face_cell->getAsCell(), new_cell_faces);
        }

        g->mesh->Destroy(face->getAsFace());
    }

    g->mesh->Destroy(edges[0]);
    g->mesh->Destroy(edges[1]);
    g->mesh->Destroy(node);
}

void resolve_edges(grid* g, Cell* c)
{
    if (c)
    {
        ElementArray<Node> nodes = c->getNodes();
        for (int i = 0; i < nodes.size(); i++)
            if (nodes[i].getEdges().size() == 2)
            {
                time_r1 += Timer() - ttt;
                ttt = Timer();
                unite_2_edges(g,nodes[i]);
                time_r2 += Timer() - ttt;
                ttt = Timer();
            }
        return;
    }

    for(Mesh::iteratorNode it = g->mesh->BeginNode(); it != g->mesh->EndNode(); it++)
    {
        ElementArray<Edge> edges = it->getEdges();
        if (edges.size() == 2)
        {
            unite_2_edges(g,it->getAsNode());
        }
    }
}

void gridCoarse(struct grid * g)
{
    int rank = g->mesh->GetProcessorRank();
    int m = 0;
    bool ready;
    cr_cells.clear();
    bool empty = true;
    cr_cells.active = true;
    double tt = Timer();
    int i = 0;
    for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); it++)
    {
        if (it->Integer(g->c_tags.level) == 0) continue;
    //    if (g->cell_should_unite(g,it->getAsCell()))
        {
            empty = false;
            it->Integer(g->c_tags.i) = i;
            cr_cells.push(it->getAsCell());
            i++;
        }
      //  else
    //        it->Integer(g->c_tags.i) = -1;
    }
//    cr_cells.print_size();
//    cout << "N time: " << Timer() - tt << endl;
    if (empty) 
    {
        cr_cells.push(g->mesh->BeginCell()->getAsCell());
        return;
    }

    time_u = 0;
    time_p = 0;
    time_r1 = 0;
    time_r2 = 0;
    ttt = Timer();
    for (int l = MAX_LEVEL-1; l >= 1; l--)
    {
        for (int i = 0; i < cr_cells.size(l); i++)
        {
            if (cr_cells.check(l,i) == false) continue;
            if (!cr_cells.get(l,i).isValid()) continue;
            if (cr_cells.get(l,i).GetStatus() == Element::Ghost) continue; 
        //    cell_check_coarse(g,cr_cells.get(l,i));
            Cell c = cr_cells.get(l,i);
            if (g->cell_should_unite(g,c))
                cellCoarse(g,c);
 
        }
    }
    /*
    for (int i = 0; i < cr_cells.size(); i++)
    {
        if (cr_cells.check(i) == false) continue;
        if (!cr_cells[i].isValid()) continue;
        if (cr_cells[i].GetStatus() == Element::Ghost) continue; 
        cell_check_coarse(g,cr_cells[i]);
    }
    */
    time_p += Timer() - ttt;
    /*
    cout << "Time Prep: " << time_p << endl;
    cout << "Time Unit: " << time_u << endl;
    cout << "Time Res1: " << time_r1 << endl;
    cout << "Time Res2: " << time_r2 << endl;
    */
    /*
    while (1)
    {
        ready = true;
        for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); it++)
        {
            if (it->GetStatus() == Element::Ghost) continue; 
            if (cell_check_coarse(g,it->getAsCell()))      
            {
                ready = false;
                break;
            }
        }

        //resolve_edges(g);
        if (ready) break;
    }
    */

    cr_cells.active = false;
}

void print_redist_tag(struct grid* g,  int rank)
{
    Tag tag_owner = g->mesh->RedistributeTag();
    for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); it++)
    {
        cout << rank << " : " <<  it->getAsCell().Integer(tag_owner) << " : ";
        print_cell_center(g, it->getAsCell());
    }
}


Face get_bottom_or_upper_face(grid* g, Cell c, bool bottom)
{
    ElementArray<Node> nodes = c.getNodes();
    if (nodes.size() == 0) TSNH

    double minmax_z = nodes[0].RealArray(g->mesh->CoordsTag())[2];
    for (int i = 1; i < nodes.size(); i++)
    {
        if (bottom && nodes[i].RealArray(g->mesh->CoordsTag())[2] < minmax_z
        ||  !bottom && nodes[i].RealArray(g->mesh->CoordsTag())[2] > minmax_z)
        {
            minmax_z = nodes[i].RealArray(g->mesh->CoordsTag())[2];
        }
    }

    ElementArray<Face> faces = c.getFaces();
    if (faces.size() < 1) TSNH
    Face result_face = faces[0];
    for (int i = 0; i < faces.size(); i++)
    {
        ElementArray<Node> f_nodes = faces[i].getNodes();
        bool done = true;
        for (int j = 0; j < f_nodes.size(); j++)
        {
            if (!equal(f_nodes[j].RealArray(g->mesh->CoordsTag())[2],minmax_z)) { done = false; break; }
        }
        if (done)
            return faces[i];
    }
    return result_face;
}


void resolve_vertical_conflicts(grid* g)
{
    Tag tag_owner = g->mesh->RedistributeTag();
    for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); it++)
    {
        Face uface = get_bottom_or_upper_face(g,it->getAsCell(),false);
        if (uface.getAdjElements(CELL).size() <= 1) 
        {
            Cell c = it->getAsCell();
            int i = 0;
            do
            {
                Face b_face = get_bottom_or_upper_face(g,c, true);
             
                ElementArray<Element> adjs = b_face.getAdjElements(CELL); 
                if (adjs.size() <= 1) break;
                if (adjs.size() != 2) TSNH
                if (adjs[0]->getAsCell() != c) c = adjs[0]->getAsCell();
                else c = adjs[1].getAsCell();

                c.Integer(tag_owner) = it->Integer(tag_owner);
            } while (i++ < 10000);
        } 
    }
}


void old_correct_brothers(struct grid* g, int size, int rank, int type)
{
    Tag tag_owner = g->mesh->RedistributeTag();

    int* m = new int[size];
    for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); it++)
    {
        if (it->Integer(g->c_tags.level) == 0) continue;
        for (int i = 0; i < size; i++) m[i] = 0;
        int chld_num = it->Integer(g->c_tags.chld_num);
       
        int center_node_num[8] = {6,5,1,2,7,4,0,3};                                              
        Node center_node = get_corner_node(g,it->getAsCell(),center_node_num[ chld_num ]);

        int i = 0;
        ElementArray<Cell> children;

        ElementArray<Element> adjs = center_node.getAdjElements(CELL); 
        if (adjs.size() != 8) { cout << rank << " " << adjs.size(); TSNH }

        int new_owner = 0;
        for (ElementArray<Element>::iterator p = adjs.begin(); p != adjs.end(); p++)
        {
            int owner = p->getAsCell().Integer(tag_owner);
            m[owner]++;
        }

        int max = 0;
        for (int i = 0; i < size; i++)
            if (max < m[i])
            {
                max = m[i];
                new_owner = i;
            }

            for (ElementArray<Element>::iterator p = adjs.begin(); p != adjs.end(); p++)
            {
                if (p->getAsCell().Integer(g->c_tags.level) == 0) continue;
                p->getAsCell().Integer(tag_owner) = new_owner;
            }
    }

    resolve_vertical_conflicts(g);
}

void octree_refine(grid* g, int rank)
{
    Tag tag_owner = g->mesh->RedistributeTag();

    for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); it++)
        it->Integer(tag_owner) = rank;
}

typedef map<int,ElementArray<Cell> > cmap;
void correct_brothers(struct grid* g, int size, int rank, int type)
{
    Tag tag_owner = g->mesh->RedistributeTag();

    cmap cells_map;
    for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); it++)
    {
        (cells_map[it->Integer(g->c_tags.base_id)]).push_back(it->getAsCell());
    }

    int* m = new int[size];
    for (cmap::iterator it = cells_map.begin(); it != cells_map.end(); it++)
    {
        for (int i = 0; i < size; i++) m[i] = 0;
        for (ElementArray<Cell>::iterator p = it->second.begin(); p != it->second.end(); p++)
        {
            int owner = p->getAsCell().Integer(tag_owner);
            m[owner]++;
        }

        int max = 0;
        int new_owner = 0;
        for (int i = 0; i < size; i++)
            if (max < m[i])
            {
                max = m[i];
                new_owner = i;
            }

        for (ElementArray<Cell>::iterator p = it->second.begin(); p != it->second.end(); p++)
            p->getAsCell().Integer(tag_owner) = new_owner;
    }
    resolve_vertical_conflicts(g);
    delete[] m;
}



void resolve_ghost(grid* g)
{
    int rank = g->mesh->GetProcessorRank();
    for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); it++)
    {
        if (IS_GHOST(it->getAsCell())) continue;
        ElementArray<Edge> edges = it->getEdges();
        for (ElementArray<Edge>::iterator e = edges.begin(); e != edges.end(); e++)
        {
            if (!IS_GHOST(e->getAsEdge())) continue;
            e->SetStatus(Element::Shared);
        }
    }
}

void remove_ghost_edges(grid* g)
{
    int rank = g->mesh->GetProcessorRank();
    for(Mesh::iteratorEdge f = g->mesh->BeginEdge(); f != g->mesh->EndEdge(); f++)
    {
        if (IS_GHOST(f->getAsEdge()))
        {
            g->mesh->Destroy(f->getAsEdge());
        }
    }
}

void ghost_to_shared(grid* g)
{
    int rank = g->mesh->GetProcessorRank();
    for(Mesh::iteratorEdge f = g->mesh->BeginEdge(); f != g->mesh->EndEdge(); f++)
    {
        if (IS_GHOST(f->getAsEdge()))
        {
            f->SetStatus(Element::Shared);
        }
    }
}

// Return true if edges in one plane
bool check_edges_plane(grid* g, ElementArray<Edge> edges)
{
    double xyz[3] = {1,1,1};
    ElementArray<Node> nodes;
    
    for (ElementArray<Edge>::iterator e = edges.begin(); e != edges.end(); e++)
    {
        nodes.push_back(e->getBeg());
        nodes.push_back(e->getEnd());
    }
    double c1[3], c2[3];
    for (int i = 0; i < nodes.size() - 1; i++)
    {
        c1[0] = nodes[i].RealArray(g->mesh->CoordsTag())[0];
        c1[1] = nodes[i].RealArray(g->mesh->CoordsTag())[1];
        c1[2] = nodes[i].RealArray(g->mesh->CoordsTag())[2];
        c2[0] = nodes[i+1].RealArray(g->mesh->CoordsTag())[0];
        c2[1] = nodes[i+1].RealArray(g->mesh->CoordsTag())[1];
        c2[2] = nodes[i+1].RealArray(g->mesh->CoordsTag())[2];

        if (!equal(c1[0],c2[0])) xyz[0] = 0;
        if (!equal(c1[1],c2[1])) xyz[1] = 0;
        if (!equal(c1[2],c2[2])) xyz[2] = 0;
    }

    return xyz[0] || xyz[1] || xyz[2];
}

void remove_border(grid* g)
{
    int rank = g->mesh->GetProcessorRank();
//    for(Mesh::iteratorNode it = g->mesh->BeginNode(); it != g->mesh->EndNode(); it++)
    for(Mesh::iteratorFace f = g->mesh->BeginFace(); f != g->mesh->EndFace(); f++)
    {
        if (!f->isValid()) continue;
        if (f->getAdjElements(CELL).size() != 1) continue;
        ElementArray<Node> nodes = f->getNodes();
        ElementArray<Node>::iterator it = nodes.begin();
        ElementArray<Element> adjs;
        bool find = false;
        for (; it != nodes.end(); it++)
        {
            if (it->getAdjElements(EDGE).size() != 4) continue;
            adjs = it->getAdjElements(FACE);
            if (adjs.size() != 4) continue;
            if (!check_edges_plane(g, it->getEdges())) continue;

            bool to_remove = true;
            for (ElementArray<Element>::iterator f = adjs.begin(); f != adjs.end(); f++)
            {
                if (f->getAsFace().getAdjElements(CELL).size() != 1) 
                {
                    to_remove = false;
                    break;
                }
            } 
            if (!to_remove) continue;
            
            find = true;
            break;
        }
        if (!find) continue;

        Cell cell = it->getCells().begin()->getAsCell();

        ElementArray<Face> faces = cell.getFaces();
        ElementArray<Face> new_faces;
        for (ElementArray<Face>::iterator face = faces.begin(); face != faces.end(); face++)
        {
            bool to_add = true;
            for (ElementArray<Element>::iterator f = adjs.begin(); f != adjs.end(); f++)
            {
                if (f->getAsFace() == face->getAsFace()) { to_add = false; break; }
            }
            if (to_add) new_faces.push_back(face->getAsFace());
        }
                
        // new_faces contains faces that should't be unite
        
        // Unite 4 faces to one
        /* OLD VERSION
        ElementArray<Node> new_nodes;
        ElementArray<Node> points_to_delete;

        ElementArray<Element> l_edges = it->getAdjElements(EDGE);
        
        for (ElementArray<Element>::iterator f = l_edges.begin(); f != l_edges.end(); f++)
        {
            if (f->getAsEdge().getBeg() == it->getAsNode()) points_to_delete.push_back(f->getAsEdge().getEnd());
            else points_to_delete.push_back(f->getAsEdge().getBeg());
        }
        for (ElementArray<Node>::iterator p = points_to_delete.begin(); p != points_to_delete.end(); p++)
        {
            if (p->getAdjElements(EDGE).size() > 3) new_nodes.push_back(p->getAsNode());
        }
        points_to_delete.push_back(it->getAsNode());

        // 4 faces
        for (ElementArray<Element>::iterator f = adjs.begin(); f != adjs.end(); f++)
        {
            ElementArray<Node> nodes = f->getNodes();
            // Nodes in face
            for (ElementArray<Node>::iterator node = nodes.begin(); node != nodes.end(); node++)
            {
                if (it->getAsNode() == node->getAsNode()) continue;
                bool to_add = true;
                for (ElementArray<Node>::iterator u = points_to_delete.begin(); u != points_to_delete.end(); u++)
                {
                    if (node->getAsNode() == u->getAsNode()) { to_add = false; break; }
                }
                if (!to_add) continue;
                new_nodes.push_back(node->getAsNode());
            }
        }
        */

        ElementArray<Node> new_nodes;
        for (ElementArray<Element>::iterator f = adjs.begin(); f != adjs.end(); f++)
        {
            ElementArray<Node> nodes = f->getNodes();
            for (ElementArray<Node>::iterator p = nodes.begin(); p != nodes.end(); p++)
                if (p->getAsNode() != it->getAsNode()) 
                {
                    int was = 0;
                    for (int i = 0; i < new_nodes.size(); i++)
                        if (new_nodes[i] == p->getAsNode()) was = 1;

                    if (was == 0)
                        new_nodes.push_back(p->getAsNode());
                }
        }
        order_nodes_in_face_a(g,&new_nodes);
        // Now we have 4 nodes for create new face
        Face new_face = g->mesh->CreateFace(new_nodes).first;
        new_faces.push_back(new_face);
        Cell new_cell = recreate_cell(g,cell,new_faces);
        g->mesh->Destroy(it->getAsNode());
    }
}

void command(grid* g)
{
    int rank = g->mesh->GetProcessorRank();
    resolve_connections(g);
	g->mesh->ResolveShared(); // Resolve duplicate nodes
    return;
}


double dist(double* a, double* b)
{
    return sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]) );
}
double dist(grid* g, Node n1, Node n2)
{
   double c1[3], c2[3];
   c1[0] = n1.RealArray(g->mesh->CoordsTag())[0];
   c1[1] = n1.RealArray(g->mesh->CoordsTag())[1];
   c1[2] = n1.RealArray(g->mesh->CoordsTag())[2];
   c2[0] = n2.RealArray(g->mesh->CoordsTag())[0];
   c2[1] = n2.RealArray(g->mesh->CoordsTag())[1];
   c2[2] = n2.RealArray(g->mesh->CoordsTag())[2];
   return dist(c1,c2);
}

double get_center(grid* g, Node n1, Node n2, double* cen)
{
   double c1[3], c2[3];
   c1[0] = n1.RealArray(g->mesh->CoordsTag())[0];
   c1[1] = n1.RealArray(g->mesh->CoordsTag())[1];
   c1[2] = n1.RealArray(g->mesh->CoordsTag())[2];
   c2[0] = n2.RealArray(g->mesh->CoordsTag())[0];
   c2[1] = n2.RealArray(g->mesh->CoordsTag())[1];
   c2[2] = n2.RealArray(g->mesh->CoordsTag())[2];
   cen[0] = (c1[0]+c2[0])/2.0;
   cen[1] = (c1[1]+c2[1])/2.0;
   cen[2] = (c1[2]+c2[2])/2.0;
}

void find_face_center(grid* g, Face f, double* res)
{
    ElementArray<Node> nodes = f->getNodes();
    double max = 0;
    
    for (int i = 0; i < nodes.size(); i++)
    for (int j = i+1; j < nodes.size(); j++)
    {
        if (dist(g,nodes[i],nodes[j]) > max)
        {
            max = dist(g,nodes[i],nodes[j]);
            get_center(g,nodes[i],nodes[j],res);
        }
    }
}

bool split_if_needed(grid* g, Cell c, Face f)
{
    int rank = g->mesh->GetProcessorRank();
    double cx = c.RealArrayDF(g->c_tags.center)[0]; 
    double cy = c.RealArrayDF(g->c_tags.center)[1];
    double cz = c.RealArrayDF(g->c_tags.center)[2];

    ElementArray<Node> nodes = f->getNodes();
    double fc[3] = {0,0,0};
    double fc_center[3] = {0,0,0};

    find_face_center(g,f,fc);
    fc_center[0] = fc[0];
    fc_center[1] = fc[1];
    fc_center[2] = fc[2];
    g->rev_transformation(fc);

    // Compute symmetrical point etwen c and fc
    fc[0] = 2*fc[0] - cx;
    fc[1] = 2*fc[1] - cy;
    fc[2] = 2*fc[2] - cz;

    c.RealArrayDF(g->c_tags.center)[0] = fc[0];
    c.RealArrayDF(g->c_tags.center)[1] = fc[1];
    c.RealArrayDF(g->c_tags.center)[2] = fc[2];

    int res = g->cell_should_split(g, c );
    c.RealArrayDF(g->c_tags.center)[0] = cx;
    c.RealArrayDF(g->c_tags.center)[1] = cy;
    c.RealArrayDF(g->c_tags.center)[2] = cz;

    if (res == 0) return false;
    // Face must split to 4
    
    ElementArray<Edge> edges = f->getEdges();
    for (ElementArray<Edge>::iterator edge = edges.begin(); edge != edges.end(); edge++)
    {
        if (!isSmallEdge(g, edge->getAsEdge(),c))
        {
            ElementArray<Edge> res = Split_Edge(g, edge->getAsEdge(), &c, &f);
        }
    }

    Node center_node = g->mesh->CreateNode(fc_center);
    ElementArray<Face> new_faces = split_face_to_4(g, f->getAsFace(), center_node, 0);
    if (new_faces.size() != 4) TSNH

    ElementArray<Face>     cell_faces = c.getFaces();
    for (ElementArray<Face>::iterator p = cell_faces.begin(); p != cell_faces.end(); p++)
        if (p->getAsFace() != f)
            new_faces.push_back(p->getAsFace()); 


    c = recreate_cell(g, c , new_faces);
    g->mesh->Destroy(f);
    return true;
}

void resolve_connections(grid* g)
{
    int rank = g->mesh->GetProcessorRank();
    for(Mesh::iteratorCell it = g->mesh->BeginCell(); it != g->mesh->EndCell(); it++)
    {
        if (!it->isValid()) continue;
        ElementArray<Face> faces = it->getFaces();
        for (ElementArray<Face>::iterator face = faces.begin(); face != faces.end(); face++)
        {
            if (!it->isValid()) break;
            if (!face->isValid()) continue;
            if (isSmallFace(g, face->getAsFace(), it->getAsCell())) continue;

            if(face->getAdjElements(CELL).size() == 1)
            {
                ElementArray<Node> nodes = face->getNodes();
                int sh_c = 0;
                for (ElementArray<Node>::iterator node = nodes.begin(); node != nodes.end(); node++)
                    if (node->isValid() && node->GetStatus() != Element::Owned) sh_c++;
                if (sh_c >= 4)
                {
                    split_if_needed(g, it->getAsCell(), face->getAsFace());
                }
            }
        }
    }
}

class IterTime
{
    vector<double> time;
    vector<string> what;
    double tt;
    double all;

    public:
    void start() 
    { 
        BARRIER 
        tt = Timer();  
    }
    void add(string _what)
    {
        double tt1 = Timer();
        double ct = tt1 - tt;
        all += ct;
        time.push_back(ct);
        what.push_back(_what);
        tt = tt1;   
    }
    void print_time(int rank)
    {
        if (rank > 0) return;
        for (int i = 0; i < time.size(); i++)
        {
            double per = (time[i]/all)*100;
            cout << setw(15) << what[i] << " time = " << setw(10) << time[i] << " : ";
            if (per > 70)
                cout <<  "\033[1;31m"  << setprecision(2) << setw(5) << (time[i]/all)*100 << " %" << "\033[0m" << endl;
            else
                cout << setprecision(2) << setw(5) << (time[i]/all)*100 << " %" << endl;
        }
    }
};

int yo = 0;
void gridAMR(struct grid * g, int action)
{
    int rank = g->mesh->GetProcessorRank(); // Get the rank of the current process
    IterTime time;
    time.start();

    remove_border(g);
    time.add("Remove border");

    resolve_edges(g);
    time.add("Resolve edges");

    if (action == 0 || action == 1 || action == 3)
        gridCoarse(g);
    time.add("Coarse");

    //if (action == 0 || action == 2 || action == 3)
    if (action != 3)
        gridRefine(g); 
    time.add("Refine");

    resolve_edges(g);
    time.add("Resolve edges");

    if (action == 4)
    {
        remove_border(g);
        resolve_edges(g);
        resolve_connections(g);
    }
	g->mesh->ResolveShared(); // Resolve duplicate nodes
    time.print_time(rank);
}
