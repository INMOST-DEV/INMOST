#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif

//g++ main.cpp rotate.cpp -L/usr/X11R6/lib -lX11 -lXi -lXmu -lGL -lglut -lGLU ../../INMOST.a -O5
// press space - explode mesh to see connection 
#include "inmost.h"
#include "rotate.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdarg.h>
#include "my_glut.h"
#include <iomanip>
#include "clipboard.h"

inline static unsigned int flip(const unsigned int * fp)
{
	unsigned int mask = -((int)(*fp >> 31)) | 0x80000000;
	return *fp ^ mask;
}
#define _0(x)	(x & 0x7FF)
#define _1(x)	(x >> 11 & 0x7FF)
#define _2(x)	(x >> 22 )

void draw_screen();

using namespace INMOST;
Mesh * mesh;
int interactive = 0;
double zoom = 1;
int width = 800, height = 800;
double sleft = 1e20, sright = -1e20, sbottom = 1e20, stop = -1e20, sfar = -1e20, snear = 1e20;
double shift[3] = {0,0,0};
bool perspective = false;
int drawedges = 0;
bool boundary = true, planecontrol = false, clipupdate = false, bndupdate = true, clipboxupdate = false, draw_volumetric = false;

Mesh::GeomParam table;

#define CLIP_NONE 0
#define CLIP_NODE 1
#define CLIP_FULL 2
#define CLIP_ENDP 3

#define CLIP_FACE_NONE      0
#define CLIP_FACE_INSIDE    1
#define CLIP_FACE_OUTSIDE   2
#define CLIP_FACE_INTERSECT 3

Storage::real p[3] = {0,0,0}, n[3] = {0,0,1};
ElementArray<Element> boundary_faces;
ElementArray<Edge> added_edges;
std::vector<double> harmonic_points, dual_harmonic_points, conormals;




static void GetBox(Element e, Storage::real min[3], Storage::real max[3])
{
  min[0] = min[1] = min[2] = 1.0e20;
  max[0] = max[1] = max[2] = -1.0e20;
	ElementArray<Node> nodes = e->getNodes();
	for (ElementArray<Node>::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		Storage::real_array c = it->Coords();
		for (int i = 0; i < (int)c.size(); i++) 
    {
			if (max[i] < c[i]) max[i] = c[i]; //max
			if (min[i] > c[i]) min[i] = c[i]; //min
		}
	}
  for(int i = 0; i < 3; ++i)
  {
    if( max[i] < min[i] )
    {
      max[i] = 0.0001;
      min[i] = -0.0001;
    }
    else if( max[i] == min[i] )
    {
      max[i] += 0.0001;
      min[i] += -0.0001;
    }
  }
}

double amplitude = 10;
double radius = 25;
char visualization_prompt[8192];
int visualization_prompt_active = 0;
Tag visualization_tag;
ElementType visualization_type;

void printtext(const char * fmt, ... )
{
	
	unsigned int i;
	char stext[4096];
	va_list ap;
	if ( fmt == NULL ) return;
	va_start(ap,fmt);
	vsprintf(stext,fmt,ap);
	va_end(ap);
	for(i=0;i<strlen(stext);i++)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, 
							stext[i]);
	}
}


class coord
{
	double p[3];
public:
	coord() { p[0] = p[1] = p[2] = 0; }
  coord(double xyz[3]) {p[0] = xyz[0]; p[1] = xyz[1]; p[2] = xyz[2];}
	coord(double x, double y, double z) {p[0] = x; p[1] = y; p[2] = z;}
	coord(const coord & other) {p[0] = other.p[0]; p[1] = other.p[1]; p[2] = other.p[2];}
	coord & operator = (coord const & other) {p[0] = other.p[0]; p[1] = other.p[1]; p[2] = other.p[2]; return *this;}
  coord & operator +=(const coord & other) {p[0] += other.p[0]; p[1] += other.p[1]; p[2] += other.p[2]; return *this;}
  coord & operator -=(const coord & other) {p[0] -= other.p[0]; p[1] -= other.p[1]; p[2] -= other.p[2]; return *this;}
  coord & operator *=(const coord & other) 
  {
    double tmp[3] = {p[1]*other.p[2] - p[2]*other.p[1], p[2]*other.p[0] - p[0]*other.p[2], p[0]*other.p[1] - p[1]*other.p[0]};
    p[0] = tmp[0]; p[1] = tmp[1]; p[2] = tmp[2];
    return *this;
  }
  coord & operator *=(double other) { p[0] *= other; p[1] *= other; p[2] *= other; return *this;}
  coord & operator /=(double other) { p[0] /= other; p[1] /= other; p[2] /= other; return *this;}
	coord operator -(const coord & other) const {return coord(p[0]-other.p[0],p[1]-other.p[1],p[2]-other.p[2]);}
	coord operator +(const coord & other) const {return coord(p[0]+other.p[0],p[1]+other.p[1],p[2]+other.p[2]);}
	coord operator *(const coord & other) const {return coord(p[1]*other.p[2] - p[2]*other.p[1],p[2]*other.p[0] - p[0]*other.p[2],p[0]*other.p[1] - p[1]*other.p[0]);}
	coord operator /(double other) const {return coord(p[0]/other,p[1]/other,p[2]/other);}
  coord operator *(double other) const {return coord(p[0]*other,p[1]*other,p[2]*other);}
	double operator ^(const coord & other) const {return p[0]*other.p[0]+p[1]*other.p[1]+p[2]*other.p[2];}
	~coord() {}
  double length() const {return sqrt((*this)^(*this));}
	double & operator [](int i) {return p[i];}
  double * data() {return p;}
};

double abs(const coord & p)
{
  return sqrt(p^p);
}


void get_matrix(const coord & a, const coord & b, double matrix[16])
{
	double d;
	coord z = (b-a)/sqrt((b-a)^(b-a));
	coord y;
	coord x;
	y = coord(z[1],-z[2],0);
	d = sqrt(y^y);
	if( d < 1e-5 )
	{
		y = coord(-z[2],0,z[0]);
		d = sqrt(y^y);
	}
	y = y / d;
	x = y*z;
	x = x / sqrt(x^x);
	y = x*z;
	matrix[0] = x[0];
	matrix[1] = x[1];
	matrix[2] = x[2];
	matrix[3] = 0;
	matrix[4] = y[0];
	matrix[5] = y[1];
	matrix[6] = y[2];
	matrix[7] = 0;
	matrix[8] = z[0];
	matrix[9] = z[1];
	matrix[10] = z[2];
	matrix[11] = 0;
	matrix[12] = 0;
	matrix[13] = 0;
	matrix[14] = 0;
	matrix[15] = 1;
}

GLUquadric * cylqs = NULL;
void drawcylinder(coord a,coord b, double width)
{
	double matrix[16];
	if( cylqs == NULL )
	{
		cylqs = gluNewQuadric();
		gluQuadricNormals(cylqs, GLU_SMOOTH);
		gluQuadricOrientation(cylqs, GLU_OUTSIDE);
		gluQuadricDrawStyle(cylqs,GLU_FILL);//GLU_SILHOUETTE
	}
	glPushMatrix();
	glTranslated(a[0],a[1],a[2]);
	get_matrix(a,b,matrix);
	glMultMatrixd(matrix);
	gluCylinder(cylqs, width,width,sqrt((b-a)^(b-a)), 4, 2);
	glPopMatrix();
}

class Octree : public ElementSet
{
  Tag save_center_tag;
  bool save_quad_tree;
  void SubConstruct(const Tag & child_tag, const Tag & center_tag, HandleType * cells, HandleType * temp, int size, bool quad_tree)
  {
    Storage::real_array center = RealArray(center_tag);
    //create 8 nodes
    Storage::real cell_center[3];
    int offsets[8], sizes[8];
    int dims = 3 - (quad_tree ? 1 : 0);
    int numchildren = (1 << dims);
    for(int k = 0; k < numchildren; ++k)
    {
      offsets[k] = 0;
      sizes[k] = 0;
    }
    for(int r = 0; r < size; ++r)
    {
      Element c = Element(GetMeshLink(),cells[r]);
      c->Centroid(cell_center);
      int child_num = 0;
      for(int k = 0; k < dims; ++k)
      {
        if( cell_center[k] > center[k] )
        {
          int m = 1<<k;
          child_num += m;
        }
      }
      c->IntegerDF(child_tag) = child_num;
      sizes[child_num]++;
    }
    for(int k = 1; k < numchildren; ++k)
    {
      offsets[k] = offsets[k-1]+sizes[k-1];
    }
    for(int k = 0; k < numchildren; ++k)
    {
      std::stringstream name;
      name << GetName() << "chld" << k;
      ElementSet child = GetMeshLink()->CreateSetUnique(name.str()).first;
      Storage::real_array child_center = child->RealArray(center_tag);
      for(int r = 0; r < dims; ++r)
      {
        int l = 1 << r;
        int m = k & l;
        child_center[r] = center[r] + ((m ? 1.0 : -1.0) * center[r+3] * 0.25);
        child_center[r+3] = center[r+3]*0.5;
      }
      int m = 0;
      for(int r = 0; r < size; ++r)
      {
        Element c = Element(GetMeshLink(),cells[r]);
        int q = c->IntegerDF(child_tag);
        if( q == k ) (temp+offsets[k])[m++] = cells[r];
      }
      AddChild(child);
      if( sizes[k] <= 16 && sizes[k] > 0 )
        child->PutElements(temp+offsets[k],sizes[k]);
    }
    // cells array is not needed anymore
    ElementSet child = GetChild();
    for(int k = 0; k < numchildren; ++k)
    {
      if( sizes[k] > 16 )
        Octree(child).SubConstruct(child_tag,center_tag,temp+offsets[k],cells+offsets[k],sizes[k], quad_tree);
      child = child->GetSibling();
    }
  }
  Cell SubFindCell(const Tag & center_tag, Storage::real pnt[3], bool quad_tree) const
  {
    if( HaveChild() )
    {
      Storage::real_array center = RealArray(center_tag);
      int child_num = 0, q;
      int dims = 3 - (quad_tree ? 1 : 0);
      for(int k = 0; k < dims; ++k)
      {
        if( pnt[k] > center[k] )
          child_num += (1 << k);
      }
      q = 0;
      ElementSet set = GetChild();
      while(q != child_num) {set = set->GetSibling(); q++;}
      return Octree(set).SubFindCell(center_tag,pnt,quad_tree);
    }
    else
    {
      HandleType * cells = getHandles();
      int ncells = (int)nbHandles();
      Node closest = InvalidNode();
      Storage::real mindist = 1.0e20, dist;
      for(int k = 0; k < ncells; ++k)
      {
        Node c = Node(GetMeshLink(),cells[k]);
        Storage::real_array cnt = c->Coords();
        dist = sqrt((cnt[0]-pnt[0])*(cnt[0]-pnt[0])+(cnt[1]-pnt[1])*(cnt[1]-pnt[1])+(cnt[2]-pnt[2])*(cnt[2]-pnt[2]));
        if( mindist > dist )
        {
          mindist = dist;
          closest = c;
        }
      }
      if( closest.isValid() )
      {
        ElementArray<Cell> cells = closest->getCells();
        for(ElementArray<Cell>::iterator c = cells.begin(); c != cells.end(); ++c)
          if( c->Inside(pnt) ) return c->self();
      }
      return InvalidCell();
    }
  }
  bool Inside(const Storage::real_array & center, Storage::real pnt[3], bool quad_tree) const
  {
    bool inside = true;
    int dims = 3 - (quad_tree ? 1 : 0);
    for(int i = 0; i < dims; ++i) 
      inside &= (pnt[i] >= center[i] - center[3+i]*0.5 && pnt[i] <= center[i] + center[3+i]*0.5);
    return inside;
  }
  Node SubFindNode(const Tag & center_tag, Storage::real pnt[3], bool quad_tree) const
  {
    if( HaveChild() )
    {
      Storage::real_array center = RealArray(center_tag);
      if( !Inside(center,pnt,quad_tree) ) return InvalidNode();
      int child_num = 0, q;
      int dims = 3 - (quad_tree ? 1 : 0);
      for(int k = 0; k < dims; ++k)
      {
        if( pnt[k] > center[k] )
          child_num += (1 << k);
      }
      q = 0;
      ElementSet set = GetChild();
      while(q != child_num) {set = set->GetSibling(); q++;}
      return Octree(set).SubFindNode(center_tag,pnt,quad_tree);
    }
    else
    {
      HandleType * cells = getHandles();
      int ncells = (int)nbHandles();
      Node closest = InvalidNode();
      Storage::real mindist = 1.0e20, dist;
      for(int k = 0; k < ncells; ++k)
      {
        Node c = Node(GetMeshLink(),cells[k]);
        Storage::real_array cnt = c->Coords();
        dist = sqrt((cnt[0]-pnt[0])*(cnt[0]-pnt[0])+(cnt[1]-pnt[1])*(cnt[1]-pnt[1])+(cnt[2]-pnt[2])*(cnt[2]-pnt[2]));
        if( mindist > dist )
        {
          mindist = dist;
          closest = c;
        }
      }
      return closest;
    }
  }
  void SubDestroy()
  {
    if( HaveChild() )
    {
      ElementSet set = GetChild(), next;
      while(set->isValid())
      {
        next = set->GetSibling();
        Octree(set).SubDestroy();
        set = next;
      }
    }
    DeleteSet();
    handle = InvalidHandle();
    handle_link = NULL;
  }
public:
  Octree() : ElementSet(InvalidElementSet()) {}
  Octree(const Octree & other) : ElementSet(other) {}
  Octree(const ElementSet & eset) : ElementSet(eset) {}
  void Construct(ElementType elem, bool quad_tree = false)
  {
    save_quad_tree = quad_tree;
    int dims = 3 - (quad_tree ? 1 : 0);
    Tag child_tag = GetMeshLink()->CreateTag("OCTREE_CHILD_NUM_"+GetName(),DATA_INTEGER,elem,NONE,1);
    save_center_tag = GetMeshLink()->CreateTag("OCTREE_CENTER_"+GetName(),DATA_REAL,ESET,ESET,6);
    Storage::real bounds[3][2];
    for(int k = 0; k < dims; ++k)
    {
      bounds[k][0] = 1.0e20;
      bounds[k][1] =-1.0e20;
    }
    //calculate bounds
    for(Mesh::iteratorNode node = GetMeshLink()->BeginNode(); node != GetMeshLink()->EndNode(); ++node)
    {
      Storage::real_array coord = node->Coords();
      for(int k = 0; k < dims; ++k)
      {
        if( coord[k] < bounds[k][0] ) bounds[k][0] = coord[k];
        if( coord[k] > bounds[k][1] ) bounds[k][1] = coord[k];
      }
    }
    Storage::real_array center_data = RealArray(save_center_tag);
    for(int k = 0; k < dims; ++k)
    {
      center_data[k] = (bounds[k][0]+bounds[k][1])*0.5; //central position
      center_data[k+3] = bounds[k][1]-bounds[k][0]; //length
    }
    //copy cells
    int size = GetMeshLink()->NumberOf(elem), k = 0;
    HandleType * cells = new HandleType[size*2];
    HandleType * temp = cells+size;
    for(Mesh::iteratorElement cell = GetMeshLink()->BeginElement(elem); cell != GetMeshLink()->EndElement(); ++cell)
      cells[k++] = *cell;
    SubConstruct(child_tag,save_center_tag,cells,temp,size,quad_tree);
    GetMeshLink()->DeleteTag(child_tag);
  }
  Cell FindCell(Storage::real pnt[3]) const
  {
    return SubFindCell(save_center_tag,pnt,save_quad_tree);
  }
  Node FindNode(Storage::real pnt[3]) const
  {
    return SubFindNode(save_center_tag,pnt,save_quad_tree);
  }
  void Destroy()
  {
    if( save_center_tag.isValid() )
      GetMeshLink()->DeleteTag(save_center_tag);
    SubDestroy();
  }
  ~Octree() { }
};



void GetVelocity(Cell c, const Tag & velocity_tag, coord pnt, coord & ret)
{
  ElementArray<Cell> adj = c->NeighbouringCells();
  adj.push_back(c);
  coord cnt;
  const Storage::real eps = 1.0e-8;
  Storage::real dist = 0;
  ret[0] = ret[1] = ret[2] = 0;
  for(ElementArray<Cell>::iterator it = adj.begin(); it != adj.end(); ++it)
  {
    it->Centroid(cnt.data());
    coord vel = coord(it->RealArray(velocity_tag).data());
    Storage::real l = (cnt-pnt).length() + eps;
    Storage::real omega = 1.0/(l*l);
    ret += vel*omega;
    dist += omega;
  }
  ret /= dist;
}

void GetVelocity(Node c, const Tag & velocity_tag, coord pnt, coord & ret)
{
  ElementArray<Cell> adj = c->getCells();
  coord cnt;
  const Storage::real eps = 1.0e-8;
  Storage::real dist = 0;
  ret[0] = ret[1] = ret[2] = 0;
  for(ElementArray<Cell>::iterator it = adj.begin(); it != adj.end(); ++it)
  {
    it->Centroid(cnt.data());
    coord vel = coord(it->RealArray(velocity_tag).data());
    Storage::real l = (cnt-pnt).length() + eps;
    Storage::real omega = 1.0/(l*l);
    ret += vel*omega;
    dist += omega;
  }
  ret /= dist;
}


Storage::real GetSize(Cell c)
{
  Storage::real bounds[3][2] = {{1.0e20,-1.0e20},{1.0e20,-1.0e20},{1.0e20,-1.0e20}};
  ElementArray<Node> nodes = c->getNodes();
  for(ElementArray<Node>::iterator n = nodes.begin(); n != nodes.end(); ++n)
  {
    Storage::real_array cnt = n->Coords();
    for(int k = 0; k < 3; ++k)
    {
      if( bounds[k][0] > cnt[k] ) bounds[k][0] = cnt[k];
      if( bounds[k][1] < cnt[k] ) bounds[k][1] = cnt[k];
    }
  }
  Storage::real ret = bounds[0][1]-bounds[0][0];
  ret = std::min(ret,bounds[1][1]-bounds[1][0]);
  ret = std::min(ret,bounds[2][1]-bounds[2][0]);
  return ret;
}


Storage::real GetSize(Node n, const Tag & size_tag)
{
  ElementArray<Cell> cells = n->getCells();
  Storage::real minsize = 1.0e+20, size;
  for(ElementArray<Cell>::iterator c = cells.begin(); c != cells.end(); ++c)
  {
    size = c->RealDF(size_tag);
    if( minsize > size ) minsize = size;
  }
  if( minsize > 1.0e+19 ) std::cout << __FILE__ << ":" << __LINE__ << " oops" << std::endl;
  return minsize;
}


class Streamline
{
private:
	std::vector<coord> points;
	std::vector<double> velarr;
public:
  Streamline() {}
	Streamline(const Octree & octsearch, coord pos, Tag velocity_tag, Tag cell_size, Storage::real velocity_min, Storage::real velocity_max, Storage::real sign, MarkerType visited)
	{
		Storage::real coef, len, size;
    coord next = pos, vel;
    Node c;
    const int maxsteps = 4000;
    points.reserve(maxsteps/2);
    velarr.reserve(maxsteps/2);
		points.push_back(pos);
		velarr.push_back(0);
	  while( points.size() < maxsteps )
		{
      c = octsearch.FindNode(next.data());
      if( !c.isValid() ) break;
      //c.SetMarker(visited);
      GetVelocity(c,velocity_tag,next,vel);
      len = vel.length();
      if( len < 1.0e-4 ) break;
      size = GetSize(c,cell_size);// c->RealDF(cell_size);
      coef = 0.35*size/len;
      next += vel*coef*sign;
      points.push_back(next);
      velarr.push_back((log(len+1.0e-25)-velocity_min)/(velocity_max-velocity_min));
		}
		//printf("%ld %ld\n",points.size(),velarr.size());
	}
	Streamline(const Streamline & other) { points = other.points; velarr = other.velarr; }
	Streamline & operator =(Streamline const & other) {points = other.points; velarr = other.velarr; return *this;}
	~Streamline() { points.clear(); velarr.clear(); }
	void Draw(int reduced)
	{
		
		if( reduced )
		{
			glBegin(GL_LINE_STRIP);
			for(unsigned int i = 0; i < points.size()-1; i++)
			{
				glColor3f(velarr[i+1]*0.65,0.65*(velarr[i+1] < 0.5 ? velarr[i] : 1.0-velarr[i]),0.65*(1-velarr[i+1]));
				glVertex3d(points[i][0],points[i][1],points[i][2]);
			}
			glEnd();
		}
		else for(unsigned int i = 0; i < points.size()-1; i++)
		{
			glColor3f(velarr[i+1]*0.65,0.65*(velarr[i+1] < 0.5 ? velarr[i] : 1.0-velarr[i]),0.65*(1-velarr[i+1]));
			drawcylinder(points[i],points[i+1],0.25*abs(points[i+1]-points[i]));
		}
		
	}
};

std::vector<Streamline> streamlines;


const int name_width = 32;
const int type_width = 14;
const int elems_width = 10;
const int sparse_width = 10;
const int length_width = 10;

void PrintTag(Tag t)
{
	std::cout << std::setw(name_width) << t.GetTagName() << std::setw(type_width) << DataTypeName(t.GetDataType());
	int num = 0;
	char elems[7] = "NEFCSM";
	std::string print = "";
	for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
	{
		if( t.isDefined(etype) )
		{
			print = print + elems[ElementNum(etype)];
			num++;
		}
	}
	std::cout << std::setw(elems_width) << print;
	print = "";
	for(ElementType etype = NODE; etype <= MESH; etype = etype << 1)
	{
		if( t.isSparse(etype) )
		{
			print = print + elems[ElementNum(etype)];
			num++;
		}
	}
	std::cout << std::setw(sparse_width) << print;
	if( t.GetSize() == ENUMUNDEF )
		std::cout << std::setw(length_width) << "VAR" << std::endl;
	else
		std::cout << std::setw(length_width) << t.GetSize() << std::endl;
}

void PrintTags(Mesh * m, ElementType etypes)
{
	std::cout << std::setw(name_width) << "Name" << std::setw(type_width) << "Type" << std::setw(elems_width) << "Element" << std::setw(sparse_width) << "Sparse" << std::setw(length_width) << "Length" << std::endl;
	for(Mesh::iteratorTag t = m->BeginTag(); t != m->EndTag(); ++t )
	{
		bool print = false;
		for(ElementType etype = NODE; etype <= MESH; etype = etype << 1) if( (etype&etypes) && t->isDefined(etype) ) {print = true; break;}
		if( print ) PrintTag(*t);
	}
}

bool write_tga( const char *filename, int w, int h, char *buffer ) 
{ 
  FILE *f = fopen( filename, "wb" ); 
  if ( !f ) 
    return false;
  putc(0,f);
  putc(0,f);
  putc(2,f);
  putc(0,f);putc(0,f);
  putc(0,f);putc(0,f);
  putc(0,f);
  putc(0,f);putc(0,f);
  putc(0,f);putc(0,f);
  putc((w & 0x00ff),f);
  putc((w & 0xff00)/256,f);
  putc((h & 0x00ff),f);
  putc((h & 0xff00)/256,f);
  putc(24,f);
  putc(0,f);
  size_t buflen = w * h * 3; 
  fwrite( buffer, 1, buflen, f ); 
  fclose(f); 
  return true; 
}

void screenshot()
{
  const int tiles = 8;
  int oldwidth = width;
  int oldheight = height;
  width *= tiles;
  height *= tiles;
  
  char * pixelbuffer = new char[width*height*3+oldwidth*oldheight*3];
  char * tempbuffer = pixelbuffer + width*height*3;

  
  
  
  for(int i = 0; i < tiles; ++i)
  {
    for(int j = 0; j < tiles; ++j )
    {
      glViewport(-oldwidth*i,-oldheight*j,width,height);
      draw_screen();
      glReadBuffer(GL_BACK);
      glReadPixels(0,0,oldwidth,oldheight,GL_BGR_EXT,GL_UNSIGNED_BYTE,tempbuffer);

      int koff = oldwidth*(i);
      int loff = oldheight*(j);
      
      for(int l = 0; l < oldheight; ++l)
      for(int k = 0; k < oldwidth; ++k)
      for(int m = 0; m < 3; ++m)
        pixelbuffer[((koff+k) + (loff+l)*width)*3+m] = tempbuffer[(k + l*oldwidth)*3+m];

      //filename[0] += i;
      //filename[1] += j;
      //write_tga(filename,width,height,pixelbuffer);
      //filename[0] = filename[1] = '0';
    }
  }
  
  /*
  glViewport(-oldwidth*3,-oldheight*3,width,height);
  draw_screen();
  glReadBuffer(GL_BACK);
  glReadPixels(0,0,width,height,GL_BGR_EXT,GL_UNSIGNED_BYTE,pixelbuffer);
  */


  
  
  
  write_tga("screenshot.tga",width,height,pixelbuffer);
  delete [] pixelbuffer;
  width = oldwidth;
  height = oldheight;
  glViewport(0,0,width,height);
}

struct color_t
{
	float c[4];
	color_t() {memset(c,0,sizeof(float)*4);}
	color_t(float r, float g, float b)
	{
		c[0] = r;
		c[1] = g;
		c[2] = b;
		c[3] = 1.0;
	}
	color_t(float r, float g, float b, float a)
	{
		c[0] = r;
		c[1] = g;
		c[2] = b;
		c[3] = a;
	}
	color_t(const color_t & other) {memcpy(c,other.c,sizeof(float)*4);}
	color_t & operator =(color_t const & other)
	{
		memmove(c,other.c,sizeof(float)*4);
		return *this;
	}
	void set_color() const {glColor4fv(c);}
	float & r() {return c[0];}
	float & g() {return c[1];}
	float & b() {return c[2];}
	float & a() {return c[3];}
	float r() const {return c[0];}
	float g() const {return c[1];}
	float b() const {return c[2];}
	float a() const {return c[3];}
	color_t operator *(float mult){return color_t(c[0]*mult,c[1]*mult,c[2]*mult,c[3]*mult);}
	color_t operator +(color_t other){return color_t(c[0]+other.c[0],c[1]+other.c[1],c[2]+other.c[2],c[3]+other.c[3]);}
	color_t operator -(color_t other){return color_t(c[0]-other.c[0],c[1]-other.c[1],c[2]-other.c[2],other.c[3]);}
};

class color_bar
{
	float min, max;
	std::vector<float> ticks; //ticks from 0 to 1 for each color
	std::vector<color_t> colors; //4 floats for each tick
	std::string comment;
	unsigned texture;
public:
	color_bar()
	{
		min = 0;
		max = 1;
		comment = "";

		/*
		ticks.push_back(0.f);
		ticks.push_back(0.2f);
		ticks.push_back(0.4f);
		ticks.push_back(0.6f);
		ticks.push_back(0.8f);
		ticks.push_back(1.f);

		 
		//colors.push_back(color_t(1,0,0));
		//colors.push_back(color_t(1,1,0));
		//colors.push_back(color_t(0,1,0));
		//colors.push_back(color_t(0,1,1));
		//colors.push_back(color_t(0,0,1));
		//colors.push_back(color_t(1,0,1));
		
		colors.push_back(color_t(1,0,1));
		colors.push_back(color_t(0,0,1));
		colors.push_back(color_t(0,1,1));
		colors.push_back(color_t(0,1,0));
		colors.push_back(color_t(1,1,0));
		colors.push_back(color_t(1,0,0));
		*/

		//inversed gnuplot color scheme
		ticks.push_back(0.f);
		ticks.push_back(0.05f);
		ticks.push_back(0.5f);
		ticks.push_back(0.75f);
		ticks.push_back(0.95f);
		ticks.push_back(1.f);

		colors.push_back(color_t(1,1,1));
		colors.push_back(color_t(1,1,0));
		colors.push_back(color_t(0.85,0,0));
		colors.push_back(color_t(0.65,0.25,0.85));
		colors.push_back(color_t(0.45,0,0.55));
		colors.push_back(color_t(0,0,0));

		float * pixel_array = new float[1026*4];

		
		
		for(int q = 0; q < 1026; ++q)
		{
			float t = 1.0f*q/1025.f;
			color_t c = pick_color(t);
			pixel_array[(q)*4+0] = c.r();
			pixel_array[(q)*4+1] = c.g();
			pixel_array[(q)*4+2] = c.b();
			pixel_array[(q)*4+3] = c.a();
		}

		pixel_array[0] = 0;
		pixel_array[1] = 1;
		pixel_array[2] = 0;
		pixel_array[3] = 1;

		pixel_array[1025*4+0] = 0;
		pixel_array[1025*4+1] = 1;
		pixel_array[1025*4+2] = 0;
		pixel_array[1025*4+3] = 1;

		
		
		glEnable(GL_TEXTURE);
		glEnable(GL_TEXTURE_1D);
		glGenTextures(1,&texture);
		glBindTexture(GL_TEXTURE_1D,texture);
		glPixelStorei(GL_UNPACK_ALIGNMENT,1);
		
		glTexImage1D(GL_TEXTURE_1D,0,4,1026,1,GL_RGBA,GL_FLOAT,pixel_array);

		std::cout << "Created texture " << texture << std::endl;
		
		glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_WRAP_S,GL_CLAMP);
		glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_1D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);

		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

		UnbindTexture();

		delete [] pixel_array;
		//colors.push_back(color_t(1,0,0));
	}
	void set_comment(std::string text) {comment = text;}
	void set_min(float newmin) { min = newmin;}
	void set_max(float newmax) { max = newmax;}
	color_t pick_color(float value) 
	{
		float t = (value-min)/(max-min);
		std::vector<float>::iterator it = std::lower_bound(ticks.begin(),ticks.end(),t);
		size_t pos = it-ticks.begin();
		if( it == ticks.end() || pos >= ticks.size() ) 
		{
			return colors.back();
		}
		if( pos == 0 ) 
		{
			return colors[0];
		}
		float interp = (t-ticks[pos-1])/(ticks[pos]-ticks[pos-1]);
		return (colors[pos]*interp+colors[pos-1]*(1-interp));
	}
	void BindTexture()
	{
		//glDisable( GL_TEXTURE_GEN_S ); 
		glDisable(GL_TEXTURE_2D);
		glEnable( GL_TEXTURE_1D );
		glBindTexture(GL_TEXTURE_1D, texture);
	}
	void UnbindTexture()
	{
		glDisable( GL_TEXTURE_1D );
	}
	double pick_texture(double value)
	{
		const double eps = 1.0/1024.0;
		return (value-min)/(max-min)*(1-2*eps) + eps;
		//return std::max(std::min((value-min)/(max-min),0.99),0.01);
	}
	void Draw()
	{
		float text_pos = -0.89;
		float left = -0.95;
		float right = -0.9;
		float bottom = -0.75;
		float top = 0.75;
		BindTexture();
		glBegin(GL_QUADS);

		glTexCoord1d(1.0/1024.0);
		glVertex2f(left,bottom);
		glVertex2f(right,bottom);
		glTexCoord1d(1.0);
		glVertex2f(right,top);
		glVertex2f(left,top);
		/*
		for(int i = 0; i < ticks.size()-1; ++i)
		{
			colors[i].set_color();
			glVertex2f(left,bottom+ticks[i]*(top-bottom));
			glVertex2f(right,bottom+ticks[i]*(top-bottom));
			colors[i+1].set_color();
			glVertex2f(right,bottom+ticks[(i+1)]*(top-bottom));
			glVertex2f(left,bottom+ticks[(i+1)]*(top-bottom));
		}
		*/
		glEnd();
		UnbindTexture();
		
		glColor4f(0,0,0,1);
		for(int i = 0; i < ticks.size(); ++i)
		{
			glRasterPos2f(text_pos,bottom+ticks[i]*(top-bottom));
			printtext("%f",min+ticks[i]*(max-min));
		}
		if( comment != "")
		{
			glRasterPos2f(left,bottom-0.04);
			printtext("%s",comment.c_str());
		}

		glBegin(GL_LINE_LOOP);
		glVertex2f(left,bottom);
		glVertex2f(right,bottom);
		glVertex2f(right,top);
		glVertex2f(left,top);
		glEnd();

		glBegin(GL_LINES);
		for(int i = 0; i < ticks.size(); ++i)
		{
			float pos = bottom+ticks[i]*(top-bottom);
			glVertex2f(left,pos);
			glVertex2f(left+(right-left)*0.2,pos);

			glVertex2f(right+(left-right)*0.25,pos);
			glVertex2f(right,pos);
		}
		glEnd();
	}
} * CommonColorBar;


class face2gl
{
	float dist;
	double c[4];
	double cnt[3];
	bool flag;
	ElementType etype;
	Storage::integer id;
	std::vector<double> verts;
	std::vector<double> texcoords;
	std::vector<color_t> colors;
	color_t cntcolor;
	double cnttexcoord;
public:
	static void radix_sort_dist(std::vector<face2gl> & set)
	{
		static std::vector<face2gl> tmp;
		tmp.resize(set.size());
		unsigned int i;
		const unsigned int kHist = 2048;
		unsigned int  b0[kHist * 3];
		unsigned int *b1 = b0 + kHist;
		unsigned int *b2 = b1 + kHist;
		memset(b0,0,sizeof(unsigned int)*kHist*3);
		for (i = 0; i < set.size(); i++) 
		{
			unsigned int fi = flip((unsigned int *)&set[i].dist);
			++b0[_0(fi)]; ++b1[_1(fi)]; ++b2[_2(fi)];
		}
		{
			unsigned int sum0 = 0, sum1 = 0, sum2 = 0;
			for (i = 0; i < kHist; i++) 
			{
				b0[kHist-1] = b0[i] + sum0; b0[i] = sum0 - 1; sum0 = b0[kHist-1];
				b1[kHist-1] = b1[i] + sum1; b1[i] = sum1 - 1; sum1 = b1[kHist-1];
				b2[kHist-1] = b2[i] + sum2; b2[i] = sum2 - 1; sum2 = b2[kHist-1];
			}
		}
		for (i = 0; i < set.size(); i++) tmp[++b0[_0(flip((unsigned int *)&set[i].dist))]] = set[i];
		for (i = 0; i < set.size(); i++) set[++b1[_1(flip((unsigned int *)&tmp[i].dist))]] = tmp[i];
		for (i = 0; i < set.size(); i++) tmp[++b2[_2(flip((unsigned int *)&set[i].dist))]] = set[i];
		for (i = 0; i < set.size(); i++) set[i] = tmp[set.size()-1-i];
	}
	face2gl():verts(),colors(),texcoords() {etype = NONE; id = 0; dist = 0; flag = false; memset(c,0,sizeof(double)*4);}
	face2gl(const face2gl & other) :verts(other.verts) 
	{
		etype = other.etype;
		id = other.id;
		dist = other.dist; 
		cnt[0] = other.cnt[0];
		cnt[1] = other.cnt[1];
		cnt[2] = other.cnt[2];
		c[0] = other.c[0];
		c[1] = other.c[1];
		c[2] = other.c[2];
		c[3] = other.c[3];
		flag = other.flag;
		colors = other.colors;
		texcoords = other.texcoords;
		cntcolor = other.cntcolor;
		cnttexcoord = other.cnttexcoord;
	}
	face2gl & operator =(face2gl const & other) 
	{ 
		etype = other.etype;
		id = other.id;
		verts = other.verts; 
		dist = other.dist; 
		cnt[0] = other.cnt[0];
		cnt[1] = other.cnt[1];
		cnt[2] = other.cnt[2];
		c[0] = other.c[0];
		c[1] = other.c[1];
		c[2] = other.c[2];
		c[3] = other.c[3];
		flag = other.flag;
		colors = other.colors;
		cntcolor = other.cntcolor;
		texcoords = other.texcoords;
		cnttexcoord = other.cnttexcoord;
		return *this;
	}
	~face2gl() {}
	void draw_colour() const
	{
		if( texcoords.empty() )
		{
			if( colors.empty() )
			{
				glColor4dv(c); 
				for(unsigned k = 0; k < verts.size(); k+=3) 
				{
					glVertex3dv(cnt);
					glVertex3dv(&verts[k]);
					glVertex3dv(&verts[(k+3)%verts.size()]);
				}
			}
			else
			{
				for(unsigned k = 0; k < verts.size(); k+=3) 
				{
					cntcolor.set_color();
					glVertex3dv(cnt);
					colors[k/3].set_color();
					glVertex3dv(&verts[k]);
					colors[(k/3+1)%colors.size()].set_color();
					glVertex3dv(&verts[(k+3)%verts.size()]);
				}
			}
		}
		else
		{
			for(unsigned k = 0; k < verts.size(); k+=3) 
			{
				glTexCoord1d(cnttexcoord);
				glVertex3dv(cnt);
				glTexCoord1d(texcoords[k/3]);
				glVertex3dv(&verts[k]);
				glTexCoord1d(texcoords[(k/3+1)%texcoords.size()]);
				glVertex3dv(&verts[(k+3)%verts.size()]);
			}
		}
	}
	void draw_colour_alpha(double alpha) const
	{
		if( texcoords.empty() )
		{
			if( colors.empty() )
			{
				//double cc[4] = {c[0],c[1],c[2],alpha};
				glColor4dv(c); 
				for(unsigned k = 0; k < verts.size(); k+=3) 
				{
					glVertex3dv(cnt);
					glVertex3dv(&verts[k]);
					glVertex3dv(&verts[(k+3)%verts.size()]);
				}
			}
			else
			{
				for(unsigned k = 0; k < verts.size(); k+=3) 
				{
					color_t t = cntcolor;
					t.a() = alpha;
					t.set_color();
					glVertex3dv(cnt);
					t = colors[k/3];
					t.a() = alpha;
					t.set_color();
					glVertex3dv(&verts[k]);
					t = colors[(k/3+1)%colors.size()];
					t.a() = alpha;
					t.set_color();
					glVertex3dv(&verts[(k+3)%verts.size()]);
				}
			}
		}
		else
		{
			for(unsigned k = 0; k < verts.size(); k+=3) 
			{
				glTexCoord1d(cnttexcoord);
				glVertex3dv(cnt);
				glTexCoord1d(texcoords[k/3]);
				glVertex3dv(&verts[k]);
				glTexCoord1d(texcoords[(k/3+1)%texcoords.size()]);
				glVertex3dv(&verts[(k+3)%verts.size()]);
			}
		}
	}
	void draw() const
	{
		for(unsigned k = 0; k < verts.size(); k+=3) 
		{
			glVertex3dv(cnt);
			glVertex3dv(&verts[k]);
			glVertex3dv(&verts[(k+3)%verts.size()]);
		}
	}
	void drawedges() const
	{
		for(unsigned k = 0; k < verts.size(); k+=3) 
		{
			glVertex3dv(&verts[k]); 
			glVertex3dv(&verts[(k+3)%verts.size()]); 
		}
	}
	bool operator <(const face2gl & other) const {return dist < other.dist;}
	void set_color(double r, double g, double b, double a) {c[0] = r; c[1] = g; c[2] = b; c[3] = a;}
	void add_vert(double x, double y, double z) {unsigned s = (unsigned)verts.size(); verts.resize(s+3); verts[s] = x; verts[s+1] = y; verts[s+2] = z;}
	void add_vert(double v[3]) {verts.insert(verts.end(),v,v+3);}
	void add_color(color_t c) {colors.push_back(c);}
	void add_texcoord(double val) {texcoords.push_back(val);}
	double * get_vert(int k) {return &verts[k*3];}
	unsigned size() {return (unsigned)verts.size()/3;}
	void set_center(double _cnt[3], color_t c = color_t(0,0,0,0))
	{
		cnt[0] = _cnt[0];
		cnt[1] = _cnt[1];
		cnt[2] = _cnt[2];
		cntcolor = c;
	}
	void get_center(float _cnt[3])
	{
		_cnt[0] = cnt[0];
		_cnt[1] = cnt[1];
		_cnt[2] = cnt[2];
	}
	void get_center(double _cnt[3])
	{
		_cnt[0] = cnt[0];
		_cnt[1] = cnt[1];
		_cnt[2] = cnt[2];
	}
	double * get_center() {return cnt;}
	void compute_center()
	{
		cnt[0] = cnt[1] = cnt[2] = 0;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < verts.size(); k+=3)
		{
			cnt[0] += verts[k+0];
			cnt[1] += verts[k+1];
			cnt[2] += verts[k+2];
		}
		cnt[0] /= (verts.size()/3)*1.0;
		cnt[1] /= (verts.size()/3)*1.0;
		cnt[2] /= (verts.size()/3)*1.0;
		compute_center_color();
		compute_center_texcoord();
	}
	void compute_center_color()
	{
		if( !colors.empty() )
		{
			cntcolor.r() = 0;
			cntcolor.g() = 0;
			cntcolor.b() = 0;
			cntcolor.a() = 0;
			for(INMOST_DATA_ENUM_TYPE k = 0; k < colors.size(); k++)
				cntcolor = cntcolor + colors[k];
			cntcolor =  cntcolor*(1.0f/static_cast<float>(colors.size()));
		}
	}
	void compute_center_texcoord()
	{
		if( !texcoords.empty() )
		{
			cnttexcoord = 0.0;
			for(INMOST_DATA_ENUM_TYPE k = 0; k < texcoords.size(); k++)
				cnttexcoord += texcoords[k];
			cnttexcoord /= static_cast<double>(texcoords.size());
		}
	}
	void compute_dist(double cam[3])
	{
		dist = sqrt((cnt[0]-cam[0])*(cnt[0]-cam[0])+(cnt[1]-cam[1])*(cnt[1]-cam[1])+(cnt[2]-cam[2])*(cnt[2]-cam[2]));
	}
	void set_flag(bool set) { flag = set;}
	bool get_flag() {return flag;}
	void set_elem(ElementType _etype, Storage::integer _id) {etype = _etype; id = _id;}
	Element get_elem(Mesh * m) {return m->ElementByLocalID(etype,id);}
};
face2gl DrawFace(Element f);

std::vector<face2gl> all_boundary;
std::vector<face2gl> added_faces;
std::vector<face2gl> clip_boundary;

void draw_faces_nc(std::vector<face2gl> & set, int highlight = -1)
{
  if( drawedges == 2 ) return;
	glColor4f(0,1,0,0.1);
	glBegin(GL_TRIANGLES);
	for(INMOST_DATA_ENUM_TYPE q = 0; q < set.size() ; q++) set[q].draw();
	glEnd();
	if( highlight != -1 )
	{
		glColor4f(1,0,0,1);
		glBegin(GL_TRIANGLES);
		set[highlight].draw();
		glEnd();
	}
}

void draw_faces(std::vector<face2gl> & set, int highlight = -1)
{
	if( drawedges == 2 ) return;
	if( visualization_tag.isValid() ) CommonColorBar->BindTexture();
	glBegin(GL_TRIANGLES);
	for(INMOST_DATA_ENUM_TYPE q = 0; q < set.size() ; q++) set[q].draw_colour();
	glEnd();
	if( visualization_tag.isValid() ) CommonColorBar->UnbindTexture();
	if( highlight != -1 )
	{
		glColor4f(1,0,0,1);
		glBegin(GL_TRIANGLES);
		set[highlight].draw();
		glEnd();
	}
}

void draw_faces_alpha(std::vector<face2gl> & set, double alpha)
{
	if( drawedges == 2 ) return;
	if( visualization_tag.isValid() ) CommonColorBar->BindTexture();
	glBegin(GL_TRIANGLES);
	for(INMOST_DATA_ENUM_TYPE q = 0; q < set.size() ; q++) set[q].draw_colour_alpha(alpha);
	glEnd();
	if( visualization_tag.isValid() ) CommonColorBar->UnbindTexture();
}

void draw_edges(std::vector<face2gl> & set, int highlight = -1)
{
  if( drawedges && drawedges != 2 )
  {
	  glBegin(GL_LINES);
	  for(INMOST_DATA_ENUM_TYPE q = 0; q < set.size() ; q++) set[q].drawedges();
	  glEnd();
	  if( highlight != -1 )
	  {
		  glColor4f(0,1,0,1);
		  glBegin(GL_LINES);
		  set[highlight].drawedges();
		  glEnd();
	  }
  }
}

void draw_faces_interactive_nc(std::vector<face2gl> & set)
{
  if( drawedges == 2 ) return;
	glColor4f(0,1,0,0.1);
	glBegin(GL_TRIANGLES);
	for(INMOST_DATA_ENUM_TYPE q = 0; q < set.size() ; q++) if( set[q].get_flag() ) set[q].draw();
	glEnd();
}

void draw_faces_interactive(std::vector<face2gl> & set)
{
  if( drawedges == 2 ) return;
  if( visualization_tag.isValid() ) CommonColorBar->BindTexture();
	glBegin(GL_TRIANGLES);
	for(INMOST_DATA_ENUM_TYPE q = 0; q < set.size() ; q++) if( set[q].get_flag() ) set[q].draw_colour();
	glEnd();
	if( visualization_tag.isValid() ) CommonColorBar->UnbindTexture();
}

void draw_faces_interactive_alpha(std::vector<face2gl> & set, double alpha)
{
  if( drawedges == 2 ) return;
  if( visualization_tag.isValid() ) CommonColorBar->BindTexture();
	glBegin(GL_TRIANGLES);
	for(INMOST_DATA_ENUM_TYPE q = 0; q < set.size() ; q++) if( set[q].get_flag() ) set[q].draw_colour_alpha(alpha);
	glEnd();
	if( visualization_tag.isValid() ) CommonColorBar->UnbindTexture();
}

void draw_edges_interactive(std::vector<face2gl> & set)
{
  if( drawedges && drawedges != 2 )
  {
	  glBegin(GL_LINES);
	  for(INMOST_DATA_ENUM_TYPE q = 0; q < set.size() ; q++) if( set[q].get_flag() ) set[q].drawedges();
	  glEnd();
  }
}

Storage::integer clip_plane_edge(double sp0[3], double sp1[3], double p[3], double n[3], double node[3])
{
	Storage::real u[3], w[3], D, N, sI;
	u[0] = sp1[0] - sp0[0]; u[1] = sp1[1] - sp0[1]; u[2] = sp1[2] - sp0[2];
	w[0] = sp0[0] - p[0];   w[1] = sp0[1] - p[1];   w[2] = sp0[2] - p[2];
	D =  (n[0]*u[0] + n[1]*u[1] + n[2]*u[2]);
	N = -(n[0]*w[0] + n[1]*w[1] + n[2]*w[2]);
	if( fabs(D) < 1.0e-9 )
	{
		if( fabs(N) < 1.0e-9 ) return CLIP_FULL;
		else return CLIP_NONE;
	}
	else
	{
		sI = N/D;
		if( sI < 0-1.0e-9 || sI > 1+1.0e-9 ) return CLIP_NONE;
		else
		{
			node[0] = sp0[0] + sI * u[0];
			node[1] = sp0[1] + sI * u[1];
			node[2] = sp0[2] + sI * u[2];
			return (sI > 1.0e-9 && sI < 1.0-1.0e-9) ? CLIP_NODE : CLIP_ENDP;
		}
	}
}

typedef struct point
{
	float coords[3];
	float diam;
	float dist;
	int id;
	/*
	point & operator = (point const & b)
	{
		coords[0] = b.coords[0];
		coords[1] = b.coords[1];
		coords[2] = b.coords[2];
		diam = b.diam;
		dist = b.dist;
		id = b.id;
	}

	point(const point & b)
	{
		coords[0] = b.coords[0];
		coords[1] = b.coords[1];
		coords[2] = b.coords[2];
		diam = b.diam;
		dist = b.dist;
		id = b.id;
	}
	*/
} point_t;

bool operator <(point_t & a, point_t & b) {return a.dist < b.dist;}




class volumetric
{
	std::vector<point_t> points;
	Mesh * m;

	void radix_sort_dist(std::vector<point_t> & set)
	{
		static std::vector<point_t> tmp;
		tmp.resize(set.size());
		unsigned int i;
		const unsigned int kHist = 2048;
		unsigned int  b0[kHist * 3];
		unsigned int *b1 = b0 + kHist;
		unsigned int *b2 = b1 + kHist;
		memset(b0,0,sizeof(unsigned int)*kHist*3);
		for (i = 0; i < set.size(); i++) 
		{
			unsigned int fi = flip((unsigned int *)&set[i].dist);
			++b0[_0(fi)]; ++b1[_1(fi)]; ++b2[_2(fi)];
		}
		{
			unsigned int sum0 = 0, sum1 = 0, sum2 = 0;
			for (i = 0; i < kHist; i++) 
			{
				b0[kHist-1] = b0[i] + sum0; b0[i] = sum0 - 1; sum0 = b0[kHist-1];
				b1[kHist-1] = b1[i] + sum1; b1[i] = sum1 - 1; sum1 = b1[kHist-1];
				b2[kHist-1] = b2[i] + sum2; b2[i] = sum2 - 1; sum2 = b2[kHist-1];
			}
		}
		for (i = 0; i < set.size(); i++) tmp[++b0[_0(flip((unsigned int *)&set[i].dist))]] = set[i];
		for (i = 0; i < set.size(); i++) set[++b1[_1(flip((unsigned int *)&tmp[i].dist))]] = tmp[i];
		for (i = 0; i < set.size(); i++) tmp[++b2[_2(flip((unsigned int *)&set[i].dist))]] = set[i];
		for (i = 0; i < set.size(); i++) set[i] = tmp[set.size()-1-i];
	}
public:
	
	volumetric(Mesh * _m)
	{
		m = _m;
		
		points.resize(m->NumberOfCells());
		int q = 0;
		for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
		{
			Storage::real cnt[3], cntf[3];
			it->Centroid(cnt);
			points[q].coords[0] = cnt[0];
			points[q].coords[1] = cnt[1];
			points[q].coords[2] = cnt[2];
			points[q].id = it->LocalID();
			points[q].diam = 0.f;
			ElementArray<Face> faces = it->getFaces();
			for(ElementArray<Face>::iterator f = faces.begin(); f != faces.end(); ++f)
			{
				f->Centroid(cntf);
				Storage::real d = sqrt((cnt[0]-cntf[0])*(cnt[0]-cntf[0])+(cnt[1]-cntf[1])*(cnt[1]-cntf[1])+(cnt[2]-cntf[2])*(cnt[2]-cntf[2]));
				if( points[q].diam < d ) points[q].diam = d;
			}
			++q;
		}
		/*
		points.reserve(m->NumberOfNodes()*20);
		int q = 0;
		for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
		{
			point_t pnt;
			Storage::real_array cnt = it->Coords();
			pnt.coords[0] = cnt[0];
			pnt.coords[1] = cnt[1];
			pnt.coords[2] = cnt[2];
			pnt.id = it->LocalID();
			points.push_back(pnt);
			
			ElementArray<Element> adj = it->getAdjElements(CELL);
			for(ElementArray<Element>::iterator jt = adj.begin(); jt != adj.end(); ++jt)
			{
				Storage::real cnt2[3];
				jt->Centroid(cnt2);
				const int ncoefs = 1;
				const Storage::real coefs[ncoefs] = {0.82};
				for(int j = 0; j < ncoefs; ++j)
				{
					pnt.coords[0] = cnt[0]*(1-coefs[j]) + cnt2[0]*coefs[j];
					pnt.coords[1] = cnt[1]*(1-coefs[j]) + cnt2[1]*coefs[j];
					pnt.coords[2] = cnt[2]*(1-coefs[j]) + cnt2[2]*coefs[j];
					pnt.id = it->LocalID();
					points.push_back(pnt);
				}
			}
			
		}
		*/
		printf("number of points %d\n",(int)points.size());
	}
	void camera(double pos[3], int interactive)
	{
		if( interactive ) return;
		float posf[3];
		posf[0] = pos[0];
		posf[1] = pos[1];
		posf[2] = pos[2];
		for(int k = 0; k < points.size(); ++k)
		{
			points[k].dist = sqrtf(
				(posf[0]-points[k].coords[0])*(posf[0]-points[k].coords[0])+
				(posf[1]-points[k].coords[1])*(posf[1]-points[k].coords[1])+
				(posf[2]-points[k].coords[2])*(posf[2]-points[k].coords[2]));
		}
		double t = Timer();
		radix_sort_dist(points);
		//std::sort(points.rbegin(),points.rend());
		printf("Time to sort %lf\n",Timer()-t);
	}
	void draw(int interactive)
	{
		double origin[3], right[3], up[3];
		GLdouble modelview[16],projection[16];
		GLint viewport[4] = {0,0,1,1};
		glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
		glGetDoublev(GL_PROJECTION_MATRIX, projection);
		GLdouble outx, outy, outz;  // Var's to save the answer in
		gluUnProject(0.5, 0.5, 0.,
               modelview, projection, viewport,
               &outx, &outy, &outz);
		origin[0] = outx;
		origin[1] = outy;
		origin[2] = outz;
		gluUnProject(1.0, 0.5, 0.,
               modelview, projection, viewport,
               &outx, &outy, &outz);
		right[0] = outx;
		right[1] = outy;
		right[2] = outz;
		gluUnProject(0.5, 1.0, 0.,
               modelview, projection, viewport,
               &outx, &outy, &outz);
		up[0] = outx;
		up[1] = outy;
		up[2] = outz;
		right[0] -= origin[0];
		right[1] -= origin[1];
		right[2] -= origin[2];
		up[0] -= origin[0];
		up[1] -= origin[1];
		up[2] -= origin[2];
    
		double l = sqrt(right[0]*right[0]+right[1]*right[1]+right[2]*right[2]);
		if( l )
		{
			right[0] /= l;
			right[1] /= l;
			right[2] /= l;
		}
		l = sqrt(up[0]*up[0]+up[1]*up[1]+up[2]*up[2]);
		if( l )
		{
			up[0] /= l;
			up[1] /= l;
			up[2] /= l;
		}
    
		const float alpha = 0.0075f;
		const float mult = 1.0f;
		const float rmult = 0.7f;
		//glPointSize(5.0);
		glColor4f(0.5f,0.5f,0.5f,alpha);
		glEnable(GL_BLEND);
		//glBegin(GL_TRIANGLES);
		//glBegin(GL_QUADS);
		if( interactive )
		{
      for(int k = 0; k < points.size(); ++k) if( k % 100 == 0 )
			{
				if( visualization_tag.isValid() )
				{
					color_t c = CommonColorBar->pick_color(m->CellByLocalID(points[k].id)->RealDF(visualization_tag));
					c.a() = alpha;
					c.set_color();
				}
				
				glBegin(GL_TRIANGLE_FAN);
				glVertex3f(points[k].coords[0],points[k].coords[1],points[k].coords[2]);
				glVertex3f(points[k].coords[0]+(right[0]*rmult+up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]*rmult+up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]*rmult+up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(right[0]*rmult*2)*points[k].diam*mult,points[k].coords[1]+(right[1]*rmult*2)*points[k].diam*mult,points[k].coords[2]+(right[2]*rmult*2)*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(right[0]*rmult-up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]*rmult-up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]*rmult-up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(-right[0]*rmult-up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]*rmult-up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]*rmult-up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(-right[0]*rmult*2)*points[k].diam*mult,points[k].coords[1]+(-right[1]*rmult*2)*points[k].diam*mult,points[k].coords[2]+(-right[2]*rmult*2)*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(-right[0]*rmult+up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]*rmult+up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]*rmult+up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(right[0]*rmult+up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]*rmult+up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]*rmult+up[2])*points[k].diam*mult);
				glEnd();
				
				/*
				glVertex3f(points[k].coords[0]+(right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]+up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(-up[0]*2)*points[k].diam*mult,points[k].coords[1]+(-up[1]*2)*points[k].diam*mult,points[k].coords[2]+(-up[2]*2)*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(-right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]+up[2])*points[k].diam*mult);
				*/
				/*
				glVertex3f(points[k].coords[0]+(right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]+up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(right[0]-up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]-up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]-up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(-right[0]-up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]-up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]-up[2])*points[k].diam*mult);
				glVertex3f(points[k].coords[0]+(-right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]+up[2])*points[k].diam*mult);
				*/
				//glVertex3f(points[k].coords[0]+right[0]*points[k].diam*mult,points[k].coords[1]+right[1]*points[k].diam*mult,points[k].coords[2]+right[2]*points[k].diam*mult);
				//glVertex3f(points[k].coords[0]+up[0]*points[k].diam*mult,points[k].coords[1]+up[1]*points[k].diam*mult,points[k].coords[2]+up[2]*points[k].diam*mult);
				//glVertex3f(points[k].coords[0]-right[0]*points[k].diam*mult,points[k].coords[1]-right[1]*points[k].diam*mult,points[k].coords[2]-right[2]*points[k].diam*mult);
				//glVertex3f(points[k].coords[0]-up[0]*points[k].diam*mult,points[k].coords[1]-up[1]*points[k].diam*mult,points[k].coords[2]-up[2]*points[k].diam*mult);
			}
		}
		else
		for(int k = 0; k < points.size(); ++k)
		{
			if( visualization_tag.isValid() )
			{
				color_t c = CommonColorBar->pick_color(m->CellByLocalID(points[k].id)->RealDF(visualization_tag));
				c.a() = alpha;
				c.set_color();
			}
			
			glBegin(GL_TRIANGLE_FAN);
			glVertex3f(points[k].coords[0],points[k].coords[1],points[k].coords[2]);
			glVertex3f(points[k].coords[0]+(right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]+up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(right[0]*2)*points[k].diam*mult,points[k].coords[1]+(right[1]*2)*points[k].diam*mult,points[k].coords[2]+(right[2]*2)*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(right[0]-up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]-up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]-up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(-right[0]-up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]-up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]-up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(-right[0]*2)*points[k].diam*mult,points[k].coords[1]+(-right[1]*2)*points[k].diam*mult,points[k].coords[2]+(-right[2]*2)*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(-right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]+up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]+up[2])*points[k].diam*mult);
			glEnd();	
			
			/*
			glVertex3f(points[k].coords[0]+(right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]+up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(-up[0]*2)*points[k].diam*mult,points[k].coords[1]+(-up[1]*2)*points[k].diam*mult,points[k].coords[2]+(-up[2]*2)*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(-right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]+up[2])*points[k].diam*mult);
			*/
			/*
			glVertex3f(points[k].coords[0]+(right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]+up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(right[0]-up[0])*points[k].diam*mult,points[k].coords[1]+(right[1]-up[1])*points[k].diam*mult,points[k].coords[2]+(right[2]-up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(-right[0]-up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]-up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]-up[2])*points[k].diam*mult);
			glVertex3f(points[k].coords[0]+(-right[0]+up[0])*points[k].diam*mult,points[k].coords[1]+(-right[1]+up[1])*points[k].diam*mult,points[k].coords[2]+(-right[2]+up[2])*points[k].diam*mult);
			*/
			//glVertex3f(points[k].coords[0]+right[0]*points[k].diam*mult,points[k].coords[1]+right[1]*points[k].diam*mult,points[k].coords[2]+right[2]*points[k].diam*mult);
			//glVertex3f(points[k].coords[0]+up[0]*points[k].diam*mult,points[k].coords[1]+up[1]*points[k].diam*mult,points[k].coords[2]+up[2]*points[k].diam*mult);
			//glVertex3f(points[k].coords[0]-right[0]*points[k].diam*mult,points[k].coords[1]-right[1]*points[k].diam*mult,points[k].coords[2]-right[2]*points[k].diam*mult);
			//glVertex3f(points[k].coords[0]-up[0]*points[k].diam*mult,points[k].coords[1]-up[1]*points[k].diam*mult,points[k].coords[2]-up[2]*points[k].diam*mult);
		}
		//glEnd();
		glDisable(GL_BLEND);
		//glPointSize(1.0);
	}
} * CommonVolumetricView;


/*
class volumetric2
{
	std::vector<face2gl> faces;
	Mesh * m;
public:
	volumetric2(Mesh * _m)
	{
		m = _m;
		
		faces.reserve(m->NumberOfFaces());
		int q = 0;
		INMOST_DATA_ENUM_TYPE pace = std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,(unsigned)m->NumberOfFaces()/100));
		for(Mesh::iteratorFace it = m->BeginFace(); it != m->EndFace(); ++it)
		{
			faces.push_back(DrawFace(it->self()));
			if( q%pace == 0 ) faces.back().set_flag(true);
			q++;
		}
		printf("number of faces %d\n",faces.size());
	}
	void camera(double pos[3], int interactive)
	{
		if( interactive ) return;
		for(int k = 0; k < faces.size(); ++k)
			faces[k].compute_dist(pos);
		//face2gl::radix_sort_dist(faces);
		std::sort(faces.rbegin(),faces.rend());
	}
	void draw(int interactive)
	{
		if( interactive ) draw_faces_interactive_alpha(faces,0.05);
		else draw_faces_alpha(faces,0.05);
	}
}* CommonVolumetricView;
*/
class kdtree
{
	int marked;
	struct entry
	{
		HandleType e;
		float xyz[3];
		struct entry & operator =(const struct entry & other)
		{
			e = other.e;
			xyz[0] = other.xyz[0];
			xyz[1] = other.xyz[1];
			xyz[2] = other.xyz[2];
			return *this;
		}
	} * set;
	Mesh * m;
	INMOST_DATA_ENUM_TYPE size;
	float bbox[6];
	kdtree * children;
	static int cmpElements0(const void * a,const void * b) 
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[0];
		float bd = eb->xyz[0];
		return (ad > bd) - (ad < bd);
	}
	static int cmpElements1(const void * a,const void * b) 
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[1];
		float bd = eb->xyz[1];
		return (ad > bd) - (ad < bd);
	}
	static int cmpElements2(const void * a,const void * b) 
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[2];
		float bd = eb->xyz[2];
		return (ad > bd) - (ad < bd);
	}
	void radix_sort(int dim, struct entry * temp)
	{
		unsigned int i;
		const unsigned int kHist = 2048;
		unsigned int  b0[kHist * 3];
		unsigned int *b1 = b0 + kHist;
		unsigned int *b2 = b1 + kHist;
		memset(b0,0,sizeof(unsigned int)*kHist*3);
		for (i = 0; i < size; i++) 
		{
			unsigned int fi = flip((unsigned int *)&set[i].xyz[dim]);
			++b0[_0(fi)]; ++b1[_1(fi)]; ++b2[_2(fi)];
		}
		{
			unsigned int sum0 = 0, sum1 = 0, sum2 = 0;
			for (i = 0; i < kHist; i++) 
			{
				b0[kHist-1] = b0[i] + sum0; b0[i] = sum0 - 1; sum0 = b0[kHist-1];
				b1[kHist-1] = b1[i] + sum1; b1[i] = sum1 - 1; sum1 = b1[kHist-1];
				b2[kHist-1] = b2[i] + sum2; b2[i] = sum2 - 1; sum2 = b2[kHist-1];
			}
		}
		for (i = 0; i < size; i++) temp[++b0[_0(flip((unsigned int *)&set[i].xyz[dim]))]] = set[i];
		for (i = 0; i < size; i++) set[++b1[_1(flip((unsigned int *)&temp[i].xyz[dim]))]] = temp[i];
		for (i = 0; i < size; i++) temp[++b2[_2(flip((unsigned int *)&set[i].xyz[dim]))]] = set[i];
		for (i = 0; i < size; i++) set[i] = temp[i];
	}
	void kdtree_build(int dim, int & done, int total, struct entry * temp)
	{
		if( size > 1 )
		{
			if( size > 128 ) radix_sort(dim,temp); else 
			switch(dim)
			{
			case 0: qsort(set,size,sizeof(entry),cmpElements0);break;
			case 1: qsort(set,size,sizeof(entry),cmpElements1);break;
			case 2: qsort(set,size,sizeof(entry),cmpElements2);break;
			}
			children = static_cast<kdtree *>(malloc(sizeof(kdtree)*2));//new kdtree[2];
			children[0].marked = 0;
			children[0].children = NULL;
			children[0].set = set;
			children[0].size = size/2;
			children[0].m = m;
			children[1].marked = 0;
			children[1].children = NULL;
			children[1].set = set+size/2;
			children[1].size = size - size/2;
			children[1].m = m;
			children[0].kdtree_build((dim+1)%3,done,total,temp);
			children[1].kdtree_build((dim+1)%3,done,total,temp);
			for(int k = 0; k < 3; k++)
			{
				bbox[0+2*k] = std::min(children[0].bbox[0+2*k],children[1].bbox[0+2*k]);
				bbox[1+2*k] = std::max(children[0].bbox[1+2*k],children[1].bbox[1+2*k]);
			}
		}
		else 
		{
			assert(size == 1);
			if( GetHandleElementType(set[0].e) == EDGE )
			{
				Storage::real_array n1 = Edge(m,set[0].e)->getBeg()->Coords();
				Storage::real_array n2 = Edge(m,set[0].e)->getEnd()->Coords();
				for(int k = 0; k < 3; k++)
				{
					bbox[0+2*k] = std::min(n1[k],n2[k]);
					bbox[1+2*k] = std::max(n1[k],n2[k]);
				}
				done++;
				if( done%150 == 0 )
				{
					printf("%3.1f%%\r",(done*100.0)/(total*1.0));
					fflush(stdout);
				}
			}
			else
			{
				ElementArray<Node> nodes = Element(m,set[0].e)->getNodes();
				bbox[0] = bbox[2] = bbox[4] = 1.0e20;
				bbox[1] = bbox[3] = bbox[5] = -1.0e20;
				for(INMOST_DATA_ENUM_TYPE k = 0; k < nodes.size(); ++k)
				{
					Storage::real_array coords = nodes[k].Coords();
					for(INMOST_DATA_ENUM_TYPE q = 0; q < 3; q++)
					{
						bbox[q*2+0] = std::min<float>(bbox[q*2+0],coords[q]);
						bbox[q*2+1] = std::max<float>(bbox[q*2+1],coords[q]);
					}
				}
			}
		}
	}
	kdtree() : marked(0), set(NULL), size(0), children(NULL) {}
	inline int plane_bbox(double p[3], double n[3])
	{
		Storage::real pv[3], nv[3];
		for(int k = 0; k < 3; ++k)
		{
			if( n[k] >= 0 ) 
			{ 
				pv[k] = bbox[1+2*k]; //max
				nv[k] = bbox[0+2*k]; //min
			} 
			else 
			{ 
				pv[k] = bbox[0+2*k]; //min
				nv[k] = bbox[1+2*k]; //max
			}
		}
		Storage::real pvD, nvD;
		pvD = n[0]*(pv[0]-p[0])+n[1]*(pv[1]-p[1])+n[2]*(pv[2]-p[2]);
		nvD = n[0]*(nv[0]-p[0])+n[1]*(nv[1]-p[1])+n[2]*(nv[2]-p[2]);
		if( nvD*pvD <= 0.0 )
			return 2;
		else if( nvD < 0.0 )
			return 1;
		else return 0;
	}
	bool sub_intersect_plane_edge(Tag clip_point, Tag clip_state, ElementArray<Cell> & cells, MarkerType mrk, double p[3], double n[3])
	{
		if( size == 1 )
		{
			assert( GetHandleElementType(set[0].e) == EDGE );
			Edge ee = Edge(m,set[0].e);
			Storage::real_array sp0 = ee->getBeg()->Coords();
			Storage::real_array sp1 = ee->getEnd()->Coords();
			Storage::integer & clip = m->IntegerDF(set[0].e,clip_state);
			clip = clip_plane_edge(&sp0[0],&sp1[0],p,n,&ee->RealArrayDF(clip_point)[0]);
			if( clip )
			{
				ElementArray<Cell> ecells = ee->getCells();
				for(INMOST_DATA_ENUM_TYPE k = 0; k < ecells.size(); ++k) if( !ecells[k].GetMarker(mrk) )
				{
					ecells[k].SetMarker(mrk);
					cells.push_back(ecells[k]);
				}
				marked = 1;
			}
		}
		else if( plane_bbox(p,n) == 2 )
		{
			bool test1 = children[0].sub_intersect_plane_edge(clip_point,clip_state,cells,mrk,p,n);
			bool test2 = children[1].sub_intersect_plane_edge(clip_point,clip_state,cells,mrk,p,n);
			if( test1 || test2 ) marked = 1;
		}
		return marked != 0;
	}
	void sub_intersect_plane_faces(Tag clip_state, double p[3], double n[3])
	{
		if( size == 1 )
		{
			Storage::integer state;
			Element ee(m,set[0].e);
			assert( ee->GetElementDimension() == 2 );
			ElementArray<Node> nodes = ee->getNodes();
			Storage::real_array coords = nodes[0].Coords();
			Storage::real dot0 = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
			if( dot0 <= 0.0 ) state = CLIP_FACE_INSIDE; else state = CLIP_FACE_OUTSIDE;
			for(INMOST_DATA_ENUM_TYPE k = 1; k < nodes.size(); k++)
			{
				coords = nodes[k].Coords();
				Storage::real dot = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
				if( dot*dot0 <= 0.0 ) 
				{
					state = CLIP_FACE_INTERSECT;
					break;
				}
			}
			m->IntegerDF(set[0].e,clip_state) = state;
		}
		else
		{
			marked = plane_bbox(p,n);
			if( marked == 0 )
			{
				for(INMOST_DATA_ENUM_TYPE k = 0; k < size; k++) m->IntegerDF(set[k].e,clip_state) = CLIP_FACE_OUTSIDE;
			}
			else if( marked == 1 )
			{
				for(INMOST_DATA_ENUM_TYPE k = 0; k < size; k++) m->IntegerDF(set[k].e,clip_state) = CLIP_FACE_INSIDE;
			}
			else
			{
				children[0].sub_intersect_plane_faces(clip_state,p,n);
				children[1].sub_intersect_plane_faces(clip_state,p,n);
			}
		}
	}
	void unmark_old_edges(Tag clip_state)
	{
		if( size == 1 )
		{
			assert(GetHandleElementType(set[0].e) == EDGE);
			marked = 0;
			if( GetHandleElementType(set[0].e) == EDGE )
				m->IntegerDF(set[0].e,clip_state) = CLIP_NONE;
			else if( GetHandleElementType(set[0].e) == FACE )
				m->IntegerDF(set[0].e,clip_state) = CLIP_FACE_NONE;
		}
		else if( children )
		{
			if(children[0].marked) {children[0].unmark_old_edges(clip_state); marked = 0;}
			if(children[1].marked) {children[1].unmark_old_edges(clip_state); marked = 0;}
		}
	}
	void clear_children() { if( children ) {children[0].clear_children(); children[1].clear_children(); free(children);}}
public:
	kdtree(Mesh * m) :  marked(0),m(m),children(NULL)
	{
		double tt;
		size = m->NumberOfEdges();
		assert(size > 1);
		set = new entry[size];
		INMOST_DATA_ENUM_TYPE k = 0;
		tt = Timer();
		printf("Prepearing edge set.\n");
		for(Mesh::iteratorEdge it = m->BeginEdge(); it != m->EndEdge(); ++it) 
		{
			set[k].e = *it;
			set[k].xyz[0] = (it->getBeg()->Coords()[0] + it->getEnd()->Coords()[0])*0.5;
			set[k].xyz[1] = (it->getBeg()->Coords()[1] + it->getEnd()->Coords()[1])*0.5;
			set[k].xyz[2] = (it->getBeg()->Coords()[2] + it->getEnd()->Coords()[2])*0.5;
			k++;
			if( k%150 == 0 ) 
			{
				printf("%3.1f%%\r",(k*100.0)/(size*1.0));
				fflush(stdout);
			}
		}
		printf("Done. Time %lg\n",Timer()-tt);
		int done = 0, total = size;
		printf("Building KD-tree.\n");
		tt = Timer();
		struct entry *  temp = new entry[size];
		kdtree_build(0,done,total,temp);
		delete [] temp;
		printf("Done. Time %lg\n",Timer()-tt);
		for(int k = 0; k < 3; k++)
		{
			bbox[0+2*k] = std::min(children[0].bbox[0+2*k],children[1].bbox[0+2*k]);
			bbox[1+2*k] = std::max(children[0].bbox[1+2*k],children[1].bbox[1+2*k]);
		}
	}
	kdtree(Mesh * m, HandleType * eset, INMOST_DATA_ENUM_TYPE size) : marked(0), m(m),size(size),children(NULL)
	{
		double tt;
		assert(size > 1);
		set = new entry[size];
		tt = Timer();
		printf("Prepearing elements set.\n");
		for(INMOST_DATA_ENUM_TYPE k = 0; k < size; k++) 
		{
			set[k].e = eset[k];
			Storage::real cnt[3];
			m->GetGeometricData(set[k].e,CENTROID,cnt);
			set[k].xyz[0] = cnt[0];
			set[k].xyz[1] = cnt[1];
			set[k].xyz[2] = cnt[2];
			if( k%150 == 0 ) 
			{
				printf("%3.1f%%\r",(k*100.0)/(size*1.0));
				fflush(stdout);
			}
		}
		printf("Done. Time %lg\n",Timer()-tt);
		int done = 0, total = size;
		printf("Building KD-tree.\n");
		tt = Timer();
		struct entry *  temp = new entry[size];
		kdtree_build(0,done,total,temp);
		delete [] temp;
		printf("Done. Time %lg\n",Timer()-tt);
		for(int k = 0; k < 3; k++)
		{
			bbox[0+2*k] = std::min(children[0].bbox[0+2*k],children[1].bbox[0+2*k]);
			bbox[1+2*k] = std::max(children[0].bbox[1+2*k],children[1].bbox[1+2*k]);
		}
	}
	void intersect_plane_edge(Tag clip_point, Tag clip_state, ElementArray<Cell> & cells, MarkerType mark_cells, double p[3], double n[3])
	{
		if( marked ) 
		{
			unmark_old_edges(clip_state);
			cells.clear();
		}
		sub_intersect_plane_edge(clip_point, clip_state, cells,mark_cells,p,n);
	}
	void intersect_plane_face(Tag clip_state, double p[3], double n[3])
	{
		sub_intersect_plane_faces(clip_state, p,n);
	}
	~kdtree()
	{
		delete [] set;
		clear_children();
	}
};

class clipper
{
	struct edge_point
	{
		double val;
		Storage::real xyz[3];
		Storage::integer edge;
		edge_point(){}
		edge_point(Storage::real _xyz[3], Storage::integer n, float v)
		{
			xyz[0] = _xyz[0];
			xyz[1] = _xyz[1];
			xyz[2] = _xyz[2];
			edge = n;
			val = v;
		}
		bool operator ==(const edge_point& b) const
		{
			Storage::real temp = 0.0;
			for(int k = 0; k < 3; k++) temp += (xyz[k]-b.xyz[k])*(xyz[k]-b.xyz[k]);
			if( temp < 1.0e-8 ) return true; else return false;
		}
		bool operator !=(const edge_point& b) const {return !(operator ==(b));}
		void print() {printf("%g %g %g e %d\n",xyz[0],xyz[1],xyz[2],edge);}
	};
	Tag clip_point, clip_state;
	kdtree * tree;
	Tag clips, clipsv;
	MarkerType marker;
	ElementArray<Cell> cells;
	Mesh * mm;
public:
	~clipper() 
	{
		delete tree; 
		for(INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); k++) cells[k]->RemMarker(marker);
		mm->ReleaseMarker(marker); 
		mm->DeleteTag(clips); 
		mm->DeleteTag(clipsv); 
		mm->DeleteTag(clip_point); 
		mm->DeleteTag(clip_state);
	}
	clipper(Mesh * m)
	{
		mm = m;
		cells.SetMeshLink(mm);
		tree = new kdtree(m);
		marker = m->CreateMarker();
		clips = m->CreateTag("CLIPS",DATA_REAL,CELL,CELL);
		clipsv = m->CreateTag("CLIPS_VAL",DATA_REAL,CELL,CELL);
		clip_point = m->CreateTag("CLIP_POINT",DATA_REAL,EDGE,NONE,3);
		clip_state = m->CreateTag("CLIP_STATE",DATA_INTEGER,EDGE,NONE,1);
	}
	double compute_value(Edge e, Storage::real * pnt)
	{
		if( visualization_tag.isValid() )
		{
			Storage::real_array c1 = e->getBeg()->Coords();
			Storage::real_array c2 = e->getEnd()->Coords();
			Storage::real d1,d2,t;
			d1 = sqrt((pnt[0]-c1[0])*(pnt[0]-c1[0])+(pnt[1]-c1[1])*(pnt[1]-c1[1])+(pnt[2]-c1[2])*(pnt[2]-c1[2]));
			d2 = sqrt((c2[0]-c1[0])*(c2[0]-c1[0])+(c2[1]-c1[1])*(c2[1]-c1[1])+(c2[2]-c1[2])*(c2[2]-c1[2]));
			t = d1/d2; //(pnt == c2, t = 1 : pnt == c1, t = 0)
			return e->getBeg()->RealDF(visualization_tag)*(1-t)+e->getEnd()->RealDF(visualization_tag)*t;
		}
		else return 0.f;
	}
	void clip_plane(Storage::real p[3], Storage::real n[3])
	{
		const bool print = false;

		for(INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); ++k) 
			if( cells[k]->GetMarker(marker) )
			{
				cells[k]->RealArray(clips).clear();
				cells[k]->RealArray(clipsv).clear();
				cells[k]->RemMarker(marker);
			}
		tree->intersect_plane_edge(clip_point,clip_state,cells,marker,p,n);
		dynarray<edge_point,128> clipcoords, loopcoords;
		std::vector<bool> closed;
		for(INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); ++k)
		{
			//assuming faces are convex we will have at most one clipping edge per polygon
			//otherwise every pair of clipping nodes forming edge should appear consequently
			//as long as we go through face's edges in ordered way
			clipcoords.clear();
			//we will gather all the pairs of nodes, then form closed loop
			ElementArray<Face> faces = cells[k]->getFaces();
			Face full_face;
			int ntotpoints = 0, ntotedges = 0;
			for(INMOST_DATA_ENUM_TYPE q = 0; q < faces.size(); ++q)
			{
				int last_edge_type = CLIP_NONE;
				int nfulledges = 0, npoints = 0, nstartedge = ntotedges;
				ElementArray<Edge> edges = faces[q].getEdges();
				for(INMOST_DATA_ENUM_TYPE r = 0; r < edges.size(); ++r)
				{
					Storage::integer state = edges[r].IntegerDF(clip_state);
					if( state == CLIP_FULL )
					{
						nfulledges++;
						edge_point n1 = edge_point(&edges[r].getBeg()->Coords()[0],ntotedges,visualization_tag.isValid() ? edges[r].getBeg()->RealDF(visualization_tag) : 0.f);
						edge_point n2 = edge_point(&edges[r].getEnd()->Coords()[0],ntotedges,visualization_tag.isValid() ? edges[r].getEnd()->RealDF(visualization_tag) : 0.f);
						if( npoints % 2 == 0 ) //all privious edges are closed, just add this one
						{
							clipcoords.push_back(n1);
							clipcoords.push_back(n2);
							npoints+=2;
							ntotedges++;
							last_edge_type = CLIP_FULL;
						} 
						else if(n1 == clipcoords.back()) //this may be prolongation of one point that hit one edge
						{
							clipcoords.push_back(n2);
							npoints++;
							ntotedges++;
							last_edge_type = CLIP_FULL;
						}
						else if( n2 == clipcoords.back() )
						{
							clipcoords.push_back(n1);
							npoints++;
							ntotedges++;
							last_edge_type = CLIP_FULL;
						}
						else printf("%s:%d strange orphan node before me\n",__FILE__,__LINE__);
					}
					else if( state == CLIP_ENDP )
					{
						edge_point n = edge_point(&edges[r].RealArrayDF(clip_point)[0],ntotedges,compute_value(edges[r],&edges[r].RealArrayDF(clip_point)[0]));
						bool add = true;
						if( last_edge_type == CLIP_ENDP )
						{
							if( n == clipcoords.back() )
								add = false;
						}
						else if( last_edge_type == CLIP_FULL )
						{
							if( n == clipcoords.back() || n == clipcoords[clipcoords.size()-2])
								add = false;
						}
						if( add ) //this one node should be prolongation of privious edge
						{
							if( print )
							{
								printf("added: ");
								n.print();
							}
							clipcoords.push_back(n);
							npoints++;
							if( npoints % 2 == 0 ) 
							{
								if( print ) printf("edge %d accepted\n",ntotedges);
								ntotedges++;
							}
							last_edge_type = CLIP_ENDP;
						}
						else if( print )
						{
							printf("ignored: ");
							n.print();
						}
					}
					else if( state == CLIP_NODE )
					{
						edge_point n = edge_point(&edges[r].RealArrayDF(clip_point)[0],ntotedges,compute_value(edges[r],&edges[r].RealArrayDF(clip_point)[0]));
						if( print )
						{
							printf("added: ");
							n.print();
						}
						clipcoords.push_back(n);
						npoints++;
						if( npoints % 2 == 0 ) 
						{
							if( print ) printf("edge %d accepted\n",ntotedges);
							ntotedges++;
						}
						last_edge_type = CLIP_NODE;
					}
				}
				if( npoints % 2 != 0 ) 
				{
					if( print ) printf("edge %d not closed - remove\n",ntotedges);
					clipcoords.pop_back();
					npoints--;
					//printf("%s:%d this should not happen!\n",__FILE__,__LINE__);
				}
				
				if( nfulledges == static_cast<int>(edges.size()) )
				{
					full_face = faces[q];
					break;
				}
				if( print )
				{
					printf("nodes on face %d\n",faces[q].LocalID());
					for(int m = nstartedge*2; m < static_cast<int>(clipcoords.size()); m++) clipcoords[m].print();
				}
				ntotpoints += npoints;
			}
			if( full_face.isValid() )
			{
				ElementArray<Node> nodes = full_face->getNodes();
				Storage::real_array cl = cells[k]->RealArray(clips);
				Storage::real_array clv = cells[k]->RealArray(clipsv);
				cl.resize(static_cast<Storage::real_array::size_type>(3*nodes.size()));
				clv.resize(static_cast<Storage::real_array::size_type>(nodes.size()));
				for(INMOST_DATA_ENUM_TYPE r = 0; r < nodes.size(); r++)
				{
					Storage::real_array p = nodes[r].Coords();
					cl[0+3*r] = p[0];
					cl[1+3*r] = p[1];
					cl[2+3*r] = p[2];
					clv[r] = visualization_tag.isValid() ? nodes[r].RealDF(visualization_tag) : 0.0;
				}
				cells[k]->SetMarker(marker);
			}
			else if( ntotedges > 2 )
			{
				if( print )
				{
					printf("coords on cell %d\n",cells[k]->LocalID());
					for(int m = 0; m < static_cast<int>(clipcoords.size()); m++) clipcoords[m].print();
				}
				//Can make this faster using hash
				closed.resize(ntotedges);
				std::fill(closed.begin(),closed.end(),false);
				loopcoords.push_back(clipcoords[0]); //this is starting point
				loopcoords.push_back(clipcoords[1]); //this is next
				closed[0] = true;
				for(int r = 0; r < ntotedges-2; ++r) //we need to add this number of points
				{
					bool hit = false;
					for(int q = 0; q < ntotedges; ++q) if( !closed[q] )
					{
						//some end of q-th edge connects to current end point - connect it
						if( clipcoords[q*2+0] == loopcoords.back() )
						{
							loopcoords.push_back(clipcoords[q*2+1]);
							closed[q] = true;
							hit = true;
							break;
						}
						else if( clipcoords[q*2+1] == loopcoords.back() )
						{
							loopcoords.push_back(clipcoords[q*2+0]);
							closed[q] = true;
							hit = true;
							break;
						}
					}
					if( !hit ) printf("%s:%d cannot find end for edge! total edges %d current loop size %ld\n",
						__FILE__,__LINE__,ntotedges,loopcoords.size());
				}
				Storage::real_array cl = cells[k]->RealArray(clips);
				Storage::real_array clv = cells[k]->RealArray(clipsv);
				cl.resize(static_cast<Storage::real_array::size_type>(3*loopcoords.size()));
				clv.resize(static_cast<Storage::real_array::size_type>(loopcoords.size()));
				for(INMOST_DATA_ENUM_TYPE r = 0; r < loopcoords.size(); ++r)
				{
					cl[r*3+0] = loopcoords[r].xyz[0];
					cl[r*3+1] = loopcoords[r].xyz[1];
					cl[r*3+2] = loopcoords[r].xyz[2];
					clv[r] = loopcoords[r].val;
				}

				loopcoords.clear();
				clipcoords.clear();

				cells[k]->SetMarker(marker);
			}
		}
	}
	void gen_clip(std::vector<face2gl> & out)
	{
		INMOST_DATA_ENUM_TYPE pace = std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,size()/100));
		for(INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); k++) if( cells[k]->GetMarker(marker))
		{
			face2gl f;
			f.set_color(0.6,0.6,0.6,1);
			Storage::real_array cl = cells[k]->RealArray(clips);
			Storage::real_array clv = cells[k]->RealArray(clipsv);
			for(INMOST_DATA_ENUM_TYPE q = 0; q < cl.size(); q+=3) f.add_vert(&cl[q]);
			if( visualization_tag.isValid() )
			{
				for(INMOST_DATA_ENUM_TYPE q = 0; q < clv.size(); q++) 
				{
					//f.add_color(CommonColorBar->pick_color(clv[q]));
					f.add_texcoord(CommonColorBar->pick_texture(clv[q]));
				}
			}
			f.compute_center();
			f.set_elem(cells[k]->GetElementType(),cells[k]->LocalID());
			if( k%pace == 0 ) f.set_flag(true);
			out.push_back(f);
		}
	}
	void draw_clip(INMOST_DATA_ENUM_TYPE pace)
	{
		if( visualization_tag.isValid() ) CommonColorBar->BindTexture();
		glBegin(GL_TRIANGLES);
		for(INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); k+=pace) if( cells[k]->GetMarker(marker))
		{
			Storage::real_array cl = cells[k]->RealArray(clips);
			Storage::real_array clv = cells[k]->RealArray(clipsv);
			Storage::real cnt[3] = {0,0,0};
			Storage::real cntv = 0.0;
			for(INMOST_DATA_ENUM_TYPE q = 0; q < cl.size(); q+=3)
			{
				cnt[0] += cl[q+0];
				cnt[1] += cl[q+1];
				cnt[2] += cl[q+2];
				cntv += clv[q/3];
			}
			cnt[0] /= static_cast<Storage::real>(cl.size()/3);
			cnt[1] /= static_cast<Storage::real>(cl.size()/3);
			cnt[2] /= static_cast<Storage::real>(cl.size()/3);
			cntv /= static_cast<Storage::real>(cl.size()/3);
			for(INMOST_DATA_ENUM_TYPE q = 0; q < cl.size(); q+=3) 
			{
				if( visualization_tag.isValid() ) 
				{
					//CommonColorBar->pick_color(cntv).set_color();
					glTexCoord1d(CommonColorBar->pick_texture(cntv));
				}
				glVertex3dv(cnt);
				if( visualization_tag.isValid() ) 
				{
					//CommonColorBar->pick_color(clv[q/3]).set_color();
					glTexCoord1d(CommonColorBar->pick_texture(clv[q/3]));
				}
				glVertex3dv(&cl[q]);
				if( visualization_tag.isValid() )
				{
					//CommonColorBar.pick_color(clv[(q/3+1)%clv.size()]).set_color();
					glTexCoord1d(CommonColorBar->pick_texture(clv[(q/3+1)%clv.size()]));
				}
				glVertex3dv(&cl[(q+3)%cl.size()]);
			}
		}
		glEnd();
		if( visualization_tag.isValid() ) CommonColorBar->UnbindTexture();
	}
	void draw_clip_edges(INMOST_DATA_ENUM_TYPE pace)
	{
    if( drawedges )
    {
		  glBegin(GL_LINES);
		  for(INMOST_DATA_ENUM_TYPE k = 0; k < cells.size(); k+=pace) if( cells[k]->GetMarker(marker))
		  {
			  Storage::real_array cl = cells[k]->RealArray(clips);
			
			  for(INMOST_DATA_ENUM_TYPE q = 0; q < cl.size(); q+=3) 
			  {
				  glVertex3dv(&cl[q]);
				  glVertex3dv(&cl[(q+3)%cl.size()]);
			  }
		  }
		  glEnd();
    }
	}
	INMOST_DATA_ENUM_TYPE size() {return (INMOST_DATA_ENUM_TYPE)cells.size();}
} * oclipper = NULL;

class bnd_clipper
{
	Tag clip_state;
	kdtree * tree;
	Mesh * mm;
	HandleType * faces;
	INMOST_DATA_ENUM_TYPE nfaces;
public:
	~bnd_clipper()
	{
		mm->DeleteTag(clip_state);
		delete tree;
		delete [ ]faces;
	}
	bnd_clipper(Mesh * m , HandleType * _faces, INMOST_DATA_ENUM_TYPE size)
	{
		mm = m;
		clip_state = m->CreateTag("CLIP_FACE_STATE",DATA_INTEGER,FACE,NONE,1);
		nfaces = size;
		faces = new HandleType[nfaces];
		for(INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k++) faces[k] = _faces[k];
		tree = new kdtree(mm,faces,nfaces);
	}
	double compute_value(Node n1, Node n2, Storage::real * c1, Storage::real * c2, Storage::real * pnt)
	{
		if( visualization_tag.isValid() )
		{
			Storage::real d1,d2,t;
			d1 = sqrt((pnt[0]-c1[0])*(pnt[0]-c1[0])+(pnt[1]-c1[1])*(pnt[1]-c1[1])+(pnt[2]-c1[2])*(pnt[2]-c1[2]));
			d2 = sqrt((c2[0]-c1[0])*(c2[0]-c1[0])+(c2[1]-c1[1])*(c2[1]-c1[1])+(c2[2]-c1[2])*(c2[2]-c1[2]));
			t = d1/d2; //(pnt == c2, t = 1 : pnt == c1, t = 0)
			return n1->RealDF(visualization_tag)*(1-t)+n2->RealDF(visualization_tag)*t;
		}
		else return 0.f;
	}
	void clip_plane(Storage::real p[3], Storage::real n[3])
	{
		tree->intersect_plane_face(clip_state,p,n);
	}
	void gen_clip(std::vector<face2gl> & out )
	{
		INMOST_DATA_ENUM_TYPE pace = std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,nfaces/100));
		for(INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k++)
		{
			int state = mm->IntegerDF(faces[k],clip_state);
			if( state == CLIP_FACE_INSIDE )
			{
				ElementArray<Node> nodes = Face(mm,faces[k])->getNodes();
				face2gl f;
				f.set_color(0.6,0.6,0.6,1);
				for(INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++) 
				{
					f.add_vert(&nodes[q].Coords()[0]);
					if( visualization_tag.isValid() ) 
					{
						//f.add_color(CommonColorBar->pick_color(nodes[q].RealDF(visualization_tag)));
						f.add_texcoord(CommonColorBar->pick_texture(nodes[q].RealDF(visualization_tag)));
					}
				}
				f.compute_center();
				f.set_elem(GetHandleElementType(faces[k]),GetHandleID(faces[k]));
				if( k%pace == 0 ) f.set_flag(true);
				out.push_back(f);
			}
			else if( state == CLIP_FACE_INTERSECT )
			{
				face2gl f;
				f.set_color(0.6,0.6,0.6,1);
				ElementArray<Node> nodes = Face(mm,faces[k])->getNodes();
				dynarray<bool,64> nodepos(nodes.size());
				dynarray<Storage::real,64> faceverts;
				Storage::real_array coords = nodes[0].Coords();
				for(INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++)
				{
					coords = nodes[q].Coords();
					Storage::real dot = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
					nodepos[q] = dot < 1.0e-10;
				}
				for(INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++)
				{
					if( nodepos[q] )
					{
						coords = nodes[q].Coords();
						f.add_vert(&coords[0]);
						if( visualization_tag.isValid() ) 
						{
							//f.add_color(CommonColorBar->pick_color(nodes[q].RealDF(visualization_tag)));
							f.add_texcoord(CommonColorBar->pick_texture(nodes[q].RealDF(visualization_tag)));
						}
					}
					if( nodepos[q] != nodepos[(q+1)%nodes.size()] )
					{
						Storage::real_array sp0 = nodes[q].Coords();
						Storage::real_array sp1 = nodes[(q+1)%nodes.size()].Coords();
						Storage::real node[3];
						if( clip_plane_edge(&sp0[0],&sp1[0],p,n,node) > CLIP_NONE) 
						{
							f.add_vert(node);
							if( visualization_tag.isValid() ) 
							{
								//f.add_color(CommonColorBar->pick_color(compute_value(nodes[q],nodes[(q+1)%nodes.size()],&sp0[0],&sp1[0],node)));
								f.add_texcoord(CommonColorBar->pick_texture(compute_value(nodes[q],nodes[(q+1)%nodes.size()],&sp0[0],&sp1[0],node)));
							}
						}
					}
				}
				f.compute_center();
				f.set_elem(GetHandleElementType(faces[k]),GetHandleID(faces[k]));
				if( k%pace == 0 ) f.set_flag(true);
				out.push_back(f);
			}
		}
	}
	void draw_clip(INMOST_DATA_ENUM_TYPE pace)
	{
		for(INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k+=pace)
		{
			int state = CLIP_FACE_NONE;
			ElementArray<Node> nodes = boundary_faces[k]->getNodes();
			Storage::real_array coords = nodes[0].Coords();
			Storage::real dot0 = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
			if( dot0 <= 0.0 ) state = CLIP_FACE_INSIDE; else state = CLIP_FACE_OUTSIDE;
			for(INMOST_DATA_ENUM_TYPE q = 1; q < nodes.size(); ++q)
			{
				coords = nodes[q].Coords();
				Storage::real dot = n[0]*(coords[0]-p[0])+n[1]*(coords[1]-p[1])+n[2]*(coords[2]-p[2]);
				if( dot*dot0 <= 0.0 ) 
				{
					state = CLIP_FACE_INTERSECT;
					break;
				}
			}
			if( state == CLIP_FACE_INSIDE )
			{
				ElementArray<Node> nodes = Face(mm,faces[k])->getNodes();
				glBegin(GL_POLYGON);
				for(INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++) 
				{
					if( visualization_tag.isValid() ) 
					{
						//CommonColorBar->pick_color(nodes[q].RealDF(visualization_tag)).set_color();
						glTexCoord1f(CommonColorBar->pick_texture(nodes[q].RealDF(visualization_tag)));
					}
					glVertex3dv(&nodes[q].Coords()[0]);
					
				}
				glEnd();
			}
			mm->IntegerDF(faces[k],clip_state) = state;
		}
	}
	void draw_clip_edges(INMOST_DATA_ENUM_TYPE pace)
	{
		for(INMOST_DATA_ENUM_TYPE k = 0; k < nfaces; k+=pace)
		{
			if( mm->IntegerDF(faces[k],clip_state) == CLIP_FACE_INSIDE )
			{
				ElementArray<Node> nodes = Face(mm,faces[k])->getNodes();
				glBegin(GL_LINE_LOOP);
				for(INMOST_DATA_ENUM_TYPE q = 0; q < nodes.size(); q++) glVertex3dv(&nodes[q].Coords()[0]);
				glEnd();
			}
		}
	}
	INMOST_DATA_ENUM_TYPE size() {return nfaces;}
} * bclipper = NULL;

class kdtree_picker
{
	struct entry
	{
		int index;
		float xyz[3];
	} * set;
	INMOST_DATA_ENUM_TYPE size;
	Storage::real bbox[6];
	kdtree_picker * children;
	static int cmpElements0(const void * a,const void * b) 
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[0];
		float bd = eb->xyz[0];
		return (ad > bd) - (ad < bd);
	}
	static int cmpElements1(const void * a,const void * b) 
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[1];
		float bd = eb->xyz[1];
		return (ad > bd) - (ad < bd);
	}
	static int cmpElements2(const void * a,const void * b) 
	{
		const entry * ea = ((const entry *)a);
		const entry * eb = ((const entry *)b);
		float ad = ea->xyz[2];
		float bd = eb->xyz[2];
		return (ad > bd) - (ad < bd);
	}
	inline static unsigned int flip(const unsigned int * fp)
	{
		unsigned int mask = -((int)(*fp >> 31)) | 0x80000000;
		return *fp ^ mask;
	}
#define _0(x)	(x & 0x7FF)
#define _1(x)	(x >> 11 & 0x7FF)
#define _2(x)	(x >> 22 )
	void radix_sort(int dim, struct entry * temp)
	{
		unsigned int i;
		const unsigned int kHist = 2048;
		unsigned int  b0[kHist * 3];
		unsigned int *b1 = b0 + kHist;
		unsigned int *b2 = b1 + kHist;
		memset(b0,0,sizeof(unsigned int)*kHist*3);
		for (i = 0; i < size; i++) 
		{
			unsigned int fi = flip((unsigned int *)&set[i].xyz[dim]);
			++b0[_0(fi)]; ++b1[_1(fi)]; ++b2[_2(fi)];
		}
		{
			unsigned int sum0 = 0, sum1 = 0, sum2 = 0;
			for (i = 0; i < kHist; i++) 
			{
				b0[kHist-1] = b0[i] + sum0; b0[i] = sum0 - 1; sum0 = b0[kHist-1];
				b1[kHist-1] = b1[i] + sum1; b1[i] = sum1 - 1; sum1 = b1[kHist-1];
				b2[kHist-1] = b2[i] + sum2; b2[i] = sum2 - 1; sum2 = b2[kHist-1];
			}
		}
		for (i = 0; i < size; i++) temp[++b0[_0(flip((unsigned int *)&set[i].xyz[dim]))]] = set[i];
		for (i = 0; i < size; i++) set[++b1[_1(flip((unsigned int *)&temp[i].xyz[dim]))]] = temp[i];
		for (i = 0; i < size; i++) temp[++b2[_2(flip((unsigned int *)&set[i].xyz[dim]))]] = set[i];
		for (i = 0; i < size; i++) set[i] = temp[i];
	}
	void kdtree_build(int dim, std::vector<face2gl> & in, struct entry * temp)
	{
		if( size > 1 )
		{
			if( size > 128 ) radix_sort(dim,temp); else 
			switch(dim)
			{
			case 0: qsort(set,size,sizeof(entry),cmpElements0);break;
			case 1: qsort(set,size,sizeof(entry),cmpElements1);break;
			case 2: qsort(set,size,sizeof(entry),cmpElements2);break;
			}
			children = static_cast<kdtree_picker *>(malloc(sizeof(kdtree_picker)*2));//new kdtree_picker[2];
			children[0].children = NULL;
			children[0].set = set;
			children[0].size = size/2;
			children[1].children = NULL;
			children[1].set = set+size/2;
			children[1].size = size - size/2;
			children[0].kdtree_build((dim+1)%3,in,temp);
			children[1].kdtree_build((dim+1)%3,in,temp);
			for(int k = 0; k < 3; k++)
			{
				bbox[0+2*k] = std::min(children[0].bbox[0+2*k],children[1].bbox[0+2*k]);
				bbox[1+2*k] = std::max(children[0].bbox[1+2*k],children[1].bbox[1+2*k]);
			}
		}
		else 
		{
			assert(size == 1);
			bbox[0] = bbox[2] = bbox[4] = 1.0e20;
			bbox[1] = bbox[3] = bbox[5] = -1.0e20;
			for(INMOST_DATA_ENUM_TYPE k = 0; k < in[set[0].index].size(); ++k)
			{
				double * coords = in[set[0].index].get_vert(k);
				for(INMOST_DATA_ENUM_TYPE q = 0; q < 3; q++)
				{
					bbox[q*2+0] = std::min(bbox[q*2+0],coords[q]);
					bbox[q*2+1] = std::max(bbox[q*2+1],coords[q]);
				}
			}
		}
	}
	void clear_children() { if( children ) {children[0].clear_children(); children[1].clear_children(); free(children);}}
	int raybox(double pos[3], double ray[3], double closest)
	{
		double tnear = -1.0e20, tfar = 1.0e20, t1,t2,c;
		for(int i = 0; i < 3; i++)
		{
			if( fabs(ray[i]) < 1.0e-15 )
			{
				if( pos[i] < bbox[i*2] || pos[i] > bbox[i*2+1] )
					return 0;
			}
			else
			{
				t1 = (bbox[i*2+0] - pos[i])/ray[i];
				t2 = (bbox[i*2+1] - pos[i])/ray[i];
				if( t1 > t2 ) {c = t1; t1 = t2; t2 = c;}
				if( t1 > tnear ) tnear = t1;
				if( t2 < tfar ) tfar = t2;
				if( tnear > closest ) return 0;
				if( tnear > tfar ) return 0;
				if( tfar < 0 ) return 0;
			}
		}
		return 1;
	}
	void sub_intersect_ray_faces(std::vector<face2gl> & in, double p[3], double dir[3], std::pair<double,int> & closest)
	{
		if( size == 1 )
		{
			face2gl & f = in[set[0].index];
			double * tri[3], btri[3][3], dot[3], prod[3][3], norm[3], d, proj[3];
			tri[0] = f.get_center();
			for(INMOST_DATA_ENUM_TYPE i = 0; i < f.size(); i++)
			{
				INMOST_DATA_ENUM_TYPE j = (i+1)%f.size();
				tri[1] = f.get_vert(i);
				tri[2] = f.get_vert(j);
				norm[0] = (tri[2][1]-tri[0][1])*(tri[1][2]-tri[0][2]) - (tri[2][2]-tri[0][2])*(tri[1][1]-tri[0][1]);
				norm[1] = (tri[2][2]-tri[0][2])*(tri[1][0]-tri[0][0]) - (tri[2][0]-tri[0][0])*(tri[1][2]-tri[0][2]);
				norm[2] = (tri[2][0]-tri[0][0])*(tri[1][1]-tri[0][1]) - (tri[2][1]-tri[0][1])*(tri[1][0]-tri[0][0]);
				d = norm[0]*(tri[0][0]-p[0])+norm[1]*(tri[0][1]-p[1])+norm[2]*(tri[0][2]-p[2]);
				d /= norm[0]*dir[0] + norm[1]*dir[1] + norm[2]*dir[2];
				proj[0] = p[0] + d*dir[0];
				proj[1] = p[1] + d*dir[1];
				proj[2] = p[2] + d*dir[2];

				if( d > closest.first ) break;

				for(int k = 0; k < 3; k++)
				{
					btri[k][0] = tri[k][0] - proj[0];
					btri[k][1] = tri[k][1] - proj[1];
					btri[k][2] = tri[k][2] - proj[2];
				}
				for(int k = 0; k < 3; k++)
				{
					int l = (k+1)%3;
					prod[k][0] = btri[k][1]*btri[l][2] - btri[k][2]*btri[l][1];
					prod[k][1] = btri[k][2]*btri[l][0] - btri[k][0]*btri[l][2];
					prod[k][2] = btri[k][0]*btri[l][1] - btri[k][1]*btri[l][0];
				}

				for(int k = 0; k < 3; k++)
				{
					int l = (k+1)%3;
					dot[k] = prod[k][0]*prod[l][0]+prod[k][1]*prod[l][1]+prod[k][2]*prod[l][2];
				}
				if( dot[0] >= 0 && dot[1] >= 0 && dot[2] >= 0 ) 
				{
					closest.first = d;
					closest.second = set[0].index;
					break; //don't expect anything better here
				}
			}
		}
		else
		{
			if( raybox(p,dir,closest.first) )
			{
				children[0].sub_intersect_ray_faces(in,p,dir,closest);
				children[1].sub_intersect_ray_faces(in,p,dir,closest);
			}
		}
	}
public:
	kdtree_picker() : set(NULL), size(0), children(NULL) {}
	~kdtree_picker() { clear_children(); delete [] set;}
	kdtree_picker(std::vector<face2gl> & in)
	{
		size = (unsigned)in.size();
		set = new struct entry[size];
		for(INMOST_DATA_ENUM_TYPE k = 0; k < in.size(); ++k)
		{
			set[k].index = k;
			in[k].get_center(set[k].xyz);
		}
		struct entry * temp = new struct entry[size];
		kdtree_build(0,in,temp);
		delete [] temp;
	}
	int ray_faces(std::vector<face2gl> & in, double p[3], double dir[3])
	{
		std::pair<double, int> closest(1.0e20,-1);
		sub_intersect_ray_faces(in,p,dir,closest);
		return closest.second;
	}
};

class picker
{
	kdtree_picker * tree;
	std::vector<face2gl> * faces;
public:
	~picker() {delete tree;}
	picker(std::vector<face2gl> & _faces)
	{
		faces = &_faces;
		tree = new kdtree_picker(*faces);
	}
	int select(double p[3], double ray[3]) {return tree->ray_faces(*faces,p,ray);};
} * current_picker = NULL;


void set_matrix3d()
{
	double aspect = (double)width/(double)height;
	double side = std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
	//double center[3] = { (sleft+sright)*0.5, (sbottom+stop)*0.5, (sfar+snear)*0.5};
	const double sc = 2;
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	//~ glOrtho(center[0]-sc*side*zoom*aspect,center[0]+sc*zoom*side*aspect,
		//~ center[1]-sc*side*zoom,center[1]+sc*zoom*side,
		//~ center[2]-sc*side*100,center[2]+sc*side*100);
	if( !perspective )
	{
		glOrtho(-sc*side*zoom*aspect,sc*side*zoom*aspect,
				-sc*side*zoom,sc*side*zoom,
				-sc*side*100,sc*side*100);
	}
	else
	{
		const double pi = 3.1415926535897932384626433832795;
		const double znear = 1;
		const double zfar  = 1000000.0;
		const double fH = znear * tan( 60.0 / 360.0 * pi);
		double fW = fH * aspect;
		glFrustum(-fW,fW,-fH,fH,znear,zfar);
	}
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity();
}

void set_matrix2d()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1,1,-1,1,-1,1);
	glMatrixMode(GL_MODELVIEW);
}

void reshape(int w, int h)
{
	width = w;
	height = h;
	set_matrix3d();
	glViewport(0, 0, w, h);
}


int actionstate  = 0;
double mymx = 0;
double mymy = 0;

void myclickmotion(int nmx, int nmy) // Mouse
{
	double lmx = 2.*(nmx/(double)width - 0.5),lmy = 2.*(0.5 - nmy/(double)height), dmx = lmx-mymx, dmy = lmy - mymy;
	if( actionstate == 1 ) //middle button
	{
		double shiftmod[3] = {0,0,0};
		shiftmod[0] += dmx*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear);
		shiftmod[1] += dmy*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear);
		if( planecontrol )
		{
			
			rotatevector_from_stack((double*)shiftmod);
			p[0] += shiftmod[0];
			p[1] += shiftmod[1];
			p[2] += shiftmod[2];
			clipupdate = true;
		}
		else
		{
			//shiftmod[0] += dmx*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
			//shiftmod[1] += dmy*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
			rotatevector((double*)shiftmod);
			shift[0] += shiftmod[0];
			shift[1] += shiftmod[1];
			shift[2] += shiftmod[2];
			bndupdate = true;
		}
		glutPostRedisplay();
		mymx = lmx;
		mymy = lmy;
	}
	else if( actionstate == 2 ) //right button
	{
		if( planecontrol )
		{
			double shiftmod[3] = {0,0,0};
			shiftmod[2] -= dmy*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear);
			rotatevector_from_stack((double*)shiftmod);
			p[0] += shiftmod[0];
			p[1] += shiftmod[1];
			p[2] += shiftmod[2];
			clipupdate = true;
		}
		else
		{
			zoom *= expf(-dmy);
			reshape(width,height);
			double shiftmod[3] = {0,0,0};
			shiftmod[2] += dmx*zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear);
			rotatevector((double*)shiftmod);
			shift[0] += shiftmod[0];
			shift[1] += shiftmod[1];
			shift[2] += shiftmod[2];
			bndupdate = true;
		}
		glutPostRedisplay();
		mymx = lmx;
		mymy = lmy;
	}
	else if( actionstate == 3 ) //left buttion
	{
		clickmotion(nmx,nmy);
		if( planecontrol )
		{
			//quatget(n);
			reverse_rotatevector_from_stack((double *)n);
			reverse_rotatevector((double *)n);
			rotatevector_from_stack((double *)n);
			clipupdate = true;
			quatinit();
		}
		else bndupdate = true;
	}
}
void mymotion(int nmx, int nmy) // Mouse
{
	motion(nmx,nmy);
	mymx = 2.*(nmx/(double)width - 0.5);
	mymy = 2.*(0.5 - nmy/(double)height);
	if( current_picker != NULL ) glutPostRedisplay();
}

void myclick(int b, int s, int nmx, int nmy) // Mouse
{
	if( b == GLUT_LEFT_BUTTON )
	{
		if( s == GLUT_DOWN )
		{
			actionstate = 3;
		}
		else
		{
			actionstate = 0;
		}
		click(b,s,nmx,nmy);
	}
	else if( b == GLUT_MIDDLE_BUTTON )
	{
		if( s == GLUT_DOWN )
		{
			actionstate = 1;
			interactive = true;
		}
		else
		{
			actionstate = 0;
			interactive = false;
		}
		mymx = 2.*(nmx/(double)width - 0.5);
		mymy = 2.*(0.5 - nmy/(double)height);
	}
	else if( b == GLUT_RIGHT_BUTTON )
	{
		if( s == GLUT_DOWN )
		{
			actionstate = 2;
			interactive = true;
		}
		else
		{
			actionstate = 0;
			interactive = false;
		}
		mymx = 2.*(nmx/(double)width - 0.5);
		mymy = 2.*(0.5 - nmy/(double)height);
	}
	
	
	glutPostRedisplay();
}


class Input
{
public:
	enum InputType {Double,Integer,String};
private:
	std::string str;
	std::string comment;
	void * input_link;
	InputType type;
	bool done;
	bool canceled;
	
public:
	Input(int * val, std::string comment) : comment(comment) {input_link = val; type = Integer; canceled = false; done = false; str = "";}
	Input(double * val, std::string comment) : comment(comment) {input_link = val; type = Double; canceled = false; done = false; str = "";}
	Input(char * val, std::string comment) : comment(comment) {input_link = val; type = String; canceled = false; done = false; str = "";}
	Input(void * link, InputType type, std::string comment) : comment(comment), input_link(link), type(type) {canceled = false; done = false; str = "";}
	Input(const Input & other):str(other.str),comment(other.comment), input_link(other.input_link),  type(other.type), done(other.done), canceled(other.canceled) {}
	Input & operator =(Input const & other) {comment = other.comment; input_link = other.input_link; str = other.str; type = other.type; canceled = other.canceled; done = other.done; return *this;}
	~Input() {}
	void KeyPress(char c) 
	{
		if( c == 13 )
		{
			
			done = true;
			if( type == Double ) *((double *)input_link) = atof(str.c_str());
			else if( type == Integer ) *((int *)input_link) = atoi(str.c_str());
			else if( type == String ) strcpy((char *)input_link,str.c_str());
			glutPostRedisplay();
		}
		else if( c == 8 )
		{
			if( !str.empty() ) str.erase(str.size()-1);
			glutPostRedisplay();
		}
		else if( c == 27 )
		{
			canceled = true;
			done = true;
			glutPostRedisplay();
		}
		else if( type == String || ( (c >= '0' && c <= '9') || ((str.empty() || tolower(*str.rbegin()) == 'e') && (c=='+' || c=='-')) || (type == Double && (c=='.' || c=='e' || c == 'E'))) )
		{
			str += c;
			glutPostRedisplay();
		}
	}
	bool Done() {return done;}
	bool Canceled() {return canceled;}
	void Draw()
	{
		float h = 24.0f/(float)height;
		
		glColor3f(1,1,1);
		glBegin(GL_QUADS);
		glVertex3f(-0.99,-0.99,1);
		glVertex3f(-0.99,-0.99+h,1);
		glVertex3f(0.99,-0.99+h,1);
		glVertex3f(0.99,-0.99,1);
		glEnd();
		glColor3f(0,0,0);
		glBegin(GL_LINE_LOOP);
		glVertex3f(-0.99,-0.99,1);
		glVertex3f(-0.99,-0.99+h,1);
		glVertex3f(0.99,-0.99+h,1);
		glVertex3f(0.99,-0.99,1);
		glEnd();
		
		glColor4f(0,0,0,1);
		glRasterPos2d(-0.985,-0.985);
		//printtext(str.c_str());
		char oldval[4096];
		if( type == Double ) sprintf(oldval,"%g",*(double*)input_link);
		else if( type == Integer ) sprintf(oldval,"%d",*(int*)input_link);
		else if( type == String ) sprintf(oldval,"%s",(char *)input_link);
		printtext("input number (%s[%s]:%s): %s",comment.c_str(),oldval,type == Integer ? "integer": (type == Double ? "double" : "string"), str.c_str());
	}
	std::string GetString() {return str;}
} * CommonInput = NULL;



void keyboard(unsigned char key, int x, int y)
{
	(void) x;
	(void) y;
	if( glutGetModifiers() & (GLUT_ACTIVE_CTRL) )
		std::cout << "pressed " << ((char)(key)) << " int " << ((int)key) << " ctrl " << (glutGetModifiers() & GLUT_ACTIVE_CTRL ? "yes" : "no") << " shift " << (glutGetModifiers() & GLUT_ACTIVE_SHIFT ? "yes" : "no") << " alt " << (glutGetModifiers() & GLUT_ACTIVE_ALT ? "yes" : "no") << std::endl;
	if( CommonInput != NULL )
	{
		if( (key == 'v' || key == 'V' || key == 22) && (glutGetModifiers() & GLUT_ACTIVE_CTRL) ) //paste
		{
			std::string paste = getTextFromPasteboard();
			std::cout << "paste: " << paste << std::endl;
			if( !paste.empty() )
			{
				for(int k = 0; k < paste.length(); ++k)
					CommonInput->KeyPress(paste[k]);
			}
		}
		else if( (key == 'c' || key == 'C' || key == 3) && (glutGetModifiers() & GLUT_ACTIVE_CTRL) ) //copy
		{
			std::string copy = CommonInput->GetString();
			std::cout << "copy: " << copy << std::endl;
			setTextToPasteboard(copy);
		}
		else CommonInput->KeyPress(key);
		return;
	}
	if( key == 27 )
	{
		if( oclipper ) delete oclipper;
		if( bclipper ) delete bclipper;
		if( current_picker ) delete current_picker;
		delete mesh;
		exit(-1);
	}
	else if( key == '=' || key == '+')
	{
		zoom /= 1.1;
		reshape(width,height);
		glutPostRedisplay();
		interactive = true;
		//reset_timer = Timer();
	}
	else if( key == '_' || key == '-')
	{
		zoom *= 1.1;	
		reshape(width,height);
		glutPostRedisplay();
		interactive = true;
		//reset_timer = Timer();
	}
  else if( key == 'e' )
  {
    drawedges = (drawedges+1)%3;
    glutPostRedisplay();
  }
	else if( key == 'w' )
	{
		double shiftmod[3] = {0,0,0};
		shiftmod[1] -= 0.03f*expf(zoom-1)*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
		rotatevector((double*)shiftmod);
		shift[0] += shiftmod[0];
		shift[1] += shiftmod[1];
		shift[2] += shiftmod[2];
		glutPostRedisplay();
		interactive = true;
		bndupdate = true;
	}
	else if( key == 's' )
	{
		if( !planecontrol )
		{
			double shiftmod[3] = {0,0,0};
			shiftmod[1] += 0.03f*expf(zoom-1)*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
			rotatevector((double*)shiftmod);
			shift[0] += shiftmod[0];
			shift[1] += shiftmod[1];
			shift[2] += shiftmod[2];
			glutPostRedisplay();
			interactive = true;
			bndupdate = true;
		}
	}
	else if( key == 'a' )
	{
		if( !planecontrol )
		{
			double shiftmod[3] = {0,0,0};
			shiftmod[0] += 0.03f*expf(zoom-1)*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
			rotatevector((double*)shiftmod);
			shift[0] += shiftmod[0];
			shift[1] += shiftmod[1];
			shift[2] += shiftmod[2];
			glutPostRedisplay();
			interactive = true;
			bndupdate = true;
		}
	}
	else if( key == 'd' )
	{
		if( !planecontrol )
		{
			double shiftmod[3] = {0,0,0};
			shiftmod[0] -= 0.03f*expf(zoom-1)*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5;
			rotatevector((double*)shiftmod);
			shift[0] += shiftmod[0];
			shift[1] += shiftmod[1];
			shift[2] += shiftmod[2];
			glutPostRedisplay();
			interactive = true;
			bndupdate = true;
		}
	}
	else if( key == 'r' )
	{
		if( planecontrol )
		{
			p[0] = (sleft+sright)*0.5;
			p[1] = (sbottom+stop)*0.5;
			p[2] = (sfar+snear)*0.5;
			n[0] = 0;
			n[1] = 0;
			n[2] = 1;
			quatinit();
			clipupdate = true;
		}
		else
		{
			shift[0] = 0.0f;
			shift[1] = 0.0f;
			shift[2] = 0.0f;
			zoom = 1;
			quatinit();
			bndupdate = true;
		}
		glutPostRedisplay();
		interactive = true;
	}
	else if( key == 'z' )
	{
		printf("Enter point of plane (x,y,z):\n");
		scanf("%lf %lf %lf",p,p+1,p+2);
		printf("Enter normal of plane (x,y,z):\n");
		scanf("%lf %lf %lf",n,n+1,n+2);
		Storage::real l = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
		if( l )
		{
			n[0] /= l;
			n[1] /= l;
			n[2] /= l;
		}
		if( oclipper ) 
		{
			oclipper->clip_plane(p,n);
			bclipper->clip_plane(p,n);
			clipboxupdate = true;
		}
			
		glutPostRedisplay();
	}
	else if( key == 'x' )
	{
		draw_volumetric = !draw_volumetric;
		glutPostRedisplay();
	}
	else if( key == 'b' )
	{
		boundary = !boundary;
		glutPostRedisplay();
	}
	else if( key == 'l' )
	{
		planecontrol = !planecontrol;
		if( planecontrol ) 
		{
			quatpush();
			quatinit();
		}
		else quatpop();
	}
	else if( key == 'n' )
	{
		if( CommonInput == NULL ) CommonInput = new Input(&amplitude, "Amplitude");
		glutPostRedisplay();
	}
	else if( key == 'm' )
	{
		if( CommonInput == NULL ) CommonInput = new Input(&radius, "Radius");
		glutPostRedisplay();
	}
	else if( key == ',' || key == '<' )
	{
		std::cout << "positive shift" << std::endl;
		Storage::real radius2 = radius*radius;
		for(Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
		{
			Storage::real_array coords = it->Coords();
			Storage::real proj[3] = {0,0,0};
			Storage::real d = -(n[0]*(coords[0]-p[0]) + n[1]*(coords[1]-p[1]) + n[2]*(coords[2]-p[2]));
			proj[0] = coords[0] + d*n[0];
			proj[1] = coords[1] + d*n[1];
			proj[2] = coords[2] + d*n[2];
			d = (n[0]*(proj[0]-p[0]) + n[1]*(proj[1]-p[1]) + n[2]*(proj[2]-p[2]));
			assert(fabs(d) < 1.0e-7); // check that it is indeed a projection
			Storage::real dist = sqrt((p[0]-proj[0])*(p[0]-proj[0])+(p[1]-proj[1])*(p[1]-proj[1])+(p[2]-proj[2])*(p[2]-proj[2]));
			Storage::real delta = exp(-dist*dist/radius2);
			coords[0] += n[0]*delta*amplitude;
			coords[1] += n[1]*delta*amplitude;
			coords[2] += n[2]*delta*amplitude;
		}
		all_boundary.clear();
		INMOST_DATA_ENUM_TYPE pace = std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,(unsigned)boundary_faces.size()/100));
		for(INMOST_DATA_ENUM_TYPE k = 0; k < boundary_faces.size(); k++) 
		{
			all_boundary.push_back(DrawFace(boundary_faces[k]));
			if( k%pace == 0 ) all_boundary.back().set_flag(true);
		}
		if( oclipper ) delete oclipper;
		if( bclipper ) delete bclipper;
		bclipper = new bnd_clipper(mesh,boundary_faces.data(),(int)boundary_faces.size());
		oclipper = new clipper(mesh);
		clipupdate = true;
		glutPostRedisplay();
	}
	else if( key == '.' || key == '>' )
	{
		std::cout << "negative shift" << std::endl;
		Storage::real radius2 = radius*radius;
		for(Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
		{
			Storage::real_array coords = it->Coords();
			Storage::real proj[3] = {0,0,0};
			Storage::real d = -(n[0]*(coords[0]-p[0]) + n[1]*(coords[1]-p[1]) + n[2]*(coords[2]-p[2]));
			proj[0] = coords[0] + d*n[0];
			proj[1] = coords[1] + d*n[1];
			proj[2] = coords[2] + d*n[2];
			d = (n[0]*(proj[0]-p[0]) + n[1]*(proj[1]-p[1]) + n[2]*(proj[2]-p[2]));
			assert(fabs(d) < 1.0e-7); // check that it is indeed a projection
			Storage::real dist = sqrt((p[0]-proj[0])*(p[0]-proj[0])+(p[1]-proj[1])*(p[1]-proj[1])+(p[2]-proj[2])*(p[2]-proj[2]));
			Storage::real delta = exp(-dist*dist/radius2);
			coords[0] -= n[0]*delta*amplitude;
			coords[1] -= n[1]*delta*amplitude;
			coords[2] -= n[2]*delta*amplitude;
		}
		all_boundary.clear();
		INMOST_DATA_ENUM_TYPE pace = std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,(unsigned)boundary_faces.size()/100));
		for(INMOST_DATA_ENUM_TYPE k = 0; k < boundary_faces.size(); k++) 
		{
			all_boundary.push_back(DrawFace(boundary_faces[k]));
			if( k%pace == 0 ) all_boundary.back().set_flag(true);
		}
		if( oclipper ) delete oclipper;
		if( bclipper ) delete bclipper;
		bclipper = new bnd_clipper(mesh,boundary_faces.data(),(unsigned)boundary_faces.size());
		oclipper = new clipper(mesh);
		clipupdate = true;
		glutPostRedisplay();
	}
	else if( key == 'p' )
	{
		perspective = !perspective;
		glutPostRedisplay();
	}
	else if( key == 'v' )
	{
		if( CommonInput == NULL ) 
		{
			PrintTags(mesh,CELL|FACE|EDGE|NODE);

			CommonInput = new Input(visualization_prompt, "Enter data for visualization as Element:Name:Component");
			visualization_prompt_active = 1;
			clipupdate = true;
			if( visualization_tag.isValid() ) visualization_tag =  mesh->DeleteTag(visualization_tag);
		}
		glutPostRedisplay();
	}
	else if( key == 'c' )
	{
		if( CommonInput == NULL ) 
		{
			CommonInput = new Input(visualization_prompt, "Enter data for color bounds as min:max");
			visualization_prompt_active = 2;
		}
		glutPostRedisplay();
	}
	else if( key == 'q' )
	{
		mesh->Save("mesh.vtk");
		mesh->Save("mesh.pmf");
		mesh->Save("mesh.xml");
	}
	else if( key == 't' )
		screenshot();
}

void keyboard2(unsigned char key, int x, int y)
{
	(void) x;
	(void) y;
	if( key == '=' || key == '+' ||  key == '_' || key == '-' || key == 'w' || key == 's' || key == 'a' || key == 'd' || key == 'r' || key == 'p' || key == 'z')
	{
		//printf("depressed\n");
		interactive = false;
		glutPostRedisplay();
	}
	
}

Tag face_center;

face2gl DrawFace(Element f)
{
	face2gl ret;
	ElementArray<Node> nodes = f->getNodes();

	if( f->nbAdjElements(CELL) == 0 ) ret.set_color(1,0,0,0.1);
	else ret.set_color(0,1,0,0.1);
	
	for(ElementArray<Node>::iterator kt = nodes.begin(); kt != nodes.end(); kt++)
		ret.add_vert(&(kt->Coords()[0]));
	ret.set_elem(f->GetElementType(),f->LocalID());
  ret.compute_center();
	return ret;
}

void whereami(double & cx, double & cy, double & cz)
{
   // Get the viewing matrix
   GLdouble modelview[16],projection[16];
   GLint viewport[4] = {0,0,1,1};
   glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
   glGetDoublev(GL_PROJECTION_MATRIX, projection);
   
   GLdouble outx, outy, outz;  // Var's to save the answer in

   gluUnProject(0.5, 0.5, 0.,
               modelview, projection, viewport,
               &outx, &outy, &outz);

   // Return the result.
   cx = outx;
   cy = outy;
   cz = outz;
}


void pick_mouse(double origin[3], double direction[3])
{
   // Get the viewing matrix
   GLdouble modelview[16],projection[16];
   GLint viewport[4] = {0,0,1,1};
   glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
   glGetDoublev(GL_PROJECTION_MATRIX, projection);
   
   GLdouble outx, outy, outz;  // Var's to save the answer in

   gluUnProject((mymx+1.0)/2.0, (mymy+1.0)/2.0, 0.,
               modelview, projection, viewport,
               &outx, &outy, &outz);

   // Return the result.
   origin[0] = outx;
   origin[1] = outy;
   origin[2] = outz;
   
   gluUnProject((mymx+1.0)/2.0, (mymy+1.0)/2.0, 1.,
               modelview, projection, viewport,
               &outx, &outy, &outz);
   direction[0] = outx - origin[0];
   direction[1] = outy - origin[1];
   direction[2] = outz - origin[2];
   double l = sqrt(direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);
   if( l ) 
   {
	   direction[0] /= l;
	   direction[1] /= l;
	   direction[2] /= l;
   }
}

void draw_screen()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	set_matrix3d();

	Storage::real mult = zoom*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5*0.1;
	if( perspective )
		glTranslated(0,0,-zoom*2*std::max(std::max( sright-sleft, stop-sbottom ), sfar-snear)*0.5);
	if( planecontrol ) 
		rotate_from_stack(); 
	else 
		rotate();

	//axis
	{
		glLineWidth(3.0);
		glBegin(GL_LINES);
		glColor3f(1,0,0);
		glVertex3d(0,0,0);
		glVertex3d(mult,0,0);
		glColor3f(0,1,0);
		glVertex3d(0,0,0);
		glVertex3d(0,mult,0);
		glColor3f(0,0,1);
		glVertex3d(0,0,0);
		glVertex3d(0,0,mult);
		glEnd();
		glLineWidth(1.0);
	}
	//glPointSize(1);
	
	//glTranslated(-(sleft+sright)*0.5+shift[0],-(sbottom+stop)*0.5 + shift[1],-(snear+sfar)*0.5 + shift[2]);
	
	glTranslated(shift[0],shift[1],shift[2]);


	//if( planecontrol )
	{
		glColor3f(0.6,0.4,0.2);
		glPointSize(5);
		glBegin(GL_POINTS);
		glVertex3dv(p);
		glEnd();
		glPointSize(1);

		glColor3f(0.2,0.4,0.6);
		glLineWidth(3.0);
		glBegin(GL_LINES);
		glVertex3dv(p);
		glVertex3d(p[0]+n[0]*mult*2,p[1]+n[1]*mult*2,p[2]+n[2]*mult*2);
		glEnd();
		glLineWidth(1.0);
	}

	

	double campos[3] = {0.5,0.5,0}, pickp[3], pickd[3];
	whereami(campos[0],campos[1],campos[2]);
	int picked = -1;

	
	//glTranslated((l+r)*0.5,(b+t)*0.5,(near+far)*0.5);

	{
		
		if( oclipper ) 
		{
			
			if( clipupdate ) 
			{
				if( current_picker != NULL ) {delete current_picker; current_picker = NULL;}
				oclipper->clip_plane(p,n);
				clipupdate = false;
				clipboxupdate = true;
			}
		
			if( !interactive && clipboxupdate )
			{
				clip_boundary.clear();
				oclipper->gen_clip(clip_boundary);
				bclipper->clip_plane(p,n);
				bclipper->gen_clip(clip_boundary);
				clipboxupdate = false;

				if( current_picker != NULL ) {delete current_picker; current_picker = NULL;}
				if( !clip_boundary.empty() ) current_picker = new picker(clip_boundary);
			}

			

			if( current_picker != NULL )
			{
				pick_mouse(pickp,pickd);
				picked = current_picker->select(pickp,pickd);
			}
	
			if( interactive && clipboxupdate )
			{
				//printf("draw1 %d %d\n",interactive, clipboxupdate);
				INMOST_DATA_ENUM_TYPE opace = !planecontrol ? std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,oclipper->size()/100)) : 1;
				INMOST_DATA_ENUM_TYPE bpace = std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,bclipper->size()/100));
				glColor4f(0.6,0.6,0.6,1);
				oclipper->draw_clip(opace);
				bclipper->draw_clip(bpace);
				glColor4f(0,0,0,1); 
				oclipper->draw_clip_edges(opace);
				bclipper->draw_clip_edges(bpace);
			}
			else
			{
				//printf("draw2 %d %d\n",interactive, clipboxupdate);
				if( interactive )
				{
					draw_faces_interactive(clip_boundary);
					glColor4f(0,0,0,1); 
					draw_edges_interactive(clip_boundary);
				}
				else
				{
					draw_faces(clip_boundary,picked);
					glColor4f(0,0,0,1); 
					draw_edges(clip_boundary,picked);
				}
			}
		}


		if( draw_volumetric )
		{
			if( bndupdate ) CommonVolumetricView->camera(campos,interactive);
			CommonVolumetricView->draw(interactive);
		}


    for(int k = 0; k < streamlines.size(); ++k)
      streamlines[k].Draw(true);//interactive);

		if( boundary )
		{
			glEnable(GL_BLEND);
			if( !interactive && bndupdate)
			{
				for(INMOST_DATA_ENUM_TYPE q = 0; q < all_boundary.size() ; q++) 
					all_boundary[q].compute_dist(campos);
				//std::sort(all_boundary.rbegin(),all_boundary.rend());
				face2gl::radix_sort_dist(all_boundary);


        for(INMOST_DATA_ENUM_TYPE q = 0; q < added_faces.size() ; q++) 
					added_faces[q].compute_dist(campos);
				//std::sort(all_boundary.rbegin(),all_boundary.rend());
				face2gl::radix_sort_dist(added_faces);
			}
			glColor4f(0,0,0,0.25); 
			if( interactive ) draw_edges_interactive(all_boundary);
			else draw_edges(all_boundary);
			

			if( interactive ) draw_faces_interactive_nc(all_boundary);
			else draw_faces_nc(all_boundary);
			
			glDisable(GL_BLEND);
		}

		if( !interactive && bndupdate) bndupdate = false;
	}





  glEnable(GL_BLEND);
  glColor4f(0,0,0,0.25); 
  draw_edges(added_faces);
  draw_faces_nc(added_faces);
  glDisable(GL_BLEND);

  
  if( !added_edges.empty() )
  {
    glColor3f(0,0,0);
    glBegin(GL_LINES);
    for(ElementArray<Edge>::iterator it = added_edges.begin(); it != added_edges.end(); ++it)
    {
      glVertex3dv(it->getBeg()->Coords().data());
      glVertex3dv(it->getEnd()->Coords().data());
    }
    glEnd();
  }

  glLineWidth(2.0);

  glBegin(GL_LINES);
  for(int k = 0; k < conormals.size(); k+=3)
  {
    glVertex3dv(&conormals[k]);
  }
  glEnd();

  glLineWidth(1.0);

  glPointSize(5.0);

  glColor3f(0,0,1);
  glBegin(GL_POINTS);
  for(int k = 0; k < harmonic_points.size(); k+=3)
    glVertex3dv(&harmonic_points[k]);
  glEnd();

  glColor3f(1,0,0);
  glBegin(GL_POINTS);
  for(int k = 0; k < dual_harmonic_points.size(); k+=9)
    glVertex3dv(&dual_harmonic_points[k+3]);
  glEnd();

  glPointSize(1.0);

  glColor3f(0.5,0.5,0.5);
  glBegin(GL_LINES);
  for(int k = 0; k < dual_harmonic_points.size(); k+=9)
  {
    glVertex3dv(&dual_harmonic_points[k+0]);
    glVertex3dv(&dual_harmonic_points[k+3]);
    glVertex3dv(&dual_harmonic_points[k+3]);
    glVertex3dv(&dual_harmonic_points[k+6]);
  }
  glEnd();




	if( picked != -1 )
	{
		glDisable(GL_DEPTH_TEST);
		glLoadIdentity();
		set_matrix2d();
		

		double top = 0.96, left = 0.25, interval = 0.04, bottom = top;
		
		Element e = clip_boundary[picked].get_elem(mesh);
		for(Mesh::iteratorTag t = mesh->BeginTag(); t != mesh->EndTag(); ++t) if( t->isDefined(e->GetElementType()) )
			if( e->HaveData(*t) )
				bottom -= interval;

		glColor3f(1,1,1);
		glBegin(GL_QUADS);
		glVertex2f(left-0.01,bottom-0.01);
		glVertex2f(left-0.01,0.99);
		glVertex2f(0.99,0.99);
		glVertex2f(0.99,bottom-0.01);
		glEnd();
		glColor3f(0,0,0);
		glBegin(GL_LINE_LOOP);
		glVertex2f(left-0.01,bottom-0.01);
		glVertex2f(left-0.01,0.99);
		glVertex2f(0.99,0.99);
		glVertex2f(0.99,bottom-0.01);
		glEnd();

		
		glColor3f(0.2,0.2,0.2);
		glRasterPos2f(left,top);
		printtext("%s %d",ElementTypeName(e->GetElementType()),e->LocalID());
		top -= interval;
		glColor3f(0.2,0.2,0.2);
		for(Mesh::iteratorTag t = mesh->BeginTag(); t != mesh->EndTag(); ++t) if( t->isDefined(e->GetElementType()) )
		{
			if( e->HaveData(*t) )
			{
				char str[4096];
				char temp[4096];
				str[0] = '\0';
				switch(t->GetDataType())
				{
				case DATA_INTEGER:
					{
						Storage::integer_array arr = e->IntegerArray(*t);
						for(INMOST_DATA_ENUM_TYPE k = 0; k < arr.size(); k++) 
						{
							sprintf(temp,"%s %d",str,arr[k]);
							strcpy(str,temp);
						}
						break;
					}
				case DATA_REAL:
					{
						Storage::real_array arr = e->RealArray(*t);
						for(INMOST_DATA_ENUM_TYPE k = 0; k < arr.size(); k++)
						{
							sprintf(temp,"%s %lf",str,arr[k]);
							strcpy(str,temp);
						}
						break;
					}
				case DATA_BULK:
					{
						Storage::bulk_array arr = e->BulkArray(*t);
						for(INMOST_DATA_ENUM_TYPE k = 0; k < arr.size(); k++) 
						{
							sprintf(temp,"%s %d",str,arr[k]);
							strcpy(str,temp);
						}
						break;
					}
				case DATA_REFERENCE:
					{
						Storage::reference_array arr = e->ReferenceArray(*t);
						for(INMOST_DATA_ENUM_TYPE k = 0; k < arr.size(); k++) 
						{
							if(arr.at(k) == InvalidHandle()) sprintf(temp,"%s NULL",str);
							else sprintf(temp,"%s %s:%d",str,ElementTypeName(arr[k]->GetElementType()),arr[k]->LocalID());
							strcpy(str,temp);
						}
						break;
					}
        case DATA_REMOTE_REFERENCE:
					{
						Storage::remote_reference_array arr = e->RemoteReferenceArray(*t);
						for(INMOST_DATA_ENUM_TYPE k = 0; k < arr.size(); k++) 
						{
              if(arr.at(k).first == NULL || arr.at(k).second == InvalidHandle()) sprintf(temp,"%s NULL",str);
              else sprintf(temp,"%s %p:%s:%d",str,arr[k]->GetMeshLink(),ElementTypeName(arr[k]->GetElementType()),arr[k]->LocalID());
							strcpy(str,temp);
						}
						break;
					}
                    case DATA_VARIABLE:
                    {
                        Storage::var_array arr = e->VariableArray(*t);
                        for(INMOST_DATA_ENUM_TYPE k = 0; k < arr.size(); k++)
                        {
                            std::stringstream stream;
                            stream << arr[k].GetValue() << " {[" << arr[k].GetRow().Size() << "] ";
                            for(INMOST_DATA_ENUM_TYPE q = 0; q < arr[k].GetRow().Size(); ++q)
                            {
                                stream << "(" << arr[k].GetRow().GetValue(q) << "," << arr[k].GetRow().GetIndex(q) << ") ";
                            }
                            stream << "}";
                            sprintf(temp,"%s %s",str,stream.str().c_str());
                            strcpy(str,temp);
                        }
                        break;
                    }
				}
				sprintf(temp,"%s %s %s",t->GetTagName().c_str(),DataTypeName(t->GetDataType()),str);
				strcpy(str,temp);
				glRasterPos2f(left,top);
				printtext(str);
				top -= interval;
			}
		}
		glEnable(GL_DEPTH_TEST);
	}

	if( CommonInput != NULL )
	{
		glDisable(GL_DEPTH_TEST);
		glLoadIdentity();
		set_matrix2d();
		CommonInput->Draw();
		if( CommonInput->Done() )
		{
			if( ! CommonInput->Canceled() )
			{
        if( visualization_prompt_active == 2 )
        {
					int k = 0, slen = (int)strlen(visualization_prompt);
					for(k = 0; k < slen; ++k) 
					{
						if( visualization_prompt[k] == ':' )
						{
							visualization_prompt[k] = '\0';
							break;
						}
					}
          if( k < slen && k+1 < slen )
          {
            double minv, maxv;
            visualization_prompt[k] = ':';
            minv = atof(visualization_prompt);
            maxv = atof(visualization_prompt+k+1);
            CommonColorBar->set_min(minv);
            CommonColorBar->set_max(maxv);
            clipupdate = true;
          }
          else printf("malformed string %s for color map bounds\n",visualization_prompt);
					visualization_prompt_active = 0;
        }
        else if( visualization_prompt_active == 1 )
				{
					char typen[1024],name[1024];
					unsigned comp;
					int k = 0,l, slen = (int)strlen(visualization_prompt);
					for(k = 0; k < slen; ++k) 
					{
						if( visualization_prompt[k] == ':' )
						{
							visualization_prompt[k] = '\0';
							break;
						}
					}
					for(l = k+1; l < slen; ++l) 
					{
						if( visualization_prompt[l] == ':' )
						{
							visualization_prompt[l] = '\0';
							break;
						}
					}
				
					if( k < slen && l < slen && l+1 < slen )
					{
					
						strcpy(typen,visualization_prompt);
						strcpy(name,visualization_prompt+k+1);
						comp = atoi(visualization_prompt+l+1);
						visualization_prompt[k] = ':';
						visualization_prompt[l] = ':';
						printf("type %s name %s comp %d\n",typen,name,comp);
						std::string stype(typen), sname(name);
						if( mesh->HaveTag(sname) )
						{
							Tag source_tag = mesh->GetTag(sname);
							if( source_tag.GetDataType() == DATA_REAL || source_tag.GetDataType() == DATA_INTEGER )
							{
								if( comp >= 0 && comp < source_tag.GetSize() )
								{
									visualization_type = NONE;
									for(size_t q = 0; q < stype.size(); ++q) 
									{
										stype[q] = tolower(stype[q]);
										typen[q] = tolower(typen[q]);
									}
									if( stype == "node" ) visualization_type = NODE;
									else if ( stype == "edge" ) visualization_type = EDGE;
									else if ( stype == "face" ) visualization_type = FACE;
									else if ( stype == "cell" ) visualization_type = CELL;

									if( visualization_type != NONE )
									{
										if( source_tag.isDefined(visualization_type) )
										{
											float min = 1.0e20, max = -1.0e20;
											printf("prepearing data for visualization\n");
											visualization_tag = mesh->CreateTag("VISUALIZATION_TAG",DATA_REAL,NODE,NONE,1);
											
											for(Mesh::iteratorNode it = mesh->BeginNode(); it != mesh->EndNode(); ++it)
											{
												ElementArray<Element> elems = it->getAdjElements(visualization_type);
												Storage::real_array coords = it->Coords();
												Storage::real cnt[3], dist, wgt;
												Storage::real val = 0.0, vol = 0.0, res;
												for(ElementArray<Element>::iterator jt = elems.begin(); jt != elems.end(); ++jt) if( jt->HaveData(source_tag) && jt->GetDataSize(source_tag) > comp )
												{
													jt->Centroid(cnt);
													dist = (cnt[0]-coords[0])*(cnt[0]-coords[0])+(cnt[1]-coords[1])*(cnt[1]-coords[1])+(cnt[2]-coords[2])*(cnt[2]-coords[2]);
													wgt = 1.0/(dist+1.0e-8);
													val += wgt * (source_tag.GetDataType() == DATA_REAL ? jt->RealArray(source_tag)[comp] : static_cast<Storage::real>(jt->IntegerArray(source_tag)[comp]) );
													vol += wgt;
												}
												res = val/vol;
												if( res < min ) min = res;
												if( res > max ) max = res;
												it->RealDF(visualization_tag) = res;
											}
											visualization_tag = mesh->CreateTag("VISUALIZATION_TAG",DATA_REAL,CELL,NONE,1);
											for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
											{
												ElementArray<Element> elems = it->getAdjElements(visualization_type);
												Storage::real coords[3];
												it->Centroid(coords);
												Storage::real cnt[3], dist, wgt;
												Storage::real val = 0.0, vol = 0.0, res;
												for(ElementArray<Element>::iterator jt = elems.begin(); jt != elems.end(); ++jt) if( jt->HaveData(source_tag) && jt->GetDataSize(source_tag) > comp )
												{
													jt->Centroid(cnt);
													dist = (cnt[0]-coords[0])*(cnt[0]-coords[0])+(cnt[1]-coords[1])*(cnt[1]-coords[1])+(cnt[2]-coords[2])*(cnt[2]-coords[2]);
													wgt = 1.0/(dist+1.0e-8);
													val += wgt * (source_tag.GetDataType() == DATA_REAL ? jt->RealArray(source_tag)[comp] : static_cast<Storage::real>(jt->IntegerArray(source_tag)[comp]) );
													vol += wgt;
												}
												res = val/vol;
												if( res < min ) min = res;
												if( res > max ) max = res;
												it->RealDF(visualization_tag) = res;
											}

											CommonColorBar->set_min(min);
											CommonColorBar->set_max(max);
											char comment[1024];
											sprintf(comment,"%s[%d] on %s, [%g:%g]",name,comp,typen,min,max);
											CommonColorBar->set_comment(comment);
											clipupdate = true;
											/*
											all_boundary.clear();
											INMOST_DATA_ENUM_TYPE pace = std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,(unsigned)boundary_faces.size()/100));
											for(INMOST_DATA_ENUM_TYPE k = 0; k < boundary_faces.size(); k++) 
											{
												all_boundary.push_back(DrawFace(boundary_faces[k]));
												if( k%pace == 0 ) all_boundary.back().set_flag(true);
											}
											*/
										}
										else printf("tag %s is not defined on element type %s\n",name, typen);
									}
									else printf("do not understand element type %s, should be: node, edge, face, cell\n",typen);
								}
								else printf("component is out of range for tag %s of size %u\n",name,source_tag.GetSize());
							}
							else printf("tag %s is not real or integer\n",name);
						}
						else printf("mesh do not have tag with name %s\n",name);
					}
					else printf("malformed string %s for visualization\n",visualization_prompt);
					visualization_prompt_active = 0;
					//visualization_prompt[0] = '\0';
				
				}

				glutPostRedisplay();
			}
			delete CommonInput;
			CommonInput = NULL;
		}
		glEnable(GL_DEPTH_TEST);
	}

	if( visualization_tag.isValid() )
	{
		glDisable(GL_DEPTH_TEST);
		glLoadIdentity();
		set_matrix2d();
		CommonColorBar->Draw();
		glEnable(GL_DEPTH_TEST);
	}
}

void draw()
{
  draw_screen();
	glutSwapBuffers();
}



int main(int argc, char ** argv)
{
	double tt;
	table[CENTROID]    = CELL | FACE;
	//table[NORMAL]      = FACE;
	//table[ORIENTATION] = FACE;
	//table[MEASURE]     = CELL|FACE;
	mesh = new Mesh();
	printf("Started loading mesh.\n");
	tt = Timer();
  //mesh->SetCommunicator(INMOST_MPI_COMM_WORLD);
  mesh->SetParallelFileStrategy(0);
  mesh->SetParallelStrategy(1);
	mesh->SetFileOption("VERBOSITY","2");
	if( argc < 2 )
	{
		printf("Usage: %s mesh_file [dims]\n",argv[0]);
		return 0;
	}
	//try
	{
		if( argc > 2 )	mesh->SetFileOption("VTK_GRID_DIMS",argv[2]);
		mesh->Load(argv[1]);
	} 
	/*
	catch(...)
	{
		printf("Error during loading %s\n",argv[1]);
		return -1;
	}
	*/
	printf("Done. Time %lg\n",Timer()-tt);
	//int fixed = 0;
	//tt = Timer();
	//delete mesh;
	//printf("Delete %lg\n",Timer()-tt);
	//return 0;

	std::map<Element::GeometricType,int> elems;

	for(Mesh::iteratorElement it = mesh->BeginElement(CELL|FACE|EDGE|NODE); it != mesh->EndElement(); ++it)
		elems[it->GetGeometricType()]++;

	for(std::map<Element::GeometricType,int>::iterator it = elems.begin(); it != elems.end(); ++it )
		std::cout << Element::GeometricTypeName(it->first) << ": " << it->second << std::endl;


	for(ElementType mask = NODE; mask <= CELL; mask = mask << 1)	
		std::cout << ElementTypeName(mask) << " " << mesh->NumberOf(mask) << std::endl;
	
	//return 0;
	//for(Mesh::iteratorFace f = mesh->BeginFace(); f != mesh->EndFace(); f++)
	//	if( Geometry::FixNormaleOrientation(&*f) ) fixed++;
	//printf("fixed: %d\n",fixed);
	printf("Computing geometric quantities\n");
	tt = Timer();
	mesh->PrepareGeometricData(table);
	printf("Done. %lg\n",Timer()-tt);
	
	printf("Calculating bounds.\n");
	tt = Timer();
	for(Mesh::iteratorNode n = mesh->BeginNode(); n != mesh->EndNode(); n++)
	{
		Storage::real_array c = n->Coords();
		if( c[0] > sright ) sright = c[0];
		if( c[0] < sleft ) sleft = c[0];
		if( c[1] > stop ) stop = c[1];
		if( c[1] < sbottom ) sbottom = c[1];
		if( c[2] > sfar ) sfar = c[2];
		if( c[2] < snear ) snear = c[2];
	}
	printf("Done. Time %lg\n",Timer()-tt);
	printf("%g:%g %g:%g %g:%g\n",sleft,sright,sbottom,stop,snear,sfar);

	printf("Gathering boundary faces.\n");
	tt = Timer();
	boundary_faces.SetMeshLink(mesh);
	for(Mesh::iteratorCell f = mesh->BeginCell(); f != mesh->EndCell(); f++) 
	{
		if( f->GetElementDimension() == 2 )
			boundary_faces.push_back(*f);
	}
	for(Mesh::iteratorFace f = mesh->BeginFace(); f != mesh->EndFace(); f++) if( f->GetElementDimension() == 2 )
	{
		if( f->Boundary() )
			boundary_faces.push_back(*f);
	}
	printf("Done. Time %lg\n",Timer()-tt);

	if( boundary_faces.empty() )
	{
		printf("Haven't found any boundary elements of the mesh. Nothing to display.\n");
		return -1;
	}

	printf("Prepearing set of boundary faces for drawing.\n");
	tt = Timer();
	INMOST_DATA_ENUM_TYPE pace = std::max<INMOST_DATA_ENUM_TYPE>(1,std::min<INMOST_DATA_ENUM_TYPE>(15,(unsigned)boundary_faces.size()/100));
	for(INMOST_DATA_ENUM_TYPE k = 0; k < boundary_faces.size(); k++) 
	{
		all_boundary.push_back(DrawFace(boundary_faces[k]));
		if( k%pace == 0 ) all_boundary.back().set_flag(true);
	}
	printf("Done. Time %g\n",Timer() - tt);

  if(mesh->HaveTag("ADDED_ELEMENTS") )
  {
    Tag add = mesh->GetTag("ADDED_ELEMENTS");
    for(Mesh::iteratorEdge it = mesh->BeginEdge(); it != mesh->EndEdge(); ++it) if( it->Integer(add) ) added_edges.push_back(it->self());
    for(Mesh::iteratorFace it = mesh->BeginFace(); it != mesh->EndFace(); ++it) if( it->Integer(add) && !it->Boundary() ) added_faces.push_back(DrawFace(it->self()));
  }
  if(mesh->HaveTag("CONORMALS"))
  {
    Tag cnrmls = mesh->GetTag("CONORMALS");
    conormals.reserve(mesh->NumberOfFaces()*2*3);
    Storage::real side = 1.0e20;
    for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
    {
      Storage::real max[3],min[3];
      GetBox(it->self(),min,max);
      side = std::min(side,max[0]-min[0]);
      side = std::min(side,max[1]-min[1]);
      side = std::min(side,max[2]-min[2]);
    }
    side *=0.25;
    for(Mesh::iteratorFace it = mesh->BeginFace(); it != mesh->EndFace(); ++it)
    {
      Storage::real_array fconormals = it->RealArray(cnrmls);
      Storage::real cnt[3];
      it->Centroid(cnt);
      conormals.push_back(cnt[0]);
      conormals.push_back(cnt[1]);
      conormals.push_back(cnt[2]);
      conormals.push_back(cnt[0]+side*fconormals[0]);
      conormals.push_back(cnt[1]+side*fconormals[1]);
      conormals.push_back(cnt[2]+side*fconormals[2]);
      conormals.push_back(cnt[0]);
      conormals.push_back(cnt[1]);
      conormals.push_back(cnt[2]);
      conormals.push_back(cnt[0]+side*fconormals[3]);
      conormals.push_back(cnt[1]+side*fconormals[4]);
      conormals.push_back(cnt[2]+side*fconormals[5]);
    }
  }
  if(mesh->HaveTag("HARMONIC_POINT"))
  {
    Tag h1 = mesh->GetTag("HARMONIC_POINT");
    harmonic_points.reserve(mesh->NumberOfFaces()*3);
    for(Mesh::iteratorFace it = mesh->BeginFace(); it != mesh->EndFace(); ++it)
    {
      Storage::real_array h1a = it->RealArray(h1);
      harmonic_points.insert(harmonic_points.end(),h1a.begin(),h1a.end());
    }
  }
  std::cout << "Harmonic points array size: " << harmonic_points.size()/3 << std::endl;

  if(mesh->HaveTag("DUAL_HARMONIC_POINT"))
  {
    Tag h2 = mesh->GetTag("DUAL_HARMONIC_POINT");
    for(Mesh::iteratorElement it = mesh->BeginElement(CELL|FACE); it != mesh->EndElement(); ++it) if( it->HaveData(h2) )
    {
      Storage::real_array h2a = it->RealArray(h2);
      dual_harmonic_points.insert(dual_harmonic_points.end(),h2a.begin(),h2a.end());
    }
  }

  std::cout << "Boundary harmonic points array size: " << dual_harmonic_points.size()/9 << std::endl;

	printf("Prepearing interactive mesh clipper.\n");
	//tt = Timer();
	oclipper = new clipper(mesh);
	bclipper = new bnd_clipper(mesh,boundary_faces.data(),(unsigned)boundary_faces.size());
	clipupdate = true;
	//printf("Done. Time %g\n",Timer() - tt);

	CommonVolumetricView = new volumetric(mesh);
	//CommonVolumetricView = new volumetric2(mesh);
	
	shift[0] = -(sleft+sright)*0.5;
	shift[1] = -(sbottom+stop)*0.5;
	shift[2] =  -(sfar+snear)*0.5;


	p[0] = (sleft+sright)*0.5;
	p[1] = (sbottom+stop)*0.5;
	p[2] = (sfar+snear)*0.5;

  if( mesh->HaveTag("FACE_CENTER") && mesh->GetTag("FACE_CENTER").isDefined(FACE) && mesh->GetTag("FACE_CENTER").GetSize() == 3 )
  {
    std::cout << "Have face centers!" << std::endl;
    face_center = mesh->GetTag("FACE_CENTER");
  }


  if( false )
  //if( mesh->HaveTag("VELOCITY") && mesh->GetTag("VELOCITY").isDefined(CELL) )
  {
    printf("preparing octree around mesh, was sets %d\n",mesh->NumberOfSets());
    Octree octsearch = Octree(mesh->CreateSet("octsearch").first);
    octsearch.Construct(NODE,true);
    printf("done, sets %d\n",mesh->NumberOfSets());
    printf("building streamlines\n");
    Tag cell_size = mesh->CreateTag("CELL_SIZES",DATA_REAL,CELL,NONE,1);
    Tag vel = mesh->GetTag("VELOCITY");
    Storage::real velmax = 0, velmin = 1.0e20, l;
    for(Mesh::iteratorCell c = mesh->BeginCell(); c != mesh->EndCell(); ++c)
    {
      coord velv(c->RealArray(vel).data());
      l = velv.length();
      if( l > velmax ) velmax = l;
      if( l < velmin ) velmin = l;
      c->RealDF(cell_size) = GetSize(c->self());
    }
    velmax = log(velmax+1.0e-25);
    velmin = log(velmin+1.0e-25);
    Tag flux = mesh->GetTag("FACE_FLUX");
    /*
    MarkerType visited = mesh->CreateMarker();
    if( false )
    for(Mesh::iteratorFace f = mesh->BeginFace(); f != mesh->EndFace(); ++f) if( f->Boundary() )
    {
      if( f->Real(flux) > 1.0e-5 )
      {
        coord cntf;
        f->Centroid(cntf.data());
        streamlines.push_back(Streamline(octsearch,cntf,vel,cell_size,velmin,velmax,1.0,visited));

        ElementArray<Edge> edges = f->getEdges();
        for(ElementArray<Edge>::iterator n = edges.begin(); n != edges.end(); ++n)
        {
          coord cntn;
          n->Centroid(cntn.data());
          if( cntn[2] == cntf[2] )
          {
            const Storage::real coef[4] = {0.4,0.8};
            for(int q = 0; q < 2; ++q)
              streamlines.push_back(Streamline(octsearch,cntf*coef[q]+cntn*(1-coef[q]),vel,cell_size,velmin,velmax,1.0,visited));
          }
        }
      }
      else if( f->Real(flux) < -1.0e-5 )
      {
        coord cntf;
        f->Centroid(cntf.data());
        streamlines.push_back(Streamline(octsearch,cntf,vel,cell_size,velmin,velmax,-1.0,visited));

        ElementArray<Edge> edges = f->getEdges();
        for(ElementArray<Edge>::iterator n = edges.begin(); n != edges.end(); ++n)
        {
          coord cntn;
          n->Centroid(cntn.data());
          if( cntn[2] == cntf[2] )
          {
            const Storage::real coef[4] = {0.4,0.8};
            for(int q = 0; q < 2; ++q)
              streamlines.push_back(Streamline(octsearch,cntf*coef[q]+cntn*(1-coef[q]),vel,cell_size,velmin,velmax,-1.0,visited));
          }
        }
      }
    }
    
    printf("done from boundary faces, total streamlines = %d\n",streamlines.size());
    
    for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
    {
      if( !it->GetMarker(visited) )
      {
        coord cntc;
        it->Centroid(cntc.data());
        if( coord(it->RealArray(vel).data()).length() > 1.0e-4 )
        {
          streamlines.push_back(Streamline(octsearch,cntc,vel,cell_size,velmin,velmax,1.0,0));
          streamlines.push_back(Streamline(octsearch,cntc,vel,cell_size,velmin,velmax,-1.0,0));
        }
      }
    }
    printf("done from unvisited cells, total streamlines = %d\n",streamlines.size());
    for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
      it->RemMarker(visited);
    mesh->ReleaseMarker(visited);
    */
    int numlines = 16;
    for(int i = 0; i < numlines; ++i)
    {
      for(int j = 0; j < numlines; ++j)
      {
        coord xyz;
        double div = static_cast<double>(numlines);
        xyz[0] = i/div + 0.5/div;
        xyz[1] = j/div + 0.5/div;
        xyz[2] = 0.5;
        streamlines.push_back(Streamline(octsearch,xyz,vel,cell_size,velmin,velmax,1.0,0));
        streamlines.push_back(Streamline(octsearch,xyz,vel,cell_size,velmin,velmax,-1.0,0));
      }
    }
    printf("done from %d by %d grid, total streamlines = %d\n",numlines,numlines,(int)streamlines.size());
    mesh->DeleteTag(cell_size);
    printf("done, total streamlines = %d\n",streamlines.size());
    printf("killing octree, was sets %d\n",mesh->NumberOfSets());
    octsearch.Destroy();
    printf("done, sets %d\n",mesh->NumberOfSets());
  }

	quatinit();
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutInitWindowPosition (100, 100);
	glutCreateWindow("Graph");
	
	glDepthFunc(GL_LEQUAL);
	glClearDepth(1.f);
	glEnable(GL_DEPTH_TEST);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);


	CommonColorBar = new color_bar;
	
	//glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);
	//glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	
	glClearColor (1.0f, 1.0f, 1.0f, 1.f);
	glutDisplayFunc(draw);
	glutReshapeFunc(reshape);
	
	glutKeyboardFunc(keyboard);
	glutKeyboardUpFunc(keyboard2);
	glutMouseFunc(myclick);
	glutMotionFunc(myclickmotion);
	glutPassiveMotionFunc(mymotion);
	
	//glutIdleFunc(idle);
	
	glutPostRedisplay();
	glutMainLoop();
}
