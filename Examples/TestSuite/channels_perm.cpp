#include "inmost.h"

using namespace INMOST;

std::string problem_name = "channels";

int main(int argc, char ** argv)
{
    if( argc < 3 )
    {
        std::cout << "Usage: " << argv[0] << " mesh mesh_out [Kx:1000 Ky:10 Kz:10 Reduction:0.001]" << std::endl;
        return 0;
    }
    
    Mesh * m = new Mesh;
    m->SetFileOption("VERBOSITY","2");
    try{m->Load(argv[1]);} catch(...) { std::cout << "Cannot load the mesh " << argv[1] << std::endl; return -1;}
    
    
    double max[3] = {-1.0e20, -1.0e20, -1.0e20}, min[3] = {1.0e20,1.0e20,1.0e20};
    Storage::real c[3] = {0,0,0};
    for(Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); ++it)
    {
        it->Centroid(c);
        if( c[0] > max[0] ) max[0] = c[0];
        if( c[1] > max[1] ) max[1] = c[1];
        if( c[2] > max[2] ) max[2] = c[2];
        if( c[0] < min[0] ) min[0] = c[0];
        if( c[1] < min[1] ) min[1] = c[1];
        if( c[2] < min[2] ) min[2] = c[2];
    }
    
    if( max[0] <= min[0] ) {std::cout << "strange X " << min[0] << ":" << max[0] << std::endl; return -1;}
    if( max[1] <= min[1] ) {std::cout << "strange Y " << min[1] << ":" << max[1] << std::endl; return -1;}
    if( max[2] <= min[2] )
    {
        //2d mesh
        if( m->GetDimensions() == 3 )
        {
            //offset from z-plane
            min[2] -= 0.0001;
            max[2] += 0.0001;
        }
        else
        {
            min[2] = -0.0001;
            max[2] = +0.0001;
        }
    }
    
    Tag material;
    if( m->HaveTag("MATERIAL") ) material = m->GetTag("MATERIAL");
    if( m->HaveTag("PERM") ) m->DeleteTag(m->GetTag("PERM"));
    
    Tag tensor = m->CreateTag("PERM",DATA_REAL,CELL,NONE,9);
    
    
    {
        Storage::bulk_array name = m->self()->BulkArray(m->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
        name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
    }
    
    Storage::real Kx0 = 1000.0, Ky0 = 100.0, Kz0 = 100.0;//, Reduction = 0.001;
    
    if( argc > 3 ) Kx0 = atof(argv[3]);
    if( argc > 4 ) Ky0 = atof(argv[4]);
    if( argc > 5 ) Kz0 = atof(argv[5]);
    //~ if( argc > 6 ) Reduction = atof(argv[6]);
    
    double params[5][6] =  //shift, channel width, amplitude1, phase1, amplitude2, phase2, f(x) = x + shift + amplitude1*sin(phase1*x) + amplitude2*sin(phase2*x)
    {
        {0.5  , 0.05, 0.01 ,20,0.05 ,30},
        {0.25 , 0.05, 0.05 ,15,0.01 ,45},
        {0.0  , 0.05, 0.05 ,35,0.01 ,60},
        {-0.25, 0.05, 0.025,20,0.025,30},
        {-0.48 , 0.05, 0.075,10,0.055,30}
    };
    
    std::cout << "Setting channels-guided permeability" << std::endl;
    
    for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
    {
        it->Centroid(c);
        // normalized coordinates in [0,1]
        Storage::real alpha = (c[0]-min[0])/(max[0]-min[0]);
        Storage::real beta  = (c[1]-min[1])/(max[1]-min[1]);
        //Storage::real gamma = (c[2]-min[2])/(max[2]-min[2]);
        
        
        alpha = floor(alpha*64.0)/64.0 + 1.0/128.0;
        beta = floor(beta*64.0)/64.0 + 1.0/128.0;
        
        //initiate channel
        int block_type = 0;
        
        
        for(int q = 0; q < 5; ++q) //go over 5 channels
        {
            double fx = alpha + params[q][0] + params[q][2]*sin(params[q][3]*alpha) + params[q][4]*sin(params[q][5]*alpha);
            if( fabs(fx-beta) < params[q][1] )
            {
                block_type = q+1;
                break;
            }
        }
        
        
        
        Storage::real Kx = Kx0, Ky = Ky0, Kz = Kz0;
        
        
        
        Storage::real_array perm = it->RealArrayDF(tensor);
        std::fill(perm.begin(),perm.end(),0.0);
        if( (beta < 0.5 - alpha) || (beta > 1.5 - alpha) )
        {
            perm[0] = 100;
            perm[1] = 0.0;
            perm[2] = 0.0;
            perm[3] = 0.0;
            perm[4] = 100;
            perm[5] = 0.0;
            perm[6] = 0.0;
            perm[7] = 0.0;
            perm[8] = 100;
        }
        else if( !block_type  )
        {
         /*
            Kx *= Reduction;
            Ky *= Reduction;
            Kz *= Reduction;
            //45 degree rotation
            perm[0] = Kx;//(Kx + Ky)*0.5;
            perm[1] = perm[3] = 0.0;//(Ky - Kx)*0.5;
            perm[2] = perm[6] = 0.0;
            perm[4] = Ky;//perm[0];
            perm[5] = 0.0;
            perm[6] = 0.0;
            perm[7] = 0.0;
            perm[8] = Kz;
           */
			
             perm[0] = 1;//(Kx + Ky)*0.5;
             perm[1] = perm[3] = 0.0;//(Ky - Kx)*0.5;
             perm[2] = perm[6] = 0.0;
             perm[4] = 1;//0.01;//perm[0];
             perm[5] = 0.0;
             perm[6] = 0.0;
             perm[7] = 0.0;
             perm[8] = Kz;
			
        }
        else
        {
            int q = block_type - 1;
            Storage::real dfx = 1.0 + params[q][2]*params[q][3]*cos(params[q][3]*alpha) + params[q][4]*params[q][5]*cos(params[q][5]*alpha);
            Storage::real angle = atan2(dfx,1.0);
            Storage::real acos = cos(angle);
            Storage::real asin = sin(angle);
            
            perm[0] = Kx * acos * acos + Ky * asin * asin;
            perm[1] = (Kx - Ky) * acos * asin;
            perm[2] = 0.0;
            perm[3] = perm[1];
            perm[4] = Ky * acos * acos + Kx * asin * asin;
            perm[5] = 0.0;
            perm[6] = 0.0;
            perm[7] = 0.0;
            perm[8] = Kz;	
        }
        
        
        if( false )//material.isValid() )
        {
            //isotropic tensor
            //handle corners with wells for meshes generated by Bradley
            if( it->Integer(material) == 1 || it->Integer(material) == 2 || it->Integer(material) < 0 )
            {
                perm[0] = 100;
                perm[1] = 0;
                perm[2] = 0;
                perm[3] = 0;
                perm[4] = 100;
                perm[5] = 0;
                perm[6] = 0;
                perm[7] = 0;
                perm[8] = 100;
            }
            else if( it->Integer(material) > 2 ) //set impermiable inclusions
                std::fill(perm.begin(),perm.end(),0.0);
        }
        
    }
    
    std::cout << "Saving output to " << argv[2] << std::endl;
    
    m->Save(argv[2]);
    
    delete m;
}


