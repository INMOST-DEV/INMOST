//g++ main.cpp -lpetsc -L/usr/X11R6/lib -lX11 -ldmumps -lmumps_common -lmpi_f77 -lscalapack -lpord -lblacs -lparmetis -lmetis -lmpi -lHYPRE -lgfortran -lblas -llapack ../../INMOST.a
//mpicxx main.cpp -lpetsc -L/usr/X11R6/lib -lX11 -ldmumps -lmumps_common -lmpi_f77 -lscalapack -lpord -lblacs -lparmetis -lmetis -lmpi -lHYPRE -lgfortran -lblas -llapack ../../INMOST.a -lzoltan -lparmetis -lmetis
// run ./a.out grids/rezultMesh.vtk MATERIALS -ksp_monitor -ksp_view -mat_type aij -vec_type standard

#include "inmost.h"
#include <stdio.h>
#ifndef M_PI
#define M_PI 3.141592653589
#endif
using namespace INMOST;

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

void make_vec(Storage::real v1[3], Storage::real v2[3], Storage::real out[3])
{
    out[0] = v1[0] - v2[0];
    out[1] = v1[1] - v2[1];
    out[2] = v1[2] - v2[2];
}

Storage::real dot_prod(Storage::real v1[3], Storage::real v2[3])
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

void tensor_prod(Storage::real K[9], Storage::real v[3], Storage::real out[3])
{
    out[0] = v[0]*K[0] + v[1]*K[1] + v[2]*K[2];
    out[1] = v[0]*K[3] + v[1]*K[4] + v[2]*K[5];
    out[2] = v[0]*K[6] + v[1]*K[7] + v[2]*K[8];
}

Storage::real transmissibility(Storage::real vec[3], Storage::real K_tensor[9], Storage::real normal_face[3])
{
    Storage::real Kn[3];
    tensor_prod(K_tensor,normal_face,Kn);
    return dot_prod(vec,Kn);
}


Storage::real tensor_K_mat_const(Storage::integer mat, int comp)
{
    if( comp == 0 ) //X
    {
        switch(mat)
        {
            case 0: return 850;
            case 1: return 950;
            case 2: return 1050;
            case 3: return 450;
            case 4: return 650;
            case 5: return 550;
            default: return 1000;
        }
    }
    else if( comp == 1 ) //Y
    {
        switch(mat)
        {
            case 0: return 750;
            case 1: return 450;
            case 2: return 950;
            case 3: return 350;
            case 4: return 650;
            case 5: return 750;
            default: return 1000;
        }
    }
    else if( comp == 2 ) //Z
    {
        switch(mat)
        {
            case 0: return 12;
            case 1: return 11;
            case 2: return 13;
            case 3: return 15;
            case 4: return 14;
            case 5: return 13;
            default: return 10;
        }
    }
    else return 0;
}

Storage::real tensor_K_mat_var(Storage::integer mat, int comp)
{
    if( comp == 0 ) //X
    {
        switch(mat)
        {
            case 0: return 150;
            case 1: return 50;
            case 2: return 80;
            case 3: return 90;
            case 4: return 100;
            case 5: return 20;
            default: return 10;
        }
    }
    else if( comp == 1 ) //Y
    {
        switch(mat)
        {
            case 0: return 150;
            case 1: return 50;
            case 2: return 80;
            case 3: return 90;
            case 4: return 100;
            case 5: return 20;
            default: return 10;
        }
    }
    else if( comp == 2 ) //Z
    {
        switch(mat)
        {
            case 0: return 2;
            case 1: return 5;
            case 2: return 3;
            case 3: return 2;
            case 4: return 1;
            case 5: return 1;
            default: return 2;
        }
    }
    else return 0;
}

Storage::real tensor_K_mat_rotOxy(Storage::integer mat)
{
    switch(mat)
    {
        case 0: return 6;
        case 1: return 3;
        case 2: return 5;
        case 3: return 4;
        case 4: return 5;
        case 5: return 7;
        default: return 4;
    }
}

int main(int argc,char ** argv)
{
    Mesh::Initialize(&argc,&argv);
    Solver::Initialize(&argc,&argv,"database.xml");
#if defined(USE_PARTITIONER)
    Partitioner::Initialize(&argc,&argv);
#endif
    /*
    {
        Mesh m;
        //m.SetFileOption("VERBOSITY","2");
        m.Load(argv[1]);
        m.Save("dump.pvtk");
    }
    Solver::Finalize();
    Partitioner::Finalize();
    Mesh::Finalize();
    return 0;
    */

    if( argc > 1 )
    {
        Tag phi, tensor_K, id, mat;
        Mesh * m = new Mesh();
        double ttt = Timer();
        bool repartition = false;
        m->SetCommunicator(INMOST_MPI_COMM_WORLD);
        //~ m->RemTopologyCheck(THROW_EXCEPTION);
        //~ m->SetTopologyCheck(PRINT_NOTIFY | MARK_ON_ERROR | DUPLICATE_EDGE | DUPLICATE_FACE | DEGENERATE_EDGE | DEGENERATE_FACE | DEGENERATE_CELL | FACE_PLANARITY | INTERLIEVED_FACES | TRIPLE_SHARED_FACE | FLATTENED_CELL | ADJACENT_DUPLICATE | ADJACENT_DIMENSION | PROHIBIT_MULTILINE | PROHIBIT_POLYGON | PROHIBIT_MULTIPOLYGON | PROHIBIT_POLYHEDRON | NEED_TEST_CLOSURE | DISABLE_2D);
        if( m->GetProcessorRank() == 0 ) std::cout << argv[0] << std::endl;


        if( m->isParallelFileFormat(argv[1]) )
        {
            m->Load(argv[1]);
            repartition = true;

        }
        else
        {
            if( m->GetProcessorRank() == 0 )
                m->Load(argv[1]);
        }
        BARRIER

        if( m->GetProcessorRank() == 0 ) std::cout << "Processors: " << m->GetProcessorsNumber() << std::endl;
        if( m->GetProcessorRank() == 0 ) std::cout << "Load(MPI_File): " << Timer()-ttt << std::endl;


        //~ {
        //~ double ttt2 = Timer();
        //~ Mesh t;
        //~ t.SetCommunicator(INMOST_MPI_COMM_WORLD);
        //~ t.SetParallelFileStrategy(0);
        //~ t.Load(argv[1]);
        //~ BARRIER
        //~ if( m->GetProcessorRank() == 0 ) std::cout << "Load(MPI_Scatter): " << Timer()-ttt2 << std::endl;
        //~ }

#if defined(USE_PARTITIONER)
        ttt = Timer();
        Partitioner * p = new Partitioner(m);
        p->SetMethod(Partitioner::Inner_RCM,Partitioner::Partition);
        p->Evaluate();
        delete p;
        BARRIER

        if( m->GetProcessorRank() == 0 ) std::cout << "Evaluate: " << Timer()-ttt << std::endl;


        ttt = Timer();
        m->Redistribute();
        m->ReorderEmpty(CELL|FACE|EDGE|NODE);
        BARRIER

        if( m->GetProcessorRank() == 0 ) std::cout << "Redistribute: " << Timer()-ttt << std::endl;
#endif

        if( !repartition && m->HaveTag("CM_INDEX") )
        {
            id = m->GetTag("CM_INDEX");
            m->ExchangeData(id,CELL,0);
        }
        else
        {
            ttt = Timer();
            m->AssignGlobalID(CELL | EDGE | FACE | NODE);
            BARRIER
            if( m->GetProcessorRank() == 0 ) std::cout << "Assign id: " << Timer()-ttt << std::endl;
            id = m->GlobalIDTag();
        }



        Tag test = m->CreateTag("test",DATA_INTEGER,CELL,CELL,1);

        phi = m->CreateTag("PHI",DATA_REAL,CELL,NONE,1);
        tensor_K = m->CreateTag("K",DATA_REAL,CELL,NONE,9);

        if( m->HaveTag("MATERIAL") ) mat = m->GetTag("MATERIAL");

        for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); ++it)
        {
            Storage::real_array K = it->RealArray(tensor_K);
            int material = mat.isValid() ? it->Integer(mat) : 1000;
            double rot = rand() / RAND_MAX * (tensor_K_mat_rotOxy(material)/180.0) * M_PI;
            double sinrot = sin(rot), cosrot = cos(rot);
            double sincosrot = sinrot*cosrot, sinrot2 = sinrot*sinrot, cosrot2 = cosrot*cosrot;
            double k1 = tensor_K_mat_const(material,0) + tensor_K_mat_var(material,0) * rand() / RAND_MAX;
            double k2 = tensor_K_mat_const(material,1) + tensor_K_mat_var(material,1) * rand() / RAND_MAX;
            double k3 = tensor_K_mat_const(material,2) + tensor_K_mat_var(material,2) * rand() / RAND_MAX;
            K[0] = K[4] = k1*cosrot2 + k2*sinrot2;
            K[1] = K[3] = (k2-k1)*sincosrot;
            K[8] = k3;

            if( rand() / RAND_MAX < 0.3 )
                it->Integer(test) = 1;
            //~ K[0] = K[4] = 100;
            //~ K[8] = 10;
        }

        ttt = Timer();
        m->ExchangeGhost(1,FACE);
        BARRIER
        if( m->GetProcessorRank() == 0 ) std::cout << "Exchange ghost: " << Timer()-ttt << std::endl;


        ttt = Timer();
        Solver S("inner_ilu2");
        Sparse::Matrix A;
        Sparse::Vector x,b;

        Mesh::GeomParam table;

        table[MEASURE] = CELL | FACE;
        table[CENTROID] = CELL | FACE;
        table[BARYCENTER] = CELL | FACE;
        table[NORMAL] = FACE;
        table[ORIENTATION] = FACE;
        m->PrepareGeometricData(table);
        //~ BARRIER
        //~ if( m->GetProcessorRank() == 0 ) std::cout << "Prepare geometric data: " << Timer()-ttt << std::endl;

        unsigned idmax = 0, idmin = UINT_MAX;
        for(Mesh::iteratorCell cell = m->BeginCell(); cell != m->EndCell(); ++cell)
            if( cell->GetStatus() != Element::Ghost )
            {
                unsigned pid = cell->Integer(id);
                if( pid < idmin ) idmin = pid;
                if( pid+1 > idmax ) idmax = pid+1;
            }

        A.SetInterval(idmin,idmax);
        x.SetInterval(idmin,idmax);
        b.SetInterval(idmin,idmax);
        //std::cout << m->GetProcessorRank() << " A,x,b interval " << idmin << ":" << idmax << " size " << idmax-idmin << std::endl;

        //solve \nabla \cdot \nabla phi = 0 equation
        for(Mesh::iteratorFace face = m->BeginFace(); face != m->EndFace(); ++face)
        {
            //~ std::cout << face->LocalID() << " / " << m->NumberOfFaces() << std::endl;
            Element::Status s1,s2;
            Cell r1 = face->BackCell();
            Cell r2 = face->FrontCell();
            if( ((!r1.isValid() || (s1 = r1->GetStatus()) == Element::Ghost)?0:1)+
                ((!r2.isValid() || (s2 = r2->GetStatus()) == Element::Ghost)?0:1) == 0) continue;
            Storage::real f_nrm[3], r1_cnt[3], r2_cnt[3], f_cnt[3], d1[3], d2[3], T1, T2, Coef;
            Storage::real f_area = face->Area();
            Storage::real vol1 = r1->Volume(), vol2;
            Storage::integer id1 = r1->Integer(id), id2;
            Storage::real_array K1 = r1->RealArray(tensor_K), K2;
            face->Normal(f_nrm);
            f_nrm[0] /= f_area;
            f_nrm[1] /= f_area;
            f_nrm[2] /= f_area;
            r1->Barycenter(r1_cnt);
            face->Barycenter(f_cnt);
            if( !r2.isValid() ) //boundary condition
            {
                if( r1->Integer(mat) == 2 ) //direchlet boundary for second material
                {
                    make_vec(f_cnt,r1_cnt,d1);
                    Coef = transmissibility(d1,&K1[0],f_nrm)/dot_prod(d1,d1);
                    //~ std::cout << Coef << " area " << face->Area() << " |norm| " << sqrt(dot_prod(f_nrm,f_nrm)) << " |d| " << sqrt(dot_prod(d1,d1)) << " vol " << vol1 << std::endl;
                    A[id1][id1] += -Coef / vol1 * f_area;
                    b[id1] += -Coef / vol1 * r1->Integer(mat);
                }
                //Neuman boundary on other
            }
            else
            {
                vol2 = r2->Volume();
                K2 = r2->RealArray(tensor_K);
                id2 = r2->Integer(id);
                r2->Barycenter(r2_cnt);
                make_vec(f_cnt,r1_cnt,d1);
                make_vec(r2_cnt,f_cnt,d2);
                T1 = transmissibility(d1,&K1[0],f_nrm)/dot_prod(d1,d1);
                T2 = transmissibility(d2,&K2[0],f_nrm)/dot_prod(d2,d2);
                /*
                if( T1 < 0 )
                {
                    Storage::real Kn[3];
                    std::cout << "T1 is negative " << T1 << " face id " << face->LocalID() << " orientation fix " << face->FixNormalOrientation() << " fix2 " << face->FixNormalOrientation() << std::endl;
                    {

                        {
                            Storage::real nf[3], cf[3],fc[3], cc[3], dot;
                            face->Normal(nf);
                            face->Centroid(fc);
                            r1->Centroid(cc);
                            make_vec(cc,fc,cf);
                            std::cout << "cc " << cc[0] << " " << cc[1] << " " << cc[2] << std::endl;
                            std::cout << "fc " << fc[0] << " " << fc[1] << " " << fc[2] << std::endl;
                            std::cout << "cf " << cf[0] << " " << cf[1] << " " << cf[2] << std::endl;
                            std::cout << "nf " << nf[0] << " " << nf[1] << " " << nf[2] << std::endl;
                            dot = dot_prod(nf,cf);
                            std::cout << "dot " << dot << std::endl;
                        }
                    }
                    adjacent<Node> nodes = face->getNodes();
                    for(unsigned i = 0; i < nodes.size(); ++i) std::cout << "node " << i << " " << nodes[i].Coords()[0] << " " << nodes[i].Coords()[1] << " " << nodes[i].Coords()[2] << " id " << nodes[i].LocalID() << std::endl;
                    std::cout << "  f_nrm " << f_nrm[0] << " " << f_nrm[1] << " " << f_nrm[2] << std::endl;

                    std::cout << "  f_cnt " << f_cnt[0] << " " << f_cnt[1] << " " << f_cnt[2] << std::endl;
                    std::cout << " r1_cnt " << r1_cnt[0] << " " << r1_cnt[1] << " " << r1_cnt[2] << std::endl;
                    std::cout << "   d1   " << d1[0] << " " << d1[1] << " " << d1[2] << std::endl;
                    std::cout << " tensor " << K1[0] << " " << K1[1] << " " << K1[2] << std::endl;
                    std::cout << "        " << K1[3] << " " << K1[4] << " " << K1[5] << std::endl;
                    std::cout << "        " << K1[6] << " " << K1[7] << " " << K1[8] << std::endl;
                    std::cout << " d1*nrm " << dot_prod(d1,f_nrm) << std::endl;
                    tensor_prod(&K1[0],f_nrm,Kn);
                    std::cout << " d1*Kn  " << dot_prod(d1,Kn) << std::endl;
                    std::cout << " T2     " << T2 << std::endl;
                    std::cout << " Coef " << 1.0 / ( 1.0 / T1 + 1.0 / T2 ) << std::endl;
                }
                if( T2 < 0 )
                {
                    Storage::real Kn[3];
                    std::cout << "T2 is negative " << T2 << " face id " << face->LocalID() << " orientation fix " << face->FixNormalOrientation() << " fix2 " << face->FixNormalOrientation() << std::endl;
                    std::cout << "  f_nrm " << f_nrm[0] << " " << f_nrm[1] << " " << f_nrm[2] << std::endl;
                    std::cout << "  f_cnt " << f_cnt[0] << " " << f_cnt[1] << " " << f_cnt[2] << std::endl;
                    std::cout << " r2_cnt " << r2_cnt[0] << " " << r2_cnt[1] << " " << r2_cnt[2] << std::endl;
                    std::cout << "   d2   " << d2[0] << " " << d2[1] << " " << d2[2] << std::endl;
                    std::cout << " tensor " << K2[0] << " " << K2[1] << " " << K2[2] << std::endl;
                    std::cout << "        " << K2[3] << " " << K2[4] << " " << K2[5] << std::endl;
                    std::cout << "        " << K2[6] << " " << K2[7] << " " << K2[8] << std::endl;
                    std::cout << " d2*nrm " << dot_prod(d2,f_nrm) << std::endl;
                    tensor_prod(&K2[0],f_nrm,Kn);
                    std::cout << " d2*Kn  " << dot_prod(d2,Kn) << std::endl;
                    std::cout << " T1     " << T1 << std::endl;
                    std::cout << " Coef   " << 1.0 / ( 1.0 / T1 + 1.0 / T2 ) << std::endl;
                }*/
                Coef = 1.0 / ( 1.0 / T1 + 1.0 / T2 ) * f_area; // Harmonic mean
                if( Coef < 0 )
                {
                    if( T1 > 0 ) Coef = T1;
                    else if( T2 > 0 ) Coef = T2;
                    else Coef = 0;//std::cout << "no positive transmissibility" << std::endl;
                }
                if( s1 != Element::Ghost )
                {
                    A[id1][id1] += -Coef/vol1;
                    A[id1][id2] += Coef/vol1;
                }
                if( s2 != Element::Ghost )
                {
                    A[id2][id1] += Coef/vol2;
                    A[id2][id2] += -Coef/vol2;
                }
            }

        }
        BARRIER
        if( m->GetProcessorRank() == 0 ) std::cout << "Matrix assemble: " << Timer()-ttt << std::endl;

        m->RemoveGeometricData(table);

        //~ ttt = Timer();
        //~ A.Save("matrix.mtx");
        //~ BARRIER
        //~ if( m->GetProcessorRank() == 0 ) std::cout << "Save matrix: " << Timer()-ttt << std::endl;

        ttt = Timer();

        S.SetMatrix(A);
        S.Solve(b,x);

        BARRIER
        if( m->GetProcessorRank() == 0 ) std::cout << "Solve system: " << Timer()-ttt << std::endl;

        ttt = Timer();

        for(Mesh::iteratorCell it = m->BeginCell(); it != m->EndCell(); it++)
            if( it->GetStatus() != Element::Ghost )
                it->Real(phi) = x[it->Integer(id)];


        BARRIER
        if( m->GetProcessorRank() == 0 ) std::cout << "Retrieve data: " << Timer()-ttt << std::endl;

        ttt = Timer();
        m->ExchangeData(phi,CELL,0);
        BARRIER
        if( m->GetProcessorRank() == 0 ) std::cout << "Exchange phi: " << Timer()-ttt << std::endl;

        for(int s = 0; s < 3; s++)
        {
            if( m->GetProcessorRank() == 0 ) std::cout << "strategy: " << s << std::endl;

            m->SetParallelStrategy(s);

            ttt = Timer();
            m->ExchangeData(tensor_K,CELL,0);
            BARRIER
            if( m->GetProcessorRank() == 0 ) std::cout << "Exchange dense fixed data: " << Timer()-ttt << std::endl;


            ttt = Timer();
            m->ExchangeData(m->ProcessorsTag(),CELL,0);
            BARRIER
            if( m->GetProcessorRank() == 0 ) std::cout << "Exchange dense variable data: " << Timer()-ttt << std::endl;


            ttt = Timer();
            m->ExchangeData(test,CELL,0);
            BARRIER
            if( m->GetProcessorRank() == 0 ) std::cout << "Exchange sparse fixed data: " << Timer()-ttt << std::endl;
        }




        ttt = Timer();
        m->RemoveGhost();
        BARRIER
        if( m->GetProcessorRank() == 0 ) std::cout << "Remove Ghost: " << Timer()-ttt << std::endl;

        ttt = Timer();
        m->ExchangeGhost(1,FACE);
        BARRIER
        if( m->GetProcessorRank() == 0 ) std::cout << "Exchange ghost: " << Timer()-ttt << std::endl;

        ttt = Timer();
        m->RemoveGhost();
        BARRIER
        if( m->GetProcessorRank() == 0 ) std::cout << "Remove Ghost: " << Timer()-ttt << std::endl;


        ttt = Timer();
        m->Save("grid.pmf");
        BARRIER
        if( m->GetProcessorRank() == 0 ) std::cout << "Save pmf(MPI_File): " << Timer()-ttt << std::endl;


        m->SetParallelFileStrategy(0);

        ttt = Timer();
        m->Save("grid2.pmf");
        BARRIER
        if( m->GetProcessorRank() == 0 ) std::cout << "Save pmf(MPI_Gather): " << Timer()-ttt << std::endl;

        ttt = Timer();
        m->Save("grid.pvtk");
        BARRIER
        if( m->GetProcessorRank() == 0 ) std::cout << "Save pvtk: " << Timer()-ttt << std::endl;

        delete m;
    }
    else
    {
        std::cout << argv[0] << " [mesh_file]" << std::endl;
    }
#if defined(USE_PARTITIONER)
    Partitioner::Finalize();
#endif
    Solver::Finalize();
    Mesh::Finalize();
    return 0;
}