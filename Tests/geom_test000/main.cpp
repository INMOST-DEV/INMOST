#include <cstdio>
#include <cmath>

#include "inmost.h"
using namespace INMOST;

typedef Storage::real_array real_array;
typedef Storage::real real;

#define ORIGINAL

__INLINE static void vec_cross_product(const Storage::real * vecin1, const Storage::real * vecin2, Storage::real * vecout)
{
	Storage::real temp[3];
	temp[0] = vecin1[1]*vecin2[2] - vecin1[2]*vecin2[1];
	temp[1] = vecin1[2]*vecin2[0] - vecin1[0]*vecin2[2];
	temp[2] = vecin1[0]*vecin2[1] - vecin1[1]*vecin2[0];
	vecout[0] = temp[0];
	vecout[1] = temp[1];
	vecout[2] = temp[2];
}


int main(int argc,char ** argv)
{
	int errors = 0;
	Mesh::Initialize(&argc,&argv);
	{
#if defined(ORIGINAL)
		std::cout << "inmost algorithms" << std::endl;
#else
		std::cout << "non-inmost algorithms" << std::endl;
#endif
		Mesh m;
		m.Load((argc>1)?argv[1]:"d4.pmf");
		rMatrix c(3,1), n(3,1), x(3,1), y(3,1), z(3,1), ct(3,1), nt(3,1);
		rMatrix e1(3,1,0.0), e2(3,1,0.0), e3(3,1,0.0);
		rMatrix x0(3,1),y0(3,1),z0(3,1);
		rMatrix Et(3,3), Ef(3,3), d(3,1), W(9,3);
		rMatrix nsum(3,1), ntsum(3,1);
		e1(0,0) = 1.0;
		e2(1,0) = 1.0;
		e3(2,0) = 1.0;
		Mesh::GeomParam t;
		t[BARYCENTER] = CELL|FACE;
		t[ORIENTATION] = FACE;
		t[NORMAL] = FACE;
		m.RemoveGeometricData(t);
		//m.PrepareGeometricData(t);
		for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
			it->FixNormalOrientation();
		int np = 0;
		for(Mesh::iteratorFace it = m.BeginFace(); it != m.EndFace(); ++it)
			if( !it->Planarity() ) np++;
		std::cout << "nonplanar: " << np << std::endl;
		
		/*
		for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
		{
			ElementArray<Face> faces = it->getFaces();
			rMatrix R(faces.size(),3), N(faces.size(),3), N2(3,faces.size()), R2(3,faces.size());
			for(unsigned k = 0; k < faces.size(); ++k)
			{
				faces[k].Barycenter(c.data());
				faces[k].OrientedNormal(it->self(),n.data());
				R(k,k+1,0,3) = c.Transpose();
				N(k,k+1,0,3) = n.Transpose();
			}
			std::cout << "inv: " << std::endl;
			(N.Transpose()*R).Invert().first.Print();
			N2 = (N.Transpose()*R).Invert().first*N.Transpose();
			R2 = (R.Transpose()*N).Invert().first*R.Transpose();
			nsum.Zero();
			for(unsigned k = 0; k < faces.size(); ++k)
				nsum += N2(0,3,k,k+1);
			std::cout << "Original N: " << std::endl;
			N.Transpose().Print();
			std::cout << "Modified N: " << std::endl;
			(N2*it->Volume()).Print();
			std::cout << "Original R: " << std::endl;
			R.Transpose().Print();
			std::cout << "Modified R: " << std::endl;
			(R2*it->Volume()).Print();
			nsum.Transpose().Print();
			(N2*R).Print();
			//((N.Transpose()*R).Invert().first*N.Transpose()*R).Print();
		}
		return 0;
		 */
		/*
		for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
		{
			ElementArray<Face> faces = it->getFaces();
			rMatrix E(3,3),NN(3,3), W(9,9),K(3,3);
			rMatrix R(faces.size(),3), N(faces.size(),3), R2(3,faces.size()),RK(faces.size(),3);
			real v = 0, vk = 0;
			
			E.Zero(); N.Zero(); W.Zero();
			for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
			{
				jt->Barycenter(c.data());
				jt->OrientedNormal(it->self(),n.data());
				E += c*n.Transpose()-(c.DotProduct(n)/3.0*rMatrix::Unit(3));
				NN += n*n.Transpose();
				for(int q = 0; q < 3; ++q)
				{
					W(q*4,0) += n(0,0)*n(0,0)/3.0;
					W(q*4,1) += n(1,0)*n(0,0)/3.0;
					W(q*4,2) += n(2,0)*n(0,0)/3.0;
					W(q*4,3) += n(0,0)*n(1,0)/3.0;
					W(q*4,4) += n(1,0)*n(1,0)/3.0;
					W(q*4,5) += n(2,0)*n(1,0)/3.0;
					W(q*4,6) += n(0,0)*n(2,0)/3.0;
					W(q*4,7) += n(1,0)*n(2,0)/3.0;
					W(q*4,8) += n(2,0)*n(2,0)/3.0;
				}
			}
			std::cout << "Cell " << it->LocalID() << " vol " << it->Volume() << std::endl;
			std::cout << "E:" << std::endl;
			E.Print();
			std::cout << "N:" << std::endl;
			N.Print();
			W += rMatrix::Unit(3).Kronecker(NN);
			std::cout << "W" << std::endl;
			W.Print();
			K = W.Solve(E.Repack(9,1),true).first.Repack(3,3);
			std::cout << "K" << std::endl;
			K.Print();
			E.Zero(); N.Zero();
			int k = 0;
			for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
			{
				jt->Barycenter(c.data());
				jt->OrientedNormal(it->self(),n.data());
				R(k,k+1,0,3) = c.Transpose();
				c -= K*n;
				E += c*n.Transpose()-c.DotProduct(n)/3.0*rMatrix::Unit(3);
				N(k,k+1,0,3) = n.Transpose();
				RK(k,k+1,0,3) = c.Transpose();
				k++;
			}
			v = 0; vk = 0;
			R2 = (R.Transpose()*N).Invert().first*R.Transpose();
			iMatrix nf(1,faces.size());
			for(int k = 0; k < faces.size(); ++k)
			{
				//c = R(0,3,k,k+1);
				c = R(k,k+1,0,3).Transpose();
				n = N(k,k+1,0,3).Transpose();
				v += c.DotProduct(n)/3.0;
				
				c = RK(k,k+1,0,3).Transpose();
				vk += c.DotProduct(n)/3.0;
				
				nf(0,k) = faces[k].LocalID();
			}
			std::cout << "v " << v << " vk " << vk << std::endl;
			std::cout << "E'" << std::endl;
			E.Print();
			std::cout << "faces:" << std::endl;
			nf.Print();
			std::cout << "R" << std::endl;
			R.Transpose().Print();
			std::cout << "R2" << std::endl;
			(R2*v).Print();
			std::cout << "||R2-R||: " << (R2*v-R).FrobeniusNorm() << std::endl;
			std::cout << "RK" << std::endl;
			RK.Transpose().Print();
			std::cout << "||RK-R||: " << (RK-R).FrobeniusNorm() << std::endl;
		}
		return 0;
		*/
		for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
		{
			ElementArray<Face> faces = it->getFaces();
			//face center of mass, face normal, face area, cell volume
			x.Zero(), y.Zero(), z.Zero();
			real vol = 0;
			for(ElementArray<Face>::iterator jt = faces.begin(); jt != faces.end(); ++jt)
			{
#if defined(ORIGINAL)
				jt->Barycenter(c.data());
				jt->OrientedNormal(it->self(),n.data());
				x += c(0,0) * n;
				y += c(1,0) * n;
				z += c(2,0) * n;
#else
				double s = jt->FaceOrientedOutside(it->self()) ? 1.0 : -1.0;
				ElementArray<Node> nodes = jt->getNodes();
				real_array v0 = nodes[0].Coords();
				real_array v1 = nodes[1].Coords();
				real_array v2 = nodes[2].Coords();
				real a = 0, at;
				n.Zero();
				c.Zero();
				real i1 = 0, i2 = 0, ss;
				Et.Zero();
				rMatrix n0(3,1);
				for(int q = 0; q < 3; ++q)
				{
					n0(q,0) = (v1[(q+1)%3]-v0[(q+1)%3])*(v2[(q+2)%3]-v0[(q+2)%3]);
					n0(q,0)-= (v1[(q+2)%3]-v0[(q+2)%3])*(v2[(q+1)%3]-v0[(q+1)%3]);
					n0(q,0)*=-1;
				}
				for(int k = 1; k < nodes.size()-1; ++k)
				{
					real_array v1 = nodes[k+0].Coords();
					real_array v2 = nodes[k+1].Coords();
					for(int q = 0; q < 3; ++q)
					{
						ct(q,0) = (v0[q]+v1[q]+v2[q])/3.0;
						nt(q,0) = (v1[(q+1)%3]-v0[(q+1)%3])*(v2[(q+2)%3]-v0[(q+2)%3]);
						nt(q,0)-= (v1[(q+2)%3]-v0[(q+2)%3])*(v2[(q+1)%3]-v0[(q+1)%3]);
					}
					ss = n0.DotProduct(nt);
					ss /= fabs(ss);
					nt *= 0.5;
					i1 += ct.DotProduct(nt) * s;
					//vol += ct.DotProduct(nt) * s;
					at = sqrt(nt.DotProduct(nt)) * ss;
					n += nt;
					c += ct*at;
					a += at;
					Et += ct * nt.Transpose()*s;
					//std::cout <<  "ct " << ct(0,0) << " " << ct(1,0) << " " << ct(2,0);
					//std::cout << " nt " << nt(0,0) << " " << nt(1,0) << " " << nt(2,0);
					//std::cout << " at " << at << std::endl;
					//x += ct(0,0) * nt * s;
					//y += ct(1,0) * nt * s;
					//z += ct(2,0) * nt * s;
				}
				c /= a;
				//std::cout <<  "c " << c(0,0) << " " << c(1,0) << " " << c(2,0);
				//std::cout << " n " << n(0,0) << " " << n(1,0) << " " << n(2,0);
				//std::cout << " a " << a << std::endl;
				//i2 == i1!
				Et -= i1/3.0*rMatrix::Unit(3);
				
				
				Ef = c * n.Transpose()*s;
				i2 = c.DotProduct(n) * s;
				Ef -= i2/3.0*rMatrix::Unit(3);
				
				if( (Et-Ef).FrobeniusNorm() > 1.0e-9 )
				{
					std::cout << "||Et-Ef||: " << (Et-Ef).FrobeniusNorm() << std::endl;
					real_array v0 = nodes[0].Coords();
					real_array v1 = nodes[1].Coords();
					real_array v2 = nodes[2].Coords();
					n.Zero();
					a = 0;
					rMatrix n0(3,1);
					real l1[3],l2[3];
					for(int q = 0; q < 3; ++q)
					{
						l1[q] = v1[q]-v0[q];
						l2[q] = v2[q]-v0[q];
					}
					vec_cross_product(l1,l2,n0.data());
					real ss;
					for(int k = 1; k < nodes.size()-1; ++k)
					{
						v1 = nodes[k+0].Coords();
						v2 = nodes[k+1].Coords();
						for(int q = 0; q < 3; ++q)
						{
							ct(q,0) = (v0[q]+v1[q]+v2[q])/3.0;
							//nt(q,0) = (v1[(q+1)%3]-v0[(q+1)%3])*(v2[(q+2)%3]-v0[(q+2)%3]);
							//nt(q,0)-= (v1[(q+2)%3]-v0[(q+2)%3])*(v2[(q+1)%3]-v0[(q+1)%3]);
							l1[q] = v1[q]-v0[q];
							l2[q] = v2[q]-v0[q];
						}
						vec_cross_product(l1,l2,nt.data());
						ss = n0.DotProduct(nt);
						ss /= fabs(ss);
						nt *= 0.5;
						n += nt;
						at = sqrt(nt.DotProduct(nt))*ss;
						a += at;
						std::cout << "ct " << ct(0,0) << " " << ct(1,0) << " " << ct(2,0) << " ";
						std::cout << "nt " << nt(0,0)/at << " " << nt(1,0)/at << " " << nt(2,0)/at << " ";
						std::cout << "at " << at << std::endl;
					}
					real dd = n.DotProduct(rMatrix(v0.data(),3,1));
					std::cout << "n " << n(0,0)/a << " " << n(1,0)/a << " " << n(2,0)/a << " ";
					std::cout << "planarity:";
					for(int k = 1; k < nodes.size(); ++k)
					{
						real_array v1 = nodes[k].Coords();
						std::cout << " " << n.DotProduct(rMatrix(v1.data(),3,1))-dd;
					}
					std::cout << std::endl;
					// c = c+K*n
					// Ef = c*n^T - c^T*n/3.0*I
					// (Kc)*n^T - (Kc)^T*n/3.0*I = Et
					// K*c*n^T - c^T*K^T*n/3.0*I = Et
					//            | k00 k01 k02 | nx
					// (cx cy cz) | k10 k11 k12 | ny
					//            | k20 k21 k22 | nz
					// d*n^T - d^T*n/3.0*I = Et
					//
					/*
					W(0,0) = n(0,0) - n(0,0)/3.0; // e_xx
					W(0,1) = - n(1,0)/3.0; // e_xx
					W(0,2) = - n(2,0)/3.0; // e_xx
					W(1,1) = n(1,0); // e_xy
					W(2,2) = n(2,0); // e_xz
					W(3,0) = n(0,0); // e_yx
					W(4,0) = -n(0,0)/3.0;
					W(4,1) = n(1,0) - n(1,0)/3.0; // e_yy
					W(4,2) = -n(2,0)/3.0;
					W(5,2) = n(2,0); // e_yz
					W(6,0) = n(0,0); // e_zx
					W(7,1) = n(1,0); // e_zy
					W(8,0) = -n(0,0)/3.0;
					W(8,1) = -n(1,0)/3.0;
					W(8,2) = n(2,0)-n(2,0)/3.0; // e_zz
					W *= s;
					d = W.Solve((Et-Ef).Repack(9,1)).first;
					c += d;
					
					Ef = c * n.Transpose()*s;
					i2 = c.DotProduct(n) * s;
					Ef -= i2/3.0*rMatrix::Unit(3);
					
					std::cout << "guess ||Et-Ef||: " << (Et-Ef).FrobeniusNorm() << std::endl;
					*/
					//std::cout << "Et: " << std::endl;
					//Et.Print();
					//std::cout << "Before iters:" << std::endl;
					//((c*n.Transpose() - (c.DotProduct(n)/3.0)*rMatrix::Unit(3))*s).Print();
					/*
					vMatrix vn(3,1), vc(3,1), vE(9,1);
					rMatrix W(9,6), E(9,1), up(6,1);
					for(int nit = 0; nit < 50; ++nit)
					{
						for(int q = 0; q < 3; ++q)
						{
							vn(q,0) = unknown(n(q,0)+(rand()/(1.0*RAND_MAX))*1.0e-9,q);
							vc(q,0) = unknown(c(q,0)+(rand()/(1.0*RAND_MAX))*1.0e-9,q);
						}
						vE = ((vc*vn.Transpose() - rMatrix::Unit(3)/3.0*vc.DotProduct(vn))*s - Et).Repack(9,1);
						W.Zero();
						for(int q = 0; q < 9; ++q)
						{
							E(q,0) = vE(q,0).GetValue();
							for(Sparse::Row::iterator qt = vE(q,0).GetRow().Begin(); qt != vE(q,0).GetRow().End(); ++qt)
								W(q,qt->first) = qt->second;
						}
						std::cout << nit << " " << E.DotProduct(E) << " ";
						E.Transpose().Print();
						up = W.Solve(E).first;
						n -= up(0,3,0,1);
						c -= up(3,6,0,1);
					}
					
					Ef = c * n.Transpose()*s;
					i2 = c.DotProduct(n) * s;
					Ef -= i2/3.0*rMatrix::Unit(3);
					
					std::cout << "iters ||Et-Ef||: " << (Et-Ef).FrobeniusNorm() << std::endl;
					*/
					
					//std::cout << "After iters:" << std::endl;
					//((c*n.Transpose() - (c.DotProduct(n)/3.0)*rMatrix::Unit(3))*s).Print();
					
				}
				
				
				
				//
				//std::cout << "area: " << a << std::endl;
				
				vol += c.DotProduct(n) * s;
				x += c(0,0) * n * s;
				y += c(1,0) * n * s;
				z += c(2,0) * n * s;
#endif
			}
			
			//std::cout << "nsum:  " << nsum(0,0) << " " << nsum(1,0) << " " << nsum(2,0) << std::endl;
			//std::cout << "ntsum: " << ntsum(0,0) << " " << ntsum(1,0) << " " << ntsum(2,0) << std::endl;
			
#if defined(ORIGINAL)
			x /= it->Volume();
			y /= it->Volume();
			z /= it->Volume();
			//std::cout << "inmost" << std::endl;
			//x0.ConcatCols(y0).ConcatCols(z0).Print();
			//std::cout << "vol: " << it->Volume() << std::endl;
#else
			vol /= 3.0;
			x /= vol;
			y /= vol;
			z /= vol;
			//std::cout << "non-inmost" << std::endl;
			//x.ConcatCols(y).ConcatCols(z).Print();
			//std::cout << "vol: " << vol << std::endl;
#endif
			if( (x-e1).FrobeniusNorm() > 1.0e-9 || x.CheckNans() )
			{
				std::cout << "On CELL:" << it->LocalID() << " dx: " << x(0,0) << " " << x(1,0) << " " << x(2,0) << std::endl;
				errors++;
			}
			if( (y-e2).FrobeniusNorm() > 1.0e-9|| y.CheckNans() )
			{
				std::cout << "On CELL:" << it->LocalID() << " dy: " << y(0,0) << " " << y(1,0) << " " << y(2,0) << std::endl;
				errors++;
			}
			if( (z-e3).FrobeniusNorm() > 1.0e-9|| z.CheckNans() )
			{
				std::cout << "On CELL:" << it->LocalID() << " dz: " << z(0,0) << " " << z(1,0) << " " << z(2,0) << std::endl;
				errors++;
			}
		}
	}
	Mesh::Finalize();
	if( errors )
	{
		std::cout << "Errors: " << errors << std::endl;
		return 1;
	}
	return 0;
}
