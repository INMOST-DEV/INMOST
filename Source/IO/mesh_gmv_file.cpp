#ifdef _MSC_VER //kill some warnings
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "inmost.h"

#if defined(USE_MESH)


namespace INMOST
{
  void Mesh::SaveGMV(std::string File)
  {
		Storage::integer keynum;
		Storage::real keyval;
		//ReorderEmpty(CELL | FACE | NODE | ESET);
		FILE * file = fopen(File.c_str(),"wb");
		char keyword[2048];
		sprintf(keyword,"gmvinput"); fwrite(keyword,1,8,file);
		sprintf(keyword,"ieee");
		if( sizeof(Storage::real) != 8 || sizeof(Storage::integer) != 8 )
			sprintf(keyword,"ieeei%ldr%ld",sizeof(Storage::integer),sizeof(Storage::real));
		fwrite(keyword,1,8,file);
		sprintf(keyword,"nodev"); fwrite(keyword,1,8,file);
		keynum = static_cast<Storage::integer>(NumberOfNodes()); fwrite(&keynum,sizeof(Storage::integer),1,file);
		for(Mesh::iteratorNode n = BeginNode(); n != EndNode(); n++)
		{
			fwrite(n->Coords().data(),sizeof(Storage::real),GetDimensions(),file);
			if( GetDimensions() < 3 )
			{
				Storage::real zero = 0.0;
				fwrite(&zero,sizeof(Storage::real),3-GetDimensions(),file);
			}
		}
		sprintf(keyword,"faces"); fwrite(keyword,1,8,file);
		keynum = static_cast<Storage::integer>(NumberOfFaces()); fwrite(&keynum,sizeof(Storage::integer),1,file);
		keynum = static_cast<Storage::integer>(NumberOfCells()); fwrite(&keynum,sizeof(Storage::integer),1,file);
		Tag set_id = CreateTag("TEMPORARY_ELEMENT_ID",DATA_INTEGER,CELL|FACE|NODE,NONE,1);
		Storage::integer cur_num = 0;
		for(Mesh::iteratorNode it = BeginNode(); it != EndNode(); ++it) it->IntegerDF(set_id) = cur_num++;
		cur_num = 0;
		for(Mesh::iteratorFace it = BeginFace(); it != EndFace(); ++it) it->IntegerDF(set_id) = cur_num++;
		cur_num = 0;
		for(Mesh::iteratorCell it = BeginCell(); it != EndCell(); ++it) it->IntegerDF(set_id) = cur_num++;
		for(Mesh::iteratorFace f = BeginFace(); f != EndFace(); f++)
		{
			ElementArray<Node> fnodes = f->getNodes();
			//sprintf(keyword,"nverts"); fwrite(keyword,1,8,file);
			keynum = static_cast<Storage::integer>(fnodes.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
			//sprintf(keyword,"vertex_ids"); fwrite(keyword,1,8,file);
			for(ElementArray<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
			{
				keynum = fn->IntegerDF(set_id)+1;
				fwrite(&keynum,sizeof(Storage::integer),1,file);
			}
			ElementArray<Cell> fcells = f->getCells();
			//sprintf(keyword,"cellno1"); fwrite(keyword,1,8,file);
			keynum = fcells.size() > 0 ? fcells[0].IntegerDF(set_id)+1 : 0;
			fwrite(&keynum,sizeof(Storage::integer),1,file);
			//sprintf(keyword,"cellno2"); fwrite(keyword,1,8,file);
			keynum = fcells.size() > 1 ? fcells[1].IntegerDF(set_id)+1 : 0;
			fwrite(&keynum,sizeof(Storage::integer),1,file);
		}
			
			
		if( HaveTag("MATERIAL") )
		{
			Tag t = GetTag("MATERIAL");
			if( t.isDefined(CELL) )
			{
				std::set<Storage::integer> mat;
				for(Mesh::iteratorCell c = BeginCell(); c != EndCell(); c++)
					mat.insert(c->Integer(t)+1);
				sprintf(keyword,"material"); fwrite(keyword,1,8,file);
				keynum = static_cast<Storage::integer>(mat.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
				keynum = 0; fwrite(&keynum,sizeof(Storage::integer),1,file);
				for(std::set<Storage::integer>::iterator it = mat.begin(); it != mat.end(); it++)
				{
					sprintf(keyword,"mat%d",*it); fwrite(keyword,1,8,file);
				}
				for(Mesh::iteratorCell c = BeginCell(); c != EndCell(); c++)
				{
					keynum = c->Integer(t)+1; fwrite(&keynum,sizeof(Storage::integer),1,file);
				}
			}
				
		}
			
		sprintf(keyword,"variable"); fwrite(keyword,1,8,file);
		for(ElementType etype = CELL; etype >= NODE; etype = etype >> 1)
		{
			if( etype == EDGE ) continue;
			for(Mesh::iteratorTag t = BeginTag(); t != EndTag(); t++) 
				if( t->isDefined(etype) && 
            t->GetSize() == 1 && 
            !t->isSparse(etype) && 
            t->GetTagName().substr(0,9) != "PROTECTED" && 
            t->GetDataType() != DATA_REFERENCE &&
            t->GetDataType() != DATA_REMOTE_REFERENCE)
				{
					sprintf(keyword,"%s",t->GetTagName().substr(0,8).c_str());
					fwrite(keyword,1,8,file);
					switch(etype)
					{
						case NODE: keynum = 1; break;
						case FACE: keynum = 2; break;
						case CELL: keynum = 0; break;
						default: throw NotImplemented;
					}
					fwrite(&keynum,sizeof(Storage::integer),1,file);
					for(Mesh::iteratorElement e = BeginElement(etype); e != EndElement(); e++)
					{
						keyval = 1e+20;
						switch(t->GetDataType())
						{
							case DATA_INTEGER: keyval = static_cast<Storage::real>(e->Integer(*t)); break;
							case DATA_REAL: keyval = e->Real(*t); break;
							case DATA_BULK: keyval = static_cast<Storage::real>(e->Bulk(*t)); break;
#if defined(USE_AUTODIFF)
              case DATA_VARIABLE: keyval = e->Variable(*t).GetValue(); break;
#endif
							default: throw NotImplemented;
						}
						fwrite(&keyval,sizeof(Storage::real),1,file);
					}
				}
		}
		sprintf(keyword,"endvars"); fwrite(keyword,1,8,file);
			
		sprintf(keyword,"subvars"); fwrite(keyword,1,8,file);
		for(ElementType etype = CELL; etype >= NODE; etype = etype >> 1)
		{
			if( etype == EDGE ) continue;
			for(Mesh::iteratorTag t = BeginTag(); t != EndTag(); t++) 
				if( t->isDefined(etype) && 
            t->GetSize() == 1 && 
            t->isSparse(etype) && 
            t->GetDataType() != DATA_REFERENCE && 
            t->GetDataType() != DATA_REMOTE_REFERENCE)
				{
					Storage::integer temp;
					keynum = 0;
					for(Mesh::iteratorElement e = BeginElement(etype); e != EndElement(); e++)
						if( e->HaveData(*t) ) 
							keynum++;
					if( keynum )
					{
						sprintf(keyword,"%s",t->GetTagName().substr(0,8).c_str());
						fwrite(keyword,1,8,file);
						switch(etype)
						{
							case NODE: temp = 1; break;
							case FACE: temp = 2; break;
							case CELL: temp = 0; break;
							default: throw NotImplemented;
						}
						fwrite(&temp,sizeof(Storage::integer),1,file);
						fwrite(&keynum,sizeof(Storage::integer),1,file);
						for(Mesh::iteratorElement e = BeginElement(etype); e != EndElement(); e++)
							if( e->HaveData(*t) ) 
							{
								keynum = e->IntegerDF(set_id)+1;
								fwrite(&keynum,sizeof(Storage::integer),1,file);
							}
						for(Mesh::iteratorElement e = BeginElement(etype); e != EndElement(); e++)
							if( e->HaveData(*t) )
							{
								switch(t->GetDataType())
								{
									case DATA_INTEGER: keyval = static_cast<Storage::real>(e->Integer(*t)); break;
									case DATA_REAL: keyval = e->Real(*t); break;
									case DATA_BULK: keyval = static_cast<Storage::real>(e->Bulk(*t)); break;
#if defined(USE_AUTODIFF)
                  case DATA_VARIABLE: keyval = e->Variable(*t).GetValue(); break;
#endif
									default: throw NotImplemented;
								}
								fwrite(&keyval,sizeof(Storage::real),1,file);
							}
					}
				}
		}
		sprintf(keyword,"endsubv"); fwrite(keyword,1,8,file);
			
		sprintf(keyword,"groups"); fwrite(keyword,1,8,file);
		for(Mesh::iteratorSet set = BeginSet(); set != EndSet(); set++)
		{

			for(ElementType etype = CELL; etype >= NODE; etype = etype >> 1)
			{
				if( etype == EDGE ) continue;
				keynum = 0;
				for(ElementSet::iterator it = set->Begin(); it != set->End(); it++)
					if( it->GetElementType() == etype ) keynum++;
						
				if( keynum )
				{
					Storage::integer temp;
					switch(etype)
					{
						case NODE: temp = 1; break;
						case FACE: temp = 2; break;
						case CELL: temp = 0; break;
						default: throw NotImplemented;
					}
					//std::cout << "set name: " << set->GetName() << " size " << set->Size() << " id " << set->LocalID() << std::endl;
					//sprintf(keyword,"set%d_%s\n",set->IntegerDF(set_id)+1,ElementTypeName(etype));
					sprintf(keyword,"%s",set->GetName().c_str());
					fwrite(keyword,1,8,file);
					fwrite(&temp,sizeof(Storage::integer),1,file);						
					fwrite(&keynum,sizeof(Storage::integer),1,file);
					for(ElementSet::iterator it = set->Begin(); it != set->End(); it++)
						if( it->GetElementType() == etype )
						{
							//std::cout << "Write element " << ElementTypeName(it->GetElementType()) << " num " << it->LocalID() << " handle " << it->GetHandle() << std::endl;
							keynum = it->IntegerDF(set_id)+1;
							fwrite(&keynum,sizeof(Storage::integer),1,file);
						} 
				}	
			}
		}
		sprintf(keyword,"endgrp"); fwrite(keyword,1,8,file);
			
		sprintf(keyword,"vectors"); fwrite(keyword,1,8,file);
		for(ElementType etype = CELL; etype >= NODE; etype = etype >> 1)
		{
			if( etype == EDGE ) continue;
			for(Mesh::iteratorTag t = BeginTag(); t != EndTag(); t++) 
				if( t->isDefined(etype) && 
            t->GetSize() != 1 && 
            t->GetSize() != ENUMUNDEF && 
            !t->isSparse(etype) && 
            t->GetTagName().substr(0,9) != "PROTECTED" && 
            t->GetDataType() != DATA_REFERENCE &&
            t->GetDataType() != DATA_REMOTE_REFERENCE
            )
				{
					sprintf(keyword,"%s",t->GetTagName().substr(0,8).c_str());
					fwrite(keyword,1,8,file);
					switch(etype)
					{
						case NODE: keynum = 1; break;
						case FACE: keynum = 2; break;
						case CELL: keynum = 0; break;
						default: throw NotImplemented;
					}
					fwrite(&keynum,sizeof(Storage::integer),1,file);
					keynum = t->GetSize();
					fwrite(&keynum,sizeof(Storage::integer),1,file);
					keynum = 0;
					fwrite(&keynum,sizeof(Storage::integer),1,file);
						
					for(INMOST_DATA_ENUM_TYPE q = 0; q < t->GetSize(); q++)
					{
						for(Mesh::iteratorElement e = BeginElement(etype); e != EndElement(); e++)
						{
							switch(t->GetDataType())
							{
								case DATA_INTEGER: keyval = static_cast<Storage::real>(e->IntegerArray(*t)[q]); break;
								case DATA_REAL: keyval = e->RealArray(*t)[q]; break;
								case DATA_BULK: keyval = static_cast<Storage::real>(e->BulkArray(*t)[q]); break;
#if defined(USE_AUTODIFF)
                case DATA_VARIABLE: keyval = e->VariableArray(*t)[q].GetValue(); break;
#endif
								default: throw NotImplemented;
							}
							fwrite(&keyval,sizeof(Storage::real),1,file);
						}
					}
				}
		}
		sprintf(keyword,"endvect"); fwrite(keyword,1,8,file);
			
		sprintf(keyword,"polygons"); fwrite(keyword,1,8,file);
		for(Mesh::iteratorFace f = BeginFace(); f != EndFace(); f++) if( f->Boundary() )
		{
			ElementArray<Node> fnodes = f->getNodes();
			keynum = 1; fwrite(&keynum,sizeof(Storage::integer),1,file);
			keynum = static_cast<Storage::integer>(fnodes.size()); fwrite(&keynum,sizeof(Storage::integer),1,file);
			for(ElementArray<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
				fwrite(&fn->Coords()[0],sizeof(Storage::real),1,file);
			for(ElementArray<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
				fwrite(&fn->Coords()[1],sizeof(Storage::real),1,file);
			if( GetDimensions() > 2 )
			{
				for(ElementArray<Node>::iterator fn = fnodes.begin(); fn != fnodes.end(); fn++)
					fwrite(&fn->Coords()[2],sizeof(Storage::real),1,file);
			}
			else 
			{
				Storage::real zero = 0.0;
				fwrite(&zero,sizeof(Storage::real),fnodes.size(),file);
			}
		}
		sprintf(keyword,"endpoly"); fwrite(keyword,1,8,file);
		sprintf(keyword,"endgmv"); fwrite(keyword,1,8,file);
		fclose(file);
	}
}

#endif