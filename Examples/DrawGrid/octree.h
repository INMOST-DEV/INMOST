#ifndef _OCTREE_H
#define _OCTREE_H

#include "inmost.h"

namespace INMOST
{

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
			for (int k = 0; k < numchildren; ++k)
			{
				offsets[k] = 0;
				sizes[k] = 0;
			}
			for (int r = 0; r < size; ++r)
			{
				Element c = Element(GetMeshLink(), cells[r]);
				c->Centroid(cell_center);
				int child_num = 0;
				for (int k = 0; k < dims; ++k)
				{
					if (cell_center[k] > center[k])
					{
						int m = 1 << k;
						child_num += m;
					}
				}
				c->IntegerDF(child_tag) = child_num;
				sizes[child_num]++;
			}
			for (int k = 1; k < numchildren; ++k)
			{
				offsets[k] = offsets[k - 1] + sizes[k - 1];
			}
			for (int k = 0; k < numchildren; ++k)
			{
				std::stringstream name;
				name << GetName() << "chld" << k;
				ElementSet child = GetMeshLink()->CreateSetUnique(name.str()).first;
				Storage::real_array child_center = child->RealArray(center_tag);
				for (int r = 0; r < dims; ++r)
				{
					int l = 1 << r;
					int m = k & l;
					child_center[r] = center[r] + ((m ? 1.0 : -1.0) * center[r + 3] * 0.25);
					child_center[r + 3] = center[r + 3] * 0.5;
				}
				int m = 0;
				for (int r = 0; r < size; ++r)
				{
					Element c = Element(GetMeshLink(), cells[r]);
					int q = c->IntegerDF(child_tag);
					if (q == k) (temp + offsets[k])[m++] = cells[r];
				}
				AddChild(child);
				if (sizes[k] <= 16 && sizes[k] > 0)
					child->PutElements(temp + offsets[k], sizes[k]);
			}
			// cells array is not needed anymore
			ElementSet child = GetChild();
			for (int k = 0; k < numchildren; ++k)
			{
				if (sizes[k] > 16)
					Octree(child).SubConstruct(child_tag, center_tag, temp + offsets[k], cells + offsets[k], sizes[k], quad_tree);
				child = child->GetSibling();
			}
		}
		Element SubFindClosestCell(const Tag & center_tag, Storage::real pnt[3], bool quad_tree) const
		{
			Storage::real_array center = RealArray(center_tag);
			if (Inside(center, pnt, quad_tree))
			{
				if (HaveChild())
				{
					int child_num = 0, q;
					int dims = 3 - (quad_tree ? 1 : 0);
					for (int k = 0; k < dims; ++k)
					{
						if (pnt[k] > center[k])
							child_num += (1 << k);
					}
					q = 0;
					ElementSet set = GetChild();
					while (q != child_num) { set = set->GetSibling(); q++; }
					return Octree(set).SubFindClosestCell(center_tag, pnt, quad_tree);
				}
				else
				{
					HandleType * cells = getHandles();
					int ncells = (int)nbHandles();
					Element closest = InvalidElement();
					Storage::real mindist = 1.0e20, dist, cnt[3];
					for (int k = 0; k < ncells; ++k)
					{
						Element c(GetMeshLink(), cells[k]);
						c->Centroid(cnt);
						dist = sqrt((cnt[0] - pnt[0])*(cnt[0] - pnt[0]) + (cnt[1] - pnt[1])*(cnt[1] - pnt[1]) + (cnt[2] - pnt[2])*(cnt[2] - pnt[2]));
						if (mindist > dist)
						{
							mindist = dist;
							closest = c;
						}
					}
					return closest;
				}
			}
			else return InvalidCell();
		}
		bool Inside(const Storage::real_array & center, Storage::real pnt[3], bool quad_tree) const
		{
			bool inside = true;
			int dims = 3 - (quad_tree ? 1 : 0);
			for (int i = 0; i < dims; ++i)
				inside &= (pnt[i] >= center[i] - center[3 + i] * 0.5 && pnt[i] <= center[i] + center[3 + i] * 0.5);
			return inside;
		}
		Node SubFindClosestNode(const Tag & center_tag, Storage::real pnt[3], bool quad_tree) const
		{
			if (HaveChild())
			{
				Storage::real_array center = RealArray(center_tag);
				if (!Inside(center, pnt, quad_tree)) return InvalidNode();
				int child_num = 0, q;
				int dims = 3 - (quad_tree ? 1 : 0);
				for (int k = 0; k < dims; ++k)
				{
					if (pnt[k] > center[k])
						child_num += (1 << k);
				}
				q = 0;
				ElementSet set = GetChild();
				while (q != child_num) { set = set->GetSibling(); q++; }
				return Octree(set).SubFindClosestNode(center_tag, pnt, quad_tree);
			}
			else
			{
				HandleType * cells = getHandles();
				int ncells = (int)nbHandles();
				Node closest = InvalidNode();
				Storage::real mindist = 1.0e20, dist;
				for (int k = 0; k < ncells; ++k)
				{
					Element c = Element(GetMeshLink(), cells[k]);
					ElementArray<Node> nodes = c->getNodes();
					for (ElementArray<Node>::iterator q = nodes.begin(); q != nodes.end(); ++q)
					{
						Storage::real_array cnt = q->Coords();
						dist = sqrt((cnt[0] - pnt[0])*(cnt[0] - pnt[0]) + (cnt[1] - pnt[1])*(cnt[1] - pnt[1]) + (cnt[2] - pnt[2])*(cnt[2] - pnt[2]));
						if (mindist > dist)
						{
							mindist = dist;
							closest = q->self();
						}
					}
				}
				return closest;
			}
		}
		void SubDestroy()
		{
			if (HaveChild())
			{
				ElementSet set = GetChild(), next;
				while (set->isValid())
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
			Tag child_tag = GetMeshLink()->CreateTag("OCTREE_CHILD_NUM_" + GetName(), DATA_INTEGER, elem, NONE, 1);
			save_center_tag = GetMeshLink()->CreateTag("OCTREE_CENTER_" + GetName(), DATA_REAL, ESET, ESET, 6);
			Storage::real bounds[3][2];
			for (int k = 0; k < dims; ++k)
			{
				bounds[k][0] = 1.0e20;
				bounds[k][1] = -1.0e20;
			}
			//calculate bounds
			for (Mesh::iteratorNode node = GetMeshLink()->BeginNode(); node != GetMeshLink()->EndNode(); ++node)
			{
				Storage::real_array coord = node->Coords();
				for (int k = 0; k < dims; ++k)
				{
					if (coord[k] < bounds[k][0]) bounds[k][0] = coord[k];
					if (coord[k] > bounds[k][1]) bounds[k][1] = coord[k];
				}
			}
			Storage::real_array center_data = RealArray(save_center_tag);
			for (int k = 0; k < dims; ++k)
			{
				center_data[k] = (bounds[k][0] + bounds[k][1])*0.5; //central position
				center_data[k + 3] = bounds[k][1] - bounds[k][0]; //length
			}
			//copy cells
			int size = GetMeshLink()->NumberOf(elem), k = 0;
			HandleType * cells = new HandleType[size * 2];
			HandleType * temp = cells + size;
			for (Mesh::iteratorElement cell = GetMeshLink()->BeginElement(elem); cell != GetMeshLink()->EndElement(); ++cell)
				cells[k++] = *cell;
			SubConstruct(child_tag, save_center_tag, cells, temp, size, quad_tree);
			GetMeshLink()->DeleteTag(child_tag);
		}
		Element FindClosestCell(Storage::real pnt[3]) const
		{
			return SubFindClosestCell(save_center_tag, pnt, save_quad_tree);
		}
		Node FindClosestNode(Storage::real pnt[3]) const
		{
			return SubFindClosestNode(save_center_tag, pnt, save_quad_tree);
		}
		void Destroy()
		{
			if (save_center_tag.isValid())
				GetMeshLink()->DeleteTag(save_center_tag);
			SubDestroy();
		}
		~Octree() { }
	};
}

#endif