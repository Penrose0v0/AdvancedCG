#include "HalfEdgeDataStructure.h"
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <unordered_map>
#include "fast_obj.h"

namespace HalfEdge {

	using namespace std;
	using namespace glm;

	int Vertex::countValence() const
	{
		int valence = 0;
		bool onBoundary = false;
		auto he = pHalfEdge;

		do
		{
			++valence;

			if (he->pPair == nullptr)
			{
				onBoundary = true;
				break;
			}

			he = he->pPair->pNext;
		} while (he != pHalfEdge);

		if (onBoundary)
		{
			he = pHalfEdge->pPrev;

			do
			{
				++valence;

				if (he->pPair == nullptr)
				{
					break;
				}

				he = he->pPair->pPrev;

			} while (he != pHalfEdge);
		}

		return valence;
	}

	bool Vertex::onBoundary() const
	{
		auto he = pHalfEdge;

		do
		{
			if (he->pPair == nullptr)
			{
				return true;
			}

			he = he->pPair->pNext;
		} while (he != pHalfEdge);

		return false;
	}

	vec3 Face::calcCentroidPosition() const
	{
		int nFaceVertices = 0;
		vec3 positionSum(0.f);
		auto he = pHalfEdge;

		do
		{
			positionSum += he->pStartVertex->position;
			nFaceVertices++;
			he = he->pNext;
		} while (he != pHalfEdge);

		return positionSum / (float)nFaceVertices;
	}

	void Mesh::deepCopy(const Mesh& rhs)
	{
		if (this == &rhs)
			return;

		clear();

		const int nVertices = (int)rhs.vertices.size();
		vertices.resize(nVertices);
		for (int vi = 0; vi < nVertices; ++vi)
			vertices[vi] = new Vertex(vi);

		const int nHalfEdges = (int)rhs.halfEdges.size();
		halfEdges.resize(nHalfEdges);
		for (int hi = 0; hi < nHalfEdges; ++hi)
			halfEdges[hi] = new HalfEdge(hi);

		const int nFaces = (int)rhs.faces.size();
		faces.resize(nFaces);
		for (int fi = 0; fi < nFaces; ++fi)
			faces[fi] = new Face(fi);

		for (int vi = 0; vi < nVertices; ++vi)
		{
			vertices[vi]->pHalfEdge = halfEdges[rhs.vertices[vi]->pHalfEdge->id];
			vertices[vi]->position = rhs.vertices[vi]->position;
		}

		for (int fi = 0; fi < nFaces; ++fi)
			faces[fi]->pHalfEdge = halfEdges[rhs.faces[fi]->pHalfEdge->id];

		for (int hi = 0; hi < nHalfEdges; ++hi)
		{
			HalfEdge* he = halfEdges[hi];
			he->pStartVertex = vertices[rhs.halfEdges[hi]->pStartVertex->id];
			if (rhs.halfEdges[hi]->pPair != nullptr)
				he->pPair = halfEdges[rhs.halfEdges[hi]->pPair->id];
			he->pNext = halfEdges[rhs.halfEdges[hi]->pNext->id];
			he->pPrev = halfEdges[rhs.halfEdges[hi]->pPrev->id];
			if (rhs.halfEdges[hi]->pFace != nullptr)
				he->pFace = faces[rhs.halfEdges[hi]->pFace->id];
		}
	}

	void Mesh::build(const PolygonMesh& mesh, bool checkConsistency /*= false*/)
	{
		if (mesh.getVertices().empty() || mesh.getFaceIndices().empty())
		{
			cerr << __FUNCTION__ << ": geometry not defined" << endl;
			return;
		}

		// prepare vertices

		const int nVertices = mesh.getNumVertices();
		clearVertices();
		vertices.resize(nVertices);

		for (int vi = 0; vi < nVertices; ++vi)
		{
			vertices[vi] = new Vertex(vi);
			vertices[vi]->position = mesh.getVertices()[vi];
		}

		// prepare faces and half-edges

		const int nFaces = (int)mesh.getFaceIndices().size();
		clearFaces();
		clearHalfEdges();
		faces.resize(nFaces);

		unordered_map<long long, HalfEdge*> halfEdgeDict;

		for (int fi = 0; fi < nFaces; ++fi)
		{
			Face* newFace = new Face(fi);
			faces[fi] = newFace;

			vector<HalfEdge*> sameFaceHalfEdges;

			for (const auto t : mesh.getFaceIndices()[fi].indices)
			{
				HalfEdge* he = new HalfEdge(halfEdges.size());
				halfEdges.push_back(he);
				sameFaceHalfEdges.push_back(he);

				Vertex* pStartVertex = vertices[t.vertexIdx];
				he->pStartVertex = pStartVertex;
				if (pStartVertex->pHalfEdge == nullptr)
					pStartVertex->pHalfEdge = he;

				he->pFace = newFace;
			}

			newFace->pHalfEdge = sameFaceHalfEdges.front();

			const int nSameFaceHE = (int)sameFaceHalfEdges.size();
			for (int si = 0; si < nSameFaceHE; ++si)
			{
				HalfEdge* curr = sameFaceHalfEdges[si];
				curr->pPrev = sameFaceHalfEdges[(si - 1 + nSameFaceHE) % nSameFaceHE];
				curr->pNext = sameFaceHalfEdges[(si + 1) % nSameFaceHE];

				const long long startVertexID = curr->pStartVertex->id;
				const long long endVertexID = curr->pNext->pStartVertex->id;
				const long long hashKey = startVertexID * (long long)nVertices + endVertexID;
				halfEdgeDict[hashKey] = curr;
			}
		}

		// setting pairs of half-edges

		for (auto pHalfEdge : halfEdges)
		{
			const long long startVertexID = pHalfEdge->pStartVertex->id;
			const long long endVertexID = pHalfEdge->pNext->pStartVertex->id;
			const long long pairHashKey = endVertexID * (long long)nVertices + startVertexID;
			auto itr = halfEdgeDict.find(pairHashKey);

			if (itr != halfEdgeDict.end())
			{
				HalfEdge* pPair = itr->second;
				pHalfEdge->pPair = pPair;
			}
		}

		if (checkConsistency)
			checkDataConsistency();
	}

	void Mesh::restore(PolygonMesh& mesh) const
	{
		// restore positions

		const int nPositions = (int)vertices.size();
		vector<glm::vec3> tmpVertices(nPositions);
		for (int pi = 0; pi < nPositions; ++pi)
			tmpVertices[pi] = vertices[pi]->position;

		mesh.setVertices(tmpVertices);

		// restore faces

		const int nFaces = (int)faces.size();
		vector<FaceIndices> tmpFaces(nFaces);
		for (int fi = 0; fi < nFaces; ++fi)
		{
			auto& targetIndices = tmpFaces[fi].indices;
			targetIndices.clear();

			HalfEdge* startHE = faces[fi]->pHalfEdge;
			HalfEdge* he = startHE;
			do
			{
				const int vertexIdx = he->pStartVertex->id;
				targetIndices.emplace_back(vertexIdx, -1, vertexIdx);

				he = he->pNext;
			} while (he != startHE);
		}

		mesh.setFaces(tmpFaces);
	}

	void Mesh::checkDataConsistency() const
	{
		int nInconsistentVertices = 0;
		for (int vi =0; vi<(int)vertices.size(); ++vi)
		{
			auto v = vertices[vi];

			bool isInconsistent = false;

			if (v->id < 0)
			{
				cerr << __FUNCTION__ << ": Warning: Vertex[" << vi << "] has a negative id" << endl;
				isInconsistent = true;
			}
			else if (v->pHalfEdge == nullptr)
			{
				cerr << __FUNCTION__ << ": Warning: " << (*v) << "'s pHalfEdge is null" << endl;
				isInconsistent = true;
			}
			else if (v->pHalfEdge != halfEdges[v->pHalfEdge->id])
			{
				cerr << __FUNCTION__ << ": Warning: " << (*v) << "'s pHalfEdge has an invalid address" << endl;
				isInconsistent = true;
			}

			nInconsistentVertices += isInconsistent;
		}

		int nInconsistentFaces = 0;
		for (int fi = 0; fi <(int)faces.size(); ++fi)
		{
			auto f = faces[fi];
			bool isInconsistent = false;

			if (f->id < 0)
			{
				cerr << __FUNCTION__ << ": Warning: Face[" << fi << "] has a negative id" << endl;
				isInconsistent = true;
			}
			else if (f->pHalfEdge == nullptr)
			{
				cerr << __FUNCTION__ << ": Warning: " << (*f) << "'s pHalfEdge is null" << endl;
				isInconsistent = true;
			}
			else if (f->pHalfEdge != halfEdges[f->pHalfEdge->id])
			{
				cerr << __FUNCTION__ << ": Warning: " << (*f) << "'s pHalfEdge has an invalid address" << endl;
				isInconsistent = true;
			}

			nInconsistentFaces += isInconsistent;
		}

		int nInconsistentHalfEdges = 0;
		for (int hi = 0; hi<(int)halfEdges.size(); ++hi)
		{
			auto he = halfEdges[hi];
			bool isInconsistent = false;

			if (he->id < 0)
			{
				cerr << __FUNCTION__ << ": Warning: HalfEdge[" << hi << "] has a negative id" << endl;
				isInconsistent = true;
			}
			else if (he->pStartVertex == nullptr)
			{
				cerr << __FUNCTION__ << ": Warning: " << (*he) << "'s pStartVertex is null" << endl;
				isInconsistent = true;
			}
			else if (he->pStartVertex != vertices[he->pStartVertex->id])
			{
				cerr << __FUNCTION__ << ": Warning: " << (*he) << "'s pStartVertex has an invalid address" << endl;
				isInconsistent = true;
			}

			if (he->pNext == nullptr)
			{
				cerr << __FUNCTION__ << ": Warning: " << (*he) << "'s pNext is null" << endl;
				isInconsistent = true;
			}
			else if (he->pNext != halfEdges[he->pNext->id])
			{
				cerr << __FUNCTION__ << ": Warning: " << (*he) << "'s pNext has an invalid address" << endl;
				isInconsistent = true;
			}

			if (he->pPrev == nullptr)
			{
				cerr << __FUNCTION__ << ": Warning: " << (*he) << "'s pPrev is null" << endl;
				isInconsistent = true;
			}
			else if (he->pPrev != halfEdges[he->pPrev->id])
			{
				cerr << __FUNCTION__ << ": Warning: " << (*he) << "'s pPrev has an invalid address" << endl;
				isInconsistent = true;
			}

			if (he->pPair != nullptr && he->pPair != halfEdges[he->pPair->id])
			{
				cerr << __FUNCTION__ << ": Warning: " << (*he) << "'s pPair has an invalid address" << endl;
				isInconsistent = true;
			}

			if (he->pFace == nullptr)
			{
				cerr << __FUNCTION__ << ": Warning: " << (*he) << "'s pFace is null" << endl;
				isInconsistent = true;
			}
			else if (he->pFace != faces[he->pFace->id])
			//if (he->pFace != nullptr && he->pFace != faces[he->pFace->id])
			{
				cerr << __FUNCTION__ << ": Warning: " << (*he) << "'s pFace has an invalid address" << endl;
				isInconsistent = true;
			}

			nInconsistentHalfEdges += isInconsistent;
		}

		cerr << __FUNCTION__ << ": inconsistent # verts = " << nInconsistentVertices << "/" << vertices.size()
			<< ", # half edges = " << nInconsistentHalfEdges << "/" << halfEdges.size()
			<< ", # faces = " << nInconsistentFaces << "/" << faces.size() << endl;
	}

	ostream& operator<<(ostream& stream, const Vertex& v)
	{
		stream << "Vertex[id=" << v.id << ",pos=" << to_string(v.position) << "]";
		return stream;
	}

	ostream& operator<<(ostream& stream, const Face& f)
	{
		stream << "Face[id=" << f.id << "]";
		return stream;
	}

	ostream& operator<<(ostream& stream, const HalfEdge& he)
	{
		stream << "HalfEdge[id=" << he.id << ",startVertex=" << he.pStartVertex->id << ",endVertex=" << he.pNext->pStartVertex->id << "]";
		return stream;
	}

}

void HalfEdge::Helper::SetPair(HalfEdge* he1, HalfEdge* he2)
{
	he1->pPair = he2;
	he2->pPair = he1;
}

void HalfEdge::Helper::SetPrevNext(HalfEdge* prevHE, HalfEdge* nextHE)
{
	prevHE->pNext = nextHE;
	nextHE->pPrev = prevHE;
}
