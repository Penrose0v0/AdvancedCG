#include "LoopSubdivision.h"
#include <unordered_map>

using namespace std;
using namespace glm;

void LoopSubdivision::subdivide(PolygonMesh& mesh, int nSubdiv)
{
	if (mesh.getVertices().empty() || mesh.getFaceIndices().empty())
	{
		std::cerr << __FUNCTION__ << ": mesh not ready" << std::endl;
		return;
	}

	mesh.triangulate();

	HalfEdge::Mesh _mesh;
	_mesh.build(mesh);

	for (int iter = 0; iter < nSubdiv; ++iter)
		_mesh = apply(_mesh);

	_mesh.restore(mesh);
	mesh.calcVertexNormals();
}

HalfEdge::Mesh LoopSubdivision::apply(HalfEdge::Mesh& mesh)
{
	// return mesh;	// TODO: delete this line

	HalfEdge::Mesh newMesh;

	const int nOldVertices = (int)mesh.vertices.size();
	const int nOldFaces = (int)mesh.faces.size();
	const int nOldHalfEdges = (int)mesh.halfEdges.size();

	// Step 0: allocate memory for even (i.e., old) vertices

	newMesh.vertices.reserve(nOldFaces + nOldHalfEdges);

	for (int vi = 0; vi < nOldVertices; ++vi)
		newMesh.addVertex();

	// Step 1: create odd (i.e., new) vertices by splitting half edges

	unordered_map<long, HalfEdge::Vertex*> newEdgeMidpointDict;
	unordered_map<long, pair<HalfEdge::HalfEdge*, HalfEdge::HalfEdge*>> newHalfEdgeDict;

	for (int hi = 0; hi < nOldHalfEdges; ++hi)
	{
		auto oldHE = mesh.halfEdges[hi];

		vec3 startVertexPosition = mesh.vertices[oldHE->pStartVertex->id]->position;
		vec3 endVertexPosition = mesh.vertices[oldHE->getEndVertex()->id]->position;

		HalfEdge::Vertex* edgeMidpoint = nullptr;

		if (oldHE->pPair == nullptr)	// on boundary
		{
			edgeMidpoint = newMesh.addVertex();
			// TODO: calculate the position of edge midpoint
			edgeMidpoint->position = 0.5f * (startVertexPosition + endVertexPosition);
		}
		else
		{
			auto edgeMidpointIter = newEdgeMidpointDict.find(oldHE->pPair->id); // check if pair has been already registered

			if (edgeMidpointIter == newEdgeMidpointDict.end())
			{
				edgeMidpoint = newMesh.addVertex();
				// TODO: calculate the position of edge midpoint
				HalfEdge::HalfEdge* pairHE = oldHE->pPair;
				vec3 upVertexPosition = mesh.vertices[oldHE->pNext->getEndVertex()->id]->position; 
				vec3 downVertexPosition = mesh.vertices[pairHE->pNext->getEndVertex()->id]->position;
				edgeMidpoint->position = (3.0f / 8.0f) * (startVertexPosition + endVertexPosition) +
					(1.0f / 8.0f) * (upVertexPosition + downVertexPosition);
			}
			else
			{
				edgeMidpoint = edgeMidpointIter->second;
			}
		}

		newEdgeMidpointDict[oldHE->id] = edgeMidpoint; // used in Step 3

		auto formerHE = newMesh.addHalfEdge();
		auto latterHE = newMesh.addHalfEdge();

		auto evenStartVertex = newMesh.vertices[oldHE->pStartVertex->id];
		auto evenEndVertex = newMesh.vertices[oldHE->getEndVertex()->id];

		formerHE->pStartVertex = evenStartVertex;
		if (evenStartVertex->pHalfEdge == nullptr)
			evenStartVertex->pHalfEdge = formerHE;

		latterHE->pStartVertex = edgeMidpoint;
		if (edgeMidpoint->pHalfEdge == nullptr)
			edgeMidpoint->pHalfEdge = latterHE;

		newHalfEdgeDict[hi] = make_pair(formerHE, latterHE);

		// register pairs

		if (oldHE->pPair != nullptr)
		{
			auto iter = newHalfEdgeDict.find(oldHE->pPair->id);

			if (iter != newHalfEdgeDict.end())
			{
				HalfEdge::HalfEdge* pairFormerHE = iter->second.first;
				HalfEdge::HalfEdge* pairLatterHE = iter->second.second;

				HalfEdge::Helper::SetPair(pairFormerHE, latterHE);
				HalfEdge::Helper::SetPair(pairLatterHE, formerHE);
			}
		}
	}

	// Step 2: update even (i.e., old) vertices

	for (int vi = 0; vi < nOldVertices; ++vi)
	{
		auto newVertex = newMesh.vertices[vi];
		const auto oldVertex = mesh.vertices[vi];
		const auto oldVertexPosition = oldVertex->position;

		// TODO: calculate the new vertex position
		// c.f., HalfEdge::Vertex::countValence() in HalfEdgeDataStructure.cpp
		vec3 sum = vec3(0.0f);
		HalfEdge::HalfEdge* startHE = oldVertex->pHalfEdge;
		HalfEdge::HalfEdge* he = startHE;
		int valence; 
		float beta;

		bool onBoundary = oldVertex->onBoundary(); 
		if (onBoundary) {
			while (he->pPair != nullptr) he = he->pPair->pNext; 
			auto p1 = he->pNext->pStartVertex; 

			he = he->pPrev;
			while (he->pPair != nullptr)  he = he->pPair->pPrev; 
			auto p2 = he->pStartVertex; 

			sum = p1->position + p2->position; 
			valence = 2; 
			beta = 1.0f / 8.0f; 
		}
		else {
			valence = oldVertex->countValence();
			beta = valence == 3 ? 3.0f / 16.0f : (3.0f / (8.0f * valence)); 
			do {
				sum += he->pPair->pStartVertex->position;
				he = he->pPair->pNext;
			} while (he != startHE);
		}

		newVertex->position = (1.0f - valence * beta) * oldVertexPosition + beta * sum;
	}

	// Step 3: create new faces

	for (int fi = 0; fi < nOldFaces; ++fi)
	{
		auto oldFace = mesh.faces[fi];

		// TODO: update the half-edge data structure within each old face
		// HINT: the number of new faces within each old face is always 4 in the case of Loop subdivision,
		//       so you can write down all the steps without using a "for" or "while" loop 
		auto newPair = newHalfEdgeDict.find(oldFace->pHalfEdge->id)->second; 
		auto he1 = oldFace->pHalfEdge; 
		auto he2 = he1->pNext; 
		auto he3 = he2->pNext; 
		auto o1 = he1->pStartVertex; 
		auto o2 = he2->pStartVertex; 
		auto o3 = he3->pStartVertex; 
		auto m1 = newEdgeMidpointDict[he1->id]; 
		auto m2 = newEdgeMidpointDict[he2->id];
		auto m3 = newEdgeMidpointDict[he3->id];

		auto newPair1 = newHalfEdgeDict.find(he1->id)->second;
		auto o1_m1 = newPair1.first, m1_o2 = newPair1.second; 
		auto newPair2 = newHalfEdgeDict.find(he2->id)->second;
		auto o2_m2 = newPair2.first, m2_o3 = newPair2.second;
		auto newPair3 = newHalfEdgeDict.find(he3->id)->second;
		auto o3_m3 = newPair3.first, m3_o1 = newPair3.second;

		// 1
		HalfEdge::Face* newFace1 = newMesh.addFace();
		HalfEdge::HalfEdge* m1_m3 = newMesh.addHalfEdge(); 
		newFace1->pHalfEdge = o1_m1; 

		//o1_m1->pStartVertex = o1; 
		HalfEdge::Helper::SetPrevNext(o1_m1, m1_m3); 
		HalfEdge::Helper::SetPrevNext(m3_o1, o1_m1);
		o1_m1->pFace = newFace1; 

		m1_m3->pStartVertex = m1;
		HalfEdge::Helper::SetPrevNext(m1_m3, m3_o1);
		HalfEdge::Helper::SetPrevNext(o1_m1, m1_m3);
		m1_m3->pFace = newFace1;

		//m3_o1->pStartVertex = m3;
		HalfEdge::Helper::SetPrevNext(m3_o1, o1_m1);
		HalfEdge::Helper::SetPrevNext(m1_m3, m3_o1);
		m3_o1->pFace = newFace1;

		// 2
		HalfEdge::Face* newFace2 = newMesh.addFace();
		HalfEdge::HalfEdge* m2_m1 = newMesh.addHalfEdge();
		newFace2->pHalfEdge = o2_m2;

		//o2_m2->pStartVertex = o2;
		HalfEdge::Helper::SetPrevNext(o2_m2, m2_m1);
		HalfEdge::Helper::SetPrevNext(m1_o2, o2_m2);
		o2_m2->pFace = newFace2;

		m2_m1->pStartVertex = m2;
		HalfEdge::Helper::SetPrevNext(m2_m1, m1_o2);
		HalfEdge::Helper::SetPrevNext(o2_m2, m2_m1);
		m2_m1->pFace = newFace2;

		//m1_o2->pStartVertex = m1;
		HalfEdge::Helper::SetPrevNext(m1_o2, o2_m2);
		HalfEdge::Helper::SetPrevNext(m2_m1, m1_o2);
		m1_o2->pFace = newFace2;

		// 3
		HalfEdge::Face* newFace3 = newMesh.addFace();
		HalfEdge::HalfEdge* m3_m2 = newMesh.addHalfEdge();
		newFace3->pHalfEdge = o3_m3;

		//o3_m3->pStartVertex = o3;
		HalfEdge::Helper::SetPrevNext(o3_m3, m3_m2);
		HalfEdge::Helper::SetPrevNext(m2_o3, o3_m3);
		o3_m3->pFace = newFace3;

		m3_m2->pStartVertex = m3;
		HalfEdge::Helper::SetPrevNext(m3_m2, m2_o3);
		HalfEdge::Helper::SetPrevNext(o3_m3, m3_m2);
		m3_m2->pFace = newFace3;

		//m2_o3->pStartVertex = m2;
		HalfEdge::Helper::SetPrevNext(m2_o3, o3_m3);
		HalfEdge::Helper::SetPrevNext(m3_m2, m2_o3);
		m2_o3->pFace = newFace3;

		// 4
		HalfEdge::Face* newFace4 = newMesh.addFace();
		HalfEdge::HalfEdge* m1_m2 = newMesh.addHalfEdge();
		HalfEdge::HalfEdge* m2_m3 = newMesh.addHalfEdge();
		HalfEdge::HalfEdge* m3_m1 = newMesh.addHalfEdge();
		newFace4->pHalfEdge = m1_m2;

		m1_m2->pStartVertex = m1;
		HalfEdge::Helper::SetPrevNext(m1_m2, m2_m3);
		HalfEdge::Helper::SetPrevNext(m3_m1, m1_m2);
		m1_m2->pFace = newFace4;

		m2_m3->pStartVertex = m2;
		HalfEdge::Helper::SetPrevNext(m2_m3, m3_m1);
		HalfEdge::Helper::SetPrevNext(m1_m2, m2_m3);
		m2_m3->pFace = newFace4;

		m3_m1->pStartVertex = m3;
		HalfEdge::Helper::SetPrevNext(m3_m1, m1_m2);
		HalfEdge::Helper::SetPrevNext(m2_m3, m3_m1);
		m3_m1->pFace = newFace4;

		HalfEdge::Helper::SetPair(m1_m2, m2_m1); 
		HalfEdge::Helper::SetPair(m2_m3, m3_m2);
		HalfEdge::Helper::SetPair(m3_m1, m1_m3);
	}

	cerr << __FUNCTION__ << ": check data consistency" << endl;
	newMesh.checkDataConsistency();

	return move(newMesh);
}
