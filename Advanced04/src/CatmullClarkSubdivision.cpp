#include "CatmullClarkSubdivision.h"
#include <unordered_map>

using namespace std;
using namespace glm;

void CatmullClarkSubdivision::subdivide(PolygonMesh& mesh, int nSubdiv)
{
	if (mesh.getVertices().empty() || mesh.getFaceIndices().empty())
	{
		std::cerr << __FUNCTION__ << ": mesh not ready" << std::endl;
		return;
	}

	HalfEdge::Mesh _mesh;
	_mesh.build(mesh);

	for (int iter = 0; iter < nSubdiv; ++iter)
		_mesh = apply(_mesh);

	_mesh.restore(mesh);
	mesh.calcVertexNormals();
}

HalfEdge::Mesh CatmullClarkSubdivision::apply(HalfEdge::Mesh& mesh)
{
	// return mesh;	// TODO: delete this line

	HalfEdge::Mesh newMesh;

	const int nOldVertices = (int)mesh.vertices.size();
	const int nOldFaces = (int)mesh.faces.size();
	const int nOldHalfEdges = (int)mesh.halfEdges.size();

	// Step 0: allocate memory for even (i.e., old) vertices

	newMesh.vertices.reserve(nOldFaces + nOldVertices + nOldHalfEdges);

	for (int vi = 0; vi < nOldVertices; ++vi)
		newMesh.addVertex();

	// Step 1: generate face centroids

	for (int fi = 0; fi < nOldFaces; ++fi)
	{
		auto oldFace = mesh.faces[fi];
		auto newFaceCentroid = newMesh.addVertex();
		// TODO: calculate the positions of face centroids
		newFaceCentroid->position = oldFace->calcCentroidPosition(); 
	}

	// Step 2: create odd (i.e., new) vertices by splitting half edges

	unordered_map<long, HalfEdge::Vertex*> newEdgeMidpointDict;
	unordered_map<long, pair<HalfEdge::HalfEdge*, HalfEdge::HalfEdge*>> newHalfEdgeDict;

	for (int hi = 0; hi < nOldHalfEdges; ++hi)
	{
		auto oldHE = mesh.halfEdges[hi];

		vec3 startVertexPosition = oldHE->getStartVertex()->position;
		vec3 endVertexPosition = oldHE->getEndVertex()->position;

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
				vec3 oppositeStart = oldHE->pPair->pFace->calcCentroidPosition();
				vec3 oppositeEnd = oldHE->pFace->calcCentroidPosition();
				edgeMidpoint->position = (startVertexPosition + endVertexPosition + oppositeStart + oppositeEnd) / 4.0f;
			}
			else
			{
				edgeMidpoint = edgeMidpointIter->second;
			}
		}

		newEdgeMidpointDict[oldHE->id] = edgeMidpoint;

		auto formerHE = newMesh.addHalfEdge();
		auto latterHE = newMesh.addHalfEdge();

		auto evenStartVertex = newMesh.vertices[oldHE->pStartVertex->id];
		auto evenEndVertex = newMesh.vertices[oldHE->pNext->pStartVertex->id];

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

	// Step 3: update even (i.e., old) vertex positions

	for (int vi = 0; vi < nOldVertices; ++vi)
	{
		auto newVertex = newMesh.vertices[vi];
		const auto oldVertex = mesh.vertices[vi];
		const auto oldVertexPosition = oldVertex->position;

		// TODO: calculate the new vertex position
		// c.f., HalfEdge::Vertex::countValence() in HalfEdgeDataStructure.cpp

		bool onBoundary = oldVertex->onBoundary();

		int n = oldVertex->countValence();

		vec3 F(0.0f);
		vec3 R(0.0f);

		HalfEdge::HalfEdge* he = oldVertex->pHalfEdge;

		if (onBoundary) {
			while (he->pPair != nullptr) he = he->pPair->pNext;
			auto p1 = he->pNext->pStartVertex;

			he = he->pPrev;
			while (he->pPair != nullptr)  he = he->pPair->pPrev;
			auto p2 = he->pStartVertex;

			newVertex->position = (3.0f / 4.0f) * oldVertex->position + (1.0f / 8.0f) * (p1->position + p2->position); 
			continue; 
		}

		do
		{
			if (he->pFace != nullptr)
			{
				F += he->pFace->calcCentroidPosition();
			}
			he = he->pPair->pNext;
		} while (he != oldVertex->pHalfEdge && he->pPair != nullptr);

		F /= static_cast<float>(n);

		he = oldVertex->pHalfEdge;
		do
		{
			if (he->pPair != nullptr)
			{
				vec3 midpoint = 0.5f * (he->getStartVertex()->position + he->getEndVertex()->position);
				R += midpoint;
			}
			he = he->pPair->pNext;
		} while (he != oldVertex->pHalfEdge && he->pPair != nullptr);

		R /= static_cast<float>(n);

		newVertex->position = (F + 2.0f * R + (n - 3.0f) * oldVertexPosition) / static_cast<float>(n);
	}

	// Step 4: set up new faces

	for (int fi = 0; fi < nOldFaces; ++fi)
	{
		auto oldFace = mesh.faces[fi];
		auto centroidVertex = newMesh.vertices[oldFace->id + nOldVertices];

		// TODO: update the half-edge data structure within each old face
		// HINT: use the following std::vector to store temporal data and process step by step
		//vector<HalfEdge::HalfEdge*> tmpToCentroidHalfEdges;
		//vector<HalfEdge::Face*> tmpNewFaces;
		vector<HalfEdge::HalfEdge*> faceHalfEdges;

		// 元の面のハーフエッジを巡回し、新しい四角形を作成
		HalfEdge::HalfEdge* cur_he = oldFace->pHalfEdge;
		HalfEdge::HalfEdge* prev_he = cur_he->pPrev;
		do {
			// previous
			auto prevPair = newHalfEdgeDict.find(prev_he->id)->second; 
			auto m0_o1 = prevPair.second;

			// current
			auto curPair = newHalfEdgeDict.find(cur_he->id)->second;
			auto o1_m1 = curPair.first;
			auto m1_o2 = curPair.second;

			// this loop
			HalfEdge::Face* newFace = newMesh.addFace();
			HalfEdge::HalfEdge* m1_c = newMesh.addHalfEdge();
			HalfEdge::HalfEdge* c_m1 = newMesh.addHalfEdge();
			newFace->pHalfEdge = o1_m1;

			HalfEdge::Helper::SetPrevNext(m0_o1, o1_m1);
			HalfEdge::Helper::SetPrevNext(o1_m1, m1_c);
			HalfEdge::Helper::SetPair(m1_c, c_m1);
			HalfEdge::Helper::SetPrevNext(c_m1, m1_o2);

			m0_o1->pFace = newFace;
			o1_m1->pFace = newFace;
			m1_c->pFace = newFace; 
			m1_c->pStartVertex = newEdgeMidpointDict[cur_he->id];
			c_m1->pStartVertex = centroidVertex; 

			// c_m0 should be set in another loop
			// HalfEdge::Helper::SetPrevNext(m1_c, m0_o1->pPrev);

			cur_he = cur_he->pNext; 
			prev_he = cur_he->pPrev; 
		} while (cur_he != oldFace->pHalfEdge); 

		do {
			auto curPair = newHalfEdgeDict.find(cur_he->id)->second;
			auto o1_m1 = curPair.first;
			auto c_m0 = o1_m1->pPrev->pPrev; 

			c_m0->pFace = o1_m1->pFace; 
			HalfEdge::Helper::SetPrevNext(o1_m1->pNext, c_m0);

			cur_he = cur_he->pNext;
		} while (cur_he != oldFace->pHalfEdge);

		auto curPair = newHalfEdgeDict.find(cur_he->id)->second;
		auto o1_m1 = curPair.first;
		centroidVertex->pHalfEdge = o1_m1->pNext->pPair; 
	}

	cerr << __FUNCTION__ << ": check data consistency" << endl;
	newMesh.checkDataConsistency();

	return move(newMesh);
}
