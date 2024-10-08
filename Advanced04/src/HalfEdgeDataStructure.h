#pragma once

#include <iostream>
#include <glm/glm.hpp>
#include "PolygonMesh.h"

namespace HalfEdge {
	struct HalfEdge;

	struct Vertex
	{
		Vertex() : pHalfEdge(nullptr), id(-1) {}
		Vertex(int _id) : pHalfEdge(nullptr), id(_id) {}

		int countValence() const;
		bool onBoundary() const;

		HalfEdge* pHalfEdge;
		glm::vec3 position;
		long id;
	};

	std::ostream& operator<< (std::ostream& stream, const Vertex &v);

	struct Face
	{
		Face() : pHalfEdge(nullptr), id(-1) {}
		Face(int _id) : pHalfEdge(nullptr), id(_id) {}

		glm::vec3 calcCentroidPosition() const;

		HalfEdge* pHalfEdge;
		long id;
	};

	std::ostream& operator<< (std::ostream& stream, const Face& f);

	struct HalfEdge
	{
		HalfEdge() : pStartVertex(nullptr), pPair(nullptr), pNext(nullptr), pPrev(nullptr), pFace(nullptr), id(-1) {}
		HalfEdge(int _id) : pStartVertex(nullptr), pPair(nullptr), pNext(nullptr), pPrev(nullptr), pFace(nullptr), id(_id) {}

		inline Vertex* getStartVertex() { return pStartVertex; }
		inline Vertex* getEndVertex() { return pNext->pStartVertex; }

		Vertex* pStartVertex;
		HalfEdge* pPair, *pNext, *pPrev;
		Face* pFace;
		long id;
	};

	std::ostream& operator<< (std::ostream& stream, const HalfEdge& he);

	struct Mesh
	{
		Mesh() {}
		~Mesh() { clear(); }
		Mesh(const Mesh& rhs) { deepCopy(rhs); }

		Mesh& operator=(const Mesh& rhs) { deepCopy(rhs); return *this; }

		void deepCopy(const Mesh& rhs);

		void clear() { clearVertices(); clearFaces(); clearHalfEdges(); }
		void clearVertices() { for (auto p : vertices) delete p; vertices.clear(); }
		void clearFaces() { for (auto p : faces) delete p; faces.clear(); }
		void clearHalfEdges() { for (auto p : halfEdges) delete p; halfEdges.clear(); }

		void build(const PolygonMesh& mesh, bool checkConsistency = false);
		void restore(PolygonMesh& mesh) const;

		void checkDataConsistency() const;

		Vertex* addVertex() { auto v = new Vertex(vertices.size()); vertices.push_back(v); return v; }
		Face* addFace() { auto f = new Face(faces.size()); faces.push_back(f); return f; }
		HalfEdge* addHalfEdge() { auto he = new HalfEdge(halfEdges.size()); halfEdges.push_back(he); return he; }

		std::vector<Vertex*> vertices;
		std::vector<Face*> faces;
		std::vector<HalfEdge*> halfEdges;
	};

	namespace Helper
	{
		void SetPair(HalfEdge* he1, HalfEdge* he2);
		void SetPrevNext(HalfEdge* prevHE, HalfEdge* nextHE);
	}
}