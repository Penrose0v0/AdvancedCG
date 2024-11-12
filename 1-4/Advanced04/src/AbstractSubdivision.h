#pragma once

#include "PolygonMesh.h"
#include "HalfEdgeDataStructure.h"
#include <iostream>

class AbstractSubdivision
{
public:
	virtual ~AbstractSubdivision() {}
	virtual void subdivide(PolygonMesh& mesh, int nSubdiv) = 0;

protected:
	virtual HalfEdge::Mesh apply(HalfEdge::Mesh &mesh) = 0;
};
