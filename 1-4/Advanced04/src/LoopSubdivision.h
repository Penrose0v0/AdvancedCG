#pragma once

#include "AbstractSubdivision.h"

class LoopSubdivision : public AbstractSubdivision
{
public:
	static AbstractSubdivision* Create() { return new LoopSubdivision(); }

	void subdivide(PolygonMesh& mesh, int nSubdiv);

private:
	HalfEdge::Mesh apply(HalfEdge::Mesh& mesh);
};