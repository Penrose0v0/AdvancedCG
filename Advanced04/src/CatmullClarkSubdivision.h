#pragma once

#include "AbstractSubdivision.h"

class CatmullClarkSubdivision : public AbstractSubdivision
{
public:
	static AbstractSubdivision* Create() { return new CatmullClarkSubdivision(); }

	void subdivide(PolygonMesh& mesh, int nSubdiv);

private:
	HalfEdge::Mesh apply(HalfEdge::Mesh& mesh);

};