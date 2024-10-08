#pragma once

#include "GLSLProgramObject.h"
#include "PolygonMesh.h"

class AbstractRenderer
{
public:
	virtual ~AbstractRenderer() {}
	virtual void render(const PolygonMesh& mesh) const = 0;
};