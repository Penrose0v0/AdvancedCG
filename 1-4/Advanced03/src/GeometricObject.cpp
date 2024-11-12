#include "GeometricObject.h"


using namespace std;

vector<GeometricObject *> GeometricObject::s_GeometricObjectCache;

void GeometricObject::RegisterObject(GeometricObject* obj)
{
	//s_PseudoColorCache.emplace_back(0.5f * s_RandDist(s_RandSrc) + 0.5f,
	//								0.5f * s_RandDist(s_RandSrc) + 0.5f,
	//								0.5f * s_RandDist(s_RandSrc) + 0.5f);
	s_GeometricObjectCache.emplace_back(obj);
}
