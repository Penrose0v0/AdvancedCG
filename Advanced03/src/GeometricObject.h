#pragma once

#include "Ray.h"
#include "Material.h"
#include "Texture.h"

#include <vector>

class GeometricObject
{
public:
	static void ClearGeometricObjectCache()
	{
		for (int i=0; i<(int)s_GeometricObjectCache.size(); i++)
			delete s_GeometricObjectCache[i];
		s_GeometricObjectCache.clear();
	}

	GeometricObject()
		: m_pMaterial(0)
	{
	}

	virtual ~GeometricObject()
	{
		// material deletion is managed in a static function of Material class
		// if (material) delete material;
	}

	virtual bool hit(const Ray &r, float tmin, float tmax, HitRecord &record) const = 0;

	virtual void drawGL() const = 0;	// for preview using OpenGL

	Material *getMaterial() const { return m_pMaterial; }
	void setMaterial(Material *m) { m_pMaterial = m; }

	virtual glm::vec3 getBoundingBoxMin() const { return glm::vec3(0,0,0); }	// to be implemented
	virtual glm::vec3 getBoundingBoxMax() const { return glm::vec3(0,0,0); }	// to be implemented

protected:
	typedef glm::vec3 vec3;
	typedef glm::vec2 vec2;
	typedef float Real;

	static std::vector<GeometricObject *> s_GeometricObjectCache;

	Material *m_pMaterial;

	static void RegisterObject(GeometricObject* obj);
};

