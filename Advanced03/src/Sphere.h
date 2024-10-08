#pragma once

#include <GL/glew.h>
#include "GeometricObject.h"

class Sphere : public GeometricObject
{
private:
	Sphere()
	{
	}

	Sphere(const Sphere &s)
		: m_Center(s.m_Center), m_Radius(s.m_Radius)
	{
		m_pMaterial = s.m_pMaterial;
	}

	Sphere(const vec3 &_center, Real r, Material *m)
		: m_Center(_center), m_Radius(r)
	{
		m_pMaterial = m;
	}

public:
	static Sphere *CreateGeometricObject()
	{
		Sphere *obj = new Sphere();
		GeometricObject::RegisterObject(obj);
		//s_GeometricObjectCache.push_back( obj );
		return obj;
	}

	static Sphere *CreateGeometricObject(const Sphere &s)
	{
		Sphere *obj = new Sphere(s);
		GeometricObject::RegisterObject(obj);
		//s_GeometricObjectCache.push_back( obj );
		return obj;
	}

	static Sphere *CreateGeometricObject(const vec3 &_center, Real r, Material *m)
	{
		Sphere *obj = new Sphere(_center, r, m);
		GeometricObject::RegisterObject(obj);
		//s_GeometricObjectCache.push_back( obj );
		return obj;
	}

	bool hit(const Ray &r, Real tmin, Real tmax, HitRecord &record) const;
	//bool shadowHit(const Ray &r, Real tmin, Real tmax) const;

	void drawGL() const;	// for preview using OpenGL

	inline vec3 getCenter() const { return m_Center; }
	inline Real getRadius() const { return m_Radius; }

	inline void setCenter(const vec3 &c) { m_Center = c; }
	inline void setRadius(Real r) { m_Radius = r; }

private:
	vec3 m_Center;
	Real m_Radius;

	static int s_NumSphereVertices;
	static GLuint s_VBO;

	static void BakeVBO();
};

