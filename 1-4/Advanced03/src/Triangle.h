#pragma once

#include "GeometricObject.h"
#include <algorithm>

class Triangle : public GeometricObject
{
public:
	static Triangle *CreateGeometricObject()
	{
		Triangle *obj = new Triangle();
		s_GeometricObjectCache.push_back( obj );
		return obj;
	}

	static Triangle *CreateGeometricObject(const Triangle &tri)
	{
		Triangle *obj = new Triangle(tri);
		s_GeometricObjectCache.push_back( obj );
		return obj;
	}

	static Triangle *CreateGeometricObject(const vec3 &_v0, const vec3 &_v1, const vec3 &_v2, Material *m)
	{
		Triangle *obj = new Triangle(_v0, _v1, _v2, m);
		s_GeometricObjectCache.push_back( obj );
		return obj;
	}

	// TriangleMesh クラスで vector に格納するので、コンストラクタを public にする必要あり 

	Triangle()
	{
		m_TexCoords[0] = vec2( 0.f );
		m_TexCoords[1] = vec2( 0.f );
	}

	Triangle(const Triangle &tri)
	{
		copy(tri);
	}

	Triangle(const vec3 &_v0, const vec3 &_v1, const vec3 &_v2, Material *m)
	{
		m_Vertices[0] = _v0;
		m_Vertices[1] = _v1;
		m_Vertices[2] = _v2;

		m_Normals[0] = m_Normals[1] = m_Normals[2] = getFaceNormal();

		m_TexCoords[0] = vec2( 0.f );
		m_TexCoords[1] = vec2( 0.f );

		m_pMaterial = m;
	}

	void copy(const Triangle& tri)
	{
		m_Vertices[0] = tri.m_Vertices[0];
		m_Vertices[1] = tri.m_Vertices[1];
		m_Vertices[2] = tri.m_Vertices[2];

		m_Normals[0] = tri.m_Normals[0];
		m_Normals[1] = tri.m_Normals[1];
		m_Normals[2] = tri.m_Normals[2];

		m_TexCoords[0] = tri.m_TexCoords[0];
		m_TexCoords[1] = tri.m_TexCoords[1];

		m_pMaterial = tri.m_pMaterial;
	}

	Triangle& operator=(const Triangle& tri)
	{
		copy(tri);
		return *this;
	}

	bool hit(const Ray &r, Real tmin, Real tmax, HitRecord &record) const;
	//bool shadowHit(const Ray &r, Real tmin, Real tmax) const;

	void drawGL() const;	// for preview using OpenGL

	vec3 getFaceNormal() const
	{ 
		//vec3 n = (m_Vertices[1] - m_Vertices[0]) ^ (m_Vertices[2] - m_Vertices[0]);
		//return glm::normalize(n);
		return glm::normalize(glm::cross(m_Vertices[1] - m_Vertices[0], m_Vertices[2] - m_Vertices[0]));
	}

	inline vec3 getVertex0() const { return m_Vertices[0]; }
	inline vec3 getVertex1() const { return m_Vertices[1]; }
	inline vec3 getVertex2() const { return m_Vertices[2]; }

	inline void setVertex0(const vec3 &_v0) { m_Vertices[0] = _v0; }
	inline void setVertex1(const vec3 &_v1) { m_Vertices[1] = _v1; }
	inline void setVertex2(const vec3 &_v2) { m_Vertices[2] = _v2; }

	inline vec3 getNormal0() const { return m_Normals[0]; }
	inline vec3 getNormal1() const { return m_Normals[1]; }
	inline vec3 getNormal2() const { return m_Normals[2]; }

	inline void setNormal0(const vec3 &_n0) { m_Normals[0] = _n0; }
	inline void setNormal1(const vec3 &_n1) { m_Normals[1] = _n1; }
	inline void setNormal2(const vec3 &_n2) { m_Normals[2] = _n2; }

	inline vec2 getTexcoord0() const { return m_TexCoords[0]; }
	inline vec2 getTexcoord1() const { return m_TexCoords[1]; }
	inline vec2 getTexcoord2() const { return m_TexCoords[2]; }

	inline void setTexCoord0(const vec2 &_t0) { m_TexCoords[0] = _t0; }
	inline void setTexcoord1(const vec2 &_t1) { m_TexCoords[1] = _t1; }
	inline void setTexcoord2(const vec2 &_t2) { m_TexCoords[2] = _t2; }

	vec3 getBoundingBoxMin() const
	{
		const float x = std::min( m_Vertices[0].x, std::min(m_Vertices[1].x, m_Vertices[2].x));
		const float y = std::min( m_Vertices[0].y, std::min(m_Vertices[1].y, m_Vertices[2].y));
		const float z = std::min( m_Vertices[0].z, std::min(m_Vertices[1].z, m_Vertices[2].z));
		return vec3(x,y,z);
	}
	
	vec3 getBoundingBoxMax() const
	{
		const float x = std::max( m_Vertices[0].x, std::max(m_Vertices[1].x, m_Vertices[2].x));
		const float y = std::max( m_Vertices[0].y, std::max(m_Vertices[1].y, m_Vertices[2].y));
		const float z = std::max( m_Vertices[0].z, std::max(m_Vertices[1].z, m_Vertices[2].z));
		return vec3(x,y,z);
	}

private:
	vec3 m_Vertices[3];
	vec3 m_Normals[3];
	vec2 m_TexCoords[3];

};

