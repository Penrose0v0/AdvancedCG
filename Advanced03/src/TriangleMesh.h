#pragma once

#include <GL/glew.h>
#include "Triangle.h"
#include <array>
#include <vector>

struct VertexTuple
{
	VertexTuple() {}
	VertexTuple(int _vertexIdx) : vertexIdx(_vertexIdx), textureIdx(_vertexIdx), normalIdx(_vertexIdx) {}
	VertexTuple(int _vertexIdx, int _textureIdx, int _normalIdx) : vertexIdx(_vertexIdx), textureIdx(_textureIdx), normalIdx(_normalIdx) {}

	void set(int _vertexIdx) { vertexIdx = _vertexIdx;	textureIdx = _vertexIdx; normalIdx = _vertexIdx; }
	void set(int _vertexIdx, int _textureIdx, int _normalIdx) { vertexIdx = _vertexIdx;	textureIdx = _textureIdx; normalIdx = _normalIdx; }

	int vertexIdx, textureIdx, normalIdx;
};

struct TriangleIndices
{
	TriangleIndices() {}

	TriangleIndices(int _i0, int _i1, int _i2) : indices({ VertexTuple(_i0), VertexTuple(_i1), VertexTuple(_i2) }) {}
	TriangleIndices(const VertexTuple& _v0, const VertexTuple& _v1, const VertexTuple& _v2) : indices({ _v0, _v1, _v2 }) {}

	void set(int _i0, int _i1, int _i2) { indices = { VertexTuple(_i0), VertexTuple(_i1), VertexTuple(_i2) }; }
	void set(const VertexTuple& _v0, const VertexTuple& _v1, const VertexTuple& _v2) { indices = { _v0, _v1, _v2 }; }

	std::array<VertexTuple, 3> indices;
};

class TriangleMesh : public GeometricObject
{
private:
	TriangleMesh() : m_ShadingType(Auto_Shading), m_VBO(0)
	{
	}

public:
	enum Shading_Type
	{
		Flat_Shading,
		Smooth_Shading,	// Phong smooth shading
		Auto_Shading	// specified by the mesh data
	};

	static TriangleMesh *CreateGeometricObject()
	{
		TriangleMesh *obj = new TriangleMesh();
		//s_GeometricObjectCache.push_back( obj );
		GeometricObject::RegisterObject(obj);
		return obj;
	}

	~TriangleMesh()
	{
		if (m_VBO) glDeleteBuffers(1, &m_VBO);
	}

	// class definitions

	bool hit(const Ray &r, Real tmin, Real tmax, HitRecord &record) const;
	//bool shadowHit(const Ray &r, Real tmin, Real tmax) const;

	void drawGL() const;	// for preview using OpenGL

	void clear() { m_Triangles.clear(); }

	int getNumTriangles() const { return m_Triangles.size(); }

	// setter/getter

	void addTriangle(const Triangle &t) { m_Triangles.push_back(t); }
	void setTriangle(int i, const Triangle &t) { m_Triangles[i] = t; }

	const std::vector<Triangle> &getTriangles() const { return m_Triangles; }
	std::vector<Triangle> &getTriangles() { return m_Triangles; }

	const Triangle &getTriangle(int i) const { return m_Triangles[i]; }
	Triangle &getTriangle(int i) { return m_Triangles[i]; }

	void computeBoundingBox();

	void setBoundingBoxMin(const vec3 &p) { m_BoundingBoxPos[0] = vec3( p ); }
	void setBoundingBoxMax(const vec3 &p) { m_BoundingBoxPos[1] = vec3( p ); }

	vec3 getBoundingBoxMin() const { return m_BoundingBoxPos[0]; }
	vec3 getBoundingBoxMax() const { return m_BoundingBoxPos[1]; }

	bool loadObj(const char* filename);

	void setShadingType(Shading_Type type) { m_ShadingType = type; }
	Shading_Type getShadingType() const { return m_ShadingType; }

	void bakeVBO();

private:
	Shading_Type m_ShadingType;

	GLuint m_VBO;

	std::vector<Triangle> m_Triangles;

	vec3 m_BoundingBoxPos[2];

	//bool rayBoundingBoxIntersectionTest(const Ray &r, Real tmin, Real tmax) const;

	void calcVertexNormals(std::vector<glm::vec3>& normals, const std::vector<glm::vec3> &vertices, const std::vector<TriangleIndices> &indices);
	void scaleAndCenterize(std::vector<glm::vec3>& vertices, float scale = 1.f);
};

