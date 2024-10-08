#pragma once

#include <GL/glew.h>
#include <vector>
#include <array>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"

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

class TriMesh
{
public:
	TriMesh() : m_TexID(0), m_VertexVBO(0), m_VertexNormalVBO(0), m_TexCoordVBO(0) {}
	~TriMesh();

	const int getNumVertices() const { return (int)m_Vertices.size(); }
	const int getNumTriangles() const { return (int)m_TriangleIndices.size(); }

	void setVertices(const std::vector<glm::vec3>& vertices) { m_Vertices = vertices; }
	void setTriangles(const std::vector<TriangleIndices>& triangles) { m_TriangleIndices = triangles; }

	const std::vector<glm::vec3>& getVertices() const { return m_Vertices; }
	const std::vector<glm::vec3>& getVertexNormals() const { return m_VertexNormals; }
	const std::vector<glm::vec2>& getTexCoords() const { return m_TexCoords; }
	const std::vector<TriangleIndices>& getTriangleIndices() const { return m_TriangleIndices; }
	const glm::mat4& getModelMatrix() const { return m_ModelMatrix; }

	const GLuint getTexID() const { return m_TexID; }
	const GLuint getVertexVBO() const { return m_VertexVBO; }
	const GLuint getVertexNormalVBO() const { return m_VertexNormalVBO; }
	const GLuint getTexCoordVBO() const { return m_TexCoordVBO; }

	bool loadObj(const char* filename);
	bool loadTexture(const char* filename);

	void calcVertexNormals();

	void bakeVBOs();

	void renderTexturedMesh() const;
	void renderMeshGeometry() const;
	void renderWireframeMesh() const;

private:
	std::vector<glm::vec3> m_Vertices, m_VertexNormals;
	std::vector<glm::vec2> m_TexCoords;
	std::vector<TriangleIndices> m_TriangleIndices;

	glm::mat4 m_ModelMatrix;

	GLuint m_TexID, m_VertexVBO, m_VertexNormalVBO, m_TexCoordVBO;

	void calcModelMatrix();
};
