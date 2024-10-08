#pragma once

#include <GL/glew.h>
#include <vector>
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

struct FaceIndices
{
	FaceIndices() {}

	FaceIndices(int _i0, int _i1, int _i2) : indices({ VertexTuple(_i0), VertexTuple(_i1), VertexTuple(_i2) }) {}
	FaceIndices(const VertexTuple& _v0, const VertexTuple& _v1, const VertexTuple& _v2) : indices({ _v0, _v1, _v2 }) {}

	void set(int _i0, int _i1, int _i2) { indices = { VertexTuple(_i0), VertexTuple(_i1), VertexTuple(_i2) }; }
	void set(const VertexTuple& _v0, const VertexTuple& _v1, const VertexTuple& _v2) { indices = { _v0, _v1, _v2 }; }

	std::vector<VertexTuple> indices;
};

class PolygonMesh
{
public:
	PolygonMesh() {}
	~PolygonMesh() {}

	const int getNumVertices() const { return (int)m_Vertices.size(); }
	const int getNumFaces() const { return (int)m_FaceIndices.size(); }

	void setVertices(const std::vector<glm::vec3>& vertices) { m_Vertices = vertices; }
	void setFaces(const std::vector<FaceIndices>& faceIndices) { m_FaceIndices = faceIndices; }

	const std::vector<glm::vec3>& getVertices() const { return m_Vertices; }
	const std::vector<glm::vec3>& getVertexNormals() const { return m_VertexNormals; }
	const std::vector<FaceIndices>& getFaceIndices() const { return m_FaceIndices; }

	void triangulate();

	bool loadObj(const char* filename);

	void calcVertexNormals();

	void renderMesh() const;
	void renderMeshWithoutNormals() const;
	void renderWireframeMesh() const;

private:
	std::vector<glm::vec3> m_Vertices, m_VertexNormals;
	std::vector<FaceIndices> m_FaceIndices;

	void calcModelMatrix();
};
