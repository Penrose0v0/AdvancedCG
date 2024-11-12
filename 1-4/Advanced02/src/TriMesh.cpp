#include "TriMesh.h"
#define FAST_OBJ_IMPLEMENTATION
#include "fast_obj.h"
#include "Image2OGLTexture.h"
#include <iostream>
#include <glm/ext.hpp>
#if defined(_WIN32) || defined(__APPLE__)
#include <algorithm>	// for std::replace in some environment
#else
#include <bits/stdc++.h>	// for std::replace with g++
#endif

using namespace std;

TriMesh::~TriMesh()
{
	if (m_TexID) glDeleteTextures(1, &m_TexID);
	if (m_VertexVBO) glDeleteBuffers(1, &m_VertexVBO);
	if (m_VertexNormalVBO) glDeleteBuffers(1, &m_VertexNormalVBO);
	if (m_TexCoordVBO) glDeleteBuffers(1, &m_TexCoordVBO);
}

bool TriMesh::loadObj(const char* filename)
{
	fastObjMesh* m = fast_obj_read(filename);

	if (!m)
	{
		cerr << __FUNCTION__ << ": loading " << filename << " failed" << endl;
		return false;
	}

	if (m->materials && m->materials->map_Kd.path)
	{
		string diffuseTexPath(m->materials->map_Kd.path);
		replace(diffuseTexPath.begin(), diffuseTexPath.end(), '\\', '/');

		int w, h;
		Image2OGLTexture(diffuseTexPath.c_str(), m_TexID, w, h);
	}

	if (m->position_count > 1)
	{
		m_Vertices.resize(m->position_count - 1);
		memcpy(&m_Vertices[0], m->positions + 3, 3 * sizeof(float) * (m->position_count - 1));
	}
	else
	{
		cerr << __FUNCTION__ << ": warning: no vertices defined" << endl;
		return false;
	}

	if (m->normal_count > 1)
	{
		m_VertexNormals.resize(m->normal_count - 1);
		memcpy(&m_VertexNormals[0], m->normals + 3, 3 * sizeof(float) * (m->normal_count - 1));
	}
	else
	{
		m_VertexNormals.clear();
	}

	if (m->texcoord_count > 1)
	{
		m_TexCoords.resize(m->texcoord_count - 1);
		memcpy(&m_TexCoords[0], m->texcoords + 2, 2 * sizeof(float) * (m->texcoord_count - 1));
	}
	else
	{
		m_TexCoords.clear();
	}

	m_TriangleIndices.clear();
	m_TriangleIndices.reserve(m->face_count);

	int nNonTriangles = 0;

	for (unsigned int gi = 0; gi < m->group_count; ++gi)
	{
		const fastObjGroup& grp = m->groups[gi];

		int nGroupFaceVerticesCount = 0;
		for (unsigned int groupFaceIndex = 0; groupFaceIndex < grp.face_count; ++groupFaceIndex)
		{
			unsigned int nFaceVertices = m->face_vertices[grp.face_offset + groupFaceIndex];

			nNonTriangles += (nFaceVertices != 3);

			fastObjIndex m0 = m->indices[grp.index_offset + nGroupFaceVerticesCount++];
			fastObjIndex m1 = m->indices[grp.index_offset + nGroupFaceVerticesCount++];

			for (int vi = 2; vi < nFaceVertices; ++vi)	// defining triangles
			{
				fastObjIndex m2 = m->indices[grp.index_offset + nGroupFaceVerticesCount++];

				if (m0.t)
				{
					if (m0.n)
						m_TriangleIndices.emplace_back(VertexTuple(m0.p - 1, m0.t - 1, m0.n - 1), VertexTuple(m1.p - 1, m1.t - 1, m1.n - 1), VertexTuple(m2.p - 1, m2.t - 1, m2.n - 1));
					else
						m_TriangleIndices.emplace_back(VertexTuple(m0.p - 1, m0.t - 1, -1), VertexTuple(m1.p - 1, m1.t - 1, -1), VertexTuple(m2.p - 1, m2.t - 1, -1));
				}
				else
				{
					if (m0.n)
						m_TriangleIndices.emplace_back(VertexTuple(m0.p - 1, -1, m0.n - 1), VertexTuple(m1.p - 1, -1, m1.n - 1), VertexTuple(m2.p - 1, -1, m2.n - 1));
					else
						m_TriangleIndices.emplace_back(VertexTuple(m0.p - 1), VertexTuple(m1.p - 1), VertexTuple(m2.p - 1));
				}

				m1 = m2;
			}
		}
	}

	fast_obj_destroy(m);

	if (nNonTriangles)
		cerr << __FUNCTION__ << ": Note: " << nNonTriangles << " non-triangles found (automatically converted to triangles)" << endl;

	cout << __FUNCTION__ << ": " << filename << " loaded" << endl;
	cout << "  # verts:\t" << m_Vertices.size() << endl
		 << "  # normals:\t" << m_VertexNormals.size() << endl
		 << "  # tex coords:\t" << m_TexCoords.size() << endl
		 << "  # triangles:\t" << m_TriangleIndices.size() << endl;

	if (m_VertexNormals.empty())
		calcVertexNormals();

	calcModelMatrix();
	bakeVBOs();

	return true;
}

bool TriMesh::loadTexture(const char* filename)
{
	int w, h;
	return Image2OGLTexture(filename, m_TexID, w, h);
}

void TriMesh::calcVertexNormals()
{
	const int nFaces = (int)m_TriangleIndices.size();
	const int nVerts = (int)m_Vertices.size();

	m_VertexNormals.resize(nVerts);

//#pragma omp parallel
	{
//#pragma omp for
		for (int ni = 0; ni < nVerts; ni++) m_VertexNormals[ni] = glm::vec3(0.0);

//#pragma omp for
		for (int fi = 0; fi < nFaces; fi++)
		{
			const int i0 = m_TriangleIndices[fi].indices[0].vertexIdx;
			const int i1 = m_TriangleIndices[fi].indices[1].vertexIdx;
			const int i2 = m_TriangleIndices[fi].indices[2].vertexIdx;

			if (i0 >= nVerts || i1 >= nVerts || i2 >= nVerts || i0 < 0 || i1 < 0 || i2 < 0)
			{
				cerr << __FUNCTION__ << ": warning: triangle index out of range: (" << i0 << "," << i1 << "," << i2 << ")" << endl;
			}

			const glm::vec3 weightedNormal = glm::cross((m_Vertices[i1] - m_Vertices[i0]), (m_Vertices[i2] - m_Vertices[i1]));

			m_VertexNormals[i0] += weightedNormal;
			m_VertexNormals[i1] += weightedNormal;
			m_VertexNormals[i2] += weightedNormal;
		}

//#pragma omp for
		for (int ni = 0; ni < nVerts; ni++)
		{
			const float norm = glm::length(m_VertexNormals[ni]);
			if (norm > 0.f)
				m_VertexNormals[ni] = m_VertexNormals[ni] / norm;
		}
	}
}

void TriMesh::renderTexturedMesh() const
{
	if (!m_TexID || !m_VertexVBO || !m_TexCoordVBO)
	{
		cerr << __FUNCTION__ << ": texID = " << m_TexID << ", vertex VBO = " << m_VertexVBO << ", tex coord VBO = " << m_TexCoordVBO << endl;
		return;
	}

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, m_TexID);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, m_VertexVBO);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, m_VertexNormalVBO);
	glNormalPointer(GL_FLOAT, 0, 0);

	glActiveTexture(GL_TEXTURE0);
	glBindBuffer(GL_ARRAY_BUFFER, m_TexCoordVBO);
	glTexCoordPointer(2, GL_FLOAT, 0, 0);

	glDrawArrays(GL_TRIANGLES, 0, (GLsizei)m_TriangleIndices.size() * 3);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);

	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);
}

void TriMesh::renderMeshGeometry() const
{
	if (!m_VertexVBO || !m_VertexNormalVBO)
	{
		//cerr << __FUNCTION__ << ": vertex or normal VBO not ready" << endl;
		return;
	}

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, m_VertexVBO);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, m_VertexNormalVBO);
	glNormalPointer(GL_FLOAT, 0, 0);

	glDrawArrays(GL_TRIANGLES, 0, (GLsizei)m_TriangleIndices.size() * 3);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

void TriMesh::renderWireframeMesh() const
{
	if (!m_VertexVBO || !m_VertexNormalVBO)
	{
		//cerr << __FUNCTION__ << ": vertex or normal VBO not ready" << endl;
		return;
	}

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1, 1);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	renderMeshGeometry();

	glDisable(GL_POLYGON_OFFSET_FILL);
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

	glColor3f(0.4f, 0.4f, 0.4f);
	glLineWidth(0.1f);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	renderMeshGeometry();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void TriMesh::calcModelMatrix()
{
	glm::vec3 bboxMin(1000000.f), bboxMax(-1000000.f);

	for (int i = 0; i < (int)m_Vertices.size(); ++i)
	{
		const auto v = m_Vertices[i];
		if (bboxMin.x > v.x) bboxMin.x = v.x;
		if (bboxMin.y > v.y) bboxMin.y = v.y;
		if (bboxMin.z > v.z) bboxMin.z = v.z;
		if (bboxMax.x < v.x) bboxMax.x = v.x;
		if (bboxMax.y < v.y) bboxMax.y = v.y;
		if (bboxMax.z < v.z) bboxMax.z = v.z;
	}

	glm::vec3 middle = 0.5f * (bboxMin + bboxMax);
	const float radius = glm::distance(bboxMax, middle);

	const float scale = 1.f / radius;
	const glm::vec3 trans = -scale * middle;

#if 0
	m_ModelMatrix[0][0] = m_ModelMatrix[1][1] = m_ModelMatrix[2][2] = scale;
	m_ModelMatrix[1][0] = m_ModelMatrix[2][0] = 0.f;
	m_ModelMatrix[0][1] = m_ModelMatrix[2][1] = 0.f;
	m_ModelMatrix[0][2] = m_ModelMatrix[1][2] = 0.f;
	m_ModelMatrix[0][3] = m_ModelMatrix[1][3] = m_ModelMatrix[2][3] = 0.f;
	m_ModelMatrix[3][0] = trans.x;
	m_ModelMatrix[3][1] = trans.y;
	m_ModelMatrix[3][2] = trans.z;
	m_ModelMatrix[3][3] = 1.f;
#else
	m_ModelMatrix = glm::mat4(1.f);

	glm::vec3 _bboxMin(1000000.f), _bboxMax(-1000000.f);

	for (int i = 0; i < (int)m_Vertices.size(); ++i)
	{
		m_Vertices[i] = scale * m_Vertices[i] + trans;

		const auto v = m_Vertices[i];
		if (_bboxMin.x > v.x) _bboxMin.x = v.x;
		if (_bboxMin.y > v.y) _bboxMin.y = v.y;
		if (_bboxMin.z > v.z) _bboxMin.z = v.z;
		if (_bboxMax.x < v.x) _bboxMax.x = v.x;
		if (_bboxMax.y < v.y) _bboxMax.y = v.y;
		if (_bboxMax.z < v.z) _bboxMax.z = v.z;
	}

	cerr << __FUNCTION__ << ": new bounding box: " << glm::to_string(_bboxMin) << "x" << glm::to_string(_bboxMax) << endl;
#endif
}

void TriMesh::bakeVBOs()
{
	const int nTriangles = (int)m_TriangleIndices.size();
	const int nVertices = (int)m_Vertices.size();
	const int nNormals = (int)m_VertexNormals.size();
	const int nTexCoords = (int)m_TexCoords.size();

	vector<glm::vec3> tmpVertices, tmpNormals;
	vector<glm::vec2> tmpTexCoords;

	tmpVertices.resize(3 * nTriangles);
	if (nNormals) tmpNormals.resize(3 * nTriangles);
	if (nTexCoords) tmpTexCoords.resize(3 * nTriangles);

//#pragma omp parallel for
	for (int ti = 0; ti < nTriangles; ++ti)
		for (int vi = 0; vi < 3; ++vi)
		{
			const auto& t = m_TriangleIndices[ti].indices[vi];

			tmpVertices[3 * ti + vi] = m_Vertices[t.vertexIdx];

			if (nNormals)
				tmpNormals[3 * ti + vi] = m_VertexNormals[t.normalIdx >= 0 ? t.normalIdx : t.vertexIdx];

			if (nTexCoords)
				tmpTexCoords[3 * ti + vi] = m_TexCoords[t.textureIdx >= 0 ? t.textureIdx : t.vertexIdx];
		}

	if (!m_VertexVBO) glGenBuffers(1, &m_VertexVBO);
	glBindBuffer(GL_ARRAY_BUFFER, m_VertexVBO);
	glBufferData(GL_ARRAY_BUFFER, tmpVertices.size() * sizeof(glm::vec3), &tmpVertices[0], GL_STATIC_DRAW);

	if (nNormals)
	{
		if (!m_VertexNormalVBO) glGenBuffers(1, &m_VertexNormalVBO);
		glBindBuffer(GL_ARRAY_BUFFER, m_VertexNormalVBO);
		glBufferData(GL_ARRAY_BUFFER, tmpNormals.size() * sizeof(glm::vec3), &tmpNormals[0], GL_STATIC_DRAW);
	}
	else
	{
		cerr << __FUNCTION__ << ": warning: vertex normals not ready" << endl;
		glDeleteBuffers(1, &m_VertexNormalVBO);
		m_VertexNormalVBO = 0;
	}

	if (nTexCoords)
	{
		if (!m_TexCoordVBO) glGenBuffers(1, &m_TexCoordVBO);
		glBindBuffer(GL_ARRAY_BUFFER, m_TexCoordVBO);
		glBufferData(GL_ARRAY_BUFFER, tmpTexCoords.size() * sizeof(glm::vec2), &tmpTexCoords[0], GL_STATIC_DRAW);
	}
	else if (m_TexCoordVBO)
	{
		glDeleteBuffers(1, &m_TexCoordVBO);
		m_TexCoordVBO = 0;
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);
}
