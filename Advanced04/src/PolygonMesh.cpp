#include "PolygonMesh.h"
#define FAST_OBJ_IMPLEMENTATION
#include "fast_obj.h"
//#include "Image2OGLTexture.h"
#include <iostream>
#include <glm/ext.hpp>
#if defined(_WIN32) || defined(__APPLE__)
#include <algorithm>	// for std::replace in some environment
#else
#include <bits/stdc++.h>	// for std::replace with g++
#endif

using namespace std;

void PolygonMesh::triangulate()
{
	if (m_FaceIndices.empty())
	{
		cerr << __FUNCTION__ << ": no faces" << endl;
		return;
	}

	const int nFaces = (int)m_FaceIndices.size();

	vector<FaceIndices> newFaceIndices;
	newFaceIndices.reserve(2 * nFaces);

	for (int fi = 0; fi < nFaces; ++fi)
	{
		const int nFaceVertices = (int)m_FaceIndices[fi].indices.size();

		if (nFaceVertices == 3)
		{
			newFaceIndices.push_back(m_FaceIndices[fi]);
		}
		else
		{
			VertexTuple v0 = m_FaceIndices[fi].indices[0];
			VertexTuple v1 = m_FaceIndices[fi].indices[1];

			for (int vi = 2; vi < nFaceVertices; ++vi)
			{
				VertexTuple v2 = m_FaceIndices[fi].indices[vi];
				newFaceIndices.emplace_back(v0, v1, v2);
				v1 = v2;
			}
		}
	}

	newFaceIndices.shrink_to_fit();

	cerr << __FUNCTION__ << ": # faces " << m_FaceIndices.size() << " -> " << newFaceIndices.size() << endl;

	m_FaceIndices = move(newFaceIndices);
}

bool PolygonMesh::loadObj(const char* filename)
{
	fastObjMesh* m = fast_obj_read(filename);

	if (!m)
	{
		cerr << __FUNCTION__ << ": loading " << filename << " failed" << endl;
		return false;
	}

	//if (m->materials && m->materials->map_Kd.path)
	//{
	//	string diffuseTexPath(m->materials->map_Kd.path);
	//	replace(diffuseTexPath.begin(), diffuseTexPath.end(), '\\', '/');

	//	int w, h;
	//	Image2OGLTexture(diffuseTexPath.c_str(), m_TexID, w, h);
	//}

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

	//if (m->texcoord_count > 1)
	//{
	//	m_TexCoords.resize(m->texcoord_count - 1);
	//	memcpy(&m_TexCoords[0], m->texcoords + 2, 2 * sizeof(float) * (m->texcoord_count - 1));
	//}
	//else
	//{
	//	m_TexCoords.clear();
	//}

	m_FaceIndices.clear();
	m_FaceIndices.reserve(m->face_count);

	for (unsigned int gi = 0; gi < m->group_count; ++gi)
	{
		const fastObjGroup& grp = m->groups[gi];

		int nGroupFaceVerticesCount = 0;
		for (unsigned int groupFaceIndex = 0; groupFaceIndex < grp.face_count; ++groupFaceIndex)
		{
			unsigned int nFaceVertices = m->face_vertices[grp.face_offset + groupFaceIndex];

			FaceIndices t;

			for (int vi = 0; vi < nFaceVertices; ++vi)
			{
				fastObjIndex m0 = m->indices[grp.index_offset + nGroupFaceVerticesCount++];

				if (m0.t)
				{
					if (m0.n)
						t.indices.emplace_back(VertexTuple(m0.p - 1, m0.t - 1, m0.n - 1));
					else
						t.indices.emplace_back(VertexTuple(m0.p - 1, m0.t - 1, -1));
				}
				else
				{
					if (m0.n)
						t.indices.emplace_back(VertexTuple(m0.p - 1, -1, m0.n - 1));
					else
						t.indices.emplace_back(VertexTuple(m0.p - 1));
				}
			}

			m_FaceIndices.emplace_back(t);
		}
	}

	fast_obj_destroy(m);

	cout << __FUNCTION__ << ": " << filename << " loaded" << endl;
	cout << "  # verts:\t" << m_Vertices.size() << endl
		 << "  # normals:\t" << m_VertexNormals.size() << endl
		 //<< "  # tex coords:\t" << m_TexCoords.size() << endl
		 << "  # triangles:\t" << m_FaceIndices.size() << endl;

	if (m_VertexNormals.empty())
		calcVertexNormals();

	calcModelMatrix();

	return true;
}

//bool TriMesh::loadTexture(const char* filename)
//{
//	int w, h;
//	return Image2OGLTexture(filename, m_TexID, w, h);
//}

void PolygonMesh::calcVertexNormals()
{
	const int nFaces = (int)m_FaceIndices.size();
	const int nVerts = (int)m_Vertices.size();

	m_VertexNormals.resize(nVerts);

//#pragma omp parallel
	{
//#pragma omp for
		for (int ni = 0; ni < nVerts; ni++) m_VertexNormals[ni] = glm::vec3(0.0);

//#pragma omp for
		for (int fi = 0; fi < nFaces; fi++)
		{
			const auto& indices = m_FaceIndices[fi].indices;
			const int nFaceVertices = (int)indices.size();

			glm::vec3 weightedFaceNormal(0.f);
			for (int tj = 0; tj < nFaceVertices; ++tj)
			{
				glm::vec3 v0 = m_Vertices[indices[tj].vertexIdx];
				glm::vec3 v1 = m_Vertices[indices[(tj + 1) % nFaceVertices].vertexIdx];
				glm::vec3 v2 = m_Vertices[indices[(tj + 2) % nFaceVertices].vertexIdx];
				weightedFaceNormal += glm::cross(v1 - v0, v2 - v1);
			}

			for (int tj = 0; tj < nFaceVertices; ++tj)
				m_VertexNormals[indices[tj].vertexIdx] += weightedFaceNormal;
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

void PolygonMesh::renderMesh() const
{
	if (m_FaceIndices.empty())
		return;

	for (int fi = 0; fi < (int)m_FaceIndices.size(); ++fi)
	{
		glBegin(GL_POLYGON);
		for (const auto t : m_FaceIndices[fi].indices)
		{
			const int normalIdx = (t.normalIdx == -1) ? t.vertexIdx : t.normalIdx;
			glNormal3fv(glm::value_ptr(m_VertexNormals[t.normalIdx]));
			glVertex3fv(glm::value_ptr(m_Vertices[t.vertexIdx]));
		}
		glEnd();
	}
}

void PolygonMesh::renderMeshWithoutNormals() const
{
	if (m_FaceIndices.empty())
		return;

	for (int fi = 0; fi < (int)m_FaceIndices.size(); ++fi)
	{
		glBegin(GL_POLYGON);
		for (const auto t : m_FaceIndices[fi].indices)
			glVertex3fv(glm::value_ptr(m_Vertices[t.vertexIdx]));
		glEnd();
	}
}

void PolygonMesh::renderWireframeMesh() const
{
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1, 1);
	renderMeshWithoutNormals();
	glDisable(GL_POLYGON_OFFSET_FILL);
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	renderMeshWithoutNormals();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void PolygonMesh::calcModelMatrix()
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
	//m_ModelMatrix = glm::mat4(1.f);

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
