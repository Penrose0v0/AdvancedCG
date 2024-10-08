#include "TriangleMesh.h"
#include <cstdlib>
#include <array>
#define FAST_OBJ_IMPLEMENTATION
#include "fast_obj.h"
//#include <GL/glut.h>

using namespace std;
//using namespace MyAlgebra;

// public methods

bool TriangleMesh::hit(const Ray &r, Real tmin, Real tmax, HitRecord &record) const
{
	bool hitTriangle = false;
	float tTmp = tmax;

	for (int i=0; i<(int)m_Triangles.size(); i++)
	{
		HitRecord tmpRec;
		const bool isHit = m_Triangles[i].hit(r, tmin, tmax, tmpRec);

		if (isHit && tTmp > tmpRec.m_ParamT)
		{
			tTmp = tmpRec.m_ParamT;
			record = tmpRec;
			hitTriangle = true;
		}
	}

	if (hitTriangle)
	{
		record.m_pMaterial = m_pMaterial;
	}

	return hitTriangle;
}

//bool TriangleMesh::shadowHit(const Ray &r, Real tmin, Real tmax) const
//{
//	//if ( ! rayBoundingBoxIntersectionTest(r,tmin,tmax) )
//	//	return false;
//
//	for (int i=0; i<(int)m_Triangles.size(); i++)
//	{
//		const bool isHit = m_Triangles[i].shadowHit(r, tmin, tmax);
//		if (isHit) return true;
//	}
//
//	return false;
//}

void TriangleMesh::drawGL() const	// for preview using OpenGL
{
	if (!m_VBO)
		return;

	glEnableClientState(GL_VERTEX_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glDrawArrays(GL_TRIANGLES, 0, (GLsizei)m_Triangles.size() * 3);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glDisableClientState(GL_VERTEX_ARRAY);

	//for (int i=0; i<(int)m_Triangles.size(); i++)
	//{
	//	m_Triangles[i].drawGL();
	//}
}

void TriangleMesh::computeBoundingBox()
{
	m_BoundingBoxPos[0] = glm::vec3(  100000.f );
	m_BoundingBoxPos[1] = glm::vec3( -100000.f );

	for (int i=0; i<(int)m_Triangles.size(); i++)
	{
		const vec3 v0 = m_Triangles[i].getVertex0();
		const vec3 v1 = m_Triangles[i].getVertex1();
		const vec3 v2 = m_Triangles[i].getVertex2();

		m_BoundingBoxPos[0].x = min(m_BoundingBoxPos[0].x, min(v0.x, min(v1.x, v2.x)));
		m_BoundingBoxPos[0].y = min(m_BoundingBoxPos[0].y, min(v0.y, min(v1.y, v2.y)));
		m_BoundingBoxPos[0].z = min(m_BoundingBoxPos[0].z, min(v0.z, min(v1.z, v2.z)));

		m_BoundingBoxPos[1].x = max(m_BoundingBoxPos[1].x, max(v0.x, max(v1.x, v2.x)));
		m_BoundingBoxPos[1].y = max(m_BoundingBoxPos[1].y, max(v0.y, max(v1.y, v2.y)));
		m_BoundingBoxPos[1].z = max(m_BoundingBoxPos[1].z, max(v0.z, max(v1.z, v2.z)));
	}	
}

bool TriangleMesh::loadObj(const char* filename)
{
	fastObjMesh* m = fast_obj_read(filename);

	if (!m)
	{
		cerr << __FUNCTION__ << ": loading " << filename << " failed" << endl;
		return false;
	}

	vector<glm::vec3> vertices, normals;
	vector<glm::vec2> texCoords;

	if (m->position_count > 1)
	{
		vertices.resize(m->position_count - 1);
		memcpy(&vertices[0], m->positions + 3, 3 * sizeof(float) * (m->position_count - 1));
	}
	else
	{
		cerr << __FUNCTION__ << ": warning: no vertices defined" << endl;
		return false;
	}

	if (m->normal_count > 1 && m_ShadingType == Auto_Shading)
	{
		normals.resize(m->normal_count - 1);
		memcpy(&normals[0], m->normals + 3, 3 * sizeof(float) * (m->normal_count - 1));
	}

	if (m->texcoord_count > 1)
	{
		texCoords.resize(m->texcoord_count - 1);
		memcpy(&texCoords[0], m->texcoords + 2, 2 * sizeof(float) * (m->texcoord_count - 1));
	}

	int nNonTriangles = 0;

	vector<TriangleIndices> indices;
	indices.reserve(m->face_count);

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
						indices.emplace_back(VertexTuple(m0.p - 1, m0.t - 1, m0.n - 1), VertexTuple(m1.p - 1, m1.t - 1, m1.n - 1), VertexTuple(m2.p - 1, m2.t - 1, m2.n - 1));
					else
						indices.emplace_back(VertexTuple(m0.p - 1, m0.t - 1, -1), VertexTuple(m1.p - 1, m1.t - 1, -1), VertexTuple(m2.p - 1, m2.t - 1, -1));
				}
				else
				{
					if (m0.n)
						indices.emplace_back(VertexTuple(m0.p - 1, -1, m0.n - 1), VertexTuple(m1.p - 1, -1, m1.n - 1), VertexTuple(m2.p - 1, -1, m2.n - 1));
					else
						indices.emplace_back(VertexTuple(m0.p - 1), VertexTuple(m1.p - 1), VertexTuple(m2.p - 1));
				}

				m1 = m2;
			}
		}
	}

	fast_obj_destroy(m);

	if (nNonTriangles)
		cerr << __FUNCTION__ << ": Note: " << nNonTriangles << " non-triangles found (automatically converted to triangles)" << endl;

	cout << __FUNCTION__ << ": " << filename << " loaded" << endl;
	cout << "  # verts:\t" << vertices.size() << endl
		 << "  # normals:\t" << normals.size() << endl
		 << "  # tex coords:\t" << texCoords.size() << endl
		 << "  # triangles:\t" << m_Triangles.size() << endl;

	if (normals.empty())
		calcVertexNormals(normals, vertices, indices);

	scaleAndCenterize(vertices, 1.f);

	const int nTriangles = (int)indices.size();

	m_Triangles.clear();
	m_Triangles.resize(nTriangles);

	for (int ti = 0; ti < nTriangles; ++ti)
	{
		const auto t = indices[ti].indices;

		Triangle& tri = m_Triangles[ti];

		tri.setVertex0(vertices[t[0].vertexIdx]);
		tri.setVertex1(vertices[t[1].vertexIdx]);
		tri.setVertex2(vertices[t[2].vertexIdx]);

		if (m_ShadingType == Flat_Shading)
		{
			tri.setNormal0(normals[ti]);
			tri.setNormal1(normals[ti]);
			tri.setNormal2(normals[ti]);
		}
		else
		{
			tri.setNormal0(normals[t[0].normalIdx]);
			tri.setNormal1(normals[t[1].normalIdx]);
			tri.setNormal2(normals[t[2].normalIdx]);
		}

		if (texCoords.size())
		{
			tri.setTexCoord0(texCoords[t[0].textureIdx]);
			tri.setTexcoord1(texCoords[t[1].textureIdx]);
			tri.setTexcoord2(texCoords[t[2].textureIdx]);
		}
	}

	bakeVBO();

	return true;
}

void TriangleMesh::calcVertexNormals(vector<glm::vec3>& normals, const vector<glm::vec3>& vertices, const vector<TriangleIndices>& indices)
{
	const int nFaces = (int)indices.size();
	const int nVerts = (int)vertices.size();
	const int nNormals = (m_ShadingType == Flat_Shading) ? nFaces : nVerts;

	normals.resize(nNormals);

	//#pragma omp parallel
	{
		//#pragma omp for
		for (int ni = 0; ni < nNormals; ni++) normals[ni] = glm::vec3(0.0);

		//#pragma omp for
		for (int fi = 0; fi < nFaces; fi++)
		{
			const int i0 = indices[fi].indices[0].vertexIdx;
			const int i1 = indices[fi].indices[1].vertexIdx;
			const int i2 = indices[fi].indices[2].vertexIdx;

			if (i0 >= nVerts || i1 >= nVerts || i2 >= nVerts || i0 < 0 || i1 < 0 || i2 < 0)
			{
				cerr << __FUNCTION__ << ": warning: triangle index out of range: (" << i0 << "," << i1 << "," << i2 << ")" << endl;
			}

			const glm::vec3 weightedNormal = glm::cross((vertices[i1] - vertices[i0]), (vertices[i2] - vertices[i1]));

			if (m_ShadingType == Flat_Shading)
			{
				normals[fi] = glm::normalize(weightedNormal);
			}
			else
			{
				normals[i0] += weightedNormal;
				normals[i1] += weightedNormal;
				normals[i2] += weightedNormal;
			}
		}

		if (m_ShadingType == Smooth_Shading)
		{
			//#pragma omp for
			for (int ni = 0; ni < nVerts; ni++)
			{
				const float norm = glm::length(normals[ni]);
				if (norm > 0.f)
					normals[ni] = normals[ni] / norm;
			}
		}
	}
}

void TriangleMesh::scaleAndCenterize(vector<glm::vec3>& vertices, float scale)
{
	glm::vec3 bboxMin(1000000.f), bboxMax(-1000000.f);

	for (int i = 0; i < (int)vertices.size(); ++i)
	{
		const auto v = vertices[i];
		if (bboxMin.x > v.x) bboxMin.x = v.x;
		if (bboxMin.y > v.y) bboxMin.y = v.y;
		if (bboxMin.z > v.z) bboxMin.z = v.z;
		if (bboxMax.x < v.x) bboxMax.x = v.x;
		if (bboxMax.y < v.y) bboxMax.y = v.y;
		if (bboxMax.z < v.z) bboxMax.z = v.z;
	}

	const glm::vec3 diff = bboxMax - bboxMin;
	const float s = scale / std::max(diff.x, std::max(diff.y, diff.z));

	for (int i = 0; i < (int)vertices.size(); ++i)
	{
		vertices[i] = s * (vertices[i] - bboxMin);
	}
}

void TriangleMesh::bakeVBO()
{
	const int nTriangles = (int)m_Triangles.size();
	vector<glm::vec3> vertices(3 * nTriangles);

	for (int ti = 0; ti < nTriangles; ++ti)
	{
		vertices[3 * ti + 0] = m_Triangles[ti].getVertex0();
		vertices[3 * ti + 1] = m_Triangles[ti].getVertex1();
		vertices[3 * ti + 2] = m_Triangles[ti].getVertex2();
	}

	if (!m_VBO) glGenBuffers(1, &m_VBO);
	glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), &vertices[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

#if 0

// private methods

bool TriangleMesh::rayBoundingBoxIntersectionTest(const Ray &r, Real tmin, Real tmax) const
{
	float intervalMin = tmin;
	float intervalMax = tmax;

	const vec3 &origin = r.getOrigin();
	const vec3 &dir = r.getUnitDir();

	// test x direction

	if (dir.x == 0.f)
	{
		if (origin.x < m_BoundingBoxPos[0].x || m_BoundingBoxPos[1].x < origin.x)
			return false;
	}
	else
	{
		const int xSign = (dir.x>0.f) ? 0 : 1;
		const float tMinXDir = (m_BoundingBoxPos[xSign].x - origin.x) / dir.x;
		const float tMaxXDir = (m_BoundingBoxPos[1-xSign].x - origin.x) / dir.x;

		if (tMinXDir > intervalMin) intervalMin = tMinXDir;
		if (tMaxXDir < intervalMax) intervalMax = tMaxXDir;

		if (intervalMin > intervalMax) return false;
	}

	// test y direction

	if (dir.y == 0.f)
	{
		if (origin.y < m_BoundingBoxPos[0].y || m_BoundingBoxPos[1].y < origin.y)
			return false;
	}
	else
	{
		const int ySign = (dir.y>0.f) ? 0 : 1;
		const float tMinYDir = (m_BoundingBoxPos[ySign].y - origin.y) / dir.y;
		const float tMaxYDir = (m_BoundingBoxPos[1-ySign].y - origin.y) / dir.y;

		if (tMinYDir > intervalMin) intervalMin = tMinYDir;
		if (tMaxYDir < intervalMax) intervalMax = tMaxYDir;

		if (intervalMin > intervalMax) return false;
	}

	// test z direction

	if (dir.z == 0.f)
	{
		if (origin.z < m_BoundingBoxPos[0].z || m_BoundingBoxPos[1].z < origin.z)
			return false;
	}
	else
	{
		const int zSign = (dir.z>0.f) ? 0 : 1;
		const float tMinZDir = (m_BoundingBoxPos[zSign].z - origin.z) / dir.z;
		const float tMaxZDir = (m_BoundingBoxPos[1-zSign].z - origin.z) / dir.z;

		if (tMinZDir > intervalMin) intervalMin = tMinZDir;
		if (tMaxZDir < intervalMax) intervalMax = tMaxZDir;

		if (intervalMin > intervalMax) return false;
	}

	return true;
}

#endif
