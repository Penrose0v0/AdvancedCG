#include <GL/glew.h>
#include "Triangle.h"
#include <cstdlib>
//#include <GL/glut.h>

using namespace std;

bool Triangle::hit(const Ray &r, Real tmin, Real tmax, HitRecord &record) const
{
	const vec3 d = r.getUnitDir();
	const vec3 e1 = m_Vertices[1] - m_Vertices[0];
	const vec3 e2 = m_Vertices[2] - m_Vertices[0];

	//const vec3 pVec = d ^ e2;
	const vec3 pVec = glm::cross(d, e2);

	const Real det = glm::dot(e1, pVec);

	if (det > -1.0e-10 && det < 1.0e-10)
		return false;

	const Real invDet = 1.f/det;

	const vec3 tVec = r.getOrigin() - m_Vertices[0];

	const Real beta = invDet * glm::dot(tVec, pVec);

	if (beta < 0.f || beta > 1.f)
		return false;

	//const vec3 qVec = tVec ^ e1;
	const vec3 qVec = glm::cross(tVec, e1);

	const Real gamma = invDet * glm::dot(d, qVec);

	if (gamma < 0.f || beta + gamma > 1.f)
		return false;

	const Real t = invDet * glm::dot(e2, qVec);

	if (t < tmin || t > tmax)
		return false;

	const float alpha = 1.f - (beta + gamma);

	const vec3 normal = alpha * m_Normals[0] + beta * m_Normals[1] + gamma * m_Normals[2];
	const vec2 texCoord = alpha * m_TexCoords[0] + beta * m_TexCoords[1] + gamma * m_TexCoords[2];

	record.m_ParamT = t;
	record.m_Normal = normal;
	record.m_HitPos = r.calculatePosition(t);
	record.m_TexCoords = vec3(texCoord.x, texCoord.y, 0.f);
	record.m_pMaterial = m_pMaterial;

	return true;
}

//bool Triangle::shadowHit(const Ray &r, Real tmin, Real tmax) const
//{
//	const vec3 d = r.getUnitDir();
//	const vec3 e1 = m_Vertices[1] - m_Vertices[0];
//	const vec3 e2 = m_Vertices[2] - m_Vertices[0];
//
//	const vec3 pVec = d ^ e2;
//
//	const Real det = glm::dot(e1, pVec);
//
//	if (det > -1.0e-10 && det < 1.0e-10)
//		return false;
//
//	const Real invDet = 1.0/det;
//
//	const vec3 tVec = r.getOrigin() - m_Vertices[0];
//
//	const Real u = invDet * dot(tVec, pVec);
//
//	if (u < 0.0 || u > 1.0)
//		return false;
//
//	const vec3 qVec = tVec ^ e1;
//
//	const Real v = invDet * dot(d, qVec);
//
//	if (v < 0.0 || u + v > 1.0)
//		return false;
//
//	const Real t = invDet * dot(e2, qVec);
//
//	if (t < tmin || t > tmax)
//		return false;
//
//	return true;
//}

// for preview using OpenGL
void Triangle::drawGL() const
{
	glBegin(GL_TRIANGLES);
	glVertex3fv(glm::value_ptr(m_Vertices[0]));
	glVertex3fv(glm::value_ptr(m_Vertices[1]));
	glVertex3fv(glm::value_ptr(m_Vertices[2]));
	glEnd();

//#if 1
//	glBegin(GL_TRIANGLES);
//	glColor3fv( 0.5f*m_Normals[0] + vec3(0.5f,0.5f,0.5f) ); glVertex3fv( m_Vertices[0] );
//	glColor3fv( 0.5f*m_Normals[1] + vec3(0.5f,0.5f,0.5f) ); glVertex3fv( m_Vertices[1] );
//	glColor3fv( 0.5f*m_Normals[2] + vec3(0.5f,0.5f,0.5f) ); glVertex3fv( m_Vertices[2] );
//	glEnd();
//#endif
//#if 0
//	glBegin(GL_TRIANGLES);
//	glNormal3fv( m_Normals[0] ); glVertex3fv( m_Vertices[0] );
//	glNormal3fv( m_Normals[1] ); glVertex3fv( m_Vertices[1] );
//	glNormal3fv( m_Normals[2] ); glVertex3fv( m_Vertices[2] );
//	glEnd();
//#endif
}
