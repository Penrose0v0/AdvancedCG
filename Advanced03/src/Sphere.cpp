#include "Sphere.h"
#include <cstdlib>
#include <cmath>
//#include <GL/glut.h>

using namespace std;
//using namespace MyAlgebra;

int Sphere::s_NumSphereVertices = 0;
GLuint Sphere::s_VBO = 0;

bool Sphere::hit(const Ray &r, Real tmin, Real tmax, HitRecord &record) const
{
	const vec3 v = r.getOrigin() - m_Center;
	const vec3 d = r.getUnitDir();

	const Real b = glm::dot(v, d);
	const Real c = glm::dot(v, v) - m_Radius*m_Radius;

	const Real D = b*b - c;

	if (D > 0.0)
	{
		const Real sqrtD = sqrtf(D);
		Real t = -b - sqrtD;

		// check if t is within a valid interval
		if (t < tmin)
			t = -b + sqrtD;
		if (t < tmin || t > tmax)
			return false;

		// valid hit

		vec3 n = glm::normalize(r.calculatePosition( t ) - m_Center);
		//n.normalize();

		record.m_ParamT = t;
		record.m_Normal = n;
		record.m_HitPos = r.calculatePosition(t);
		//record.m_TexCoords.set( 0.f );
		record.m_TexCoords = glm::vec3(0.f);
		record.m_pMaterial = m_pMaterial;

		return true;
	}

	return false;
}

//bool Sphere::shadowHit(const Ray &r, Real tmin, Real tmax) const
//{
//	const vec3 v = r.getOrigin() - center;
//	const vec3 d = r.getUnitDir();
//
//	const Real b = dot(v, d);
//	const Real c = v.sq_norm() - radius*radius;
//
//	const Real D = b*b - c;
//
//	if (D > 0.0)
//	{
//		const Real sqrtD = sqrtf(D);
//		Real t = -b - sqrtD;
//
//		// check if t is within a valid interval
//		if (t < tmin)
//			t = -b + sqrtD;
//		if (t < tmin || t > tmax)
//			return false;
//
//		// valid hit
//
//		return true;
//	}
//
//	return false;
//}

// for preview using OpenGL
void Sphere::drawGL() const
{
	if (!s_VBO)
		BakeVBO();

	glPushMatrix();
	glTranslatef( m_Center.x, m_Center.y, m_Center.z );
	glScalef(m_Radius, m_Radius, m_Radius);

	glEnableClientState(GL_VERTEX_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, s_VBO);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	//glDrawArrays(GL_TRIANGLES, 0, s_NumSphereVertices);
	glDrawArrays(GL_QUADS, 0, s_NumSphereVertices);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glDisableClientState(GL_VERTEX_ARRAY);

	glPopMatrix();
}

void Sphere::BakeVBO()
{
	const int nDivSphereLatitude = 16;
	const int nDivSphereLongitude = 32;

	const float radius = 1.f;
	const float dRadLat = glm::pi<float>() / (float)nDivSphereLatitude;
	const float dRadLon = 2 * glm::pi<float>() / (float)nDivSphereLongitude;

	float radLatitude = 0.f;

	vector<glm::vec3> vertices;

	for (int vi = 0; vi < nDivSphereLatitude; vi++)
	{
		const float lat0 = radLatitude;
		const float lat1 = radLatitude + dRadLat;

		float radLongitude = 0.f;

		for (int ui = 0; ui < nDivSphereLongitude; ui++)
		{
			const float lon0 = radLongitude;
			const float lon1 = radLongitude + dRadLon;

			//vertices.emplace_back(radius * sinf(lat0) * cosf(lon1), radius * cosf(lat0), radius * sinf(lat0) * sinf(lon1));
			//vertices.emplace_back(radius * sinf(lat0) * cosf(lon0), radius * cosf(lat0), radius * sinf(lat0) * sinf(lon0));
			//vertices.emplace_back(radius * sinf(lat1) * cosf(lon0), radius * cosf(lat1), radius * sinf(lat1) * sinf(lon0));
			//vertices.emplace_back(radius * sinf(lat0) * cosf(lon1), radius * cosf(lat0), radius * sinf(lat0) * sinf(lon1));
			//vertices.emplace_back(radius * sinf(lat1) * cosf(lon0), radius * cosf(lat1), radius * sinf(lat1) * sinf(lon0));
			//vertices.emplace_back(radius * sinf(lat1) * cosf(lon1), radius * cosf(lat1), radius * sinf(lat1) * sinf(lon1));
			vertices.emplace_back(radius * sinf(lat0) * cosf(lon1), radius * cosf(lat0), radius * sinf(lat0) * sinf(lon1));
			vertices.emplace_back(radius * sinf(lat0) * cosf(lon0), radius * cosf(lat0), radius * sinf(lat0) * sinf(lon0));
			vertices.emplace_back(radius * sinf(lat1) * cosf(lon0), radius * cosf(lat1), radius * sinf(lat1) * sinf(lon0));
			vertices.emplace_back(radius * sinf(lat1) * cosf(lon1), radius * cosf(lat1), radius * sinf(lat1) * sinf(lon1));

			radLongitude = lon1;
		}

		radLatitude = lat1;
	}

	if (!s_VBO) glGenBuffers(1, &s_VBO);
	glBindBuffer(GL_ARRAY_BUFFER, s_VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * vertices.size(), &vertices[0].x, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	s_NumSphereVertices = (int)vertices.size();
}
