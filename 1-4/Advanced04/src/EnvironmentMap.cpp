#include "EnvironmentMap.h"
#include "glm/glm.hpp"
#include "arcball_camera.h"
#include "ImageRect.h"
#include <vector>
#undef _UNICODE
#include <IL/il.h>

using namespace std;
using namespace glm;

extern ArcballCamera g_Camera;

void EnvironmentMap::genStripTexture(int size, float bandWidth, vec3 color0, vec3 color1)
{
	m_Texture.allocate(size, size);
	const int bandWidthInt = size * bandWidth;

	for (int yi = 0; yi < size; ++yi)
	{
		const vec3 color = ((yi % bandWidthInt) / float(bandWidthInt) < 0.5f) ? color0 : color1;

		for (int xi = 0; xi < size; ++xi)
		{
			m_Texture(xi, yi) = color;
		}
	}

	if (!m_TexID) glGenTextures(1, &m_TexID);
	glBindTexture(GL_TEXTURE_2D, m_TexID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, size, size, 0, GL_RGB, GL_FLOAT, m_Texture.getData());
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);	// deprecated from OpenGL 3.0
	glBindTexture(GL_TEXTURE_2D, 0);

	bakeVBO();
}

bool EnvironmentMap::load(const char* filename)
{
	ILuint imgName;
	ilGenImages(1, &imgName);
	ilBindImage(imgName);

	ilEnable(IL_ORIGIN_SET);
	ilOriginFunc(IL_ORIGIN_LOWER_LEFT);

	if (!ilLoadImage(filename))
	{
		cerr << __FUNCTION__ << ": Error: cannot open " << filename << endl;
		ilDeleteImages(1, &imgName);
		return false;
	}

	const int width = ilGetInteger(IL_IMAGE_WIDTH);
	const int height = ilGetInteger(IL_IMAGE_HEIGHT);
	m_Texture.allocate(width, height);
	ilCopyPixels(0, 0, 0, width, height, 1, IL_RGB, IL_FLOAT, m_Texture.getData());

	if (!m_TexID) glGenTextures(1, &m_TexID);
	glBindTexture(GL_TEXTURE_2D, m_TexID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, width, height, 0, GL_RGB, GL_FLOAT, m_Texture.getData());
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);	// deprecated from OpenGL 3.0
	glBindTexture(GL_TEXTURE_2D, 0);

	ilDeleteImages(1, &imgName);

	cerr << __FUNCTION__ << ": file loaded: " << filename << " (" << width << "x" << height << ")" << endl;

	bakeVBO();

	return true;
}

void EnvironmentMap::drawGL() const
{
	if (!m_TexID || !m_VBO)
		return;

	mat4 M = g_Camera.transform();
	M[3][0] = M[3][1] = M[3][2] = 0.f;

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadMatrixf(value_ptr(M));

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
	glVertexPointer(3, GL_FLOAT, 0, 0);
	glTexCoordPointer(2, GL_FLOAT, 0, (GLvoid *)(sizeof(vec3) * m_NumSphereVertices));

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, m_TexID);

	glDrawArrays(GL_QUADS, 0, m_NumSphereVertices);

	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);

	glPopMatrix();
}

void EnvironmentMap::bakeVBO()
{
	const int nDivSphereLatitude = 64;
	const int nDivSphereLongitude = 128;

	const float radius = 100.f;
	const float dRadLat = glm::pi<float>() / (float)nDivSphereLatitude;
	const float dRadLon = 2 * glm::pi<float>() / (float)nDivSphereLongitude;

	float radLatitude = 0.f;

	vector<glm::vec3> vertices;
	vector<glm::vec2> texCoords;

	for (int vi = 0; vi < nDivSphereLatitude; vi++)
	{
		const float lat0 = radLatitude;
		const float lat1 = radLatitude + dRadLat;

		float radLongitude = 0.f;

		for (int ui = 0; ui < nDivSphereLongitude; ui++)
		{
			const float lon0 = radLongitude;
			const float lon1 = radLongitude + dRadLon;

			vertices.emplace_back(radius * sinf(lat0) * cosf(lon1), radius * cosf(lat0), radius * sinf(lat0) * sinf(lon1));
			vertices.emplace_back(radius * sinf(lat0) * cosf(lon0), radius * cosf(lat0), radius * sinf(lat0) * sinf(lon0));
			vertices.emplace_back(radius * sinf(lat1) * cosf(lon0), radius * cosf(lat1), radius * sinf(lat1) * sinf(lon0));
			vertices.emplace_back(radius * sinf(lat1) * cosf(lon1), radius * cosf(lat1), radius * sinf(lat1) * sinf(lon1));

			texCoords.emplace_back(lon1 / (2 * pi<float>()), 1.f - lat0 / pi<float>());
			texCoords.emplace_back(lon0 / (2 * pi<float>()), 1.f - lat0 / pi<float>());
			texCoords.emplace_back(lon0 / (2 * pi<float>()), 1.f - lat1 / pi<float>());
			texCoords.emplace_back(lon1 / (2 * pi<float>()), 1.f - lat1 / pi<float>());

			radLongitude = lon1;
		}

		radLatitude = lat1;
	}

	size_t vertexSize = sizeof(glm::vec3) * vertices.size();
	size_t texCoordSize = sizeof(glm::vec2) * texCoords.size();

	if (!m_VBO) glGenBuffers(1, &m_VBO);
	glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
	glBufferData(GL_ARRAY_BUFFER, vertexSize + texCoordSize, NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, vertexSize, &vertices[0].x);
	glBufferSubData(GL_ARRAY_BUFFER, vertexSize, texCoordSize, &texCoords[0].x);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	m_NumSphereVertices = (int)vertices.size();
}
