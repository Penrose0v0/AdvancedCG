#pragma once

#include <GL/glew.h>
#include "ImageRect.h"
#include "Ray.h"

class EnvironmentMap
{
public:
	typedef ImageRect<glm::vec3> ImageRGBf;

	EnvironmentMap() : m_NumSphereVertices(0), m_VBO(0), m_TexID(0) {}
	~EnvironmentMap()
	{
		if (m_VBO) glDeleteBuffers(1, &m_VBO);
		if (m_TexID) glDeleteTextures(1, &m_TexID);
	}

	glm::vec3 fetchColor(const Ray &ray) const;

	bool load(const char* filename);

	void drawGL() const;

private:
	ImageRGBf m_Texture;
	int m_NumSphereVertices;
	GLuint m_VBO, m_TexID;

	void bakeVBO();
};