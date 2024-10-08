#pragma once

#include <GL/glew.h>
#include "ImageRect.h"
#include <glm/glm.hpp>

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

	inline GLuint getTexID() const { return m_TexID; }

	void genStripTexture(int size = 1024, float bandWidth = 0.01f, glm::vec3 color0 = glm::vec3(0.f), glm::vec3 color1 = glm::vec3(1,1,0));

	bool load(const char* filename);

	void drawGL() const;

private:
	ImageRGBf m_Texture;
	int m_NumSphereVertices;
	GLuint m_VBO, m_TexID;

	void bakeVBO();
};