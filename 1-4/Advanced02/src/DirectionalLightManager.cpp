#include "DirectionalLightManager.h"
#include "PathFinder.h"
#include "arcball_camera.h"
#include <cstdlib>
#include <random>
#include <iostream>

using namespace std;
using namespace glm;

TriMesh DirectionalLightManager::s_ArrowMesh;
GLSLProgramObject* DirectionalLightManager::s_pShader = 0;
GLuint DirectionalLightManager::s_VAO = 0;
std::vector<glm::vec3> DirectionalLightManager::s_Colors;
std::vector<glm::vec2> DirectionalLightManager::s_Dirs;

extern ArcballCamera g_Camera;
extern glm::mat4 g_ProjMatrix;

void DirectionalLightManager::Init()
{
	if (s_pShader) delete s_pShader;
	s_pShader = new GLSLProgramObject();

	PathFinder finder;
#ifdef __APPLE__
	finder.addSearchPath("GLSL_Mac");
	finder.addSearchPath("../GLSL_Mac");
	finder.addSearchPath("../../GLSL_Mac");
#else
	finder.addSearchPath("GLSL");
	finder.addSearchPath("../GLSL");
	finder.addSearchPath("../../GLSL");
#endif

	s_pShader->attachShaderSourceFile(finder.find("arrow.vert").c_str(), GL_VERTEX_SHADER);
	s_pShader->attachShaderSourceFile(finder.find("arrow.frag").c_str(), GL_FRAGMENT_SHADER);

	s_pShader->link();

	if (!s_pShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": shader link failed" << endl;
		s_pShader->printProgramLog();
		return;
	}

	BuildArrowMesh();

	if (!s_ArrowMesh.getVertexVBO() || !s_ArrowMesh.getVertexNormalVBO())
	{
		cerr << __FUNCTION__ << ": mesh VBOs not ready" << endl;
		return;
	}

	GenRandomLights(1);
	UpdateLightVAO();
}

void DirectionalLightManager::GenRandomLights(int nLights)
{
	s_Colors.clear();
	s_Dirs.clear();

	mt19937 randSrc(12345);
	uniform_real_distribution<float> randDist(0, 1);

	const float colorScale = 1.f;

	for (int li = 0; li < nLights; ++li)
	{
		s_Colors.emplace_back(colorScale * (0.5f * randDist(randSrc) + 0.5f),
							  colorScale * (0.5f * randDist(randSrc) + 0.5f),
							  colorScale * (0.5f * randDist(randSrc) + 0.5f));
		s_Dirs.emplace_back(2 * glm::pi<float>() * randDist(randSrc), 
							glm::pi<float>() / 2.25f * randDist(randSrc));
	}
}

void DirectionalLightManager::UpdateLightVAO()
{
	GLuint vertexColorVBO, sphericalCoordVBO;
	glGenBuffers(1, &vertexColorVBO);
	glBindBuffer(GL_ARRAY_BUFFER, vertexColorVBO);
	glBufferData(GL_ARRAY_BUFFER, s_Colors.size() * sizeof(glm::vec3), &s_Colors[0], GL_STATIC_DRAW);

	glGenBuffers(1, &sphericalCoordVBO);
	glBindBuffer(GL_ARRAY_BUFFER, sphericalCoordVBO);
	glBufferData(GL_ARRAY_BUFFER, s_Dirs.size() * sizeof(glm::vec2), &s_Dirs[0], GL_STATIC_DRAW);

	if (!s_VAO) glGenVertexArrays(1, &s_VAO);
	glBindVertexArray(s_VAO);

	glBindBuffer(GL_ARRAY_BUFFER, s_ArrowMesh.getVertexVBO());
	glEnableVertexAttribArray(s_pShader->getAttributeLocation("vertexPosition"));
	glVertexAttribPointer(s_pShader->getAttributeLocation("vertexPosition"), 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	glBindBuffer(GL_ARRAY_BUFFER, s_ArrowMesh.getVertexNormalVBO());
	glEnableVertexAttribArray(s_pShader->getAttributeLocation("vertexNormal"));
	glVertexAttribPointer(s_pShader->getAttributeLocation("vertexNormal"), 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	glBindBuffer(GL_ARRAY_BUFFER, vertexColorVBO);
	glEnableVertexAttribArray(s_pShader->getAttributeLocation("vertexColor"));
	glVertexAttribPointer(s_pShader->getAttributeLocation("vertexColor"), 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	glBindBuffer(GL_ARRAY_BUFFER, sphericalCoordVBO);
	glEnableVertexAttribArray(s_pShader->getAttributeLocation("sphericalCoord"));
	glVertexAttribPointer(s_pShader->getAttributeLocation("sphericalCoord"), 2, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	glVertexAttribDivisor(s_pShader->getAttributeLocation("vertexPosition"), 0);	// no division
	glVertexAttribDivisor(s_pShader->getAttributeLocation("vertexNormal"), 0);	// no division
	glVertexAttribDivisor(s_pShader->getAttributeLocation("vertexColor"), 1);	// instance-wise
	glVertexAttribDivisor(s_pShader->getAttributeLocation("sphericalCoord"), 1);	// instance-wise

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glDeleteBuffers(1, &vertexColorVBO);
	glDeleteBuffers(1, &sphericalCoordVBO);
}

void DirectionalLightManager::Draw(int nLights /*= -1*/)
{
	if (!s_pShader || !s_pShader->linkSucceeded())
		return;

	//glPushAttrib(GL_ENABLE_BIT);	// deprecated from OpenGL 3.0
	glEnable(GL_DEPTH_TEST);

	auto modelViewMatrix = g_Camera.transform();

	float arrowPhi = glm::pi<float>() / 3.f;
	float arrowTheta = glm::pi<float>() / 4.f;

	s_pShader->use();
	s_pShader->sendUniformMatrix4fv("projMatrix", glm::value_ptr(g_ProjMatrix));
	s_pShader->sendUniformMatrix4fv("modelViewMatrix", glm::value_ptr(modelViewMatrix));
	s_pShader->sendUniform1f("scale", 0.3f);
	s_pShader->sendUniform1f("radius", 3.f);

	if (nLights == -1)
		nLights = s_Colors.size();

	glBindVertexArray(s_VAO);
	//glDrawArrays(GL_TRIANGLES, 0, 3 * s_ArrowMesh.getNumTriangles());
	glDrawArraysInstanced(GL_TRIANGLES, 0, 3 * s_ArrowMesh.getNumTriangles(), nLights);
	glBindVertexArray(0);

	s_pShader->disable();

	//glPopAttrib();	// deprecated from OpenGL 3.0
}

void DirectionalLightManager::BuildArrowMesh()
{
	const float coneHeight = 0.4f;
	const float cylinderHeight = 0.6f;
	const float coneRadius = 0.4f;
	const float cylinderRadius = 0.2f;

	const int nDiv = 64;
	const float dRad = 2.f * pi<float>() / nDiv;

	vector<vec3> vertices;
	vector<TriangleIndices> indices;

	// top cone

	int partStartIndex = 0;
	int prevIndex = 0;
	vertices.emplace_back(coneRadius * cosf(0.f), cylinderHeight, coneRadius * sinf(0.f));

	for (int i = 1; i <= nDiv; ++i)
	{
		if (i == nDiv)
		{
			vertices.emplace_back(0.f, cylinderHeight + coneHeight, 0.f);
			indices.emplace_back(prevIndex, partStartIndex, prevIndex + 1);
		}
		else
		{
			const float rad = i * dRad;
			vertices.emplace_back(0.f, cylinderHeight + coneHeight, 0.f);
			vertices.emplace_back(coneRadius * cosf(rad), cylinderHeight, coneRadius * sinf(rad));
			indices.emplace_back(prevIndex, prevIndex + 2, prevIndex + 1);
			prevIndex = (int)vertices.size() - 1;
		}
	}

	// cone bottom

	{
		vertices.emplace_back(0.f, cylinderHeight, 0.f);
		const int coneBottomCenterIndex = (int)vertices.size() - 1;

		vertices.emplace_back(coneRadius * cosf(0.f), cylinderHeight, coneRadius * sinf(0.f));
		prevIndex = (int)vertices.size() - 1;
		partStartIndex = prevIndex;

		for (int i = 1; i <= nDiv; ++i)
		{
			if (i == nDiv)
			{
				indices.emplace_back(prevIndex, coneBottomCenterIndex, partStartIndex);
			}
			else
			{
				const float rad = i * dRad;
				vertices.emplace_back(coneRadius * cosf(rad), cylinderHeight, coneRadius * sinf(rad));
				indices.emplace_back(prevIndex, coneBottomCenterIndex, prevIndex + 1);
				prevIndex = (int)vertices.size() - 1;
			}
		}
	}

	// cylinder body

	vertices.emplace_back(cylinderRadius * cosf(0.f), cylinderHeight, cylinderRadius * sinf(0.f));
	vertices.emplace_back(cylinderRadius * cosf(0.f), 0.f, cylinderRadius * sinf(0.f));
	prevIndex = (int)vertices.size() - 2;
	partStartIndex = prevIndex;

	for (int i = 1; i <= nDiv; ++i)
	{
		if (i == nDiv)
		{
			indices.emplace_back(prevIndex, prevIndex + 1, partStartIndex);
			indices.emplace_back(prevIndex + 1, partStartIndex + 1, partStartIndex);
		}
		else
		{
			const float rad = i * dRad;
			vertices.emplace_back(cylinderRadius * cosf(rad), cylinderHeight, cylinderRadius * sinf(rad));
			vertices.emplace_back(cylinderRadius * cosf(rad), 0.f, cylinderRadius * sinf(rad));
			indices.emplace_back(prevIndex, prevIndex + 1, prevIndex + 2);
			indices.emplace_back(prevIndex + 1, prevIndex + 3, prevIndex + 2);
			prevIndex = (int)vertices.size() - 2;
		}
	}

	// cylinder bottom

	{
		vertices.emplace_back(0.f, 0.f, 0.f);
		const int cylinderBottomCenterIndex = (int)vertices.size() - 1;

		vertices.emplace_back(cylinderRadius * cosf(0.f), 0.f, cylinderRadius * sinf(0.f));
		prevIndex = (int)vertices.size() - 1;
		partStartIndex = prevIndex;

		for (int i = 1; i <= nDiv; ++i)
		{
			if (i == nDiv)
			{
				indices.emplace_back(prevIndex, cylinderBottomCenterIndex, partStartIndex);
			}
			else
			{
				const float rad = i * dRad;
				vertices.emplace_back(cylinderRadius * cosf(rad), 0.f, cylinderRadius * sinf(rad));
				indices.emplace_back(prevIndex, cylinderBottomCenterIndex, prevIndex + 1);
				prevIndex = (int)vertices.size() - 1;
			}
		}
	}

	s_ArrowMesh.setVertices(vertices);
	s_ArrowMesh.setTriangles(indices);
	s_ArrowMesh.calcVertexNormals();
	s_ArrowMesh.bakeVBOs();
}
