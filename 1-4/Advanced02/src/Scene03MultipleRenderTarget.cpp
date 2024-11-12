#include "Scene03MultipleRenderTarget.h"
#include <iostream>
#include <cmath>
#include "imgui.h"
#include "PathFinder.h"
#include "TriMesh.h"
#include "arcball_camera.h"
#include "DirectionalLightManager.h"
#include "MaterialParameters.h"
#include "CheckGLError.h"

using namespace std;

GLSLProgramObject* Scene03MultipleRenderTarget::s_pTexShader = 0;
GLSLProgramObject* Scene03MultipleRenderTarget::s_pShader = 0;

GLuint Scene03MultipleRenderTarget::s_MeshVAO = 0;
GLuint Scene03MultipleRenderTarget::s_PlaneVAO = 0;
GLuint Scene03MultipleRenderTarget::s_QuadVAO = 0;
GLuint Scene03MultipleRenderTarget::s_MRTFBO = 0;
GLuint Scene03MultipleRenderTarget::s_DepthRB = 0;
GLuint* Scene03MultipleRenderTarget::s_pTexIDs = 0;

int Scene03MultipleRenderTarget::s_CachedWidth = 0;
int Scene03MultipleRenderTarget::s_CachedHeight = 0;

// global variables defined in main.cpp

extern TriMesh g_TriMesh, g_Plane;
extern MaterialParameters g_Material;
extern ArcballCamera g_Camera;
extern glm::mat4 g_ProjMatrix;

void Scene03MultipleRenderTarget::Init()
{
	ReloadShaders();
}

void Scene03MultipleRenderTarget::ReloadShaders()
{
	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

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

	const GLuint vertexPositionLocation = 0;
	const GLuint vertexNormalLocation = 1;
	const GLuint inTexCoordLocation = 1;

	if (!s_pTexShader) delete s_pTexShader;
	s_pTexShader = new GLSLProgramObject();
	s_pTexShader->attachShaderSourceFile(finder.find("tex.vert").c_str(), GL_VERTEX_SHADER);
	s_pTexShader->attachShaderSourceFile(finder.find("tex.frag").c_str(), GL_FRAGMENT_SHADER);
	s_pTexShader->setAttributeLocation("vertexPosition", vertexPositionLocation);
	s_pTexShader->setAttributeLocation("inTexCoord", inTexCoordLocation);
	s_pTexShader->link();

	if (!s_pTexShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": tex shader link failed" << endl;
		s_pTexShader->printProgramLog();
	}

	if (!s_pShader) delete s_pShader;
	s_pShader = new GLSLProgramObject();
	s_pShader->attachShaderSourceFile(finder.find("mrt.vert").c_str(), GL_VERTEX_SHADER);
	s_pShader->attachShaderSourceFile(finder.find("mrt.frag").c_str(), GL_FRAGMENT_SHADER);
	s_pShader->setAttributeLocation("vertexPosition", vertexPositionLocation);
	s_pShader->setAttributeLocation("vertexNormal", vertexNormalLocation);
	s_pShader->link();

	if (!s_pShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": mrt link failed" << endl;
		s_pShader->printProgramLog();
	}

	if (!s_MeshVAO) glGenVertexArrays(1, &s_MeshVAO);
	glBindVertexArray(s_MeshVAO);

	glBindBuffer(GL_ARRAY_BUFFER, g_TriMesh.getVertexVBO());
	glEnableVertexAttribArray(vertexPositionLocation);
	glVertexAttribPointer(vertexPositionLocation, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	glBindBuffer(GL_ARRAY_BUFFER, g_TriMesh.getVertexNormalVBO());
	glEnableVertexAttribArray(vertexNormalLocation);
	glVertexAttribPointer(vertexNormalLocation, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	if (!s_PlaneVAO) glGenVertexArrays(1, &s_PlaneVAO);
	glBindVertexArray(s_PlaneVAO);

	glBindBuffer(GL_ARRAY_BUFFER, g_Plane.getVertexVBO());
	glEnableVertexAttribArray(vertexPositionLocation);
	glVertexAttribPointer(vertexPositionLocation, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	glBindBuffer(GL_ARRAY_BUFFER, g_Plane.getVertexNormalVBO());
	glEnableVertexAttribArray(vertexNormalLocation);
	glVertexAttribPointer(vertexNormalLocation, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	float quadVertices[] = { -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1 };
	float quadTexCoords[] = { 0,  0, 1,  0, 1, 1,  0,  0, 1, 1,  0, 1 };

	if (!s_QuadVAO) glGenVertexArrays(1, &s_QuadVAO);
	glBindVertexArray(s_QuadVAO);

	GLuint quadVBO;
	glGenBuffers(1, &quadVBO);
	glBindBuffer(GL_ARRAY_BUFFER, quadVBO);

	glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices) + sizeof(quadTexCoords), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(quadVertices), quadVertices);
	glBufferSubData(GL_ARRAY_BUFFER, sizeof(quadVertices), sizeof(quadTexCoords), quadTexCoords);

	glEnableVertexAttribArray(vertexPositionLocation);
	glVertexAttribPointer(vertexPositionLocation, 2, GL_FLOAT, GL_FALSE, 0, (const void*)0); // vertices

	glEnableVertexAttribArray(inTexCoordLocation);
	glVertexAttribPointer(inTexCoordLocation, 2, GL_FLOAT, GL_FALSE, 0, (const void*)sizeof(quadVertices)); // texture coordinates
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);
}

void Scene03MultipleRenderTarget::Draw()
{
	if (!s_PlaneVAO || !s_MeshVAO)
		return;

	if (!s_pTexIDs || s_CachedWidth != s_WindowWidth)
		AllocateFBO(s_WindowWidth, s_WindowHeight);

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

	auto viewMatrix = g_Camera.transform();
	auto modelViewMatrix = viewMatrix * g_TriMesh.getModelMatrix();

	const GLenum colorAttachements[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };

	//glPushAttrib(GL_ENABLE_BIT);	// deprecated from OpenGL 3.0
	glEnable(GL_DEPTH_TEST);

	glBindFramebuffer(GL_FRAMEBUFFER, s_MRTFBO);
	glDrawBuffers(3, colorAttachements);
	glViewport(0, 0, s_WindowWidth, s_WindowHeight);
	glClearColor(0, 0, 0, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

	s_pShader->use();
	s_pShader->sendUniformMatrix4fv("projMatrix", glm::value_ptr(g_ProjMatrix));
	s_pShader->sendUniformMatrix4fv("modelViewMatrix", glm::value_ptr(modelViewMatrix));

	// 0v0
	glm::mat3 modelViewInvTransposed = glm::transpose(glm::inverse(glm::mat3(modelViewMatrix)));
	s_pShader->sendUniformMatrix3fv("modelViewInvTransposed", glm::value_ptr(modelViewInvTransposed));

	glBindVertexArray(s_MeshVAO);
	glDrawArrays(GL_TRIANGLES, 0, 3 * g_TriMesh.getNumTriangles());
	glBindVertexArray(0);

	s_pShader->sendUniformMatrix4fv("modelViewMatrix", glm::value_ptr(viewMatrix));

	glBindVertexArray(s_PlaneVAO);
	glDrawArrays(GL_TRIANGLES, 0, 3 * g_Plane.getNumTriangles());
	glBindVertexArray(0);

	s_pShader->disable();

	//glPopAttrib();		// deprecated from OpenGL 3.0
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

	// display textures

	//glMatrixMode(GL_PROJECTION);
	//glPushMatrix();
	//glLoadIdentity();
	//glMatrixMode(GL_MODELVIEW);
	//glPushMatrix();
	//glLoadIdentity();
	//
	//glEnable(GL_TEXTURE_2D);
	//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	//glm::vec2 offsets[3] = { glm::vec2(-1,-1), glm::vec2(-1, 0), glm::vec2(0, 0) };
	glm::vec2 offsets[3] = { glm::vec2(-0.5f, -0.5f), glm::vec2(-0.5f, 0.5f), glm::vec2(0.5f, 0.5f) };

	for (int ti = 0; ti < 3; ++ti)
	{
		const auto p = offsets[ti];

		glBindTexture(GL_TEXTURE_2D, s_pTexIDs[ti]);

		glm::mat4 M(0.f);
		M[0][0] = M[1][1] = 0.5f;
		M[2][2] = M[3][3] = 1.f;
		M[3][0] = offsets[ti].x;
		M[3][1] = offsets[ti].y;

		s_pTexShader->use();
		s_pTexShader->sendUniform1i("tex", 0);
		s_pTexShader->sendUniformMatrix4fv("projModelViewMatrix", glm::value_ptr(M));

		glBindVertexArray(s_QuadVAO);
		glDrawArrays(GL_TRIANGLES, 0, 6);
		glBindVertexArray(0);

		s_pTexShader->disable();

		glBindTexture(GL_TEXTURE_2D, 0);

		//glBegin(GL_QUADS);
		//glTexCoord2f(0.f, 0.f);		glVertex2f(p.x, p.y);
		//glTexCoord2f(1.f, 0.f);		glVertex2f(p.x + 1.f, p.y);
		//glTexCoord2f(1.f, 1.f);		glVertex2f(p.x + 1.f, p.y + 1.f);
		//glTexCoord2f(0.f, 1.f);		glVertex2f(p.x, p.y + 1.f);
		//glEnd();
	}

	glBindTexture(GL_TEXTURE_2D, 0);

	//glMatrixMode(GL_PROJECTION);
	//glPopMatrix();
	//glMatrixMode(GL_MODELVIEW);
	//glPopMatrix();

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);
}

void Scene03MultipleRenderTarget::Resize(GLFWwindow* window, int w, int h)
{
	AbstractScene::Resize(window, w, h);
	AllocateFBO(w, h);
}

void Scene03MultipleRenderTarget::ImGui()
{
}

void Scene03MultipleRenderTarget::Destroy()
{
	if (s_pTexShader) delete s_pTexShader;
	if (s_pShader) delete s_pShader;

	if (s_pTexIDs)
	{
		glDeleteTextures(3, s_pTexIDs);
		delete[] s_pTexIDs;
	}

	if (s_MRTFBO) glDeleteFramebuffers(1, &s_MRTFBO);
	if (s_DepthRB) glDeleteRenderbuffers(1, &s_DepthRB);

	if (s_MeshVAO) glDeleteVertexArrays(1, &s_MeshVAO);
	if (s_PlaneVAO) glDeleteVertexArrays(1, &s_PlaneVAO);
}

void Scene03MultipleRenderTarget::AllocateFBO(int w, int h)
{
	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

	s_CachedWidth = w;
	s_CachedHeight = h;

	if (!s_MRTFBO) glGenFramebuffers(1, &s_MRTFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, s_MRTFBO);

	if (!s_DepthRB) glGenRenderbuffers(1, &s_DepthRB);
	glBindRenderbuffer(GL_RENDERBUFFER, s_DepthRB);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32, w, h);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, s_DepthRB);

	const GLenum colorAttachements[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

	if (!s_pTexIDs)
	{
		s_pTexIDs = new GLuint[3];
		glGenTextures(3, s_pTexIDs);
	}

	for (int ti = 0; ti < 3; ti++)
	{
		glBindTexture(GL_TEXTURE_2D, s_pTexIDs[ti]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F_ARB, w, h, 0, GL_RGB, GL_FLOAT, 0);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glFramebufferTexture2D(GL_FRAMEBUFFER, colorAttachements[ti], GL_TEXTURE_2D, s_pTexIDs[ti], 0);
	}

	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
	{
		cerr << __FUNCTION__ << ": something wrong with FBO" << endl;
		return;
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glBindRenderbuffer(GL_RENDERBUFFER, 0);
	glBindTexture(GL_TEXTURE_2D, 0);

	//cerr << __FUNCTION__ << ": fbo = " << s_MRTFBO << ", rb = " << s_DepthRB << ", tex[0] = " << s_pTexIDs[0]
	//	<< ", tex[1] = " << s_pTexIDs[1] << ", tex[2] = " << s_pTexIDs[2] << endl;

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);
}
