#include "Scene02ImageSmoothing.h"
#include <iostream>
#include "PathFinder.h"
#include "Image2OGLTexture.h"

using namespace std;

string Scene02ImageSmoothing::s_VertexShaderFilename = "scene02_image_smoothing.vert";
string Scene02ImageSmoothing::s_FragmentShaderFilename = "scene02_image_smoothing.frag";
GLSLProgramObject* Scene02ImageSmoothing::s_pShader = 0;

int Scene02ImageSmoothing::s_TexWidth = 0;
int Scene02ImageSmoothing::s_TexHeight = 0;
GLuint Scene02ImageSmoothing::s_TexID = 0;

int Scene02ImageSmoothing::s_HalfKernelSize = 1;

GLuint Scene02ImageSmoothing::s_VBO = 0;
GLuint Scene02ImageSmoothing::s_VAO = 0;

void Scene02ImageSmoothing::Init()
{
	ReloadShaders();

	PathFinder finder;
	finder.addSearchPath("Resources");
	finder.addSearchPath("../Resources");
	finder.addSearchPath("../../Resources");

	Image2OGLTexture(finder.find("kitten.jpg").c_str(), s_TexID, s_TexWidth, s_TexHeight, GL_TEXTURE_2D, GL_LINEAR, GL_MIRRORED_REPEAT);
}

void Scene02ImageSmoothing::ReloadShaders()
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

	const GLuint vertexPositionLocation = 0;
	const GLuint inTexCoordLocation = 1;

	s_pShader->attachShaderSourceFile(finder.find(s_VertexShaderFilename).c_str(), GL_VERTEX_SHADER);
	s_pShader->attachShaderSourceFile(finder.find(s_FragmentShaderFilename).c_str(), GL_FRAGMENT_SHADER);
	s_pShader->setAttributeLocation("vertexPosition", vertexPositionLocation);
	s_pShader->setAttributeLocation("inTexCoord", inTexCoordLocation);
	s_pShader->link();

	if (!s_pShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": shader link failed" << endl;
		s_pShader->printProgramLog();
		return;
	}

	// build a screen-sized quad

	float quadVertices[] = { -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1};
	float quadTexCoords[] = { 0,  0, 1,  0, 1, 1,  0,  0, 1, 1,  0, 1};

	if (!s_VAO) glGenVertexArrays(1, &s_VAO);
	glBindVertexArray(s_VAO);

	if (!s_VBO) glGenBuffers(1, &s_VBO);
	glBindBuffer(GL_ARRAY_BUFFER, s_VBO);

	glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices) + sizeof(quadTexCoords), NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(quadVertices), quadVertices);
	glBufferSubData(GL_ARRAY_BUFFER, sizeof(quadVertices), sizeof(quadTexCoords), quadTexCoords);

	glEnableVertexAttribArray(vertexPositionLocation);
	glVertexAttribPointer(vertexPositionLocation, 2, GL_FLOAT, GL_FALSE, 0, (const void*)0); // vertices

	glEnableVertexAttribArray(inTexCoordLocation);
	glVertexAttribPointer(inTexCoordLocation, 2, GL_FLOAT, GL_FALSE, 0, (const void*)sizeof(quadVertices)); // texture coordinates

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
}

void Scene02ImageSmoothing::Draw()
{
	if (!s_pShader || !s_pShader->linkSucceeded())
		return;

	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_DEPTH_TEST);

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, s_TexID);

	s_pShader->use();
	// TODO: uncomment these lines
	//s_pShader->sendUniform1i("tex", 0);
	//s_pShader->sendUniform1i("halfKernelSize", s_HalfKernelSize);
	//s_pShader->sendUniform1f("uScale", 1.f / float(s_TexWidth));
	//s_pShader->sendUniform1f("vScale", 1.f / float(s_TexHeight));

	glBindVertexArray(s_VAO);
	glDrawArrays(GL_TRIANGLES, 0, 6);
	glBindVertexArray(0);

	s_pShader->disable();

	glDisable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);

	glPopAttrib();
}

void Scene02ImageSmoothing::Resize(GLFWwindow* window, int w, int h)
{
	AbstractScene::Resize(window, w, h);
}

void Scene02ImageSmoothing::ImGui()
{
	ImGui::Text("Scene02ImageSmoothing Menu:");

	ImGui::SliderInt("Half Kernel Size", &s_HalfKernelSize, 0, 100);

	if (ImGui::Button("Reload Shaders"))
	{
		ReloadShaders();
	}
}

void Scene02ImageSmoothing::Destroy()
{
	if (s_pShader) delete s_pShader;
	if (s_TexID) glDeleteTextures(1, &s_TexID);
	if (s_VAO) glDeleteVertexArrays(1, &s_VAO);
	if (s_VBO) glDeleteBuffers(1, &s_VBO);
}
