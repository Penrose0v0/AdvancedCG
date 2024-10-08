#include "Scene01Checker2D.h"
#include "PathFinder.h"
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#include <iostream>

using namespace std;

string Scene01Checker2D::s_VertexShaderFilename = "scene01_checker2D.vert";
string Scene01Checker2D::s_FragmentShaderFilename = "scene01_checker2D.frag";
GLSLProgramObject* Scene01Checker2D::s_pShader = 0;

ImVec4 Scene01Checker2D::s_CheckerColor0 = ImVec4(0, 0, 0, 0);
ImVec4 Scene01Checker2D::s_CheckerColor1 =
ImVec4(79.f / 255.f, 172.f / 255.f, 135.f / 255.f, 0);
ImVec2 Scene01Checker2D::s_CheckerScale = ImVec2(0.2f, 0.2f);

GLuint Scene01Checker2D::s_VBO = 0;
GLuint Scene01Checker2D::s_VAO = 0;

void Scene01Checker2D::Init() { ReloadShaders(); }

void Scene01Checker2D::ReloadShaders() {
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

	s_pShader->attachShaderSourceFile(finder.find(s_VertexShaderFilename).c_str(), GL_VERTEX_SHADER_ARB);
	s_pShader->attachShaderSourceFile(finder.find(s_FragmentShaderFilename).c_str(), GL_FRAGMENT_SHADER_ARB);
	s_pShader->setAttributeLocation("vertexPosition", vertexPositionLocation);
	s_pShader->setAttributeLocation("inTexCoord", inTexCoordLocation);
	s_pShader->link();

	if (!s_pShader->linkSucceeded()) {
		cerr << __FUNCTION__ << ": shader link failed" << endl;
		s_pShader->printProgramLog();
		return;
	}

	// build a screen-sized quad

	float quadVertices[] = { -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1 };
	float quadTexCoords[] = { 0,  0, 1,  0, 1, 1,  0,  0, 1, 1,  0, 1 };

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

void Scene01Checker2D::Draw() {
	if (!s_pShader || !s_pShader->linkSucceeded())
		return;

	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_DEPTH_TEST);

	s_pShader->use();
	// TODO: uncomment these lines
	 s_pShader->sendUniform4f("checkerColor0", s_CheckerColor0.x,
	 s_CheckerColor0.y, s_CheckerColor0.z, s_CheckerColor0.w);
	 s_pShader->sendUniform4f("checkerColor1", s_CheckerColor1.x,
	 s_CheckerColor1.y, s_CheckerColor1.z, s_CheckerColor1.w);
	 s_pShader->sendUniform2f("checkerScale", s_CheckerScale.x,
	 s_CheckerScale.y);

	glBindVertexArray(s_VAO);
	glDrawArrays(GL_TRIANGLES, 0, 6);
	glBindVertexArray(0);

	s_pShader->disable();

	glPopAttrib();
}

void Scene01Checker2D::Resize(GLFWwindow* window, int w, int h) {
	AbstractScene::Resize(window, w, h);
}

void Scene01Checker2D::ImGui() {
	ImGui::Text("Scene01Checker2D Menu:");

	ImGui::ColorEdit3("Checker Color #0", (float*)& s_CheckerColor0.x);
	ImGui::ColorEdit3("Checker Color #1", (float*)& s_CheckerColor1.x);

	ImGui::SliderFloat("Horizontal Scale", &s_CheckerScale.x, 0.05f, 1.0f);
	ImGui::SliderFloat("Vertical Scale", &s_CheckerScale.y, 0.05f, 1.0f);

	if (ImGui::Button("Reload Shaders")) {
		ReloadShaders();
	}
}

void Scene01Checker2D::Destroy() {
	if (s_pShader) delete s_pShader;
	if (s_VAO) glDeleteVertexArrays(1, &s_VAO);
	if (s_VBO) glDeleteBuffers(1, &s_VBO);
}
