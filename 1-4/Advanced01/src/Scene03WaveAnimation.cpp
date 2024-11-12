#include "Scene03WaveAnimation.h"
#include <iostream>
#include "PathFinder.h"
#include "Image2OGLTexture.h"

using namespace std;

string Scene03WaveAnimation::s_VertexShaderFilename = "scene03_wave_animation.vert";
string Scene03WaveAnimation::s_FragmentShaderFilename = "scene03_wave_animation.frag";
GLSLProgramObject* Scene03WaveAnimation::s_pShader = 0;

int Scene03WaveAnimation::s_NumQuadLineVertices = 0;

glm::vec2 Scene03WaveAnimation::s_PrevMouse = glm::vec2(0,0);
ArcballCamera Scene03WaveAnimation::s_Camera(glm::vec3(2.f), glm::vec3(0.f), glm::vec3(0, 1, 0));

GLuint Scene03WaveAnimation::s_VBO = 0;
GLuint Scene03WaveAnimation::s_VAO = 0;

extern glm::mat4 g_ProjMatrix;

void Scene03WaveAnimation::Init()
{
	ReloadShaders();
}

void Scene03WaveAnimation::ReloadShaders()
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

	s_pShader->attachShaderSourceFile(finder.find(s_VertexShaderFilename).c_str(), GL_VERTEX_SHADER);
	s_pShader->attachShaderSourceFile(finder.find(s_FragmentShaderFilename).c_str(), GL_FRAGMENT_SHADER);
	s_pShader->setAttributeLocation("vertexPosition", vertexPositionLocation);
	s_pShader->link();

	if (!s_pShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": shader link failed" << endl;
		s_pShader->printProgramLog();
		return;
	}

	// build a screen-sized quad

	vector<glm::vec3> quads;
	const int nTilesHalf = 500;
	const float tileSize = 0.05f;

	for (int yi = -nTilesHalf; yi <= nTilesHalf; ++yi)
	{
		const float y0 = yi * tileSize;
		const float y1 = y0 + tileSize;

		for (int xi = -nTilesHalf; xi <= nTilesHalf; ++xi)
		{
			const float x0 = xi * tileSize;
			const float x1 = x0 + tileSize;

			quads.emplace_back(y0, 0.f, x0);
			quads.emplace_back(y0, 0.f, x1);

			quads.emplace_back(y0, 0.f, x1);
			quads.emplace_back(y1, 0.f, x1);

			quads.emplace_back(y1, 0.f, x1);
			quads.emplace_back(y1, 0.f, x0);

			quads.emplace_back(y1, 0.f, x0);
			quads.emplace_back(y0, 0.f, x0);
		}
	}

	if (!s_VAO) glGenVertexArrays(1, &s_VAO);
	glBindVertexArray(s_VAO);

	if (!s_VBO) glGenBuffers(1, &s_VBO);
	glBindBuffer(GL_ARRAY_BUFFER, s_VBO);

	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * quads.size(), &quads[0].x, GL_STATIC_DRAW);

	glEnableVertexAttribArray(vertexPositionLocation);
	glVertexAttribPointer(vertexPositionLocation, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0); // vertices

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	s_NumQuadLineVertices = (int)quads.size();
}

void Scene03WaveAnimation::Draw()
{
	if (!s_pShader || !s_pShader->linkSucceeded())
		return;

	//const float aspect = s_WindowWidth / float(s_WindowHeight);
	//auto projMatrix = glm::perspective(45.f, aspect, 0.01f, 100.f);
	auto modelViewMatrix = s_Camera.transform();
	//auto projModelViewMatrix = projMatrix * modelViewMatrix;
	auto projModelViewMatrix = g_ProjMatrix * modelViewMatrix;
	
	s_pShader->use();
	s_pShader->sendUniformMatrix4fv("projModelViewMatrix", glm::value_ptr(projModelViewMatrix));
	// TODO: uncomment these lines
	s_pShader->sendUniform1f("temporalSignal", ImGui::GetTime() / 5.f);
	s_pShader->sendUniform4f("lineColor", 1, 1, 1, 1);

	glBindVertexArray(s_VAO);
	glDrawArrays(GL_LINES, 0, s_NumQuadLineVertices);
	glBindVertexArray(0);

	s_pShader->disable();
}

void Scene03WaveAnimation::Cursor(GLFWwindow* window, double xpos, double ypos)
{
	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE &&
		glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_RELEASE &&
		glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_RELEASE)
	{
		return;
	}

	glm::vec2 currPos(xpos / s_WindowWidth, (s_WindowHeight - ypos) / s_WindowHeight);

	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
	{
		s_Camera.rotate(s_PrevMouse, currPos);
	}
	else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS)
	{
		s_Camera.zoom(currPos.y - s_PrevMouse.y);
	}
	else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
	{
		s_Camera.pan(currPos - s_PrevMouse);
	}

	s_PrevMouse = currPos;
}

void Scene03WaveAnimation::Mouse(GLFWwindow* window, int button, int action, int mods)
{
	if (action == GLFW_PRESS)
	{
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);
		s_PrevMouse = glm::vec2(xpos / s_WindowWidth, (s_WindowHeight - ypos) / s_WindowHeight);
	}
}

void Scene03WaveAnimation::Resize(GLFWwindow* window, int w, int h)
{
	AbstractScene::Resize(window, w, h);
}

void Scene03WaveAnimation::ImGui()
{
	ImGui::Text("Scene03WaveAnimation Menu:");

	if (ImGui::Button("Reload Shaders"))
	{
		ReloadShaders();
	}
}

void Scene03WaveAnimation::Destroy()
{
	if (s_pShader) delete s_pShader;
	if (s_VAO) glDeleteVertexArrays(1, &s_VAO);
	if (s_VBO) glDeleteBuffers(1, &s_VBO);
}
