#include "Scene04PseudoNormal.h"
#include <iostream>
#include "PathFinder.h"
#include "Image2OGLTexture.h"

using namespace std;

string Scene04PseudoNormal::s_VertexShaderFilename = "scene04_pseudo_normal.vert";
string Scene04PseudoNormal::s_FragmentShaderFilename = "scene04_pseudo_normal.frag";
GLSLProgramObject* Scene04PseudoNormal::s_pTexShader = 0;
GLSLProgramObject* Scene04PseudoNormal::s_pColorShader = 0;
GLSLProgramObject* Scene04PseudoNormal::s_pShader = 0;

bool Scene04PseudoNormal::s_RenderWireframe = false;
int Scene04PseudoNormal::s_RenderType = Render_Pseudo_Normal;

TriMesh Scene04PseudoNormal::s_TriMesh;

glm::vec2 Scene04PseudoNormal::s_PrevMouse = glm::vec2(0,0);
ArcballCamera Scene04PseudoNormal::s_Camera(glm::vec3(2.f), glm::vec3(0.f), glm::vec3(0, 1, 0));

GLuint Scene04PseudoNormal::s_TexVAO = 0;
GLuint Scene04PseudoNormal::s_NormalVAO = 0;

extern glm::mat4 g_ProjMatrix;

void Scene04PseudoNormal::Init()
{
	PathFinder finder;
	finder.addSearchPath("Resources");
	finder.addSearchPath("../Resources");
	finder.addSearchPath("../../Resources");

	s_TriMesh.loadObj(finder.find("duck.obj").c_str());
	s_TriMesh.loadTexture(finder.find("duckCM.jpg").c_str());

	ReloadShaders();
}

void Scene04PseudoNormal::ReloadShaders()
{
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
	const GLuint vertexNormalLocation = 1;

	if (s_pTexShader) delete s_pTexShader;
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
		return;
	}

	if (s_pColorShader) delete s_pColorShader;
	s_pColorShader = new GLSLProgramObject();
	s_pColorShader->attachShaderSourceFile(finder.find("color.vert").c_str(), GL_VERTEX_SHADER);
	s_pColorShader->attachShaderSourceFile(finder.find("color.frag").c_str(), GL_FRAGMENT_SHADER);
	s_pColorShader->setAttributeLocation("vertexPosition", vertexPositionLocation);
	s_pColorShader->link();

	if (!s_pColorShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": color shader link failed" << endl;
		s_pColorShader->printProgramLog();
		return;
	}

	if (s_pShader) delete s_pShader;
	s_pShader = new GLSLProgramObject();
	s_pShader->attachShaderSourceFile(finder.find(s_VertexShaderFilename).c_str(), GL_VERTEX_SHADER);
	s_pShader->attachShaderSourceFile(finder.find(s_FragmentShaderFilename).c_str(), GL_FRAGMENT_SHADER);
	s_pShader->setAttributeLocation("vertexPosition", vertexPositionLocation);
	s_pShader->setAttributeLocation("vertexNormal", vertexNormalLocation);
	s_pShader->link();

	if (!s_pShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": shader link failed" << endl;
		s_pShader->printProgramLog();
		return;
	}

	if (!s_TriMesh.getVertexVBO() || !s_TriMesh.getVertexNormalVBO())
	{
		cerr << __FUNCTION__ << ": mesh VBOs not ready" << endl;
		return;
	}

	if (!s_TexVAO) glGenVertexArrays(1, &s_TexVAO);
	glBindVertexArray(s_TexVAO);

	glBindBuffer(GL_ARRAY_BUFFER, s_TriMesh.getVertexVBO());
	glEnableVertexAttribArray(vertexPositionLocation);
	glVertexAttribPointer(vertexPositionLocation, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0); // vertices

	glBindBuffer(GL_ARRAY_BUFFER, s_TriMesh.getTexCoordVBO());
	glEnableVertexAttribArray(inTexCoordLocation);
	glVertexAttribPointer(inTexCoordLocation, 2, GL_FLOAT, GL_FALSE, 0, (const void*)0); // vertices

	if (!s_NormalVAO) glGenVertexArrays(1, &s_NormalVAO);
	glBindVertexArray(s_NormalVAO);

	glBindBuffer(GL_ARRAY_BUFFER, s_TriMesh.getVertexVBO());
	glEnableVertexAttribArray(vertexPositionLocation);
	glVertexAttribPointer(vertexPositionLocation, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0); // vertices

	glBindBuffer(GL_ARRAY_BUFFER, s_TriMesh.getVertexNormalVBO());
	glEnableVertexAttribArray(vertexNormalLocation);
	glVertexAttribPointer(vertexNormalLocation, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0); // vertices

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void Scene04PseudoNormal::Draw()
{
	if (!s_pShader || !s_pShader->linkSucceeded())
		return;

	glPushAttrib(GL_ENABLE_BIT);
	glEnable(GL_DEPTH_TEST);

	//const float aspect = s_WindowWidth / float(s_WindowHeight);
	//auto projMatrix = glm::perspective(45.f, aspect, 0.01f, 100.f);
	auto modelViewMatrix = s_Camera.transform() * s_TriMesh.getModelMatrix();

	if (s_RenderWireframe)
	{
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0f, 1.0f);
	}

	if (s_RenderType == Render_Texture)
	{
		glBindTexture(GL_TEXTURE_2D, s_TriMesh.getTexID());
		
		s_pTexShader->use();
		s_pTexShader->sendUniform1i("tex", 0);
		//s_pTexShader->sendUniformMatrix4fv("projModelViewMatrix", glm::value_ptr(projMatrix * modelViewMatrix));
		s_pTexShader->sendUniformMatrix4fv("projModelViewMatrix", glm::value_ptr(g_ProjMatrix * modelViewMatrix));

		glBindVertexArray(s_TexVAO);
		glDrawArrays(GL_TRIANGLES, 0, 3 * s_TriMesh.getNumTriangles());
		glBindVertexArray(0);

		s_pTexShader->disable();

		glBindTexture(GL_TEXTURE_2D, 0);

		//glMatrixMode(GL_PROJECTION);
		//glLoadMatrixf(glm::value_ptr(projMatrix));
		//glMatrixMode(GL_MODELVIEW);
		//glLoadMatrixf(glm::value_ptr(modelViewMatrix));
		//
		//s_TriMesh.renderTexturedMesh();
	}
	else if (s_RenderType == Render_Pseudo_Normal)
	{
		glm::mat3 modelViewInvTransposed = glm::transpose(glm::inverse(glm::mat3(modelViewMatrix)));

		s_pShader->use();
		s_pShader->sendUniformMatrix4fv("modelViewMatrix", glm::value_ptr(modelViewMatrix));
		//s_pShader->sendUniformMatrix4fv("projMatrix", glm::value_ptr(projMatrix));
		s_pShader->sendUniformMatrix4fv("projMatrix", glm::value_ptr(g_ProjMatrix));
		// TODO: uncomment these lines
		//s_pShader->sendUniformMatrix3fv("modelViewInvTransposed", glm::value_ptr(modelViewInvTransposed));

		glBindVertexArray(s_NormalVAO);
		glDrawArrays(GL_TRIANGLES, 0, 3 * s_TriMesh.getNumTriangles());
		glBindVertexArray(0);

		s_pShader->disable();
	}

	if (s_RenderWireframe)
	{
		glDisable(GL_POLYGON_OFFSET_FILL);

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

		s_pColorShader->use();
		//s_pColorShader->sendUniformMatrix4fv("projModelViewMatrix", glm::value_ptr(projMatrix * modelViewMatrix));
		s_pColorShader->sendUniformMatrix4fv("projModelViewMatrix", glm::value_ptr(g_ProjMatrix * modelViewMatrix));
		s_pColorShader->sendUniform3f("color", 0.7f, 0.7f, 0.7f);

		glBindVertexArray(s_TexVAO);
		glDrawArrays(GL_TRIANGLES, 0, 3 * s_TriMesh.getNumTriangles());
		glBindVertexArray(0);

		s_pColorShader->disable();

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		//glMatrixMode(GL_PROJECTION);
		//glLoadMatrixf(glm::value_ptr(projMatrix));
		//glMatrixMode(GL_MODELVIEW);
		//glLoadMatrixf(glm::value_ptr(modelViewMatrix));
		//
		//s_TriMesh.renderWireframeMesh();
	}

	glPopAttrib();
}

void Scene04PseudoNormal::Cursor(GLFWwindow* window, double xpos, double ypos)
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

void Scene04PseudoNormal::Mouse(GLFWwindow* window, int button, int action, int mods)
{
	if (action == GLFW_PRESS)
	{
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);
		s_PrevMouse = glm::vec2(xpos / s_WindowWidth, (s_WindowHeight - ypos) / s_WindowHeight);
	}
}

void Scene04PseudoNormal::Resize(GLFWwindow* window, int w, int h)
{
	AbstractScene::Resize(window, w, h);
}

void Scene04PseudoNormal::ImGui()
{
	ImGui::Text("Scene04PseudoNormal Menu:");

	ImGui::Checkbox("Wireframe", &s_RenderWireframe);

	const char* renderTypes[] = { "(None)", "Texture", "Pseudo Normal" };
	ImGui::ListBox("Render Type", &s_RenderType, renderTypes, sizeof(renderTypes)/sizeof(const char*));

	if (ImGui::Button("Reload Shaders"))
	{
		ReloadShaders();
	}
}

void Scene04PseudoNormal::Destroy()
{
	if (s_pTexShader) delete s_pTexShader;
	if (s_pColorShader) delete s_pColorShader;
	if (s_pShader) delete s_pShader;
	if (s_TexVAO) glDeleteVertexArrays(1, &s_TexVAO);
	if (s_NormalVAO) glDeleteVertexArrays(1, &s_NormalVAO);
}
