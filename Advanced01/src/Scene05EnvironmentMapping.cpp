#include "Scene05EnvironmentMapping.h"
#include <iostream>
#include <cmath>
#include "PathFinder.h"
#include "Image2OGLTexture.h"
#include "CheckGLError.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

string Scene05EnvironmentMapping::s_VertexShaderFilename = "scene05_envmap.vert";
string Scene05EnvironmentMapping::s_FragmentShaderFilename = "scene05_envmap.frag";
GLSLProgramObject* Scene05EnvironmentMapping::s_pTexShader = 0;
GLSLProgramObject* Scene05EnvironmentMapping::s_pShader = 0;

TriMesh Scene05EnvironmentMapping::s_TriMesh;

glm::vec2 Scene05EnvironmentMapping::s_PrevMouse = glm::vec2(0,0);
ArcballCamera Scene05EnvironmentMapping::s_Camera(glm::vec3(2.f), glm::vec3(0.f), glm::vec3(0, 1, 0));

int Scene05EnvironmentMapping::s_NumSkyDomeVertices = 0;

GLuint Scene05EnvironmentMapping::s_VAO = 0;
GLuint Scene05EnvironmentMapping::s_EnvMapTexID = 0;
GLuint Scene05EnvironmentMapping::s_SkyDomeVAO = 0;

extern glm::mat4 g_ProjMatrix;

void Scene05EnvironmentMapping::Init()
{
	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

	PathFinder finder;
	finder.addSearchPath("Resources");
	finder.addSearchPath("../Resources");
	finder.addSearchPath("../../Resources");

	s_TriMesh.loadObj(finder.find("sphere_642verts.obj").c_str());
	//s_TriMesh.loadObj(finder.find("bunny10k.obj").c_str());

	{
		int w, h;
		Image2OGLTexture(finder.find("venice_sunrise_2k.jpg").c_str(), s_EnvMapTexID, w, h, GL_TEXTURE_2D, GL_LINEAR, GL_MIRRORED_REPEAT);
		//Image2OGLTexture(finder.find("kasuga_envmap.jpg").c_str(), s_EnvMapTexID, w, h, GL_TEXTURE_2D, GL_LINEAR, GL_MIRRORED_REPEAT);
	}

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

	BuildSkyDome();

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

	ReloadShaders();

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);
}

void Scene05EnvironmentMapping::ReloadShaders()
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

	if (!s_VAO) glGenVertexArrays(1, &s_VAO);
	glBindVertexArray(s_VAO);

	glBindBuffer(GL_ARRAY_BUFFER, s_TriMesh.getVertexVBO());
	glEnableVertexAttribArray(vertexPositionLocation);
	glVertexAttribPointer(vertexPositionLocation, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	glBindBuffer(GL_ARRAY_BUFFER, s_TriMesh.getVertexNormalVBO());
	glEnableVertexAttribArray(vertexNormalLocation);
	glVertexAttribPointer(vertexNormalLocation, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void Scene05EnvironmentMapping::Draw()
{
	if (!s_pShader || !s_pShader->linkSucceeded())
		return;

	//glPushAttrib(GL_ENABLE_BIT);	// deprecated from OpenGL 3.0
	glEnable(GL_DEPTH_TEST);

	//const float aspect = s_WindowWidth / float(s_WindowHeight);
	//auto projMatrix = glm::perspective(45.f, aspect, 0.01f, 100.f);
	auto modelViewMatrix = s_Camera.transform() * s_TriMesh.getModelMatrix();
	//auto projModelViewMatrix = projMatrix * modelViewMatrix;
	auto projModelViewMatrix = g_ProjMatrix * modelViewMatrix;
	
	{
		//auto skyDomeMatrix = s_Camera.transform();
		//skyDomeMatrix[3][0] = skyDomeMatrix[3][1] = skyDomeMatrix[3][2] = 0.f;
		//
		//glMatrixMode(GL_PROJECTION);
		//glLoadMatrixf(glm::value_ptr(projMatrix));
		//glMatrixMode(GL_MODELVIEW);
		//glLoadMatrixf(glm::value_ptr(skyDomeMatrix));

		RenderSkyDome();
	}

	glm::vec3 eye = s_Camera.eye();

	//glEnable(GL_TEXTURE_2D);	// deprecated from OpenGL 3.0
	glBindTexture(GL_TEXTURE_2D, s_EnvMapTexID);

	s_pShader->use();
	s_pShader->sendUniformMatrix4fv("projModelViewMatrix", glm::value_ptr(projModelViewMatrix));
	// TODO: uncomment these lines
	//s_pShader->sendUniform3fv("eye", glm::value_ptr(eye));
	//s_pShader->sendUniform1ui("envmap", 0);

	glBindVertexArray(s_VAO);
	glDrawArrays(GL_TRIANGLES, 0, 3 * s_TriMesh.getNumTriangles());
	glBindVertexArray(0);

	s_pShader->disable();

	//glDisable(GL_TEXTURE_2D);		// deprecated from OpenGL 3.0
	glBindTexture(GL_TEXTURE_2D, 0);

	//glPopAttrib();	// deprecated from OpenGL 3.0
}

void Scene05EnvironmentMapping::Cursor(GLFWwindow* window, double xpos, double ypos)
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

void Scene05EnvironmentMapping::Mouse(GLFWwindow* window, int button, int action, int mods)
{
	if (action == GLFW_PRESS)
	{
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);
		s_PrevMouse = glm::vec2(xpos / s_WindowWidth, (s_WindowHeight - ypos) / s_WindowHeight);
	}
}

void Scene05EnvironmentMapping::Resize(GLFWwindow* window, int w, int h)
{
	AbstractScene::Resize(window, w, h);
}

void Scene05EnvironmentMapping::ImGui()
{
	ImGui::Text("Scene05EnvironmentMapping Menu:");

	if (ImGui::Button("Reload Shaders"))
	{
		ReloadShaders();
	}
}

void Scene05EnvironmentMapping::Destroy()
{
	if (s_pShader) delete s_pShader;
	if (s_VAO) glDeleteVertexArrays(1, &s_VAO);
	if (s_EnvMapTexID) glDeleteTextures(1, &s_EnvMapTexID);
	if (s_SkyDomeVAO) glDeleteVertexArrays(1, &s_SkyDomeVAO);
}

void Scene05EnvironmentMapping::BuildSkyDome()
{
	const int nDivSphereLatitude = 50;
	const int nDivSphereLongitude = 50;

	const float radius = 50.f;
	const float dRadLat = M_PI / (float)nDivSphereLatitude;
	const float dRadLon = 2 * M_PI / (float)nDivSphereLongitude;

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
			vertices.emplace_back(radius * sinf(lat0) * cosf(lon1), radius * cosf(lat0), radius * sinf(lat0) * sinf(lon1));
			vertices.emplace_back(radius * sinf(lat1) * cosf(lon0), radius * cosf(lat1), radius * sinf(lat1) * sinf(lon0));
			vertices.emplace_back(radius * sinf(lat1) * cosf(lon1), radius * cosf(lat1), radius * sinf(lat1) * sinf(lon1));

			texCoords.emplace_back(lon1 / (2 * M_PI), 1.f - lat0 / M_PI);
			texCoords.emplace_back(lon0 / (2 * M_PI), 1.f - lat0 / M_PI);
			texCoords.emplace_back(lon0 / (2 * M_PI), 1.f - lat1 / M_PI);
			texCoords.emplace_back(lon1 / (2 * M_PI), 1.f - lat0 / M_PI);
			texCoords.emplace_back(lon0 / (2 * M_PI), 1.f - lat1 / M_PI);
			texCoords.emplace_back(lon1 / (2 * M_PI), 1.f - lat1 / M_PI);

			radLongitude = lon1;
		}

		radLatitude = lat1;
	}

	if (!s_SkyDomeVAO) glGenVertexArrays(1, &s_SkyDomeVAO);
	glBindVertexArray(s_SkyDomeVAO);

	GLuint skyDomeVBO;
	glGenBuffers(1, &skyDomeVBO);
	glBindBuffer(GL_ARRAY_BUFFER, skyDomeVBO);

	size_t vertexSize = sizeof(glm::vec3) * vertices.size();
	size_t texCoordSize = sizeof(glm::vec2) * texCoords.size();

	glBufferData(GL_ARRAY_BUFFER, vertexSize + texCoordSize, NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, vertexSize, &vertices[0].x);
	glBufferSubData(GL_ARRAY_BUFFER, vertexSize, texCoordSize, &texCoords[0].x);

	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0); // vertices

	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, (const void*)vertexSize); // texture coordinates

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glDeleteBuffers(1, &skyDomeVBO);

	s_NumSkyDomeVertices = (int)vertices.size();
}

void Scene05EnvironmentMapping::RenderSkyDome()
{
	if (!s_SkyDomeVAO || !s_EnvMapTexID || !s_pTexShader)
	{
		//cerr << __FUNCTION__ << ": vbo = " << s_SkyDomeVBO << ", tex ID = " << s_EnvMapTexID << endl;
		return;
	}

	auto skyDomeMatrix = s_Camera.transform();
	skyDomeMatrix[3][0] = skyDomeMatrix[3][1] = skyDomeMatrix[3][2] = 0.f;

	glBindTexture(GL_TEXTURE_2D, s_EnvMapTexID);

	s_pTexShader->use();
	s_pTexShader->sendUniform1i("tex", 0);
	s_pTexShader->sendUniformMatrix4fv("projModelViewMatrix", glm::value_ptr(g_ProjMatrix * skyDomeMatrix));

	glBindVertexArray(s_SkyDomeVAO);
	glDrawArrays(GL_TRIANGLES, 0, s_NumSkyDomeVertices);
	glBindVertexArray(0);

	s_pTexShader->disable();

	glBindTexture(GL_TEXTURE_2D, 0);

#if 0
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, s_EnvMapTexID);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, s_SkyDomeVAO);
	glVertexPointer(3, GL_FLOAT, 0, 0);
	glActiveTexture(GL_TEXTURE0);
	glTexCoordPointer(2, GL_FLOAT, 0, (GLvoid*)(sizeof(glm::vec3) * s_NumSkyDomeVertices));

	glDrawArrays(GL_TRIANGLES, 0, s_NumSkyDomeVertices);

	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);

	glDisable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);
#endif
}
