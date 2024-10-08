#include "Scene01ShadingExamples.h"
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

vector<GLSLProgramObject*> Scene01ShadingExamples::s_Shaders;
vector<string> Scene01ShadingExamples::s_VertexShaderFilenames;
vector<string> Scene01ShadingExamples::s_FragmentShaderFilenames;

int Scene01ShadingExamples::s_SelectedShaderIndex = Scene01ShadingExamples::Lambert_Shader;

GLuint Scene01ShadingExamples::s_MeshVAO = 0;
GLuint Scene01ShadingExamples::s_PlaneVAO = 0;

// global variables defined in main.cpp

extern TriMesh g_TriMesh, g_Plane;
extern MaterialParameters g_Material;
extern ArcballCamera g_Camera;
extern glm::mat4 g_ProjMatrix;

void Scene01ShadingExamples::Init()
{
	ReloadShaders();
}

void Scene01ShadingExamples::ReloadShaders()
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

	s_VertexShaderFilenames = {
		"lambert.vert",
		"phong.vert",
		"blinn_phong.vert",
		"cook_torrance.vert"
	};
	s_FragmentShaderFilenames = {
		"lambert.frag",
		"phong.frag",
		"blinn_phong.frag",
		"cook_torrance.frag"
	};

	const GLuint vertexPositionLocation = 0;
	const GLuint vertexNormalLocation = 1;

	for (auto p : s_Shaders) delete p;
	s_Shaders.clear();
	s_Shaders.resize(Num_Shaders);

	for (int si = 0; si < Num_Shaders; ++si)
	{
		s_Shaders[si] = new GLSLProgramObject();
		s_Shaders[si]->attachShaderSourceFile(finder.find(s_VertexShaderFilenames[si]).c_str(), GL_VERTEX_SHADER);
		s_Shaders[si]->attachShaderSourceFile(finder.find(s_FragmentShaderFilenames[si]).c_str(), GL_FRAGMENT_SHADER);
		s_Shaders[si]->setAttributeLocation("vertexPosition", vertexPositionLocation);
		s_Shaders[si]->setAttributeLocation("vertexNormal", vertexNormalLocation);
		s_Shaders[si]->link();

		if (!s_Shaders[si]->linkSucceeded())
		{
			cerr << __FUNCTION__ << ": shader[" << si << "] link failed" << endl;
			s_Shaders[si]->printProgramLog();
		}
	}

	if (!g_TriMesh.getVertexVBO() || !g_TriMesh.getVertexNormalVBO())
	{
		cerr << __FUNCTION__ << ": mesh VBOs not ready" << endl;
		return;
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

	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);
}

void Scene01ShadingExamples::Draw()
{
	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

	auto pShader = s_Shaders[s_SelectedShaderIndex];

	if (!pShader || !pShader->linkSucceeded())
		return;

	//glPushAttrib(GL_ENABLE_BIT);	// deprecated from OpenGL 3.0
	glEnable(GL_DEPTH_TEST);

	auto viewMatrix = g_Camera.transform();
	auto modelViewMatrix = viewMatrix * g_TriMesh.getModelMatrix();
	auto projModelViewMatrix = g_ProjMatrix * modelViewMatrix;
	glm::mat3 modelViewInverseTransposed = glm::transpose(glm::inverse(glm::mat3(modelViewMatrix)));

	auto& color = DirectionalLightManager::GetLightColors()[0];
	auto& dir = DirectionalLightManager::GetLightDirs()[0];

	const float cosPhi = cos(dir.x), sinPhi = sin(dir.x);
	const float cosTheta = cos(dir.y), sinTheta = sin(dir.y);
	const auto eLightDir = glm::mat3(viewMatrix) * glm::vec3(sinPhi * sinTheta, cosTheta, cosPhi * sinTheta);

	if (s_SelectedShaderIndex == Lambert_Shader)
	{
		pShader->use();
		pShader->sendUniformMatrix4fv("projModelViewMatrix", glm::value_ptr(projModelViewMatrix));
		pShader->sendUniformMatrix3fv("modelViewInverseTransposed", glm::value_ptr(modelViewInverseTransposed));
		pShader->sendUniform3fv("eLightDir", glm::value_ptr(eLightDir));
		pShader->sendUniform3fv("lightColor", glm::value_ptr(color));
		pShader->sendUniform3fv("diffuseCoeff", glm::value_ptr(g_Material.diffuseCoeff));
		pShader->sendUniform3fv("ambient", glm::value_ptr(g_Material.ambient));

		glBindVertexArray(s_MeshVAO);
		glDrawArrays(GL_TRIANGLES, 0, 3 * g_TriMesh.getNumTriangles());
		glBindVertexArray(0);

		pShader->sendUniformMatrix4fv("projModelViewMatrix", glm::value_ptr(g_ProjMatrix * viewMatrix));
		pShader->sendUniformMatrix3fv("modelViewInverseTransposed", glm::value_ptr(glm::mat3(viewMatrix)));

		glBindVertexArray(s_PlaneVAO);
		glDrawArrays(GL_TRIANGLES, 0, 3 * g_Plane.getNumTriangles());
		glBindVertexArray(0);

		pShader->disable();
	}
	else if (s_SelectedShaderIndex == Phong_Shader || s_SelectedShaderIndex == Blinn_Phong_Shader)
	{
		pShader->use();
		pShader->sendUniformMatrix4fv("projMatrix", glm::value_ptr(g_ProjMatrix));
		pShader->sendUniformMatrix4fv("modelViewMatrix", glm::value_ptr(modelViewMatrix));
		// TODO: uncomment these lines
		//pShader->sendUniformMatrix3fv("modelViewInverseTransposed", glm::value_ptr(modelViewInverseTransposed));
		//pShader->sendUniform3fv("eLightDir", glm::value_ptr(eLightDir));
		//pShader->sendUniform3fv("lightColor", glm::value_ptr(color));
		//pShader->sendUniform1f("shininess", g_Material.shininess);
		//pShader->sendUniform3fv("diffuseCoeff", glm::value_ptr(g_Material.diffuseCoeff));
		//pShader->sendUniform3fv("ambient", glm::value_ptr(g_Material.ambient));

		glBindVertexArray(s_MeshVAO);
		glDrawArrays(GL_TRIANGLES, 0, 3 * g_TriMesh.getNumTriangles());
		glBindVertexArray(0);

		// TODO: uncomment these lines
		//pShader->sendUniformMatrix4fv("modelViewMatrix", glm::value_ptr(viewMatrix));
		//pShader->sendUniformMatrix3fv("modelViewInverseTransposed", glm::value_ptr(glm::mat3(viewMatrix)));

		glBindVertexArray(s_PlaneVAO);
		glDrawArrays(GL_TRIANGLES, 0, 3 * g_Plane.getNumTriangles());
		glBindVertexArray(0);

		pShader->disable();
	}
	else if (s_SelectedShaderIndex == Cook_Torrance_Shader)
	{
		pShader->use();
		pShader->sendUniformMatrix4fv("projMatrix", glm::value_ptr(g_ProjMatrix));
		pShader->sendUniformMatrix4fv("modelViewMatrix", glm::value_ptr(modelViewMatrix));
		pShader->sendUniformMatrix3fv("modelViewInverseTransposed", glm::value_ptr(modelViewInverseTransposed));
		pShader->sendUniform3fv("eLightDir", glm::value_ptr(eLightDir));
		pShader->sendUniform3fv("lightColor", glm::value_ptr(color));
		pShader->sendUniform1f("roughness", 1.f / g_Material.shininess);
		pShader->sendUniform3fv("fresnel0", glm::value_ptr(g_Material.fresnel0));
		pShader->sendUniform3fv("diffuseCoeff", glm::value_ptr(g_Material.diffuseCoeff));
		pShader->sendUniform3fv("ambient", glm::value_ptr(g_Material.ambient));

		glBindVertexArray(s_MeshVAO);
		glDrawArrays(GL_TRIANGLES, 0, 3 * g_TriMesh.getNumTriangles());
		glBindVertexArray(0);

		pShader->sendUniformMatrix4fv("modelViewMatrix", glm::value_ptr(viewMatrix));
		pShader->sendUniformMatrix3fv("modelViewInverseTransposed", glm::value_ptr(glm::mat3(viewMatrix)));

		glBindVertexArray(s_PlaneVAO);
		glDrawArrays(GL_TRIANGLES, 0, 3 * g_Plane.getNumTriangles());
		glBindVertexArray(0);

		pShader->disable();
	}

	DirectionalLightManager::Draw(1);

	//glPopAttrib();	// deprecated from OpenGL 3.0

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);
}

void Scene01ShadingExamples::Resize(GLFWwindow* window, int w, int h)
{
	AbstractScene::Resize(window, w, h);
}

void Scene01ShadingExamples::ImGui()
{
	ImGui::Text("Scene01ShadingExamples Menu:");

	auto &color = DirectionalLightManager::GetLightColors()[0];
	auto &dir = DirectionalLightManager::GetLightDirs()[0];

	bool lightUpdated = false;

	lightUpdated |= ImGui::ColorEdit3("Light Color", (float*)& color.x);
	lightUpdated |= ImGui::SliderFloat("Light Phi", &dir.x, 0.0, 2.f * glm::pi<float>());
	lightUpdated |= ImGui::SliderFloat("Light Theta", &dir.y, 0.f, glm::pi<float>() / 2.25f);

	if (lightUpdated)
		DirectionalLightManager::UpdateLightVAO();

	const char* reflectionModels[] = { "Lambert", "Phong", "Blinn-Phong", "Cook-Torrance" };
	ImGui::ListBox("Reflection Model", &s_SelectedShaderIndex, reflectionModels, sizeof(reflectionModels) / sizeof(const char*));

	if (s_SelectedShaderIndex == Cook_Torrance_Shader || s_SelectedShaderIndex == Phong_Shader || s_SelectedShaderIndex == Blinn_Phong_Shader)
	{
		if (s_SelectedShaderIndex == Cook_Torrance_Shader)
		{
			ImGui::ColorEdit3("Fresnel 0", (float*)& g_Material.fresnel0.x);
		}

		ImGui::SliderFloat("Shininess", &g_Material.shininess, 0.001, 100.f);
	}

	ImGui::ColorEdit3("Diffuse Coeff", (float*)&g_Material.diffuseCoeff.x);
	ImGui::ColorEdit3("Ambient Color", (float*)&g_Material.ambient.x);

	if (ImGui::Button("Reload Shaders"))
	{
		ReloadShaders();
	}
}

void Scene01ShadingExamples::Destroy()
{
	for (auto p : s_Shaders) delete p;
	s_Shaders.clear();
	if (s_MeshVAO) glDeleteVertexArrays(1, &s_MeshVAO);
	if (s_PlaneVAO) glDeleteVertexArrays(1, &s_PlaneVAO);
}
