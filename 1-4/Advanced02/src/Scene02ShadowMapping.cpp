#include "Scene02ShadowMapping.h"
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

GLSLProgramObject* Scene02ShadowMapping::s_pTexShader = 0;
GLSLProgramObject* Scene02ShadowMapping::s_pShadow1stPassShader = 0;
GLSLProgramObject* Scene02ShadowMapping::s_pShadow2ndPassShader = 0;
GLSLProgramObject* Scene02ShadowMapping::s_pSoftShadow2ndPassShader = 0;

vector<GLSLProgramObject*> Scene02ShadowMapping::s_Shaders;
vector<string> Scene02ShadowMapping::s_VertexShaderFilenames;
vector<string> Scene02ShadowMapping::s_FragmentShaderFilenames;

bool Scene02ShadowMapping::s_UseSoftShadow = false;
bool Scene02ShadowMapping::s_DisplayShadowMap = false;

int Scene02ShadowMapping::s_ShadowTexSize = 256;

GLuint Scene02ShadowMapping::s_MeshVAO = 0;
GLuint Scene02ShadowMapping::s_PlaneVAO = 0;
GLuint Scene02ShadowMapping::s_QuadVAO = 0;

GLuint Scene02ShadowMapping::s_ShadowFBO = 0;
GLuint Scene02ShadowMapping::s_ShadowTexID = 0;

// global variables defined in main.cpp

extern TriMesh g_TriMesh, g_Plane;
extern MaterialParameters g_Material;
extern ArcballCamera g_Camera;
extern glm::mat4 g_ProjMatrix;

void Scene02ShadowMapping::Init()
{
	ReloadShaders();
	GenShadowMap();
}

void Scene02ShadowMapping::ReloadShaders()
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
	s_pTexShader->attachShaderSourceFile(finder.find("gray_tex.frag").c_str(), GL_FRAGMENT_SHADER);
	s_pTexShader->setAttributeLocation("vertexPosition", vertexPositionLocation);
	s_pTexShader->setAttributeLocation("inTexCoord", inTexCoordLocation);
	s_pTexShader->link();

	if (!s_pTexShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": tex shader link failed" << endl;
		s_pTexShader->printProgramLog();
	}

	if (!s_pShadow1stPassShader) delete s_pShadow1stPassShader;
	s_pShadow1stPassShader = new GLSLProgramObject();
	s_pShadow1stPassShader->attachShaderSourceFile(finder.find("shadow_1st_pass.vert").c_str(), GL_VERTEX_SHADER);
	s_pShadow1stPassShader->attachShaderSourceFile(finder.find("shadow_1st_pass.frag").c_str(), GL_FRAGMENT_SHADER);
	s_pShadow1stPassShader->setAttributeLocation("vertexPosition", vertexPositionLocation);
	s_pShadow1stPassShader->setAttributeLocation("vertexNormal", vertexNormalLocation);
	s_pShadow1stPassShader->link();

	if (!s_pShadow1stPassShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": shadow_1st_pass link failed" << endl;
		s_pShadow1stPassShader->printProgramLog();
	}

	if (!s_pShadow2ndPassShader) delete s_pShadow2ndPassShader;
	s_pShadow2ndPassShader = new GLSLProgramObject();
	s_pShadow2ndPassShader->attachShaderSourceFile(finder.find("shadow_blinn_phong.vert").c_str(), GL_VERTEX_SHADER);
	s_pShadow2ndPassShader->attachShaderSourceFile(finder.find("shadow_blinn_phong.frag").c_str(), GL_FRAGMENT_SHADER);
	s_pShadow2ndPassShader->setAttributeLocation("vertexPosition", vertexPositionLocation);
	s_pShadow2ndPassShader->setAttributeLocation("vertexNormal", vertexNormalLocation);
	s_pShadow2ndPassShader->link();

	if (!s_pShadow2ndPassShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": shadow_blinn_phong link failed" << endl;
		s_pShadow2ndPassShader->printProgramLog();
	}

	if (!s_pSoftShadow2ndPassShader) delete s_pSoftShadow2ndPassShader;
	s_pSoftShadow2ndPassShader = new GLSLProgramObject();
	s_pSoftShadow2ndPassShader->attachShaderSourceFile(finder.find("pcf_shadow_blinn_phong.vert").c_str(), GL_VERTEX_SHADER);
	s_pSoftShadow2ndPassShader->attachShaderSourceFile(finder.find("pcf_shadow_blinn_phong.frag").c_str(), GL_FRAGMENT_SHADER);
	s_pSoftShadow2ndPassShader->setAttributeLocation("vertexPosition", vertexPositionLocation);
	s_pSoftShadow2ndPassShader->setAttributeLocation("vertexNormal", vertexNormalLocation);
	s_pSoftShadow2ndPassShader->link();

	if (!s_pSoftShadow2ndPassShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": pcf_shadow_blinn_phong link failed" << endl;
		s_pSoftShadow2ndPassShader->printProgramLog();
	}

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

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

	if (!s_PlaneVAO) glGenVertexArrays(1, &s_PlaneVAO);
	glBindVertexArray(s_PlaneVAO);

	glBindBuffer(GL_ARRAY_BUFFER, g_Plane.getVertexVBO());
	glEnableVertexAttribArray(vertexPositionLocation);
	glVertexAttribPointer(vertexPositionLocation, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	glBindBuffer(GL_ARRAY_BUFFER, g_Plane.getVertexNormalVBO());
	glEnableVertexAttribArray(vertexNormalLocation);
	glVertexAttribPointer(vertexNormalLocation, 3, GL_FLOAT, GL_FALSE, 0, (const void*)0);

	// build a screen-sized quad

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

void Scene02ShadowMapping::Draw()
{
	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

	if (!s_PlaneVAO || !s_MeshVAO || !s_ShadowFBO || !s_ShadowTexID)
		return;

	auto viewMatrix = g_Camera.transform();
	auto modelViewMatrix = viewMatrix * g_TriMesh.getModelMatrix();
	auto projModelViewMatrix = g_ProjMatrix * modelViewMatrix;
	glm::mat3 modelViewInverseTransposed = glm::transpose(glm::inverse(glm::mat3(modelViewMatrix)));

	auto& color = DirectionalLightManager::GetLightColors()[0];
	auto& dir = DirectionalLightManager::GetLightDirs()[0];

	const float cosPhi = cos(dir.x), sinPhi = sin(dir.x);
	const float cosTheta = cos(dir.y), sinTheta = sin(dir.y);
	const auto wLightDir = glm::vec3(sinPhi * sinTheta, cosTheta, cosPhi * sinTheta);
	const auto eLightDir = glm::mat3(viewMatrix) * wLightDir;

	const float boundingSphereRadius = 1.f / g_TriMesh.getModelMatrix()[0][0];
	glm::mat4 shadowProjModelView;
	CalcLightViewTransform(wLightDir, boundingSphereRadius, shadowProjModelView);

	// 1st pass: render geometry from the light's viewpoint

	//glPushAttrib(GL_ENABLE_BIT | GL_VIEWPORT_BIT);	// deprecated from OpenGL 3.0
	glEnable(GL_DEPTH_TEST);

	glBindFramebuffer(GL_FRAMEBUFFER, s_ShadowFBO);
	glDrawBuffer(GL_NONE);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(6.0f, 6.0f);
	glViewport(0, 0, s_ShadowTexSize, s_ShadowTexSize);

	glClear(GL_DEPTH_BUFFER_BIT);

	s_pShadow1stPassShader->use();
	s_pShadow1stPassShader->sendUniformMatrix4fv("projModelViewMatrix", glm::value_ptr(shadowProjModelView * g_TriMesh.getModelMatrix()));

	glBindVertexArray(s_MeshVAO);
	glDrawArrays(GL_TRIANGLES, 0, 3 * g_TriMesh.getNumTriangles());
	glBindVertexArray(0);

	s_pShadow1stPassShader->sendUniformMatrix4fv("projModelViewMatrix", glm::value_ptr(shadowProjModelView));

	glBindVertexArray(s_PlaneVAO);
	glDrawArrays(GL_TRIANGLES, 0, 3 * g_Plane.getNumTriangles());
	glBindVertexArray(0);

	s_pShadow1stPassShader->disable();

	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	//glPopAttrib();	// deprecated from OpenGL 3.0

	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	glDisable(GL_POLYGON_OFFSET_FILL);
	glDrawBuffer(GL_BACK);

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);


	// 2nd pass

	glViewport(0, 0, s_WindowWidth, s_WindowHeight);

	glm::mat4 biasMatrix(0.5f, 0.f, 0.f, 0.f,
		0.f, 0.5f, 0.f, 0.f,
		0.f, 0.f, 0.5f, 0.f,
		0.5f, 0.5f, 0.5f, 1.f);

	//glPushAttrib(GL_ENABLE_BIT);		// deprecated from OpenGL 3.0

	//glActiveTexture(GL_TEXTURE0);
	//glEnable(GL_TEXTURE_2D);		// deprecated from OpenGL 3.0
	glBindTexture(GL_TEXTURE_2D, s_ShadowTexID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);

	if (s_UseSoftShadow)
	{
		s_pSoftShadow2ndPassShader->use();
		s_pSoftShadow2ndPassShader->sendUniformMatrix4fv("projMatrix", glm::value_ptr(g_ProjMatrix));
		s_pSoftShadow2ndPassShader->sendUniformMatrix4fv("modelViewMatrix", glm::value_ptr(modelViewMatrix));
		// TODO: uncomment these lines
		s_pSoftShadow2ndPassShader->sendUniformMatrix3fv("modelViewInverseTransposed", glm::value_ptr(modelViewInverseTransposed));
		s_pSoftShadow2ndPassShader->sendUniform1i("shadowTex", 0);
		s_pSoftShadow2ndPassShader->sendUniform2f("texMapScale", 1.f / float(s_ShadowTexSize), 1.f / float(s_ShadowTexSize));
		s_pSoftShadow2ndPassShader->sendUniformMatrix4fv("biasedShadowProjModelView", glm::value_ptr(biasMatrix * shadowProjModelView));
		s_pSoftShadow2ndPassShader->sendUniform3fv("eLightDir", glm::value_ptr(eLightDir));
		s_pSoftShadow2ndPassShader->sendUniform3fv("lightColor", glm::value_ptr(color));
		s_pSoftShadow2ndPassShader->sendUniform1f("shininess", g_Material.shininess);
		s_pSoftShadow2ndPassShader->sendUniform3fv("diffuseCoeff", glm::value_ptr(g_Material.diffuseCoeff));
		s_pSoftShadow2ndPassShader->sendUniform3fv("ambient", glm::value_ptr(g_Material.ambient));
	}
	else
	{
		s_pShadow2ndPassShader->use();
		s_pShadow2ndPassShader->sendUniformMatrix4fv("projMatrix", glm::value_ptr(g_ProjMatrix));
		s_pShadow2ndPassShader->sendUniformMatrix4fv("modelViewMatrix", glm::value_ptr(modelViewMatrix));
		// TODO: uncomment these lines
		s_pShadow2ndPassShader->sendUniformMatrix3fv("modelViewInverseTransposed", glm::value_ptr(modelViewInverseTransposed));
		s_pShadow2ndPassShader->sendUniform1i("shadowTex", 0);
		s_pShadow2ndPassShader->sendUniformMatrix4fv("biasedShadowProjModelView", glm::value_ptr(biasMatrix* shadowProjModelView));
		s_pShadow2ndPassShader->sendUniform3fv("eLightDir", glm::value_ptr(eLightDir));
		s_pShadow2ndPassShader->sendUniform3fv("lightColor", glm::value_ptr(color));
		s_pShadow2ndPassShader->sendUniform1f("shininess", g_Material.shininess);
		s_pShadow2ndPassShader->sendUniform3fv("diffuseCoeff", glm::value_ptr(g_Material.diffuseCoeff));
		s_pShadow2ndPassShader->sendUniform3fv("ambient", glm::value_ptr(g_Material.ambient));
	}

	glBindVertexArray(s_MeshVAO);
	glDrawArrays(GL_TRIANGLES, 0, 3 * g_TriMesh.getNumTriangles());
	glBindVertexArray(0);

	// TODO: uncomment these lines
	if (s_UseSoftShadow)
		s_pSoftShadow2ndPassShader->sendUniformMatrix4fv("modelViewMatrix", glm::value_ptr(viewMatrix));
	else
		s_pShadow2ndPassShader->sendUniformMatrix4fv("modelViewMatrix", glm::value_ptr(viewMatrix));

	glBindVertexArray(s_PlaneVAO);
	glDrawArrays(GL_TRIANGLES, 0, 3 * g_Plane.getNumTriangles());
	glBindVertexArray(0);

	if (s_UseSoftShadow)
		s_pSoftShadow2ndPassShader->disable();
	else
		s_pShadow2ndPassShader->disable();

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
	glBindTexture(GL_TEXTURE_2D, 0);
	//glDisable(GL_TEXTURE_2D);		// deprecated from OpenGL 3.0

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

	DirectionalLightManager::Draw(1);

	// render depth texture

	if (s_DisplayShadowMap)
	{
		glDisable(GL_DEPTH_TEST);

		const int texDisplaySize = 256;
		const int offset = 10;

		glm::mat4 M(0.f);
		M[0][0] = texDisplaySize / float(s_WindowWidth);
		M[1][1] = texDisplaySize / float(s_WindowHeight);
		M[2][2] = M[3][3] = 1.f;
		M[3][0] = ((s_WindowHeight - offset) + (s_WindowHeight - offset - texDisplaySize)) / float(s_WindowWidth) - 1.f;
		M[3][1] = ((s_WindowHeight - offset) + (s_WindowHeight - offset - texDisplaySize)) / float(s_WindowHeight) - 1.f;

		glBindTexture(GL_TEXTURE_2D, s_ShadowTexID);

		s_pTexShader->use();
		s_pTexShader->sendUniform1i("tex", 0);
		s_pTexShader->sendUniformMatrix4fv("projModelViewMatrix", glm::value_ptr(M));

		glBindVertexArray(s_QuadVAO);
		glDrawArrays(GL_TRIANGLES, 0, 6);
		glBindVertexArray(0);

		s_pTexShader->disable();

		glBindTexture(GL_TEXTURE_2D, 0);

		//glm::mat4 ortho = glm::ortho(0.f, float(s_WindowWidth), 0.f, float(s_WindowHeight));
		//
		//glMatrixMode(GL_PROJECTION);
		//glPushMatrix();
		//glLoadMatrixf(glm::value_ptr(ortho));
		//glMatrixMode(GL_MODELVIEW);
		//glPushMatrix();
		//glLoadIdentity();
		//
		////glEnable(GL_TEXTURE_2D);	// deprecated from OpenGL 3.0
		//glBindTexture(GL_TEXTURE_2D, s_ShadowTexID);
		//
		//glColor3f(1, 1, 1);
		//glBegin(GL_QUADS);
		//glTexCoord2f(0.f, 0.f);		glVertex2f(s_WindowWidth - 10 - texDisplaySize, s_WindowHeight - 10 - texDisplaySize);
		//glTexCoord2f(1.f, 0.f);		glVertex2f(s_WindowWidth - 10, s_WindowHeight - 10 - texDisplaySize);
		//glTexCoord2f(1.f, 1.f);		glVertex2f(s_WindowWidth - 10, s_WindowHeight - 10);
		//glTexCoord2f(0.f, 1.f);		glVertex2f(s_WindowWidth - 10 - texDisplaySize, s_WindowHeight - 10);
		//glEnd();
		//
		//glBindTexture(GL_TEXTURE_2D, 0);
		//
		//glMatrixMode(GL_PROJECTION);
		//glPopMatrix();
		//glMatrixMode(GL_MODELVIEW);
		//glPopMatrix();
	}

	//glPopAttrib();	// deprecated from OpenGL 3.0

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);
}

void Scene02ShadowMapping::Resize(GLFWwindow* window, int w, int h)
{
	if (h == 0) h = 1;

	AbstractScene::Resize(window, w, h);
}

void Scene02ShadowMapping::ImGui()
{
	ImGui::Text("Scene02ShadowMapping Menu:");

	auto& color = DirectionalLightManager::GetLightColors()[0];
	auto& dir = DirectionalLightManager::GetLightDirs()[0];

	bool lightUpdated = false;

	lightUpdated |= ImGui::ColorEdit3("Light Color", (float*)& color.x);
	lightUpdated |= ImGui::SliderFloat("Light Phi", &dir.x, 0.0, 2.f * glm::pi<float>());
	lightUpdated |= ImGui::SliderFloat("Light Theta", &dir.y, 0.f, glm::pi<float>() / 2.25f);

	if (lightUpdated)
		DirectionalLightManager::UpdateLightVAO();

	if (ImGui::SliderInt("Shadow Map Size", &s_ShadowTexSize, 32, 1024))
		GenShadowMap();

	ImGui::Checkbox("Use Soft Shadow", &s_UseSoftShadow);
	ImGui::Checkbox("Display Shadow Map", &s_DisplayShadowMap);

	if (ImGui::Button("Reload Shaders"))
	{
		ReloadShaders();
	}
}

void Scene02ShadowMapping::Destroy()
{
	for (auto p : s_Shaders) delete p;
	s_Shaders.clear();
	if (s_MeshVAO) glDeleteVertexArrays(1, &s_MeshVAO);
	if (s_PlaneVAO) glDeleteVertexArrays(1, &s_PlaneVAO);
}

void Scene02ShadowMapping::GenShadowMap()
{
	CheckGLError(__FUNCTION__, __FILE__, __LINE__);

	if (!s_ShadowFBO) glGenFramebuffers(1, &s_ShadowFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, s_ShadowFBO);

	if (!s_ShadowTexID) glGenTextures(1, &s_ShadowTexID);
	glBindTexture(GL_TEXTURE_2D, s_ShadowTexID);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, s_ShadowTexSize, s_ShadowTexSize, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

	GLfloat border_color[4] = { 1, 1, 1, 1 };
	glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, border_color);

	glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, s_ShadowTexID, 0);

	// Always check that our framebuffer is ok
	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
	{
		cerr << __FUNCTION__ << ": depth fbo not set" << endl;
		return;
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glBindTexture(GL_TEXTURE_2D, 0);

	CheckGLError(__FUNCTION__, __FILE__, __LINE__);
}

void Scene02ShadowMapping::CalcLightViewTransform(const glm::vec3& lightDir, const float boundingSphereRadius, glm::mat4& shadowProjModelViewMatrix)
{
	// adjust the up vector if lightDir.y is close to one
	glm::vec3 up = (lightDir.y > 0.8f) ? glm::vec3(0, 0, 1) : glm::vec3(0, 1, 0);
	glm::mat4 modelView = glm::lookAt(lightDir, glm::vec3(0.f), up);

	glm::mat4 proj = glm::ortho(-1.f, 1.f, -1.f, 1.f, 0.01f, 2.f);

	shadowProjModelViewMatrix = proj * modelView;
}
