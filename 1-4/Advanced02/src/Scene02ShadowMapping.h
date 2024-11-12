#pragma once

#include "AbstractScene.h"
#include <string>
#include <vector>
#include "imgui.h"

class Scene02ShadowMapping : public AbstractScene
{
public:
	static void Init();
	static void ReloadShaders();
	static void Draw();
	static void Resize(GLFWwindow* window, int w, int h);
	static void ImGui();
	static void Destroy();

private:
	static GLSLProgramObject* s_pTexShader, * s_pShadow1stPassShader, *s_pShadow2ndPassShader, * s_pSoftShadow2ndPassShader;
	static std::vector<GLSLProgramObject*> s_Shaders;
	static std::vector<std::string> s_VertexShaderFilenames, s_FragmentShaderFilenames;

	static bool s_UseSoftShadow, s_DisplayShadowMap;

	static int s_ShadowTexSize;

	static GLuint s_MeshVAO, s_PlaneVAO, s_QuadVAO, s_ShadowFBO, s_ShadowTexID;

	static void GenShadowMap();
	static void CalcLightViewTransform(const glm::vec3 &lightDir, const float boundingSphereRadius, glm::mat4 &shadowProjModelViewMatrix);
};
