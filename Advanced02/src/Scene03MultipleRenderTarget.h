#pragma once

#include "AbstractScene.h"
#include <string>
#include <vector>
#include "imgui.h"

class Scene03MultipleRenderTarget : public AbstractScene
{
public:
	static void Init();
	static void ReloadShaders();
	static void Draw();
	static void Resize(GLFWwindow* window, int w, int h);
	static void ImGui();
	static void Destroy();

private:
	static GLSLProgramObject* s_pTexShader, * s_pShader;

	static GLuint s_MeshVAO, s_PlaneVAO, s_QuadVAO, s_MRTFBO, s_DepthRB;
	static GLuint *s_pTexIDs;

	static int s_CachedWidth, s_CachedHeight;

	static void AllocateFBO(int w, int h);
};
