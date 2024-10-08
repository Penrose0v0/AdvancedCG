#pragma once

#include "AbstractScene.h"
#include <string>
#include <vector>
#include "imgui.h"
#include "arcball_camera.h"
#include "TriMesh.h"

class Scene04PseudoNormal : public AbstractScene
{
public:
	enum Render_Type
	{
		Render_None,
		Render_Texture,
		Render_Pseudo_Normal
	};

	static void Init();
	static void ReloadShaders();
	static void Draw();
	static void Cursor(GLFWwindow* window, double xpos, double ypos);
	static void Mouse(GLFWwindow* window, int button, int action, int mods);
	static void Resize(GLFWwindow* window, int w, int h);
	static void ImGui();
	static void Destroy();

private:
	static GLSLProgramObject* s_pTexShader, * s_pColorShader, * s_pShader;
	static std::string s_VertexShaderFilename, s_FragmentShaderFilename;

	static bool s_RenderWireframe;
	static int s_RenderType;

	static TriMesh s_TriMesh;

	static glm::vec2 s_PrevMouse;
	static ArcballCamera s_Camera;

	static GLuint s_TexVAO, s_NormalVAO;
};
