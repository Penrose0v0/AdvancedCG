#pragma once

#include "AbstractScene.h"
#include <string>
#include "imgui.h"

class Scene01Checker2D : public AbstractScene
{
public:
	static void Init();
	static void ReloadShaders();
	static void Draw();
	static void Cursor(GLFWwindow* window, double xpos, double ypos) {}
	static void Mouse(GLFWwindow* window, int button, int action, int mods) {}
	static void Resize(GLFWwindow* window, int w, int h);
	static void ImGui();
	static void Destroy();

private:
	static GLSLProgramObject* s_pShader;
	static std::string s_VertexShaderFilename, s_FragmentShaderFilename;

	static ImVec4 s_CheckerColor0, s_CheckerColor1;
	static ImVec2 s_CheckerScale;

	static GLuint s_VBO, s_VAO;
};
