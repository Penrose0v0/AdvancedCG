#pragma once

#include "AbstractScene.h"
#include <string>
#include <vector>
#include "imgui.h"
#include "arcball_camera.h"

class Scene03WaveAnimation : public AbstractScene
{
public:
	static void Init();
	static void ReloadShaders();
	static void Draw();
	static void Cursor(GLFWwindow* window, double xpos, double ypos);
	static void Mouse(GLFWwindow* window, int button, int action, int mods);
	static void Resize(GLFWwindow* window, int w, int h);
	static void ImGui();
	static void Destroy();

private:
	static GLSLProgramObject* s_pShader;
	static std::string s_VertexShaderFilename, s_FragmentShaderFilename;

	static int s_NumQuadLineVertices;

	static glm::vec2 s_PrevMouse;
	static ArcballCamera s_Camera;

	static GLuint s_VBO, s_VAO;
};
