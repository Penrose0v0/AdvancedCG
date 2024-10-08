#pragma once

#include "AbstractScene.h"
#include <string>
#include <vector>
#include "imgui.h"

class Scene01ShadingExamples : public AbstractScene
{
public:
	enum Shader_Types
	{
		Lambert_Shader,
		Phong_Shader,
		Blinn_Phong_Shader,
		Cook_Torrance_Shader,
		Num_Shaders
	};

	static void Init();
	static void ReloadShaders();
	static void Draw();
	static void Resize(GLFWwindow* window, int w, int h);
	static void ImGui();
	static void Destroy();

private:
	static std::vector<GLSLProgramObject*> s_Shaders;
	static std::vector<std::string> s_VertexShaderFilenames, s_FragmentShaderFilenames;

	static int s_SelectedShaderIndex;

	static GLuint s_MeshVAO, s_PlaneVAO;
};
