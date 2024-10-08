#pragma once

#include "GLSLProgramObject.h"
#include "TriMesh.h"

class DirectionalLightManager
{
public:
	static void Init();
	static void GenRandomLights(int nLights);
	static void UpdateLightVAO();
	static void Draw(int nLights = -1);

	static int GetNumLights() { return (int)s_Colors.size(); }
	static std::vector<glm::vec3> &GetLightColors() { return s_Colors; }
	static std::vector<glm::vec2> &GetLightDirs() { return s_Dirs; }

private:
	static GLSLProgramObject* s_pShader;
	static TriMesh s_ArrowMesh;
	static GLuint s_VAO;
	static std::vector<glm::vec3> s_Colors;
	static std::vector<glm::vec2> s_Dirs;

	static void BuildArrowMesh();
};