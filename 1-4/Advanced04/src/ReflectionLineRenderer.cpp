#include "ReflectionLineRenderer.h"
#include "EnvironmentMap.h"
#include "imgui.h"
#include "arcball_camera.h"
#include "CheckGLError.h"

using namespace std;

GLSLProgramObject* ReflectionLineRenderer::s_pShader = nullptr;
EnvironmentMap ReflectionLineRenderer::s_EnvironmentMap;
bool ReflectionLineRenderer::s_DisplaySkyDome = false;

extern ArcballCamera g_Camera;

void ReflectionLineRenderer::Init()
{
	CHECK_GL_ERROR;

	if (!s_pShader) s_pShader = new GLSLProgramObject();

	s_pShader->attachShaderCodeString(
		R"(#version 120

			varying vec3 vWorldEyeDir;
			varying vec3 vWorldNormal;

			uniform vec3 eye;

			void main(void)
			{
				vWorldNormal = gl_Normal;
				vWorldEyeDir = gl_Vertex.xyz - eye;
				gl_Position = gl_ProjectionMatrix * gl_ModelViewMatrix * gl_Vertex;
			}
		)", GL_VERTEX_SHADER);

	s_pShader->attachShaderCodeString(
		R"(#version 120

			#define PI 3.141592653589793

			varying vec3 vWorldEyeDir;
			varying vec3 vWorldNormal;

			uniform sampler2D envmap;

			float atan2(in float y, in float x)
			{
				return x == 0.0 ? sign(y) * PI / 2 : atan(y, x);
			}

			void main()
			{
				vec3 reflectDir = normalize(reflect(vWorldEyeDir, normalize(vWorldNormal)));
				vec2 texCoord = vec2(0.5 * atan2(-reflectDir.z, -reflectDir.x) / PI + 0.5, asin(reflectDir.y) / PI + 0.5);
				gl_FragColor = texture2D(envmap, texCoord);
			}
		)", GL_FRAGMENT_SHADER);

	s_pShader->link();

	if (!s_pShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": reflection line shader" << endl;
		s_pShader->printProgramLog();
	}

	CHECK_GL_ERROR;

	s_EnvironmentMap.genStripTexture(1024, 0.05f);

	CHECK_GL_ERROR;
}

void ReflectionLineRenderer::ImGui()
{
	ImGui::Text("Reflection Line Renderer Menu:");

	ImGui::Checkbox("Display Sky Dome", &s_DisplaySkyDome);
}

void ReflectionLineRenderer::render(const PolygonMesh& mesh) const
{
	if (!s_pShader || mesh.getVertices().empty() || mesh.getVertexNormals().empty())
		return;

	CHECK_GL_ERROR;

	if (s_DisplaySkyDome)
		s_EnvironmentMap.drawGL();

	glm::vec3 eye = g_Camera.eye();

	glBindTexture(GL_TEXTURE_2D, s_EnvironmentMap.getTexID());

	s_pShader->use();
	s_pShader->sendUniform1i("envmap", 0);
	s_pShader->sendUniform3fv("eye", glm::value_ptr(eye));

	mesh.renderMesh();

	s_pShader->disable();

	glBindTexture(GL_TEXTURE_2D, 0);

	CHECK_GL_ERROR;
}
