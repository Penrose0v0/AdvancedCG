#include "BlinnPhongRenderer.h"
#include <iostream>
#include "imgui.h"

using namespace std;

GLSLProgramObject* BlinnPhongRenderer::s_pShader = nullptr;
float BlinnPhongRenderer::s_Shininess = 64.f;
glm::vec3 BlinnPhongRenderer::s_DiffuseCoeff(0.4f, 0.8f, 0.4f);
glm::vec3 BlinnPhongRenderer::s_AmbientColor(0.1f, 0.1f, 0.1f);

void BlinnPhongRenderer::Init()
{
	if (!s_pShader) s_pShader = new GLSLProgramObject();

	s_pShader->attachShaderCodeString(
		R"(#version 120

			varying vec3 ePosition;
			varying vec3 eNormal;

			void main(void)
			{
				eNormal = mat3(gl_ModelViewMatrixInverseTranspose) * gl_Normal;
				ePosition = (gl_ModelViewMatrix * gl_Vertex).xyz;
				gl_Position = gl_ProjectionMatrix * vec4(ePosition, 1);
			}
		)", GL_VERTEX_SHADER);

	s_pShader->attachShaderCodeString(
		R"(#version 120

			varying vec3 ePosition;
			varying vec3 eNormal;

			uniform float shininess;
			//uniform vec3 eLightDir;
			uniform vec3 lightColor;
			uniform vec3 diffuseCoeff;
			uniform vec3 ambient;

			void main(void)
			{
				vec3 viewingDir = normalize(-ePosition);
				vec3 eLightDir = viewingDir;

				vec3 eNormalNormalized = normalize(eNormal);
				float dotDiffuse = max(dot(eNormalNormalized, eLightDir), 0.0);
				vec3 diffuseColor = dotDiffuse * lightColor * diffuseCoeff;

				vec3 specularColor = vec3(0, 0, 0);
				if (dotDiffuse > 0.0)
				{
					vec3 halfVec = normalize(eLightDir + viewingDir);
					float dotSpecular = max(dot(eNormalNormalized, halfVec), 0.0);
					specularColor = pow(dotSpecular, shininess) * lightColor;
				}

				gl_FragColor = vec4(specularColor + diffuseColor + ambient, 1);
			}
		)", GL_FRAGMENT_SHADER);

	s_pShader->link();

	if (!s_pShader->linkSucceeded())
	{
		cerr << __FUNCTION__ << ": blinn-phong shader" << endl;
		s_pShader->printProgramLog();
	}
}

void BlinnPhongRenderer::ImGui()
{
	ImGui::Text("Blinn-Phong Renderer Menu:");

	ImGui::SliderFloat("Shininess", &s_Shininess, 0.001, 100.f);
	ImGui::ColorEdit3("Diffuse Coeff", glm::value_ptr(s_DiffuseCoeff));
	ImGui::ColorEdit3("Ambient Color", glm::value_ptr(s_AmbientColor));
}

void BlinnPhongRenderer::render(const PolygonMesh& mesh) const
{
	if (!s_pShader || mesh.getVertices().empty() || mesh.getVertexNormals().empty())
		return;

	s_pShader->use();
	s_pShader->sendUniform3fv("lightColor", glm::value_ptr(glm::vec3(1.f)));
	s_pShader->sendUniform1f("shininess", s_Shininess);
	s_pShader->sendUniform3fv("diffuseCoeff", glm::value_ptr(s_DiffuseCoeff));
	s_pShader->sendUniform3fv("ambient", glm::value_ptr(s_AmbientColor));

	mesh.renderMesh();

	s_pShader->disable();
}

