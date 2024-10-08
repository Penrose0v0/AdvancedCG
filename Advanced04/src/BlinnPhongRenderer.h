#pragma once

#include "AbstractRenderer.h"

class BlinnPhongRenderer : public AbstractRenderer
{
public:
	static void Init();
	static AbstractRenderer* Create() { return new BlinnPhongRenderer(); }
	static void ImGui();

	void render(const PolygonMesh& mesh) const;

private:
	static GLSLProgramObject* s_pShader;

	static float s_Shininess;
	static glm::vec3 s_DiffuseCoeff;
	static glm::vec3 s_AmbientColor;
};