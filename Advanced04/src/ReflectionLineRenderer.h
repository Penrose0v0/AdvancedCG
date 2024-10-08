#pragma once

#include "AbstractRenderer.h"
#include "EnvironmentMap.h"

class ReflectionLineRenderer : public AbstractRenderer
{
public:
	static void Init();
	static AbstractRenderer* Create() { return new ReflectionLineRenderer(); }
	static void ImGui();

	void render(const PolygonMesh& mesh) const;

private:
	static GLSLProgramObject* s_pShader;
	static EnvironmentMap s_EnvironmentMap;
	static bool s_DisplaySkyDome;
};