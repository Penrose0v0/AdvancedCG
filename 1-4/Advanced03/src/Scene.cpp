#include "Scene.h"
#include <GL/glew.h>

using namespace std;

mt19937 Scene::s_RandSrc(12345);
uniform_real_distribution<float> Scene::s_RandDist(0, 1);

bool Scene::loadEnvironmentMap(const char* filename)
{
	auto* pEnv = new EnvironmentMap();

	if (pEnv->load(filename))
	{
		if (m_pEnvironmentMap) delete m_pEnvironmentMap;
		m_pEnvironmentMap = pEnv;
		return true;
	}

	delete m_pEnvironmentMap;
	m_pEnvironmentMap = 0;

	return false;
}

void Scene::drawGL() const
{
	const int nObjects = (int)m_Objects.size();
	if (!nObjects) return;

	if (m_pEnvironmentMap)
		m_pEnvironmentMap->drawGL();

	// colorize with pseudo colors

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0f, 1.0f);
	for (int oi = 0; oi < nObjects; ++oi)
	{
		glColor3fv(glm::value_ptr(m_PseudoColors[oi]));
		m_Objects[oi]->drawGL();
	}

	// draw wireframes

	glDisable(GL_POLYGON_OFFSET_FILL);
	glColor3f(0.3f, 0.3f, 0.3f);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	for (int oi = 0; oi < nObjects; ++oi)
		m_Objects[oi]->drawGL();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}
