#pragma once

#include <vector>
#include <random>
#include "GeometricObject.h"
#include "Ray.h"
#include "EnvironmentMap.h"
//#include "LightSource.h"
//#include "EnvironmentMap.h"

class Scene
{
public:
	Scene()
		: m_pEnvironmentMap(0), m_BackgroundColor(1.f)
	{
	}

	~Scene()
	{
		if (m_pEnvironmentMap) delete m_pEnvironmentMap;
	}

	// objects

	int getNumObjects() const { return (int)m_Objects.size(); }

	const std::vector<GeometricObject*>& getObjects() const { return m_Objects; }
	std::vector<GeometricObject*>& getObjects() { return m_Objects; }

	GeometricObject* getObject(int i) { return m_Objects[i]; }
	void addObject(GeometricObject* o)
	{
		m_Objects.push_back(o);
		m_PseudoColors.emplace_back(0.5f * s_RandDist(s_RandSrc) + 0.5f,
									0.5f * s_RandDist(s_RandSrc) + 0.5f,
									0.5f * s_RandDist(s_RandSrc) + 0.5f);
	}

	bool loadEnvironmentMap(const char* filename);

	glm::vec3 getBackgroundColor(const Ray& r) const
	{
		return (!m_pEnvironmentMap) ? m_BackgroundColor : m_pEnvironmentMap->fetchColor(r);
	}

	void drawGL() const;

	//void clearScene()
	//{
	//	m_Objects.clear();
	//	m_PseudoColors.clear();
	//
	//	//for (int i = 0; i < (int)m_Lights.size(); i++)
	//	//	delete m_Lights[i];
	//	//m_Lights.clear();
	//
	//	//if (m_pEnvironmentMap)
	//	//{
	//	//	delete m_pEnvironmentMap;
	//	//	m_pEnvironmentMap = 0;
	//	//}
	//}

private:
	//std::vector<LightSource*> m_Lights;
	std::vector<GeometricObject*> m_Objects;
	std::vector<glm::vec3> m_PseudoColors;

	EnvironmentMap* m_pEnvironmentMap;
	glm::vec3 m_BackgroundColor;

	static std::mt19937 s_RandSrc;
	static std::uniform_real_distribution<float> s_RandDist;

	//// returns value in [-1, 1]
	//inline float frand() { return 2.f * (((float)rand() / (float)RAND_MAX) - 0.5f); }
};

