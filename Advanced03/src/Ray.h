#pragma once

#include <ostream>
#include "glm/glm.hpp"
#include "glm/ext.hpp"	// for glm::to_string()

class Ray
{
public:
	Ray()
	{
	}

	Ray(const Ray& r)
		: m_Origin(r.m_Origin), m_UnitDir(r.m_UnitDir)
	{
	}

	Ray(const glm::vec3 &o, const glm::vec3 &d)
		: m_Origin(o), m_UnitDir(d)
	{
	}

	~Ray()
	{
	}

	inline glm::vec3 getOrigin() const { return m_Origin; }
	inline glm::vec3 getUnitDir() const { return m_UnitDir; }

	inline glm::vec3 calculatePosition(float t) const { return m_Origin + t*m_UnitDir; }

	inline void setOrigin(const glm::vec3 &o) { m_Origin = o; }
	inline void setUnitDir(const glm::vec3 &d) { m_UnitDir = d; }
	inline void setDirWithNormalization(const glm::vec3 &d)
	{
		m_UnitDir = glm::normalize(d);
	}

	friend std::ostream &operator<<(std::ostream &os, const Ray &r)
	{
		os << glm::to_string(r.getOrigin()) << " + t " << glm::to_string(r.getUnitDir());
		return os;
	}


private:
	glm::vec3 m_Origin;
	glm::vec3 m_UnitDir;	// should be a unit vector
};

