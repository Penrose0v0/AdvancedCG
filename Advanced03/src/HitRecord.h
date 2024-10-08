#pragma once

#include "glm/glm.hpp"

class Material;

struct HitRecord
{
	float m_ParamT;	// ray parameter (used like o + t * d)
	glm::vec3 m_Normal;	// surface normal
	glm::vec3 m_HitPos;	// hit position p (= o + t + d)
	glm::vec3 m_TexCoords;	// texture coordinate
	Material *m_pMaterial;
};

