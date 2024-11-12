#pragma once

#include "glm/glm.hpp"

struct MaterialParameters
{
	MaterialParameters() : shininess(32.f), fresnel0(1,1,1), diffuseCoeff(0.7f), ambient(0.f) {}

	float shininess;
	glm::vec3 fresnel0;
	glm::vec3 diffuseCoeff;
	glm::vec3 ambient;
};