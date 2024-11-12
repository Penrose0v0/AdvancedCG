#pragma once

#include "Material.h"

class PseudoNormalColorMaterial : public Material
{
private:
	typedef glm::vec3 vec3;
	typedef float Real;

	PseudoNormalColorMaterial()
	{
	}

	PseudoNormalColorMaterial(const PseudoNormalColorMaterial &m)
	{
	}

public:
	static PseudoNormalColorMaterial *CreateMaterial()
	{
		PseudoNormalColorMaterial *m = new PseudoNormalColorMaterial();
		material_cache.push_back( m );
		return m;
	}

	static PseudoNormalColorMaterial *CreateMaterial(const PseudoNormalColorMaterial *_m)
	{
		PseudoNormalColorMaterial *m = new PseudoNormalColorMaterial(*_m);
		material_cache.push_back( m );
		return m;
	}

	~PseudoNormalColorMaterial()
	{
	}

	Material_Type getMaterialType() const { return Pseudo_Normal_Color_Type; }

	void clone(Material **m) const { *m = new PseudoNormalColorMaterial(*this); }

};
